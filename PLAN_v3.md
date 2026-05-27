# Protein Design Studio — 最终实施方案 v3

> **一句话：加载蛋白质 → 看到3D结构 → 点一个位置 → 工具告诉你突变这里的置信度是多少、类似实验做过没有、怎么做实验验证。**
> 所有计算在你的RTX 3060上本地完成。AI只做"解释"，不做"计算"。

---

## 零、硬件验证

```
✅ RTX 3060 6GB  →  ESM-2 650M (~2.5GB)  ← 仅有的GPU负载
✅ 8GB RAM       →  Web服务 ~200MB + Python ~500MB = 宽裕
✅ 20GB磁盘      →  模型2.5GB + 数据库<100MB + PDB缓存<500MB = 够用
✅ WSL2          →  CUDA直通已验证（nvidia-smi可用）
```

**用户启动方式：**
```bash
git clone <repo>
cd protein-studio
pip install -r requirements.txt
python run.py
# 浏览器打开 http://localhost:8899
```

不需要Docker，不需要数据库，不需要配置。自动检测GPU，没有GPU就用CPU（慢但能跑）。

---

## 一、架构

```
┌─────────────────────────────────────────────────────────┐
│  浏览器 (3Dmol.js)                                       │
│  ┌───────────────────────┐ ┌──────────────────────────┐ │
│  │ 3D查看器              │ │ 右侧面板                  │ │
│  │ • 点残基→选中         │ │ • 残基信息卡片            │ │
│  │ • 多着色方案          │ │ • 突变设计面板            │ │
│  │ • 测量距离            │ │ • AI导师对话              │ │
│  │                       │ │ • 置信度指示器            │ │
│  └───────────────────────┘ └──────────────────────────┘ │
└─────────────────────┬───────────────────────────────────┘
                      │ HTTP/SSE
┌─────────────────────┴───────────────────────────────────┐
│  FastAPI (Python)                    localhost:8899      │
│                                                          │
│  ┌──────────┐  ┌───────────┐  ┌────────────────────┐   │
│  │ PDB解析  │  │ ESM打分    │  │ 置信度引擎          │   │
│  │(Biopython)│  │(本地GPU)  │  │(sklearn,本地CPU)   │   │
│  │          │  │           │  │                    │   │
│  │ 结构特征 │  │ log_likeli│  │ ProTherm校准        │   │
│  │ SASA     │  │ hood差异  │  │ Bootstrap CI       │   │
│  │ 二级结构 │  │           │  │ 相似案例检索        │   │
│  │ 氢键网络 │  │           │  │                    │   │
│  └──────────┘  └───────────┘  └────────────────────┘   │
│                                                          │
│  ┌──────────────────────────────────────────────────┐   │
│  │ AI解释层 (DeepSeek API — 仅在需要时调用)           │   │
│  │ • 残基角色解释                                    │   │
│  │ • 突变效应自然语言解读                             │   │
│  │ • 实验方案生成                                    │   │
│  │ • 回答用户问题                                    │   │
│  └──────────────────────────────────────────────────┘   │
└──────────────────────────────────────────────────────────┘
```

---

## 二、数据流（用户一次点击发生了什么）

```
用户点击残基GLU222
  │
  ├─→ 1. 前端发请求: GET /api/residue/1EMA/A/222
  │
  ├─→ 2. 后端并行计算 (全部本地，<2秒):
  │     ├─ PDB解析: 提取残基属性 (SASA, SS, 氢键, 邻居列表)
  │     ├─ ESM打分: 对该位置所有20种氨基酸打分
  │     └─ 数据库查询: 在ProTherm中搜索类似突变案例
  │
  ├─→ 3. 置信度引擎:
  │     ├─ 查校准表: "表面暴露位, Glu→Ala, 此场景下ESM的RMSE=X"
  │     ├─ 计算95% CI: "ΔΔG = -0.8 ± 1.2 kcal/mol"
  │     └─ 查相似案例: "ProTherm中有23条类似突变, 80%实验验证为稳定化"
  │
  ├─→ 4. AI解释层 (仅在需要详细解释时调用):
  │     ├─ 输入: 计算得分 + 置信度 + 相似案例 + 结构上下文
  │     ├─ 输出: 自然语言解释 + 教学引导
  │     └─ 成本: 约3K-8K tokens/次
  │
  └─→ 5. 返回前端: JSON + 渲染
```

**关键设计：步骤1-3不消耗任何API token。只有步骤4在用户需要解释时才调用DS API。**

---

## 三、精确Token消耗（修正后）

### 用确定性代码替代LLM的部分（零token）

| 原来用LLM做的 | 改为 | 效果 |
|-------------|------|------|
| ProTherm 30K条解析验证 | Python正则+Biopython cross-check | 零token, 10秒, 零幻觉 |
| 突变表示标准化 | 字典映射(3字母↔1字母) | 零token, 零错误 |
| PDB ID有效性验证 | HTTP HEAD请求+格式校验 | 零token |
| 结构环境分类 | DSSP算法+距离计算 | 零token, 100%确定 |
| 重复检测 | hash+集合运算 | 零token |
| 置信区间计算 | sklearn bootstrap | 零token, 真实统计 |

### LLM真正做的事（token用在刀刃上）

| 任务 | Token | 为什么必须用LLM |
|------|-------|----------------|
| 数据库条目自然语言摘要 | ~15M | 30K条×500 tokens, 把"PDB:1LZ1, S26G, ΔΔG=-1.2"变成人能读的段落 |
| GFP家族知识综合 | ~3M | 一篇完整综述+关键突变图谱 |
| 实验协议生成(10套) | ~2M | 引物设计/PCR/表达/纯化/表征全流程 |
| 校准报告解读 | ~2M | 把统计结果翻译成用户能理解的语言 |
| AI导师对话(持续) | ~50M | 日常使用，每次5K-15K |
| **10天内总计** | **~20M** | (协议+GFP+摘要部分+校准报告) |
| **90天总计** | **~70M** | (加上持续AI导师使用) |

**从700M→70M，但这是诚实的数字。幻觉风险降到零，统计严谨性拉满。**

---

## 四、10天可执行冲刺计划

### Day 1-2：骨架 + 确定性数据库管线

```
目标: 工具能启动，能加载PDB，能显示3D结构

上午:
  [ ] 项目结构搭建 (FastAPI + 静态文件服务)
  [ ] PDB解析模块 (Biopython: 提取原子坐标/残基信息/序列)
  [ ] ESM-2加载测试 (验证3060能跑, 记录显存占用)

下午:
  [ ] 3Dmol.js前端骨架 (加载1EMA, cartoon+sticks渲染)
  [ ] 点击残基交互 (选中→高亮→发送请求→显示残基卡片)
  [ ] ProTherm下载+解析脚本 (纯Python, 不用LLM)

晚上:
  [ ] ProTherm解析后台运行 (30K条, 预计10-30分钟CPU)
  [ ] 验证: 解析成功率>95%, 无幻觉
```

### Day 3-4：ESM打分 + 置信度引擎

```
目标: 点击残基后能看到ESM打分+置信区间

上午:
  [ ] ESM打分管线 (输入残基位置→20种aa批量打分→返回log_likelihood)
  [ ] 得分转ΔΔG近似 (pseudo-ΔΔG = -ln(P_mut/P_wt), 后续用ProTherm校准)

下午:
  [ ] 校准系统:
      - 在ProTherM子集上计算ESM预测vs实验ΔΔG的相关性
      - 按结构环境分层(核心/表面/界面)
      - 按突变类型分层(大→小/小→大/极性变化/电荷变化)
      - 每层计算RMSE + Pearson r
      - 存储校准参数到JSON

晚上:
  [ ] Bootstrap置信区间:
      - 对每个预测, 在该场景下重采样1000次
      - 返回95% CI
  [ ] 相似案例检索:
      - 在ProTherm中找"相同蛋白"或"相同结构环境+相同突变类型"
      - 返回统计: 多少例/成功率/平均ΔΔG
```

### Day 5-6：AI解释层

```
目标: AI能解释计算结果，生成实验方案

Day 5:
  [ ] DeepSeek API客户端 (支持OpenRouter + AI-POOL双线路)
  [ ] AI导师系统提示词设计:
      - 角色: 蛋白质设计导师
      - 知识: 当前蛋白信息+选中残基+计算得分+置信度+相似案例
      - 行为: 苏格拉底式提问, 先解释再提问
  [ ] 残基解释生成 (点击残基→AI自动给出一段解释)
  [ ] 前端聊天面板集成

Day 6:
  [ ] 实验协议生成系统:
      - 输入: 突变列表 + 蛋白信息
      - 输出: 引物序列(含Tm/GC%), PCR条件, 表达纯化方案
  [ ] GFP家族知识综合 (LLM生成GFP综述+突变图谱)
      - 这是token消耗的主要部分
      - ~3M tokens, 可分批次后台生成
```

### Day 7-8：前端完善 + 突变设计面板

```
目标: 工具完整可用，能做完整的"选残基→看打分→设计突变→生成方案"流程

Day 7:
  [ ] 突变设计面板:
      - 选残基→选目标氨基酸→显示预测ΔΔG+置信度
      - 颜色编码 (绿: 高置信度有益, 黄: 不确定, 红: 高置信度有害)
      - 排序+筛选
  [ ] 多着色方案 (疏水性/B因子/电荷/保守性)
  [ ] 距离测量工具

Day 8:
  [ ] 结果导出:
      - 突变列表→CSV
      - 实验方案→Markdown/PDF
      - 结构截图
  [ ] 错误处理+加载状态
  [ ] 移动端响应式适配 (基础)
```

### Day 9-10：打磨 + Demo

```
Day 9:
  [ ] 端到端测试 (加载1EMA→点击GLU222→ESM打分→置信度→AI解释→设计突变→生成引物)
  [ ] Bug修复
  [ ] 性能优化 (ESM推理缓存, 减少重复计算)

Day 10:
  [ ] README + 使用文档 (中英文)
  [ ] 录制Demo (从零启动到完成一个设计)
  [ ] 开源准备 (LICENSE, CONTRIBUTING)
```

---

## 五、项目文件结构

```
protein-studio/                      # 仓库根目录
├── README.md                        # 项目说明+快速开始
├── run.py                           # 一键启动
├── requirements.txt
│
├── backend/
│   ├── __init__.py
│   ├── main.py                      # FastAPI app
│   ├── config.py                    # 配置(模型路径/API密钥/端口)
│   │
│   ├── structure/
│   │   ├── __init__.py
│   │   ├── pdb_parser.py           # PDB解析 (Biopython)
│   │   └── features.py             # 结构特征提取 (SASA/SS/氢键)
│   │
│   ├── scoring/
│   │   ├── __init__.py
│   │   ├── esm_scorer.py           # ESM-2打分
│   │   └── consensus.py            # 多方法共识打分(未来扩展)
│   │
│   ├── confidence/
│   │   ├── __init__.py
│   │   ├── calibrator.py           # 校准模型 (sklearn)
│   │   ├── database.py             # ProTherm/SKEMPI查询
│   │   └── bootstrap.py            # Bootstrap置信区间
│   │
│   ├── ai/
│   │   ├── __init__.py
│   │   ├── client.py               # DeepSeek API客户端
│   │   ├── tutor.py                # AI导师对话管理
│   │   ├── protocol.py             # 实验协议生成
│   │   └── prompts.py              # 提示词模板
│   │
│   └── data/                        # 运行时数据
│       ├── protherm_raw.txt         # 原始ProTherm数据
│       ├── protherm_parsed.json     # 解析后
│       ├── calibration.json         # 校准参数
│       └── pdbs/                    # 缓存的PDB文件
│
├── frontend/
│   ├── index.html
│   ├── css/
│   │   └── style.css
│   └── js/
│       ├── app.js                   # 主逻辑
│       ├── viewer.js               # 3Dmol.js查看器
│       ├── panel.js                # 残基信息+突变面板
│       └── chat.js                 # AI对话
│
└── docs/
    ├── ARCHITECTURE.md
    └── API.md
```

---

## 六、API设计

```
GET  /api/health                     → {"status": "ok", "gpu": true}

# 结构
GET  /api/structure/{pdb_id}         → 结构摘要 (序列/链/残基列表)
GET  /api/structure/{pdb_id}/load    → 加载PDB(下载+解析+缓存)
GET  /api/residue/{pdb_id}/{chain}/{resnum}
                                     → 单残基详情 (属性+环境+邻居)

# 打分
POST /api/score/{pdb_id}/{chain}/{resnum}
  Body: {mutations: ["ALA","ARG",...]}
  Response: {
    wild_type: "GLU",
    scores: [
      {aa: "ALA", pseudo_ddg: -0.8, confidence_95ci: [-2.1, 0.5]},
      ...
    ],
    similar_cases: [
      {pdb: "1LZ1", mutation: "E26A", ddg_exp: -0.3, method: "CD"},
      ...
    ],
    calibration: {
      scenario: "surface_exposed_charged_to_hydrophobic",
      pearson_r: 0.68,
      rmse: 1.25,
      n_samples: 342
    }
  }

# AI
POST /api/chat
  Body: {message: "...", context: {pdb_id, chain, resnum}}
  Response: SSE流 (AI回复)

POST /api/protocol
  Body: {mutations: [...], protein_info: {...}}
  Response: {primers: [...], pcr_conditions: {...}, ...}

# 导出
GET  /api/export/csv                 → 突变列表CSV
GET  /api/export/protocol            → 实验方案Markdown
```

---

## 七、关键技术细节

### ESM-2加载策略

```python
# 启动时加载一次，常驻显存
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.cuda().eval()  # ~2.5GB VRAM

# 批量打分 (一次API调用打分20种aa)
def score_all_mutations(sequence, position):
    """对某位置所有20种氨基酸打分"""
    variants = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        if aa == sequence[position]:
            continue
        mut_seq = sequence[:position] + aa + sequence[position+1:]
        variants.append(mut_seq)
    
    # 批量推理 (~1秒)
    with torch.no_grad():
        results = model(batch_tokens)
    
    return {aa: score for aa, score in zip(aas, scores)}
```

### 置信度计算

```python
# 不是LLM说"可靠"就算可靠，是真实的统计
def predict_with_confidence(esm_score, structural_context, mutation_type):
    scenario = f"{structural_context}_{mutation_type}"
    calib = calibration_db[scenario]  # 从ProTherm预计算的
    
    # 回归校准 (ESM的log_likelihood不是物理ΔΔG)
    calibrated_ddg = calib.slope * esm_score + calib.intercept
    
    # Bootstrap 95% CI (真实统计，不是LLM编的)
    ci_low, ci_high = bootstrap_ci(
        esm_score, calib.residuals, n_samples=1000
    )
    
    # 相似案例
    similar = protherm_db.query(
        struct_context=structural_context,
        mutation_type=mutation_type,
        limit=20
    )
    
    return {
        "predicted_ddg": calibrated_ddg,
        "ci_95": [ci_low, ci_high],
        "reliability": "high" if calib.pearson_r > 0.7 else "medium" if calib.pearson_r > 0.4 else "low",
        "similar_cases_n": len(similar),
        "similar_cases_success_rate": sum(1 for s in similar if s.ddg < 0) / len(similar),
        "calibration_pearson_r": calib.pearson_r,
        "calibration_rmse": calib.rmse,
        "calibration_n": calib.n_samples,
    }
```

### 颜色编码（前端）

```
置信度显示:
  🟢 高置信度 (r>0.7, CI窄)  → "这个预测比较可靠"
  🟡 中置信度 (r 0.4-0.7)   → "有一定参考价值，但需要实验验证"
  🔴 低置信度 (r<0.3, CI宽)  → "预测不太可靠，建议谨慎"
  ⚪ 无校准数据               → "此场景暂无实验数据校准"

预测ΔΔG:
  🔵 ΔΔG < -1.5  → 强稳定化预测
  🟢 -1.5 to -0.5 → 弱稳定化预测
  ⚪ -0.5 to 0.5  → 中性
  🟡 0.5 to 1.5   → 弱不稳定预测
  🔴 > 1.5         → 强不稳定预测
```

---

## 八、Token消耗时间线（10天+后续）

```
Day 1-4:  零token (全部本地计算+数据解析)
Day 5:    ~3M tokens (GFP家族知识综合)
Day 6:    ~2M tokens (实验协议库)
Day 7:    ~1M tokens (校准报告解读)
Day 8:    ~1M tokens (测试+微调)
Day 9-10: ~1M tokens (文档+测试)
─────────────────────────────
10天小计: ~8M tokens ← 极度精简，全是刀刃

Day 11-30:  ~15M tokens (扩展到50个蛋白家族, AI导师日常使用)
Day 31-90: ~30M tokens (教材引擎逐步完成, AI导师持续使用)
─────────────────────────────
90天总计: ~53M tokens

这不是700M，但这是诚实的数字。
```
