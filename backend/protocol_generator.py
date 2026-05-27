"""Protein Design Studio — 实验方案自动生成器

基于模板引擎 + AI润色，生成规范化的湿实验方案。
覆盖: 定点突变、蛋白表达、纯化、表征。
"""
from typing import Optional


# ── 内置模板 ──────────────────────────────────────────────

SDM_TEMPLATE = """## 定点突变实验方案 (QuikChange)

### 蛋白信息
- PDB: {pdb_id}
- 突变: {mutation_desc}
- 引物长度: {primer_length} bp

### 引物设计
{primers}

### PCR 反应体系 (50 μL)
| 组分 | 体积 |
|------|------|
| 模板 DNA (10-50 ng/μL) | 1 μL |
| 正向引物 (10 μM) | 1.25 μL |
| 反向引物 (10 μM) | 1.25 μL |
| dNTP mix (2.5 mM each) | 4 μL |
| 10× PCR buffer | 5 μL |
| DNA 聚合酶 (2.5 U/μL) | 0.5 μL |
| ddH₂O | 补至 50 μL |

### PCR 程序
| 步骤 | 温度 | 时间 | 循环 |
|------|------|------|------|
| 预变性 | 95°C | 3 min | 1× |
| 变性 | 95°C | 30 s | |
| 退火 | 55-65°C (梯度) | 30 s | 18× |
| 延伸 | 68°C | {extension_time} | |
| 终延伸 | 68°C | 10 min | 1× |
| 保存 | 4°C | ∞ | |

### DpnI 消化
- 加 1 μL DpnI (20 U/μL) 到 PCR 产物
- 37°C 孵育 1-2 小时
- 可选: 65°C 20 min 灭活

### 转化
1. 取 5-10 μL DpnI 消化产物加入 50 μL 感受态细胞
2. 冰浴 30 min
3. 42°C 热激 45-90 s (根据菌株)
4. 冰浴 2 min
5. 加 500 μL SOC/LB，37°C 振荡 1 h
6. 涂布含抗生素的 LB 平板
7. 37°C 过夜培养

### 筛选
- 挑取 3-5 个单克隆
- 接种 5 mL LB + 抗生素
- 37°C 振荡过夜
- 提取质粒 → Sanger 测序验证

### 注意事项
- 模板用量宜少 (过多增加背景)
- 使用高保真聚合酶 (如 Phusion, Q5)
- 延伸时间根据质粒长度: 2 kb/min (常规聚合酶) 或 30 sec/kb (高保真)
"""

EXPRESSION_TEMPLATE = """## 蛋白表达方案 (E. coli 标准)

### 菌株选择
{strain_info}

### 转化与预培养
1. 将表达质粒转化至表达菌株
2. 挑取单克隆接种至 5 mL LB + 抗生素
3. 37°C, 220 rpm 过夜培养

### 放大培养
1. 按 1:100 接种至 {culture_volume} mL LB + 抗生素
2. 37°C, 220 rpm 培养至 OD₆₀₀ = 0.6-0.8 (~2-3 h)
3. 取样 1 mL 作为未诱导对照

### 诱导
- 加入 IPTG 至终浓度 {iptg_conc} mM
- {induction_temp}°C, 220 rpm, {induction_time} h
- 可选: 低温诱导 (15-20°C, 16-20 h) 提高可溶性

### 收菌
1. 4°C, 4000-6000 g 离心 15-20 min
2. 弃上清
3. 菌体沉淀可用 PBS 洗一次
4. -20°C 或 -80°C 保存，或直接裂解

### 裂解
1. 重悬菌体于裂解缓冲液 (PBS + 蛋白酶抑制剂)
2. 超声破碎: 振幅 30-40%, 超声 3s/间隔 5s, 总计 5-10 min
3. 4°C, 12000-15000 g 离心 20-30 min
4. 收集上清 (可溶性部分) 和沉淀 (包涵体)

### 注意事项
- OD₆₀₀ 过高诱导 → 包涵体风险
- 低温诱导 → 提高可溶性
- 添加 2-5% 乙醇或 0.5 M 山梨醇可辅助折叠
"""

PURIFICATION_TEMPLATE = """## 蛋白纯化方案

### 样品准备
- 上清经 0.45 μm 滤膜过滤
- 所有缓冲液预冷至 4°C

### Ni-NTA 亲和层析 (His-tag)
{his_tag_section}

### 缓冲液配方
| 缓冲液 | 成分 |
|--------|------|
| 结合缓冲液 | 50 mM Tris-HCl (pH 8.0), 300 mM NaCl, 20 mM 咪唑 |
| 洗涤缓冲液 | 50 mM Tris-HCl (pH 8.0), 300 mM NaCl, 50 mM 咪唑 |
| 洗脱缓冲液 | 50 mM Tris-HCl (pH 8.0), 300 mM NaCl, 250 mM 咪唑 |

### 层析步骤
1. 平衡柱子: 5-10 CV 结合缓冲液
2. 上样: 流速 0.5-1 mL/min
3. 洗涤: 10-20 CV 洗涤缓冲液 (去除杂蛋白)
4. 洗脱: 5-10 CV 洗脱缓冲液, 收集 1-2 mL/管
5. 再生: 5 CV 500 mM 咪唑 → 10 CV 结合缓冲液

### 后续处理
1. SDS-PAGE 验证纯度
2. Bradford 或 BCA 法定量
3. 透析/脱盐去除咪唑
4. 浓缩至所需浓度 (超滤管, MWCO 适当)

### 注意事项
- 咪唑影响后续实验 → 必须透析去除
- 低浓度咪唑 (20 mM) 可减少非特异性结合
- 柱子在 4°C 使用以减少蛋白酶活性
"""

CHARACTERIZATION_TEMPLATE = """## 蛋白表征方案

### SDS-PAGE
- 分离胶: 12-15% (根据蛋白大小)
- 上样量: 5-10 μg/孔
- 染色: 考马斯亮蓝 R-250

### 浓度测定
- 方法: Bradford 或 BCA
- 标准品: BSA (0.1-2.0 mg/mL)

### 圆二色光谱 (CD)
- 远紫外 (190-260 nm): 二级结构分析
- 近紫外 (260-320 nm): 三级结构指纹
- 蛋白浓度: 0.2-0.5 mg/mL
- 光程: 0.1 cm (远紫外) 或 1 cm (近紫外)
- 缓冲液: 低盐 PBS 或 10 mM 磷酸盐 (避免高Cl⁻)

### 热稳定性
- 方法: CD 或 DSF (Differential Scanning Fluorimetry)
- CD: 监测 222 nm (α-螺旋) 或 218 nm (β-折叠) 随温度变化
- DSF: SYPRO Orange 染料, 温度梯度 1°C/min
- 温度范围: 20-95°C
- 计算 Tm (中点变性温度)

### 酶活测定 (如适用)
- 底物浓度梯度 (0.1-10× Km)
- 酶浓度: 确定线性范围
- 反应时间: 测定初速度 (通常 1-10 min)
- 计算: kcat, Km, kcat/Km
"""


def generate_protocol(
    protocol_type: str,
    pdb_id: str = "",
    mutation_desc: str = "",
    primers: str = "",
    primer_length: int = 30,
    extension_time: str = "7 min",
    strain_info: str = "BL21(DE3) (常规表达)",
    culture_volume: str = "500",
    iptg_conc: float = 0.5,
    induction_temp: int = 25,
    induction_time: int = 16,
    his_tag_section: str = "- 使用 Ni-NTA 琼脂糖树脂 (如 Qiagen, GE Healthcare)\n- 柱体积: 1-2 mL (分析级)\n- 可选: 加 5 mM β-巯基乙醇防止氧化",
    ai_description: str = "",
) -> str:
    """生成实验方案

    Args:
        protocol_type: sdm | expression | purification | characterization | full
        其他参数: 模板变量

    Returns:
        Markdown 格式实验方案
    """
    if protocol_type == "sdm":
        template = SDM_TEMPLATE
        vars = {
            "pdb_id": pdb_id,
            "mutation_desc": mutation_desc,
            "primers": primers,
            "primer_length": primer_length,
            "extension_time": extension_time,
        }
    elif protocol_type == "expression":
        template = EXPRESSION_TEMPLATE
        vars = {
            "strain_info": strain_info,
            "culture_volume": culture_volume,
            "iptg_conc": iptg_conc,
            "induction_temp": induction_temp,
            "induction_time": induction_time,
        }
    elif protocol_type == "purification":
        template = PURIFICATION_TEMPLATE
        vars = {"his_tag_section": his_tag_section}
    elif protocol_type == "characterization":
        template = CHARACTERIZATION_TEMPLATE
        vars = {}
    elif protocol_type == "full":
        parts = [
            generate_protocol("sdm", pdb_id=pdb_id, mutation_desc=mutation_desc,
                            primers=primers, extension_time=extension_time),
            "\n---\n",
            generate_protocol("expression", strain_info=strain_info,
                            culture_volume=culture_volume, iptg_conc=iptg_conc,
                            induction_temp=induction_temp, induction_time=induction_time),
            "\n---\n",
            generate_protocol("purification", his_tag_section=his_tag_section),
            "\n---\n",
            generate_protocol("characterization"),
        ]
        if ai_description:
            parts.insert(0, f"> 🤖 AI备注: {ai_description}\n\n")
        return "\n".join(parts)
    else:
        raise ValueError(f"Unknown protocol type: {protocol_type}")

    result = template.format(**vars)

    if ai_description:
        result = f"> 🤖 AI备注: {ai_description}\n\n" + result

    return result


async def generate_ai_protocol(
    pdb_id: str,
    mutations: list[dict],
    ai_client_func,
    server_base: str = "https://api.deepseek.com/v1",
) -> str:
    """使用AI生成针对特定蛋白的定制化方案

    Args:
        pdb_id: PDB ID
        mutations: [{\"resnum\": 222, \"wt\": \"E\", \"mut\": \"A\"}, ...]
        ai_client_func: ai_client.chat_stream 或类似函数
    """
    mut_desc = ", ".join(f"{m['wt']}{m['resnum']}{m['mut']}" for m in mutations)
    prompt = f"""请为以下蛋白突变设计生成完整的湿实验方案:

蛋白: PDB {pdb_id}
突变: {mut_desc}

请包含:
1. 定点突变引物设计要点 (Tm, GC%, 长度建议)
2. 推荐表达系统 (E.coli/酵母/哺乳动物)
3. 纯化策略 (标签选择, 层析步骤)
4. 表征方法 (哪些实验最能验证突变效应)
5. 特殊注意事项 (该蛋白已知的稳定性/溶解度问题)

用中文, 格式为 Markdown。"""

    # 流式收集回复
    full_response = ""
    async for token in ai_client_func(prompt, context=f"蛋白: PDB {pdb_id}", model="deepseek-v4-pro"):
        full_response += token

    return generate_protocol(
        "full",
        pdb_id=pdb_id,
        mutation_desc=mut_desc,
        primers="(由AI生成, 请用 primer_designer API 获取精确序列)",
        ai_description=full_response[:2000],  # 限制长度
    )
