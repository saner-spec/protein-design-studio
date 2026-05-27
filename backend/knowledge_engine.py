"""Protein Design Studio — 文献知识引擎

对加载的蛋白质结构，自动调用AI生成文献综述和设计策略。
结果缓存在本地，避免重复API调用。
"""
import json
import hashlib
from pathlib import Path
from typing import Optional

from .config import DATA_DIR, DS_API_KEY
from .ai_client import chat_stream

KNOWLEDGE_DIR = DATA_DIR / "knowledge"
KNOWLEDGE_DIR.mkdir(parents=True, exist_ok=True)

KNOWLEDGE_PROMPT = """你是一位蛋白质工程专家。请为以下蛋白质撰写一份结构化的知识摘要。

蛋白信息:
{pdb_info}

请包含以下内容 (每项2-5句话，中文):

1. **蛋白名称与功能**: 这个蛋白是什么，在生物体内做什么
2. **结构特征**: 折叠类型、结构域、活性位点、关键残基
3. **已知突变**: 文献中报道过的重要突变及其效应（如果有）
4. **工程化历史**: 这个蛋白被改造过吗？有哪些经典的成功/失败案例
5. **设计建议**: 如果要对这个蛋白做突变设计，应该关注哪些区域？为什么？
6. **实验注意事项**: 表达、纯化、表征时的特殊考虑
7. **参考文献**: 列出3-5篇关键文献（格式: 作者, 年份, 标题, 期刊）

格式要求:
- Markdown格式
- 每项200-400字
- 不确定的内容标注"[推测]"
- 不要虚构文献，只列出你确定存在的
"""


def _get_cache_key(pdb_id: str) -> str:
    return hashlib.md5(pdb_id.encode()).hexdigest()[:12]


def _cache_path(pdb_id: str) -> Path:
    return KNOWLEDGE_DIR / f"{_get_cache_key(pdb_id)}.json"


def get_cached_knowledge(pdb_id: str) -> Optional[dict]:
    """获取缓存的知识页面"""
    path = _cache_path(pdb_id)
    if path.exists():
        with open(path, encoding="utf-8-sig") as f:
            return json.load(f)
    return None


def save_knowledge(pdb_id: str, content: str, pdb_info: dict) -> dict:
    """保存知识页面到缓存"""
    data = {
        "pdb_id": pdb_id,
        "content": content,
        "residue_count": pdb_info.get("residue_count", 0),
        "chains": pdb_info.get("chains", []),
    }
    with open(_cache_path(pdb_id), "w") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    return data


def build_pdb_info_for_prompt(structure: dict) -> str:
    """从结构数据构建提示词用的蛋白信息摘要"""
    residues = structure.get("residues", [])
    ss = {"H": 0, "E": 0, "C": 0}
    for r in residues:
        s = r.get("secondary_structure", "C")
        ss[s] = ss.get(s, 0) + 1
    total = len(residues)
    ss_summary = ", ".join(
        f"{k}={v}({v/total*100:.0f}%)" for k, v in ss.items() if v > 0
    ) if total > 0 else "未知"

    # 找出配体
    ligands = [l["resname"] for l in structure.get("ligands", [])
               if l["resname"] not in ("HOH", "WAT")]

    # 表面/埋藏比例
    surface_count = sum(1 for r in residues if r.get("rel_sasa", 0) > 0.25)
    buried_count = total - surface_count

    return f"""PDB ID: {structure.get('pdb_id', '')}
残基数: {total}
链: {', '.join(structure.get('chains', []))}
二级结构: {ss_summary}
配体/辅因子: {', '.join(ligands[:10]) if ligands else '无'}
表面残基: {surface_count} ({surface_count/total*100:.0f}%)
埋藏残基: {buried_count} ({buried_count/total*100:.0f}%)
"""


async def generate_knowledge(pdb_id: str, structure: dict, refresh: bool = False) -> dict:
    """生成或获取蛋白质知识页面

    优先从缓存读取，若无缓存或refresh=True则调用AI生成。
    """
    # 检查缓存
    if not refresh:
        cached = get_cached_knowledge(pdb_id)
        if cached:
            return {"cached": True, **cached}

    # 构建提示词
    pdb_info = build_pdb_info_for_prompt(structure)
    prompt = KNOWLEDGE_PROMPT.format(pdb_info=pdb_info)

    if not DS_API_KEY:
        return {"error": "API密钥未配置", "cached": False}

    # 流式收集AI回复
    full_text = ""
    try:
        async for token in chat_stream(prompt, context="", model="deepseek-v4-pro"):
            full_text += token
    except Exception as e:
        return {"error": str(e), "cached": False}

    if not full_text.strip():
        return {"error": "AI未返回内容", "cached": False}

    # 保存缓存
    result = save_knowledge(pdb_id, full_text, structure)
    return {"cached": False, **result}


async def generate_knowledge_brief(pdb_id: str, structure: dict) -> str:
    """生成超简知识摘要（点击残基时注入AI上下文）"""
    residues = structure.get("residues", [])
    ss = {"H": 0, "E": 0, "C": 0}
    for r in residues:
        ss[r.get("secondary_structure", "C")] = ss.get(r.get("secondary_structure", "C"), 0) + 1
    total = len(residues)

    ligands = [l["resname"] for l in structure.get("ligands", [])
               if l["resname"] not in ("HOH", "WAT")]

    brief = f"{structure['pdb_id']}: {total}aa, β{ss.get('E',0)}({ss.get('E',0)/total*100:.0f}%)/α{ss.get('H',0)}({ss.get('H',0)/total*100:.0f}%)"
    if ligands:
        brief += f", 配体: {','.join(ligands[:3])}"
    return brief
