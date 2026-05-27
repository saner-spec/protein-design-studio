"""Protein Design Studio — 设计报告生成器

一键汇总所有分析结果，生成完整Markdown报告。
"""
from datetime import datetime

from .dssp_exporter import get_ss_summary
from .knowledge_engine import get_cached_knowledge


def generate_report(
    structure: dict,
    scores: dict | None = None,
    calibration: dict | None = None,
    protocol: str | None = None,
    mutations: list[dict] | None = None,
) -> str:
    """生成综合设计报告

    Args:
        structure: parse_pdb() 返回的结构
        scores: ESM打分结果 (可选)
        calibration: 校准信息 (可选)
        protocol: 实验方案Markdown (可选)
        mutations: 选中的突变列表 (可选)
    """
    pdb_id = structure["pdb_id"]
    ss = get_ss_summary(structure)
    knowledge = get_cached_knowledge(pdb_id)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    lines = []
    lines.append(f"# Protein Design Studio — 设计报告")
    lines.append(f"**PDB:** {pdb_id} | **生成时间:** {now}")
    lines.append("")

    # ── 蛋白概览 ──
    lines.append("## 📊 蛋白概览")
    lines.append("")
    lines.append(f"- PDB ID: {pdb_id}")
    lines.append(f"- 残基数: {ss['total_residues']}")
    lines.append(f"- 链: {', '.join(structure.get('chains', []))}")
    lines.append(f"- α螺旋: {ss['helix_count']} ({ss['helix_pct']}%)")
    lines.append(f"- β折叠: {ss['strand_count']} ({ss['strand_pct']}%)")
    lines.append(f"- 环区: {ss['coil_count']} ({ss['coil_pct']}%)")

    # 配体
    ligands = [l for l in structure.get("ligands", []) if l["resname"] not in ("HOH", "WAT")]
    if ligands:
        lines.append(f"- 配体: {', '.join(l['resname'] for l in ligands[:10])}")
    lines.append("")

    # ── 知识综述 ──
    if knowledge:
        lines.append("## 📚 文献知识")
        lines.append("")
        lines.append(knowledge.get("content", "")[:3000])
        lines.append("")

    # ── ESM 打分 ──
    if scores and scores.get("scores"):
        lines.append("## 🔬 ESM突变打分")
        lines.append("")
        lines.append(f"**残基:** {scores.get('wt_aa','?')}{scores.get('resnum','?')} (链{scores.get('chain','?')})")
        lines.append("")
        lines.append("| 排名 | 氨基酸 | ESM分数 | 性质 | 评估 |")
        lines.append("|------|--------|---------|------|------|")
        for i, s in enumerate(scores["scores"][:10]):
            rank = i + 1
            color = "🟢" if s["color"] == "good" else "🔴" if s["color"] == "bad" else "⚪"
            lines.append(f"| {rank} | {s['aa']} | {s['score']} | {s.get('property','')} | {color} |")
        lines.append("")
        lines.append(f"> {scores.get('note', '')}")
        lines.append("")

    # ── 校准 ──
    if calibration:
        lines.append("## 📐 实验校准")
        lines.append("")
        lines.append(f"- 数据集: {calibration.get('name', '')}")
        lines.append(f"- 样本量: {calibration.get('n_samples', 0)}")
        lines.append(f"- R²: {calibration.get('r_squared', 'N/A')}")
        lines.append(f"- 斜率: {calibration.get('slope', 'N/A')}")
        lines.append(f"- 截距: {calibration.get('intercept', 'N/A')}")
        lines.append("")

    # ── 突变设计 ──
    if mutations:
        lines.append("## 🧬 突变设计方案")
        lines.append("")
        for m in mutations:
            lines.append(f"- **{m.get('wt','?')}{m.get('position','?')}{m.get('mut','?')}** (ESM: {m.get('score','?')})")
        lines.append("")

    # ── 实验方案 ──
    if protocol:
        lines.append("## 📋 实验方案")
        lines.append("")
        lines.append(protocol)
        lines.append("")

    # ── 免责声明 ──
    lines.append("---")
    lines.append("")
    lines.append("*本报告由 Protein Design Studio 自动生成。ESM打分仅供参考（相关性约0.4-0.6），实验方案为模板生成，实际实验需导师审核。*")

    return "\n".join(lines)
