"""Protein Design Studio — DSSP格式导出器

生成标准DSSP格式的二级结构分配输出。
参考: https://swift.cmbi.umcn.nl/gv/dssp/
"""
from typing import TextIO


def export_dssp(structure: dict, output_path: str | None = None) -> str:
    """将结构数据导出为DSSP格式文本

    Args:
        structure: parse_pdb() 返回的结构字典
        output_path: 可选，保存到文件

    Returns:
        DSSP格式文本
    """
    residues = structure["residues"]
    pdb_id = structure["pdb_id"]

    lines = []
    lines.append(f"==== Secondary Structure Definition for {pdb_id} ====")
    lines.append(f"  #  RESIDUE AA STRUCTURE BP1 BP2  ACC    PHI    PSI")
    lines.append(f"  " + "-" * 55)

    # 映射简化二级结构到DSSP标准8类
    SS_DSSP_MAP = {
        "H": "H",  # α-helix
        "G": "G",  # 3-10 helix
        "I": "I",  # π-helix
        "E": "E",  # extended strand (β-sheet)
        "B": "B",  # β-bridge
        "T": "T",  # turn
        "S": "S",  # bend
        "C": "C",  # coil/loop
    }

    for i, r in enumerate(residues):
        if not r.get("aa"):
            continue

        ss = r.get("secondary_structure", "C")
        dssp_code = SS_DSSP_MAP.get(ss, "C")

        phi = r.get("phi")
        psi = r.get("psi")
        sasa = r.get("sasa", 0)

        phi_str = f"{phi:7.1f}" if phi is not None else "       "
        psi_str = f"{psi:7.1f}" if psi is not None else "       "

        # 计算桥接伙伴 (简化版 — 非标准DSSP主链氢键配对算法)
        # 标准DSSP的BP基于主链NH-CO氢键模式区分平行/反平行折叠。
        # 此处简化为: 同链邻居中同为β-strand(E)的残基编号。
        bp1 = 0
        bp2 = 0
        for nb in r.get("neighbors", []):
            # 找同链邻居中同为E的残基
            for r2 in residues:
                if r2["resnum"] == nb and r2["chain"] == r["chain"]:
                    if r2.get("secondary_structure") == "E":
                        if bp1 == 0:
                            bp1 = nb
                        elif bp2 == 0:
                            bp2 = nb
                    break

        line = (
            f"{i+1:4d}"
            f"{i+1:5d}"
            f" {r['resname']:>3s}"
            f" {r['aa']}"
            f" {dssp_code}"
            f" {bp1:4d}{bp2:4d}"
            f" {int(sasa):4d}"
            f" {phi_str} {psi_str}"
        )
        lines.append(line)

    # 统计
    ss_counts = {}
    for r in residues:
        ss = r.get("secondary_structure", "C")
        ss_counts[ss] = ss_counts.get(ss, 0) + 1
    total = len(residues)
    lines.append("")
    lines.append("Summary:")
    for ss_type in ["H", "E", "C"]:
        count = ss_counts.get(ss_type, 0)
        pct = count / total * 100 if total > 0 else 0
        ss_name = {"H": "alpha-helix", "E": "beta-strand", "C": "coil/loop"}.get(ss_type, ss_type)
        lines.append(f"  {ss_type} ({ss_name}): {count:4d} ({pct:5.1f}%)")

    output = "\n".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(output)

    return output


def get_ss_summary(structure: dict) -> dict:
    """获取二级结构统计摘要"""
    residues = structure["residues"]
    ss_counts = {}
    phi_values = []
    psi_values = []

    for r in residues:
        ss = r.get("secondary_structure", "C")
        ss_counts[ss] = ss_counts.get(ss, 0) + 1
        if r.get("phi") is not None:
            phi_values.append(r["phi"])
        if r.get("psi") is not None:
            psi_values.append(r["psi"])

    total = len(residues)
    return {
        "pdb_id": structure["pdb_id"],
        "total_residues": total,
        "helix_count": ss_counts.get("H", 0),
        "helix_pct": round(ss_counts.get("H", 0) / total * 100, 1) if total > 0 else 0,
        "strand_count": ss_counts.get("E", 0),
        "strand_pct": round(ss_counts.get("E", 0) / total * 100, 1) if total > 0 else 0,
        "coil_count": ss_counts.get("C", 0),
        "coil_pct": round(ss_counts.get("C", 0) / total * 100, 1) if total > 0 else 0,
        "phi_mean": round(sum(phi_values) / len(phi_values), 1) if phi_values else None,
        "psi_mean": round(sum(psi_values) / len(psi_values), 1) if psi_values else None,
    }
