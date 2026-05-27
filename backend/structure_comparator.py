"""Protein Design Studio — 多PDB结构对比引擎

支持:
- 序列比对 (Biopython PairwiseAligner)
- 结构叠合 RMSD (Cα坐标)
- 表面性质对比 (SASA profile)
- 二级结构分布对比
- 保守性分析
"""
import math
from typing import Optional

from Bio.Align import PairwiseAligner

from .config import PDB_CACHE


def compare_structures(
    structures: list[dict],
    pdb_ids: list[str],
    chain: str = "A",
) -> dict:
    """对比2个PDB结构

    Args:
        structures: parse_pdb() 返回的结构列表
        pdb_ids: PDB ID列表
        chain: 对比的链

    Returns:
        {
            alignment: {score, pdb1_seq, pdb2_seq, match_line},
            rmsd: float (Å),
            sasa_comparison: [{resnum, sasa1, sasa2, delta, ...}],
            ss_comparison: {pdb1: {H:N, E:N, C:N}, pdb2: {H:N, E:N, C:N}},
            identity: float (序列一致性%),
            conserved_regions: [{start, end, identity}],
        }
    """
    if len(structures) < 2:
        return {"error": "需要至少2个结构进行对比"}

    s1, s2 = structures[0], structures[1]

    # 获取指定链的残基
    residues1 = [r for r in s1["residues"] if r["chain"] == chain]
    residues2 = [r for r in s2["residues"] if r["chain"] == chain]

    seq1 = "".join(r["aa"] or "X" for r in residues1)
    seq2 = "".join(r["aa"] or "X" for r in residues2)

    # ── 序列比对 ────────────────────────────────────────
    alignment = _align_sequences(seq1, seq2, pdb_ids)

    # ── RMSD ──────────────────────────────────────────
    rmsd = _calc_rmsd(residues1, residues2) if alignment["aligned_pairs"] else None

    # ── SASA对比 ──────────────────────────────────────
    sasa_comparison = _compare_sasa(residues1, residues2, alignment["aligned_pairs"])

    # ── 二级结构对比 ─────────────────────────────────
    ss_comparison = _compare_secondary_structure(residues1, residues2, pdb_ids)

    # ── 保守区域 ──────────────────────────────────────
    conserved = _find_conserved_regions(alignment["match_line"]) if alignment["match_line"] else []

    return {
        "pdb_ids": pdb_ids,
        "chain": chain,
        "alignment": alignment,
        "rmsd": round(rmsd, 2) if rmsd else None,
        "sasa_comparison": sasa_comparison[:50],  # 限制返回量
        "ss_comparison": ss_comparison,
        "identity": alignment.get("identity", 0),
        "conserved_regions": conserved[:20],
    }


def _align_sequences(seq1: str, seq2: str, pdb_ids: list[str]) -> dict:
    """序列比对，返回对齐信息"""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return {"score": 0, "identity": 0, "aligned_pairs": []}

    best = alignments[0]
    aligned1, aligned2 = best.aligned

    # 用Biopython对齐块构建display字符串
    display1_full = ""
    display2_full = ""
    match_full = ""
    pairs = []
    pos1 = pos2 = 0

    for (s1, e1), (s2, e2) in zip(aligned1, aligned2):
        # 前置gap: seq1在s1之前未对齐的部分 → seq2用gap填充
        while pos1 < s1 and pos2 < s2:
            gap_len = min(s1 - pos1, s2 - pos2)
            # 两者都在gap中 — 这种情况不常见
            for _ in range(gap_len):
                display1_full += seq1[pos1]
                display2_full += "-"
                match_full += " "
                pos1 += 1

        # seq1的gap区域 (seq2有, seq1无)
        while pos1 < s1:
            display1_full += seq1[pos1] if pos1 < len(seq1) else "-"
            display2_full += "-"
            match_full += " "
            pos1 += 1

        # seq2的gap区域 (seq1有, seq2无)
        while pos2 < s2:
            display1_full += "-"
            display2_full += seq2[pos2] if pos2 < len(seq2) else "-"
            match_full += " "
            pos2 += 1

        # 对齐块: seq1[s1:e1] ↔ seq2[s2:e2]
        for i in range(e1 - s1):
            a1 = seq1[s1 + i] if s1 + i < len(seq1) else "-"
            a2 = seq2[s2 + i] if s2 + i < len(seq2) else "-"
            display1_full += a1
            display2_full += a2
            if a1 == a2:
                match_full += "|"
            elif _is_similar(a1, a2):
                match_full += ":"
            else:
                match_full += "."
            pairs.append((pos1, pos2))
            pos1 += 1
            pos2 += 1

    # 尾部: 剩余未对齐的残基
    while pos1 < len(seq1) or pos2 < len(seq2):
        a1 = seq1[pos1] if pos1 < len(seq1) else "-"
        a2 = seq2[pos2] if pos2 < len(seq2) else "-"
        display1_full += a1
        display2_full += a2
        match_full += " " if a1 == "-" or a2 == "-" else "."
        if a1 != "-": pos1 += 1
        if a2 != "-": pos2 += 1

    # 序列一致性
    matches = sum(1 for c in match_full if c == "|" or c == ":")
    total = max(len(seq1), len(seq2))
    identity = round(matches / total * 100, 1) if total > 0 else 0

    return {
        "score": round(best.score, 1),
        "identity": identity,
        f"{pdb_ids[0]}_seq": display1_full,
        f"{pdb_ids[1]}_seq": display2_full,
        "match_line": match_full,
        "aligned_pairs": pairs,
    }


def _is_similar(aa1: str, aa2: str) -> bool:
    """判断两个氨基酸是否化学性质相似"""
    similar_groups = [
        set("AVILMFWY"),   # 疏水/芳香
        set("STNQ"),       # 极性小/酰胺
        set("DE"),         # 酸性
        set("KR"),         # 碱性
    ]
    for group in similar_groups:
        if aa1 in group and aa2 in group:
            return True
    return False


def _calc_rmsd(residues1: list, residues2: list) -> Optional[float]:
    """计算Cα RMSD (Å) — 需预先对齐"""
    # 取前N个残基的Cα做简单RMSD
    ca1 = []
    ca2 = []
    for r in residues1:
        for a in r.get("atoms", []):
            if a.get("name") == "CA":
                ca1.append((a["x"], a["y"], a["z"]))
                break
    for r in residues2:
        for a in r.get("atoms", []):
            if a.get("name") == "CA":
                ca2.append((a["x"], a["y"], a["z"]))
                break

    n = min(len(ca1), len(ca2))
    if n < 3:
        return None

    # 简单叠加RMSD (未做最优旋转, 两个结构可能在空间中有不同取向)
    # 使用前50个Cα做近似RMSD
    n_use = min(n, 50)
    sum_sq = 0.0
    for i in range(n_use):
        dx = ca1[i][0] - ca2[i][0]
        dy = ca1[i][1] - ca2[i][1]
        dz = ca1[i][2] - ca2[i][2]
        sum_sq += dx*dx + dy*dy + dz*dz

    return math.sqrt(sum_sq / n_use)


def _compare_sasa(residues1: list, residues2: list, aligned_pairs: list) -> list:
    """对比对齐残基的SASA"""
    result = []
    for idx1, idx2 in aligned_pairs[:50]:  # 限制50对
        if idx1 >= len(residues1) or idx2 >= len(residues2):
            continue
        r1, r2 = residues1[idx1], residues2[idx2]
        sasa1 = r1.get("sasa", 0)
        sasa2 = r2.get("sasa", 0)
        result.append({
            "resnum1": r1.get("resnum"),
            "resnum2": r2.get("resnum"),
            "aa1": r1.get("aa"),
            "aa2": r2.get("aa"),
            "sasa1": round(sasa1, 1),
            "sasa2": round(sasa2, 1),
            "delta": round(sasa2 - sasa1, 1),
            "surface1": sasa1 > 30,
            "surface2": sasa2 > 30,
        })
    return result


def _compare_secondary_structure(residues1: list, residues2: list, pdb_ids: list) -> dict:
    """对比二级结构分布"""
    def count_ss(residues):
        counts = {"H": 0, "E": 0, "C": 0, "other": 0}
        for r in residues:
            ss = r.get("secondary_structure", "C")
            if ss in counts:
                counts[ss] += 1
            else:
                counts["other"] += 1
        total = len(residues)
        return {k: {"count": v, "pct": round(v/total*100, 1) if total>0 else 0}
                for k, v in counts.items()}

    return {
        pdb_ids[0]: count_ss(residues1),
        pdb_ids[1]: count_ss(residues2),
    }


def _find_conserved_regions(match_line: str, min_length: int = 5) -> list:
    """从匹配行找出保守区域"""
    regions = []
    start = None
    for i, c in enumerate(match_line):
        if c in "|:":
            if start is None:
                start = i
        else:
            if start is not None:
                length = i - start
                if length >= min_length:
                    # 计算该区域的实际identity
                    region_matches = sum(1 for j in range(start, i) if match_line[j] in "|:")
                    identity = round(region_matches / length * 100, 1)
                    regions.append({
                        "start": start,
                        "end": i,
                        "length": length,
                        "identity": identity,
                    })
                start = None
    # 尾部
    if start is not None:
        length = len(match_line) - start
        if length >= min_length:
            region_matches = sum(1 for j in range(start, len(match_line)) if match_line[j] in "|:")
            identity = round(region_matches / length * 100, 1)
            regions.append({
                "start": start,
                "end": len(match_line),
                "length": length,
                "identity": identity,
            })

    return sorted(regions, key=lambda r: r["end"] - r["start"], reverse=True)
