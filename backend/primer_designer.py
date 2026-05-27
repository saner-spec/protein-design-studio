"""Protein Design Studio — PCR引物设计与序列设计工作流

Phase 2 新增:
  - 定点突变引物设计 (QuikChange-style)
  - Tm/GC含量/发夹/二聚体评估
  - 密码子优化与退火温度建议
"""

import re
import math
from typing import Optional
from .residue_map import AA3TO1

# ── 密码子表 (标准遗传密码) ──────────────────────────────
CODON_TABLE = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
}

# 偏好密码子 (大肠杆菌高频)
PREFERRED_CODONS = {
    "A": "GCG", "C": "TGC", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGC", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAC", "P": "CCG", "Q": "CAG", "R": "CGT",
    "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAT",
}

# ── Nearest-neighbor 热力学参数 (SantaLucia, 1998) ──────
# ΔH (kcal/mol) 和 ΔS (cal/mol·K) 用于 1M NaCl
NN_DH = {
    "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2,
    "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
    "CG": -10.6, "GC": -10.6, "GG": -8.0, "CC": -8.0,
}
NN_DS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3,
    "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -27.2, "GG": -19.9, "CC": -19.9,
}

# 末端校正
TERMINAL_DH = {"G": 0.1, "C": 0.1, "A": 2.3, "T": 2.3}
TERMINAL_DS = {"G": -2.8, "C": -2.8, "A": -4.1, "T": -4.1}

# 盐浓度校正常数 (Na+ 浓度 M)
SALT_CONC = 0.050  # 50 mM Na+
MG_CONC = 0.002    # 2 mM Mg2+

# ── 蛋白质→DNA 反向翻译 ──────────────────────────────────

def aa_to_dna(sequence: str, use_preferred: bool = True) -> str:
    """将蛋白质序列反向翻译为DNA序列
    
    Args:
        sequence: 单字母氨基酸序列
        use_preferred: 是否使用偏好密码子
        
    Returns:
        DNA序列 (5'→3')
    """
    codons = PREFERRED_CODONS if use_preferred else {}
    dna = []
    for aa in sequence.upper():
        if aa in codons:
            dna.append(codons[aa])
        elif aa in CODON_TABLE:
            dna.append(CODON_TABLE[aa][0])  # 使用第一个密码子
        else:
            dna.append("NNN")
    return "".join(dna)


def complement(dna: str) -> str:
    """DNA反向互补"""
    comp_map = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp_map.get(b, "N") for b in dna.upper())


def reverse_complement(dna: str) -> str:
    """反向互补 (5'→3')"""
    return complement(dna)[::-1]


# ── Tm 计算 ──────────────────────────────────────────────

def calc_tm_nearest_neighbor(seq: str, na_conc: float = SALT_CONC, mg_conc: float = MG_CONC) -> Optional[float]:
    """使用 Nearest-Neighbor 模型计算 Tm (°C)
    
    适用于短引物 (< 50 bp)。
    """
    seq = seq.upper().strip()
    if len(seq) < 4:
        return None
    if not re.match(r"^[ACGT]+$", seq):
        return None

    total_dh = 0.0
    total_ds = 0.0

    for i in range(len(seq) - 1):
        dimer = seq[i:i + 2]
        if dimer in NN_DH:
            total_dh += NN_DH[dimer]
            total_ds += NN_DS[dimer]
        else:
            # 互补对
            comp_dimer = complement(dimer)
            if comp_dimer in NN_DH:
                total_dh += NN_DH[comp_dimer]
                total_ds += NN_DS[comp_dimer]

    # 末端 GC 校正
    first = seq[0]
    last = seq[-1]
    if first in TERMINAL_DH:
        total_dh += TERMINAL_DH[first]
        total_ds += TERMINAL_DS[first]
    if last in TERMINAL_DH:
        total_dh += TERMINAL_DH[last]
        total_ds += TERMINAL_DS[last]

    # 盐浓度校正
    # [Na+] 等效: [Na+] + 0.12 * [Mg2+]^0.5 / 0.12
    if mg_conc > 0:
        na_eff = na_conc + mg_conc * 0.92
    else:
        na_eff = na_conc

    # 热力学常数 R = 1.987 cal/(mol·K)
    R = 1.987
    # 引物浓度 ~ 0.5 µM
    c_prime = 5e-7

    # Returns-Breslauer Tm 公式
    total_ds_corr = total_ds + 0.368 * (len(seq) - 1) * math.log(na_eff / 1.0, math.e)
    if total_ds_corr <= 0:
        return None

    tm = (1000 * total_dh) / (total_ds_corr + R * math.log(c_prime / 4, math.e)) - 273.15
    return round(tm, 1)


def calc_tm_basic(seq: str) -> float:
    """基础 Tm 计算公式 (Wallace rule): 2°C×(A+T) + 4°C×(G+C)"""
    seq = seq.upper().strip()
    if not seq:
        return 0.0
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return round(2.0 * at + 4.0 * gc, 1)


def gc_content(seq: str) -> float:
    """计算 GC 含量 (%)"""
    seq = seq.upper().strip()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 1)


def calc_hairpin_tm(seq: str) -> float:
    """估算发夹结构的 Tm (简化为回文检测)"""
    seq = seq.upper().strip()
    n = len(seq)
    max_stem = 0
    for i in range(n // 2):
        if seq[i] == complement(seq[n - 1 - i]):
            max_stem += 1
        else:
            break
    # 每个配对的碱基对贡献 ~2°C
    return round(max_stem * 2.0, 1)


def calc_self_dimer_dg(seq: str) -> float:
    """估算自二聚体 ΔG (简化为连续互补碱基对数量)"""
    seq = seq.upper().strip()
    rc = reverse_complement(seq)
    max_match = 0
    # 错位比对
    for offset in range(1, len(seq)):
        match = 0
        for i in range(len(seq) - offset):
            if seq[i] == rc[i + offset]:
                match += 1
            elif match > 0:
                break
        max_match = max(max_match, match)
        match = 0
        for i in range(len(seq) - offset):
            if seq[i + offset] == rc[i]:
                match += 1
            elif match > 0:
                break
        max_match = max(max_match, match)
    # 近似: 每个匹配碱基 ~ -1 kcal/mol
    return round(-max_match * 1.0, 1)


    # 引物质量评估 ─────────────────────────────────────────

def assess_primer(seq: str) -> dict:
    """综合评估引物质量"""
    seq = seq.upper().strip()
    n = len(seq)
    gc = gc_content(seq)
    tm_wallace = calc_tm_basic(seq)
    tm_nn = calc_tm_nearest_neighbor(seq) or tm_wallace  # NN失败时用Wallace兜底
    hairpin = calc_hairpin_tm(seq)
    dimer_dg = calc_self_dimer_dg(seq)

    # 评估标准
    issues = []
    if n < 18:
        issues.append("引物过短 (<18 nt)，特异性可能不足")
    elif n > 40:
        issues.append("引物过长 (>40 nt)，合成成本高")

    if gc < 40:
        issues.append(f"GC含量偏低 ({gc}%)，建议 40-60%")
    elif gc > 60:
        issues.append(f"GC含量偏高 ({gc}%)，建议 40-60%")

    if hairpin > 5:
        issues.append(f"发夹结构风险 (Tm={hairpin}°C)")

    if dimer_dg < -6:
        issues.append(f"自二聚体风险 (ΔG≈{dimer_dg} kcal/mol)")

    quality = "good"
    if len(issues) > 2:
        quality = "poor"
    elif len(issues) > 0:
        quality = "fair"

    # 退火温度: 优先NN, 失败则Wallace-5
    ann_temp = round(tm_nn - 5, 1) if tm_nn else round(tm_wallace - 5, 1)

    return {
        "sequence": seq,
        "length": n,
        "gc_percent": gc,
        "tm_nearest_neighbor": round(tm_nn, 1),
        "tm_wallace": tm_wallace,
        "hairpin_tm": hairpin,
        "self_dimer_dg": dimer_dg,
        "quality": quality,
        "issues": issues,
        "recommended_annealing_temp": ann_temp,
    }


# ── 定点突变引物设计 ────────────────────────────────────

def design_mutagenesis_primers(
    wt_sequence: str,
    mutation_position: int,
    target_aa: str,
    wt_aa: str = "",
    primer_length: int = 30,
    use_preferred_codons: bool = True,
    pdb_resnum: int | None = None,
) -> dict:
    """设计 QuikChange 定点突变引物
    
    Args:
        wt_sequence: 野生型蛋白质序列 (单字母)
        mutation_position: 突变位置 (0-based, 在蛋白质序列中的seqres_idx)
        target_aa: 目标氨基酸 (单字母)
        wt_aa: 野生型氨基酸 (自动检测)
        primer_length: 引物总长度 (bp)
        use_preferred_codons: 是否使用偏好密码子
        pdb_resnum: PDB残基编号 (用于突变标签显示，默认使用mutation_position+1)
        
    Returns:
        包含正向/反向引物和评估结果的字典
    """
    wt_seq = wt_sequence.upper()
    if not wt_aa:
        wt_aa = wt_seq[mutation_position] if mutation_position < len(wt_seq) else "X"

    # 反向翻译
    wt_dna = aa_to_dna(wt_seq, use_preferred=use_preferred_codons)
    target_codon_only = aa_to_dna(target_aa, use_preferred=use_preferred_codons)  # 仅是目标AA的单个密码子 (3bp)

    # 突变位点在DNA中的位置 (每个aa对应3个碱基)
    dna_pos = mutation_position * 3

    # 构建突变DNA序列
    mut_dna = wt_dna[:dna_pos] + target_codon_only + wt_dna[dna_pos + 3:]

    # 设计引物: 以突变位点为中心，向两侧延伸
    half = (primer_length - 3) // 2  # 突变两侧的长度

    # 正向引物 (sense strand, 与模板相同, 包含突变)
    start_f = max(0, dna_pos - half)
    end_f = min(len(wt_dna), dna_pos + 3 + half)

    # 调整长度
    while end_f - start_f < primer_length and start_f > 0:
        start_f -= 1
    while end_f - start_f < primer_length and end_f < len(wt_dna):
        end_f += 1

    forward_primer = mut_dna[start_f:end_f].upper()

    # 反向引物 = 正向引物的反向互补
    reverse_primer = reverse_complement(forward_primer)

    # 评估
    fwd_assessment = assess_primer(forward_primer)
    rev_assessment = assess_primer(reverse_primer)

    # 突变密码子信息
    wt_codon = wt_dna[dna_pos:dna_pos + 3].upper() if dna_pos + 3 <= len(wt_dna) else "NNN"
    target_codon = target_codon_only[0:3].upper() if len(target_codon_only) >= 3 else "NNN"

    # 突变区域在引物中的位置
    mut_in_primer_f = dna_pos - start_f
    mut_in_primer_r = len(forward_primer) - mut_in_primer_f - 3

    # 突变标签使用PDB残基编号（更直观）
    display_resnum = pdb_resnum if pdb_resnum is not None else (mutation_position + 1)

    return {
        "wt_aa": wt_aa,
        "target_aa": target_aa,
        "mutation": f"{wt_aa}{display_resnum}{target_aa}",
        "position": mutation_position,
        "wt_codon": wt_codon,
        "target_codon": target_codon,
        "codon_change": f"{wt_codon}→{target_codon}",
        "forward_primer": fwd_assessment,
        "reverse_primer": rev_assessment,
        "mut_region_in_fwd": (mut_in_primer_f, mut_in_primer_f + 3),
        "mut_region_in_rev": (mut_in_primer_r, mut_in_primer_r + 3),
        "primer_tm_diff": abs(fwd_assessment["tm_nearest_neighbor"] - rev_assessment["tm_nearest_neighbor"])
            if fwd_assessment["tm_nearest_neighbor"] and rev_assessment["tm_nearest_neighbor"] else None,
        "recommended_annealing_temp": round(
            (fwd_assessment["recommended_annealing_temp"] + rev_assessment["recommended_annealing_temp"]) / 2, 1
        ) if fwd_assessment["recommended_annealing_temp"] and rev_assessment["recommended_annealing_temp"] else None,
    }


def design_sequencing_primer(
    sequence: str,
    position: int,
    direction: str = "forward",
    primer_length: int = 20,
) -> dict:
    """设计测序引物 (用于验证突变)
    
    Args:
        sequence: DNA 序列
        position: 目标位置 (0-based, DNA坐标)
        direction: "forward" 或 "reverse"
        primer_length: 引物长度
        
    Returns:
        引物评估结果
    """
    seq = sequence.upper().strip()
    if direction == "forward":
        # 在突变位点上游设计正向引物
        end = max(primer_length, position)
        start = end - primer_length
        primer_seq = seq[start:end]
    else:
        # 在突变位点下游设计反向引物
        start = position + 3  # 跳过突变密码子
        end = start + primer_length
        if end > len(seq):
            end = len(seq)
            start = end - primer_length
        primer_seq = reverse_complement(seq[start:end])

    return assess_primer(primer_seq)


# ── 设计工作流摘要 ──────────────────────────────────────

def design_workflow_summary(
    pdb_id: str,
    chain: str,
    resnum: int,
    wt_aa: str,
    target_aa: str,
    wt_sequence: str,
    esm_score: Optional[float] = None,
) -> dict:
    """生成完整的序列设计工作流摘要"""
    # 定位突变位置
    mutation_position = -1
    # 从序列中找位置 — 需要考虑SEQRES对齐
    # 简单实现: 找到第N个指定氨基酸(但实际应使用seqres_idx)
    # 这里我们假设调用者提供的wt_sequence中mutation_position已正确

    primers = design_mutagenesis_primers(
        wt_sequence=wt_sequence,
        mutation_position=0,  # 占位, 实际由调用者填入
        target_aa=target_aa,
        wt_aa=wt_aa,
    )

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "resnum": resnum,
        "mutation": primers["mutation"],
        "wt_aa": wt_aa,
        "target_aa": target_aa,
        "esm_score": esm_score,
        "primers": {
            "forward": primers["forward_primer"],
            "reverse": primers["reverse_primer"],
        },
        "recommended_annealing_temp": primers["recommended_annealing_temp"],
        "primer_tm_diff": primers["primer_tm_diff"],
    }
