"""Protein Design Studio — 非标准残基映射表

PDB中常见非标准残基 → 标准20种氨基酸的映射。
ESM tokenizer 只认识标准氨基酸，PDB中的修饰残基必须先转换。
"""

# 三字母 → 一字母映射
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 非标准 → 标准残基映射
# 格式: "PDB_CODE": "STANDARD_3LETTER"
NONSTANDARD_MAP = {
    # 硒代甲硫氨酸 (用于晶体学相位) → 甲硫氨酸
    "MSE": "MET",
    # 组氨酸质子化态变体 → 组氨酸
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS",
    # 半胱氨酸修饰
    "CYM": "CYS", "CYX": "CYS", "CSO": "CYS", "CSD": "CYS", "CSX": "CYS",
    "CME": "CYS", "SMC": "CYS",
    # 天冬氨酸/谷氨酸质子化态
    "ASH": "ASP", "GLH": "GLU",
    # 赖氨酸修饰
    "LYN": "LYS", "MLY": "LYS", "MLZ": "LYS", "LLP": "LYS",
    # 磷酸化残基
    "SEP": "SER", "TPO": "THR", "PTR": "TYR",
    # 乙酰化
    "ALY": "LYS",
    # N/C端修饰
    "ACE": None,  # 乙酰化N端 — 跳过
    "NH2": None,  # 酰胺化C端 — 跳过
    "NME": None,  # N-甲基 — 跳过
    # 其他常见
    "CSW": "CYS", "CSS": "CYS", "OCS": "CYS",
    "PCA": "GLN",  # 焦谷氨酸 (N端GLN环化)
    "CGU": "GLU",  # γ-羧基谷氨酸
    "M3L": "LYS",  # 三甲基赖氨酸
    "HYP": "PRO",  # 羟脯氨酸
    "ORN": "ARG",  # 鸟氨酸 → 近似为精氨酸
    "DAL": "ALA", "DLE": "LEU", "DSG": "ASN",  # D-氨基酸 → L型近似
    "TYS": "TYR",  # 酪氨酸磺酸化
    "MLU": "LEU",  # 甲基亮氨酸
    "FME": "MET",  # N-甲酰甲硫氨酸
    "KCX": "LYS",  # 赖氨酸羧化
    "UNK": None,  # 未知残基 → 跳过, 不映射为GLY
    "ASX": "ASN",  # 天冬酰胺/天冬氨酸模糊 → 保守映射为ASN
    "GLX": "GLN",  # 谷氨酰胺/谷氨酸模糊 → 保守映射为GLN
    "SEC": "CYS",  # 硒代半胱氨酸(第21种氨基酸) → 近似为CYS
}

# 需要跳过的残基类型 (水/离子/配体/溶剂)
SKIP_RESNAMES = {
    "HOH", "WAT", "DOD",  # 水
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CO", "NI", "CU",  # 离子
    "SO4", "PO4", "EDO", "GOL", "ACT", "EPE", "PEG", "MPD", "TRS",  # 缓冲液
    "CRO",  # GFP发色团 — 保留但标记为配体
}


def normalize_residue(resname: str) -> str | None:
    """将PDB残基名映射为标准三字母码，或返回None表示应跳过"""
    resname = resname.strip().upper()
    if resname in SKIP_RESNAMES:
        return None
    if resname in AA3TO1:
        return resname
    if resname in NONSTANDARD_MAP:
        return NONSTANDARD_MAP[resname]
    return None  # 未知残基，保守跳过


def get_aa1(resname: str) -> str | None:
    """三字母 → 单字母，或返回None"""
    return AA3TO1.get(resname)


# 20种标准氨基酸 (单字母)
STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")
