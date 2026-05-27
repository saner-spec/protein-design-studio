"""Protein Design Studio — PDB解析与结构特征提取

CC审核关键决策: 序列从SEQRES提取（完整），ATOM仅用于3D坐标。
维护三层映射: seqres_idx → pdb_atom_resnum → esm_token_idx

v0.2.0 Phase 2 新增:
  - 二面角 (phi/psi) 计算
  - DSSP 二级结构指派 (基于主链氢键模式)
  - 序列设计工作流支持
"""
import math
import httpx
from pathlib import Path
from collections import defaultdict

from .config import PDB_CACHE, HBOND_DIST_CUTOFF, DEFAULT_SASA_THRESHOLD
from .residue_map import normalize_residue, get_aa1

# 原子共价半径 (Å)
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
    "S": 1.05, "P": 1.07,
}
SOLVENT_RADIUS = 1.4


async def fetch_pdb(pdb_id: str) -> Path:
    """从RCSB下载PDB文件（异步），缓存到本地"""
    pdb_id = pdb_id.lower()
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    cache_path = PDB_CACHE / f"{pdb_id}.pdb"

    if not cache_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.get(url)
            resp.raise_for_status()
            cache_path.write_bytes(resp.content)

    return cache_path


def parse_pdb(pdb_path: str | Path) -> dict:
    """解析PDB文件，返回结构数据"""
    pdb_path = Path(pdb_path)
    pdb_id = pdb_path.stem.upper()

    with open(pdb_path) as f:
        lines = f.readlines()

    # ── SEQRES 完整序列 ──────────────────────────────────
    seqres_records = defaultdict(list)
    for line in lines:
        if line.startswith("SEQRES"):
            chain = line[11]
            resnames = line[19:].split()
            seqres_records[chain].extend(resnames)

    seqres_sequences = {}
    for chain, resnames in seqres_records.items():
        seq = []
        for rn in resnames:
            aa = get_aa1(normalize_residue(rn) or "")
            if aa:
                seq.append(aa)
        seqres_sequences[chain] = "".join(seq)

    # ── ATOM / HETATM 解析 ──────────────────────────────
    atoms = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atoms.append({
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "altloc": line[16],
                "resname": line[17:20].strip(),
                "chain": line[21],
                "resnum": int(line[22:26]),
                "inscode": line[26],
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "occupancy": float(line[54:60]) if line[54:60].strip() else 1.0,
                "element": line[76:78].strip() or line[12:16].strip()[0],
                "is_hetatm": line.startswith("HETATM"),
            })

    protein_atoms, ligand_atoms = [], []
    for a in atoms:
        norm = normalize_residue(a["resname"])
        if a["is_hetatm"] and norm is None:
            ligand_atoms.append(a)
        elif norm:
            protein_atoms.append(a)

    # ── 按残基分组 (蛋白部分) ──────────────────────────
    residues_dict = {}
    for a in protein_atoms:
        key = (a["chain"], a["resnum"], a["inscode"])
        if key not in residues_dict:
            aa = get_aa1(normalize_residue(a["resname"]) or "")
            residues_dict[key] = {
                "chain": a["chain"],
                "resnum": a["resnum"],
                "inscode": a["inscode"],
                "resname": normalize_residue(a["resname"]),
                "aa": aa,
                "atoms": [],
            }
        residues_dict[key]["atoms"].append(a)

    residues = sorted(residues_dict.values(), key=lambda r: (r["chain"], r["resnum"]))

    # ── SEQRES对齐: 局部序列比对 ────────────────────────
    for r in residues:
        r["seqres_idx"] = -1

    for chain in set(r["chain"] for r in residues):
        seq = seqres_sequences.get(chain, "")
        if not seq:
            continue
        chain_res = [r for r in residues if r["chain"] == chain]
        atom_seq = "".join(r["aa"] or "X" for r in chain_res)

        # 使用简单但更鲁棒的策略: 在SEQRES中找ATOM序列的每个残基
        if atom_seq in seq:
            idx = seq.find(atom_seq)
            for i, r in enumerate(chain_res):
                r["seqres_idx"] = idx + i
        else:
            # fallback: 逐个残基在SEQRES中匹配 (处理缺失残基)
            seq_idx = 0
            atom_idx = 0
            while seq_idx < len(seq) and atom_idx < len(chain_res):
                if seq[seq_idx] == chain_res[atom_idx]["aa"]:
                    chain_res[atom_idx]["seqres_idx"] = seq_idx
                    seq_idx += 1
                    atom_idx += 1
                else:
                    seq_idx += 1  # SEQRES中有但ATOM中无 → 缺失残基, 跳过
            # 剩余残基回退到残基编号近似
            for r in chain_res:
                if r["seqres_idx"] < 0:
                    r["seqres_idx"] = r["resnum"] - 1

    # ── 结构特征 ─────────────────────────────────────────
    # Cα坐标: (chain, resnum) → (x, y, z)
    ca_atoms = {}
    # 主链原子: N, CA, C (用于二面角)
    backbone_atoms = {}
    for r in residues:
        key = (r["chain"], r["resnum"])
        bb = {}
        for a in r["atoms"]:
            if a["name"] == "CA":
                ca_atoms[key] = (a["x"], a["y"], a["z"])
                bb["CA"] = (a["x"], a["y"], a["z"])
            elif a["name"] == "N":
                bb["N"] = (a["x"], a["y"], a["z"])
            elif a["name"] == "C":
                bb["C"] = (a["x"], a["y"], a["z"])
        backbone_atoms[key] = bb

    # SASA: (chain, resnum) → area
    sasa_values = _calc_sasa(residues)

    SASA_REF = {
        "A": 113, "R": 241, "N": 158, "D": 151, "C": 140,
        "Q": 189, "E": 183, "G": 85, "H": 194, "I": 182,
        "L": 180, "K": 211, "M": 204, "F": 218, "P": 143,
        "S": 122, "T": 146, "W": 259, "Y": 229, "V": 160,
    }

    for i, r in enumerate(residues):
        rkey = (r["chain"], r["resnum"])

        # SASA
        r["sasa"] = sasa_values.get(rkey, 0.0)
        ref_sasa = SASA_REF.get(r["aa"] or "", 150)
        r["rel_sasa"] = r["sasa"] / ref_sasa if ref_sasa > 0 else 0.0

        # 二面角 phi / psi (Phase 2)
        phi, psi = _calc_torsion_angles(i, residues, backbone_atoms)
        r["phi"] = phi
        r["psi"] = psi

        # 二级结构 — Phase 2: 基于phi/psi的DSSP风格指派
        r["secondary_structure"] = _assign_ss_from_angles(phi, psi, r["aa"])

        # 邻居 (<5Å Cα距离)
        neighbors = []
        ca_i = ca_atoms.get(rkey)
        if ca_i:
            for r2 in residues:
                if r["chain"] == r2["chain"] and r["resnum"] == r2["resnum"]:
                    continue
                ca_j = ca_atoms.get((r2["chain"], r2["resnum"]))
                if ca_j:
                    d = math.sqrt(sum((a - b) ** 2 for a, b in zip(ca_i, ca_j)))
                    if d < 5.0:
                        neighbors.append(r2["resnum"])
        r["neighbors"] = neighbors

        # 氢键 (N-O < 3.5Å)
        hbonds = []
        for a1 in r["atoms"]:
            if a1["element"] not in ("N", "O"):
                continue
            for r2 in residues:
                if r["chain"] == r2["chain"] and r["resnum"] == r2["resnum"]:
                    continue
                if r["chain"] == r2["chain"] and abs(r["resnum"] - r2["resnum"]) < 2:
                    continue
                for a2 in r2["atoms"]:
                    if a2["element"] not in ("N", "O") or a1["element"] == a2["element"]:
                        continue
                    d = math.sqrt(
                        (a1["x"] - a2["x"]) ** 2
                        + (a1["y"] - a2["y"]) ** 2
                        + (a1["z"] - a2["z"]) ** 2
                    )
                    if d < HBOND_DIST_CUTOFF:
                        donor = f"{r['resname']}{r['resnum']} {a1['name']}"
                        acceptor = f"{r2['resname']}{r2['resnum']} {a2['name']}"
                        hbonds.append({"donor": donor, "acceptor": acceptor, "dist": round(d, 2)})
        r["hbonds"] = hbonds[:20]

    # ── Phase 2: DSSP最终优化 — 基于氢键模式的二级结构重指派 ──
    _refine_ss_from_hbonds(residues, ca_atoms, backbone_atoms)

    # ── 配体 ────────────────────────────────────────────
    ligands = []
    lig_by_res = defaultdict(list)
    for a in ligand_atoms:
        key = (a["chain"], a["resnum"])
        lig_by_res[key].append(a)
    for key, lig_atoms in lig_by_res.items():
        ligands.append({
            "chain": key[0],
            "resnum": key[1],
            "resname": lig_atoms[0]["resname"],
            "atoms": lig_atoms,
        })

    # ── 主序列: 取第一个蛋白链 ──────────────────────────
    main_chain = residues[0]["chain"] if residues else "A"
    main_seq = seqres_sequences.get(main_chain, "")

    return {
        "pdb_id": pdb_id,
        "seqres_sequence": main_seq,
        "residue_count": len(residues),         # 实际解析到的残基数量 (有3D坐标的)
        "seqres_length": len(main_seq),         # SEQRES序列全长 (包含未解析残基)
        "chains": list(seqres_sequences.keys()),
        "residues": residues,
        "ligands": ligands,
    }


# ── Phase 2: 二面角计算 ─────────────────────────────────────

def _calc_torsion_angles(
    idx: int,
    residues: list,
    backbone_atoms: dict,
) -> tuple:
    """计算残基的phi和psi二面角

    phi: C(i-1) - N(i) - CA(i) - C(i)
    psi: N(i) - CA(i) - C(i) - N(i+1)

    Returns:
        (phi_in_degrees, psi_in_degrees) 或 (None, None) 如无法计算
    """
    r = residues[idx]
    rkey = (r["chain"], r["resnum"])
    bb = backbone_atoms.get(rkey, {})

    # ---- phi: C(i-1) - N(i) - CA(i) - C(i) ----
    phi = None
    if "N" in bb and "CA" in bb and "C" in bb:
        if idx > 0:
            prev = residues[idx - 1]
            prev_key = (prev["chain"], prev["resnum"])
            prev_bb = backbone_atoms.get(prev_key, {})
            if prev["chain"] == r["chain"] and "C" in prev_bb:
                phi = _dihedral_angle(
                    prev_bb["C"], bb["N"], bb["CA"], bb["C"]
                )

    # ---- psi: N(i) - CA(i) - C(i) - N(i+1) ----
    psi = None
    if "N" in bb and "CA" in bb and "C" in bb:
        if idx < len(residues) - 1:
            nxt = residues[idx + 1]
            nxt_key = (nxt["chain"], nxt["resnum"])
            nxt_bb = backbone_atoms.get(nxt_key, {})
            if nxt["chain"] == r["chain"] and "N" in nxt_bb:
                psi = _dihedral_angle(
                    bb["N"], bb["CA"], bb["C"], nxt_bb["N"]
                )

    return phi, psi


def _dihedral_angle(p1, p2, p3, p4) -> float | None:
    """计算四个点定义的二面角 (返回度数)"""
    if None in p1 or None in p2 or None in p3 or None in p4:
        return None

    # 向量
    b1 = tuple(p2[i] - p1[i] for i in range(3))
    b2 = tuple(p3[i] - p2[i] for i in range(3))
    b3 = tuple(p4[i] - p3[i] for i in range(3))

    # 法向量
    def cross(v1, v2):
        return (
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        )

    def dot(v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    def norm(v):
        return math.sqrt(dot(v, v))

    n1 = cross(b1, b2)
    n2 = cross(b2, b3)

    n1_norm = norm(n1)
    n2_norm = norm(n2)
    if n1_norm < 1e-10 or n2_norm < 1e-10:
        return None

    cos_angle = dot(n1, n2) / (n1_norm * n2_norm)
    cos_angle = max(-1.0, min(1.0, cos_angle))
    angle = math.degrees(math.acos(cos_angle))

    # 符号: 通过b2与(n1 × n2)的点积判断
    sign = dot(cross(n1, n2), b2)
    if sign < 0:
        angle = -angle

    return round(angle, 1)


# ── Phase 2: DSSP风格二级结构指派 ─────────────────────────

def _assign_ss_from_angles(phi: float | None, psi: float | None, aa: str | None) -> str:
    """基于phi/psi二面角的简单二级结构指派

    使用Ramachandran图规则:
      H (α-螺旋):  phi ~ -57°, psi ~ -47°
      E (β-折叠):  phi ~ -130°, psi ~ +130°
      C (环区/卷曲): 其他

    Returns:
        "H", "E", 或 "C"
    """
    if phi is None or psi is None:
        return "C"

    # 允许±40°的容差
    # α-螺旋区域
    if (-97 <= phi <= -17) and (-87 <= psi <= -7):
        return "H"

    # β-折叠区域 (平行 + 反平行)
    # 反平行: phi ~ -139°, psi ~ +135°
    # 平行:   phi ~ -119°, psi ~ +113°
    if (-170 <= phi <= -90) and (90 <= psi <= 180):
        return "E"
    if (-170 <= phi <= -90) and (-180 <= psi <= -160):
        return "E"

    # 3-10螺旋 (phi ~ -49°, psi ~ -26°)
    if (-89 <= phi <= -9) and (-66 <= psi <= 14):
        return "H"

    return "C"


def _refine_ss_from_hbonds(residues: list, ca_atoms: dict, backbone_atoms: dict):
    """通过主链氢键模式优化二级结构指派 (DSSP核心逻辑)

    DSSP规则:
      α-螺旋: 残基i的CO与残基i+4的NH形成氢键
      3-10螺旋: 残基i的CO与残基i+3的NH形成氢键
      β-折叠: 链间/链内氢键对
    """
    # 构建残基索引映射: {chain: [(resnum, idx_in_list), ...]}
    chain_residues = defaultdict(list)
    for i, r in enumerate(residues):
        chain_residues[r["chain"]].append((r, i))

    # 对每条链计算氢键模式
    for chain, res_list in chain_residues.items():
        n = len(res_list)
        if n < 5:
            continue

        # 检测α-螺旋 (i, i+4 氢键)
        helix_scores = [0] * n
        for i in range(n - 4):
            r_i = res_list[i][0]
            r_i4 = res_list[i + 4][0]
            # 检查C=O(i)与N-H(i+4)之间是否存在氢键
            if _has_backbone_hbond(r_i, r_i4, "C", "O", "N"):
                helix_scores[i] += 1
                helix_scores[i + 4] += 1

        # 检测3-10螺旋 (i, i+3 氢键)
        for i in range(n - 3):
            r_i = res_list[i][0]
            r_i3 = res_list[i + 3][0]
            if _has_backbone_hbond(r_i, r_i3, "C", "O", "N"):
                helix_scores[i] += 1
                helix_scores[i + 3] += 1

        # 检测β-折叠: 寻找i→j和j→i的平行/反平行配对
        sheet_candidates = set()
        for i in range(n - 2):
            for j in range(i + 3, n):
                r_i = res_list[i][0]
                r_j = res_list[j][0]
                # 反平行: i→j 且 j→i
                hb1 = _has_backbone_hbond(r_i, r_j, "C", "O", "N")
                hb2 = _has_backbone_hbond(r_j, r_i, "C", "O", "N")
                # 平行: i→j 且 i→j-... 更复杂的模式
                hb3 = _has_backbone_hbond(r_i, r_j, "C", "O", "N")
                hb4 = _has_backbone_hbond(r_j, r_i, "N", "H", "O")
                if (hb1 and hb2) or (hb3 and hb4):
                    sheet_candidates.add(i)
                    sheet_candidates.add(j)

        # 指派: 优先α-螺旋, 其次β-折叠, 最后环区
        for i in range(n):
            r = res_list[i][0]
            if helix_scores[i] >= 2:
                r["secondary_structure"] = "H"
            elif i in sheet_candidates:
                r["secondary_structure"] = "E"
            # 如果原phi/psi指派为H但得分不足, 保持环区C


def _has_backbone_hbond(res_i: dict, res_j: dict, atom_i_name: str, atom_i_elem: str, atom_j_name: str) -> bool:
    """检查两个残基之间是否存在主链氢键"""
    # 找供体原子 (res_i中的atom_i_name)
    donor = None
    for a in res_i["atoms"]:
        if a["name"] == atom_i_name:
            donor = (a["x"], a["y"], a["z"])
            break
    if not donor:
        return False

    # 找受体原子 (res_j中的含atom_j_elem的元素)
    for a in res_j["atoms"]:
        if a["name"] == atom_j_name or (a["element"] == atom_j_name and a["name"] in ("N", "O")):
            acceptor = (a["x"], a["y"], a["z"])
            d = math.sqrt(
                (donor[0] - acceptor[0]) ** 2
                + (donor[1] - acceptor[1]) ** 2
                + (donor[2] - acceptor[2]) ** 2
            )
            if d < HBOND_DIST_CUTOFF:
                return True

    return False


# ── SASA 计算 ──────────────────────────────────────────────

def _calc_sasa(residues: list) -> dict:
    """简化的Shrake-Rupley SASA计算 — key为(chain, resnum)"""
    all_atoms = []
    for r in residues:
        for a in r["atoms"]:
            radius = COVALENT_RADII.get(a["element"], 1.7)
            all_atoms.append({
                "key": (r["chain"], r["resnum"]),
                "x": a["x"], "y": a["y"], "z": a["z"],
                "radius": radius,
            })

    sasa = {}
    n_points = 10  # 平衡计算速度和精度 (原20点, 对多链结构减半)
    for i, a in enumerate(all_atoms):
        r = a["radius"] + SOLVENT_RADIUS
        accessible = 0
        for k in range(n_points):
            phi = math.acos(2.0 * ((k + 0.5) / n_points) - 1.0)
            theta = math.pi * (1.0 + 5.0 ** 0.5) * k
            dx = r * math.sin(phi) * math.cos(theta)
            dy = r * math.sin(phi) * math.sin(theta)
            dz = r * math.cos(phi)
            px, py, pz = a["x"] + dx, a["y"] + dy, a["z"] + dz

            blocked = False
            for j, b in enumerate(all_atoms):
                if i == j:
                    continue
                d2 = (px - b["x"]) ** 2 + (py - b["y"]) ** 2 + (pz - b["z"]) ** 2
                if d2 < (b["radius"] + SOLVENT_RADIUS) ** 2:
                    blocked = True
                    break
            if not blocked:
                accessible += 1

        area = 4.0 * math.pi * r**2 * accessible / n_points
        key = a["key"]
        sasa[key] = sasa.get(key, 0) + area

    return sasa


def get_residue_info(structure: dict, chain: str, resnum: int) -> dict | None:
    for r in structure["residues"]:
        if r["chain"] == chain and r["resnum"] == resnum:
            return r
    return None


def is_surface_residue(residue: dict, threshold: float = DEFAULT_SASA_THRESHOLD) -> bool:
    return residue.get("rel_sasa", 0) > threshold
