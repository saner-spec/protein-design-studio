#!/usr/bin/env python3
"""下载 ProTherm / S669 完整基准数据集

尝试多个数据源，将完整的 ProTherm 导出数据处理为验证引擎可用的 JSON 格式。

用法:
    cd E:\AI_Agents\protein_designer
    python scripts/download_benchmarks.py           # 下载到 data/benchmarks/
    python scripts/download_benchmarks.py --merge   # 合并多个数据源
    python scripts/download_benchmarks.py --source prothermdb

数据来源（按优先级）:
    1. ProThermDB (web.iitm.ac.in) — 32,000+ 条，含实验条件
    2. S669 公开副本 — 669 条精选，多 GitHub 仓库有镜像
    3. ThermoMutDB — 备选
"""

import json
import sys
import os
from pathlib import Path
from collections import defaultdict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_DIR = PROJECT_ROOT / "data" / "benchmarks"
BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)

# 标准氨基酸 3-letter → 1-letter
AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",  # 硒代甲硫氨酸 → Met
}


def try_download_s669() -> list[dict] | None:
    """从多个 GitHub 镜像尝试下载 S669 数据集"""
    import urllib.request

    urls = [
        "https://raw.githubusercontent.com/KULL-Centre/_Data/main/s669.csv",
        "https://raw.githubusercontent.com/haskellinger/pub-data/main/s669.csv",
        "https://raw.githubusercontent.com/nicolaGit/AI-protein-design/main/data/s669.csv",
    ]

    for url in urls:
        try:
            print(f"  尝试: {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "ProteinDesignStudio/1.0"})
            with urllib.request.urlopen(req, timeout=20) as f:
                text = f.read().decode("utf-8")
            if text.strip().startswith("pdb") or text.strip().startswith("PDB"):
                entries = _parse_s669_csv(text)
                print(f"  ✓ 成功下载 {len(entries)} 条 S669 数据")
                return entries
        except Exception as e:
            print(f"  ✗ 失败: {e}")
            continue

    return None


def _parse_s669_csv(text: str) -> list[dict]:
    """解析 S669 CSV 格式"""
    lines = text.strip().split("\n")
    header = lines[0].strip().lower().split(",")

    # 常见列名映射
    col_map = {}
    for i, col in enumerate(header):
        col = col.strip().lower().replace('"', "")
        if "pdb" in col:
            col_map["pdb"] = i
        elif "chain" in col:
            col_map["chain"] = i
        elif col in ("position", "pos", "resi", "resnum"):
            col_map["position"] = i
        elif col in ("wt", "wild_type", "wildtype", "wild"):
            col_map["wt"] = i
        elif col in ("mut", "mutation", "mutant", "mutant_type"):
            col_map["mut"] = i
        elif col in ("ddg", "exp_ddg", "experimental_ddg", "ddg_exp"):
            col_map["ddg"] = i
        elif col in ("ph"):
            col_map["ph"] = i
        elif col in ("temp", "temperature", "t"):
            col_map["temp"] = i

    entries = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = [p.strip().strip('"') for p in line.split(",")]
        try:
            entry = {"pdb_id": "", "chain": "A", "position": 0, "wt": "", "mut": "", "exp_ddg": 0.0}
            for key, idx in col_map.items():
                if idx >= len(parts):
                    continue
                val = parts[idx]
                if key == "pdb":
                    entry["pdb_id"] = val.lower()
                elif key == "chain":
                    entry["chain"] = val or "A"
                elif key == "position":
                    entry["position"] = int(float(val))
                elif key == "wt":
                    entry["wt"] = _normalize_aa(val)
                elif key == "mut":
                    entry["mut"] = _normalize_aa(val)
                elif key == "ddg":
                    entry["exp_ddg"] = float(val)
                elif key == "ph":
                    try:
                        entry["ph"] = float(val)
                    except ValueError:
                        entry["ph"] = 7.0
                elif key == "temp":
                    try:
                        entry["temp"] = float(val)
                    except ValueError:
                        entry["temp"] = 25.0

            if entry["pdb_id"] and entry["wt"] and entry["mut"] and entry["position"] > 0:
                entries.append(entry)
        except (ValueError, IndexError):
            continue

    return entries


def _normalize_aa(aa: str) -> str:
    aa = aa.strip().upper()
    if len(aa) == 3:
        return AA3TO1.get(aa, "X")
    if len(aa) == 1 and aa in "ACDEFGHIKLMNPQRSTVWY":
        return aa
    return "X"


def try_download_prothermdb_tsv() -> list[dict] | None:
    """从 ProThermDB 下载 TSV 导出文件（如果有直接链接）"""
    import urllib.request

    urls = [
        "https://web.iitm.ac.in/bioinfo2/prothermdb/data/prothermdb.tsv",
        "https://web.iitm.ac.in/bioinfo2/prothermdb/data/prothermdb_all.tsv",
    ]

    for url in urls:
        try:
            print(f"  尝试: {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "ProteinDesignStudio/1.0"})
            with urllib.request.urlopen(req, timeout=30) as f:
                text = f.read().decode("utf-8")
            if "PDB" in text[:200] or "UniProt" in text[:200]:
                entries = _parse_prothermdb_tsv(text)
                print(f"  ✓ 成功下载 {len(entries)} 条 ProThermDB 数据")
                return entries
        except Exception as e:
            print(f"  ✗ 失败: {e}")
            continue

    return None


def _parse_prothermdb_tsv(text: str) -> list[dict]:
    """解析 ProThermDB TSV 导出格式"""
    entries = []
    lines = text.strip().split("\n")
    header = lines[0].strip().split("\t")

    # ProThermDB 常见列
    col_idx = {}
    for i, col in enumerate(header):
        col_lower = col.lower().replace(" ", "_")
        if "pdb" in col_lower and "id" in col_lower:
            col_idx["pdb"] = i
        elif col_lower in ("chain",):
            col_idx["chain"] = i
        elif col_lower in ("position", "pos", "mutation_site"):
            col_idx["pos"] = i
        elif col_lower in ("wild_type", "wt", "wild"):
            col_idx["wt"] = i
        elif col_lower in ("mutation", "mutant", "mut"):
            col_idx["mut"] = i
        elif col_lower in ("ddg", "ddg_exp", "exp_ddG"):
            col_idx["ddg"] = i
        elif col_lower in ("ph"):
            col_idx["ph"] = i
        elif col_lower in ("temperature", "temp", "t"):
            col_idx["temp"] = i

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        try:
            entry = {"pdb_id": "", "chain": "A", "position": 0, "wt": "", "mut": "", "exp_ddg": 0.0}
            for key, idx in col_idx.items():
                if idx >= len(parts):
                    continue
                val = parts[idx].strip()
                if key == "pdb":
                    entry["pdb_id"] = val.lower()
                elif key == "chain":
                    entry["chain"] = val or "A"
                elif key == "pos":
                    try:
                        entry["position"] = int(float(val))
                    except ValueError:
                        continue
                elif key == "wt":
                    entry["wt"] = _normalize_aa(val)
                elif key == "mut":
                    entry["mut"] = _normalize_aa(val)
                elif key == "ddg":
                    try:
                        entry["exp_ddg"] = float(val)
                    except ValueError:
                        continue
                elif key == "ph":
                    try:
                        entry["ph"] = float(val)
                    except ValueError:
                        entry["ph"] = 7.0
                elif key == "temp":
                    try:
                        entry["temp"] = float(val)
                    except ValueError:
                        entry["temp"] = 25.0

            if entry["pdb_id"] and entry["wt"] and entry["mut"] and entry["position"] > 0:
                entries.append(entry)
        except (ValueError, IndexError):
            continue

    return entries


def clean_entries(entries: list[dict]) -> list[dict]:
    """清洗数据集：去重、去反向突变、去非标准氨基酸、去明显错误"""
    # 去重（同一 PDB + chain + position + mut）
    seen = set()
    cleaned = []
    for e in entries:
        wt = _normalize_aa(e["wt"])
        mut = _normalize_aa(e["mut"])
        if wt == "X" or mut == "X":
            continue
        if wt == mut:
            continue  # 跳过"突变"到自身的条目
        if abs(e.get("exp_ddg", 0)) > 15:  # 极端值可能为错误
            continue

        key = (e["pdb_id"], e.get("chain", "A"), e["position"], mut)
        if key not in seen:
            seen.add(key)
            e["wt"] = wt
            e["mut"] = mut
            cleaned.append(e)

    print(f"  清洗后: {len(entries)} → {len(cleaned)} 条")
    return cleaned


def merge_all_sources() -> list[dict]:
    """合并所有可用数据源"""
    all_entries = []

    # 1. 尝试 S669
    s669 = try_download_s669()
    if s669:
        all_entries.extend(s669)

    # 2. 尝试 ProThermDB
    protherm = try_download_prothermdb_tsv()
    if protherm:
        all_entries.extend(protherm)

    if all_entries:
        all_entries = clean_entries(all_entries)
    return all_entries


def save_dataset(name: str, description: str, entries: list[dict]) -> Path:
    path = BENCHMARK_DIR / f"{name}.json"
    data = {
        "name": name,
        "description": description,
        "entries": entries,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return path


def main():
    import argparse
    parser = argparse.ArgumentParser(description="下载蛋白质稳定性基准数据集")
    parser.add_argument("--merge", action="store_true", help="合并所有可用数据源")
    parser.add_argument("--source", choices=["s669", "prothermdb", "all"],
                        default="all", help="指定数据源")
    parser.add_argument("--output", default="protherm_full", help="输出文件名（不含扩展名）")
    args = parser.parse_args()

    print("Protein Design Studio — 基准数据下载器\n")

    if args.merge or args.source == "all":
        entries = merge_all_sources()
    elif args.source == "s669":
        entries = try_download_s669() or []
    elif args.source == "prothermdb":
        entries = try_download_prothermdb_tsv() or []
    else:
        entries = []

    if not entries:
        print("\n  ✗ 所有数据源均无法访问。")
        print("    请手动从 ProThermDB 下载: https://web.iitm.ac.in/bioinfo2/prothermdb/Downloads.html")
        print("    填入表单后下载 TSV，然后用以下命令转换:")
        print(f"    python {__file__} --source prothermdb")
        sys.exit(1)

    path = save_dataset(args.output, "ProTherm 完整数据集（自动下载）", entries)
    print(f"\n  ✓ 已保存 {len(entries)} 条记录到 {path}")
    print(f"  │  后续在 Protein Design Studio 中运行验证: ")
    print(f"  │  GET /api/validation/{args.output}")
    print()


if __name__ == "__main__":
    main()
