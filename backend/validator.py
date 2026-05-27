"""Protein Design Studio — ESM 验证引擎

对基准数据集（S669 / ProTherm 衍生）运行 ESM-2 评分，
计算 Spearman ρ、Pearson r、方向准确率，提供验证面板数据。

设计原则:
  - 只读验证，不校准预测输出
  - 数据集打包在仓库中，不依赖外部服务器
  - 首次运行预计算 ESM 分数并缓存，后续即时返回
"""

import json
import math
import time
from pathlib import Path
from collections import defaultdict
from typing import Optional

from .config import DATA_DIR
from .pdb_parser import parse_pdb, fetch_pdb
from .esm_scorer import score_single_position, load_model
from .residue_map import STANDARD_AAS


BENCHMARK_DIR = DATA_DIR / "benchmarks"
RESULT_CACHE_DIR = DATA_DIR / "validation_cache"


def _load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def _save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


# ─────────────────── 统计指标 ───────────────────


def spearman_rho(xs: list[float], ys: list[float]) -> float:
    """Spearman 秩相关系数"""
    n = len(xs)
    if n < 3:
        return 0.0

    def rank(vals):
        sorted_vals = sorted((v, i) for i, v in enumerate(vals))
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n and abs(sorted_vals[j][0] - sorted_vals[i][0]) < 1e-12:
                j += 1
            avg_rank = (i + j - 1) / 2.0 + 1.0
            for k in range(i, j):
                ranks[sorted_vals[k][1]] = avg_rank
            i = j
        return ranks

    x_ranks = rank(xs)
    y_ranks = rank(ys)
    mean_xr = sum(x_ranks) / n
    mean_yr = sum(y_ranks) / n

    cov = sum((x_ranks[i] - mean_xr) * (y_ranks[i] - mean_yr) for i in range(n))
    sx = math.sqrt(sum((r - mean_xr) ** 2 for r in x_ranks))
    sy = math.sqrt(sum((r - mean_yr) ** 2 for r in y_ranks))

    if sx * sy < 1e-12:
        return 0.0
    return cov / (sx * sy)


def pearson_r(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx * sy < 1e-12:
        return 0.0
    return cov / (sx * sy)


# ─────────────────── 基准数据集 ───────────────────


class BenchmarkDataset:
    """基准验证数据集"""

    def __init__(self, name: str, description: str = ""):
        self.name = name
        self.description = description
        self.entries: list[dict] = []

    @classmethod
    def from_json(cls, path: Path) -> "BenchmarkDataset":
        data = _load_json(path)
        ds = cls(data["name"], data.get("description", ""))
        ds.entries = data.get("entries", [])
        return ds

    def to_json(self, path: Path) -> None:
        _save_json(path, {
            "name": self.name,
            "description": self.description,
            "entries": self.entries,
        })

    def unique_positions(self) -> list[tuple[str, str, int]]:
        """返回唯一的 (pdb_id, chain, position) 三元组"""
        seen = set()
        result = []
        for e in self.entries:
            key = (e["pdb_id"].lower(), e.get("chain", "A"), e["position"])
            if key not in seen:
                seen.add(key)
                result.append(key)
        return result


# ─────────────────── ESM 验证运行器 ───────────────────


class ValidationRunner:
    """对基准数据集运行 ESM 评分并计算验证指标"""

    def __init__(self, dataset: BenchmarkDataset):
        self.dataset = dataset
        self.esm_scores: dict[tuple, dict] = {}  # (pdb, chain, pos) -> {aa: score}
        self.results: list[dict] = []  # per-entry results
        self._cache_path = RESULT_CACHE_DIR / f"{dataset.name}_esm_scores.json"

    def load_cache(self) -> bool:
        if self._cache_path.exists():
            self.esm_scores = {
                tuple(k): v
                for k, v in _load_json(self._cache_path).items()
            }
            return len(self.esm_scores) > 0
        return False

    def save_cache(self) -> None:
        data = {str(list(k)): v for k, v in self.esm_scores.items()}
        _save_json(self._cache_path, data)

    async def run(self, progress_callback=None) -> dict:
        """运行验证——对每个唯一位点跑 ESM 打分，然后计算指标"""
        # 尝试加载缓存
        if self.load_cache():
            pass  # 使用已有缓存
        else:
            # 需要加载 ESM 模型
            load_model()

            all_positions = self.dataset.unique_positions()
            total = len(all_positions)

            for i, (pdb_id, chain, pos) in enumerate(all_positions):
                key = (pdb_id, chain, pos)
                if key in self.esm_scores:
                    continue

                if progress_callback:
                    progress_callback(i, total, f"{pdb_id}:{chain}:{pos}")

                try:
                    scores = await self._score_position(pdb_id, chain, pos)
                    if scores:
                        self.esm_scores[key] = scores
                except Exception:
                    continue  # 个别PDB失败不阻塞整体流程

            self.save_cache()

            if progress_callback:
                progress_callback(total, total, "计算指标中...")

        # 构建逐条结果
        self._build_results()
        return self.compute_metrics()

    async def _score_position(self, pdb_id: str, chain: str, position: int) -> Optional[dict]:
        """对单个位点跑 ESM 打分（先确保PDB有本地缓存）"""
        pdb_id_upper = pdb_id.upper()
        local_path = Path(DATA_DIR / ".." / "static" / "pdbs" / f"{pdb_id_upper}.pdb")

        if not local_path.exists():
            # 尝试从 RCSB 下载
            from .pdb_parser import PDB_CACHE
            cached = PDB_CACHE / f"{pdb_id.lower()}.pdb"
            if not cached.exists():
                try:
                    cached = await fetch_pdb(pdb_id)
                except Exception:
                    return None  # 下载失败，跳过此条
            local_path = cached

        struct = parse_pdb(local_path)

        if chain not in struct.get("chains", []):
            return None

        sequence = struct.get("seqres_sequence", "")
        if not sequence:
            return None

        # PDB编号 → 0-based序列索引 (通过 residues 列表中的 seqres_idx)
        esm_idx = None
        for residue in struct.get("residues", []):
            if residue.get("chain") == chain and residue.get("resnum") == position:
                esm_idx = residue.get("seqres_idx")
                break

        if esm_idx is None or esm_idx >= len(sequence):
            return None

        try:
            scores = score_single_position(sequence, esm_idx)
            return {aa: round(s, 4) for aa, s in scores.items()}
        except Exception:
            return None

    def _build_results(self) -> None:
        """构建逐条验证结果。自动跳过 PDB 残基与数据集 WT 不匹配的条目。"""
        self.results = []
        self.skipped: list[dict] = []
        for entry in self.dataset.entries:
            key = (entry["pdb_id"].lower(), entry.get("chain", "A"), entry["position"])
            scores = self.esm_scores.get(key)
            if scores is None:
                self.skipped.append({**entry, "reason": "PDB未加载或无ESM缓存"})
                continue

            wt_entry = entry["wt"].upper()
            mut = entry["mut"].upper()

            # 找出PDB中的实际残基（ESM分数为0.0的那个）
            actual_aa_in_pdb = None
            for aa, s in scores.items():
                if abs(s) < 1e-10:
                    actual_aa_in_pdb = aa
                    break

            if actual_aa_in_pdb and actual_aa_in_pdb != wt_entry:
                self.skipped.append({
                    **entry,
                    "reason": f"PDB残基({actual_aa_in_pdb})与数据集WT({wt_entry})不匹配",
                    "pdb_actual": actual_aa_in_pdb,
                })
                continue

            if wt_entry not in scores or mut not in scores:
                self.skipped.append({**entry, "reason": f"氨基酸{wt_entry}或{mut}不在ESM词汇表中"})
                continue

            esm_diff = scores[mut]  # logP(mut) - logP(wt)
            exp_ddg = entry["exp_ddg"]

            self.results.append({
                **entry,
                "esm_wt": scores.get(wt_entry, 0.0),
                "esm_mut": scores[mut],
                "esm_diff": esm_diff,
                "direction_match": (esm_diff < 0) == (exp_ddg < 0),
            })

    def compute_metrics(self) -> dict:
        """计算完整验证指标"""
        if len(self.results) < 5:
            return {"error": "数据点不足（需要≥5个有效条目）", "n": len(self.results)}

        esm_diffs = [r["esm_diff"] for r in self.results]
        exp_ddgs = [r["exp_ddg"] for r in self.results]

        # 整体指标
        matches = sum(1 for r in self.results if r["direction_match"])
        direction_accuracy = matches / len(self.results)

        sp = spearman_rho(esm_diffs, exp_ddgs)
        pr = pearson_r(esm_diffs, exp_ddgs)

        # 分类指标
        tp = sum(1 for r in self.results if r["esm_diff"] < 0 and r["exp_ddg"] < 0)  # 正确预测稳定
        tn = sum(1 for r in self.results if r["esm_diff"] > 0 and r["exp_ddg"] > 0)  # 正确预测不稳定
        fp = sum(1 for r in self.results if r["esm_diff"] < 0 and r["exp_ddg"] >= 0)  # 假稳定
        fn = sum(1 for r in self.results if r["esm_diff"] >= 0 and r["exp_ddg"] < 0)  # 假不稳定

        total = tp + tn + fp + fn
        stabilizing_accuracy = tp / (tp + fn) if (tp + fn) > 0 else 0
        destabilizing_accuracy = tn / (tn + fp) if (tn + fp) > 0 else 0

        # 按残基类型拆分
        by_residue_type = self._breakdown_by_residue_type()

        # 按 ESMBin 拆分（用于散点图分组趋势）
        by_esm_bin = self._breakdown_by_esm_bin()

        # 散点图数据
        scatter = [
            {
                "pdb_id": r["pdb_id"],
                "position": r["position"],
                "wt": r["wt"],
                "mut": r["mut"],
                "esm_diff": round(r["esm_diff"], 4),
                "exp_ddg": round(r["exp_ddg"], 3),
                "correct": r["direction_match"],
            }
            for r in self.results
        ]

        return {
            "dataset": self.dataset.name,
            "description": self.dataset.description,
            "n_entries": len(self.dataset.entries),
            "n_valid": len(self.results),
            "n_skipped": len(self.skipped),
            "n_unique_positions": len(self.esm_scores),
            "skipped_summary": self._summarize_skipped(),
            "metrics": {
                "spearman_rho": round(sp, 3),
                "pearson_r": round(pr, 3),
                "direction_accuracy": round(direction_accuracy, 3),
                "stabilizing_sensitivity": round(stabilizing_accuracy, 3),
                "destabilizing_specificity": round(destabilizing_accuracy, 3),
                "confusion_matrix": {
                    "tp": tp, "tn": tn, "fp": fp, "fn": fn,
                    "total": total,
                },
            },
            "by_residue_type": by_residue_type,
            "by_esm_bin": by_esm_bin,
            "scatter": scatter,
            "interpretation": {
                "spearman": _interpret_spearman(sp),
                "direction": f"ESM-2 在 {direction_accuracy*100:.0f}% 的突变中正确预测了方向（稳定化/去稳定化）",
                "caveat": "ESM-2 的 log-likelihood 与实验 ΔΔG 的相关性为中等水平（ρ≈0.4-0.6），"
                          "不能替代实验测量。验证数据来自 ProTherm 数据库，实验条件各异。",
            },
        }

    def _breakdown_by_residue_type(self) -> dict:
        """按野生型残基类型拆分方向准确率"""
        groups = defaultdict(list)
        for r in self.results:
            groups[r["wt"]].append(r["direction_match"])

        result = {}
        aa_groups = {
            "疏水": "AILMFWYV",
            "极性": "STNQ",
            "带电": "DEKRH",
            "特殊": "CGP",
        }
        for group_name, aas in aa_groups.items():
            matches = []
            for aa in aas:
                matches.extend(groups.get(aa, []))
            if matches:
                result[group_name] = {
                    "n": len(matches),
                    "accuracy": round(sum(matches) / len(matches), 3),
                }

        # 单个残基
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa in groups and len(groups[aa]) >= 3:
                result[aa] = {
                    "n": len(groups[aa]),
                    "accuracy": round(sum(groups[aa]) / len(groups[aa]), 3),
                }

        return result

    def _summarize_skipped(self) -> dict:
        """汇总跳过原因"""
        reasons = defaultdict(lambda: {"count": 0, "examples": []})
        for s in self.skipped:
            reason = s.get("reason", "未知")
            reasons[reason]["count"] += 1
            if len(reasons[reason]["examples"]) < 3:
                reasons[reason]["examples"].append(
                    f"{s['pdb_id']} {s.get('chain','A')}:{s['position']} {s['wt']}->{s['mut']}"
                )
        return {k: v for k, v in sorted(reasons.items(), key=lambda x: -x[1]["count"])}

    def _breakdown_by_esm_bin(self) -> list[dict]:
        """按 ESM 分数区间统计"""
        bins = [(-100, -3), (-3, -1), (-1, 0), (0, 1), (1, 3), (3, 100)]
        result = []
        for lo, hi in bins:
            in_bin = [r for r in self.results if lo <= r["esm_diff"] < hi]
            if not in_bin:
                continue
            avg_exp = sum(r["exp_ddg"] for r in in_bin) / len(in_bin)
            correct = sum(1 for r in in_bin if r["direction_match"])
            result.append({
                "esm_range": f"{lo} ~ {hi}",
                "n": len(in_bin),
                "avg_exp_ddg": round(avg_exp, 3),
                "accuracy": round(correct / len(in_bin), 3),
            })
        return result


def _interpret_spearman(rho: float) -> str:
    if abs(rho) >= 0.7:
        return "强相关 — ESM 评分与实验数据高度一致"
    elif abs(rho) >= 0.4:
        return "中等相关 — ESM 评分可作为筛选参考，但需实验验证"
    elif abs(rho) >= 0.2:
        return "弱相关 — ESM 评分方向大致可用，单点预测不可靠"
    else:
        return "极弱相关 — 当前模型/数据集上 ESM 评分不具预测力"


# ─────────────────── 数据集管理 ───────────────────


def list_benchmarks() -> list[dict]:
    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    benchmarks = []
    for path in sorted(BENCHMARK_DIR.glob("*.json")):
        try:
            data = _load_json(path)
            benchmarks.append({
                "name": data.get("name", path.stem),
                "file": path.stem,
                "n_entries": len(data.get("entries", [])),
                "description": data.get("description", ""),
            })
        except Exception:
            pass
    return benchmarks


def get_benchmark(name: str) -> Optional[BenchmarkDataset]:
    path = BENCHMARK_DIR / f"{name}.json"
    if path.exists():
        return BenchmarkDataset.from_json(path)
    return None


async def run_validation(benchmark_name: str = "s669_mini", force_rerun: bool = False) -> dict:
    """便捷接口——加载数据集 + 运行验证"""
    ds = get_benchmark(benchmark_name)
    if ds is None:
        available = [b["file"] for b in list_benchmarks()]
        return {"error": f"数据集 '{benchmark_name}' 不存在", "available": available}

    runner = ValidationRunner(ds)

    if force_rerun:
        runner._cache_path.unlink(missing_ok=True)
        runner.esm_scores = {}

    return await runner.run()
