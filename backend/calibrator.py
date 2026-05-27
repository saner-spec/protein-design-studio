"""Protein Design Studio — 私有校准引擎

用户输入自己的实验数据（ΔΔG测量值），系统用这些数据校准ESM预测分数。
核心: 简单线性回归 + per-residue-type校准 + 持久化。

Phase 2 核心功能 — 越用越准。
"""
import json
import math
from pathlib import Path
from typing import Optional

from .config import DATA_DIR


class CalibrationDataset:
    """校准数据集"""

    def __init__(self, name: str):
        self.name = name
        self.entries: list[dict] = []  # [{position, wt, mut, exp_ddg, esm_score}]
        self.slope: Optional[float] = None
        self.intercept: Optional[float] = None
        self.r_squared: Optional[float] = None
        self.mae: Optional[float] = None
        self.n_samples: int = 0

    def add(self, position: int, wt: str, mut: str, exp_ddg: float, esm_score: float):
        self.entries.append({
            "position": position,
            "wt": wt,
            "mut": mut,
            "exp_ddg": exp_ddg,
            "esm_score": esm_score,
        })
        self.n_samples = len(self.entries)
        # 重置校准参数（需要重新拟合）
        self.slope = None
        self.intercept = None
        self.r_squared = None
        self.mae = None

    def fit(self) -> dict:
        """线性回归: exp_ddg = slope * esm_score + intercept"""
        if self.n_samples < 2:
            return {"error": "需要至少2个数据点", "n": self.n_samples}

        xs = [e["esm_score"] for e in self.entries]
        ys = [e["exp_ddg"] for e in self.entries]

        n = len(xs)
        mean_x = sum(xs) / n
        mean_y = sum(ys) / n

        # 最小二乘法
        ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
        ss_xx = sum((x - mean_x) ** 2 for x in xs)

        if abs(ss_xx) < 1e-10:
            self.slope = 0.0
            self.intercept = mean_y
        else:
            self.slope = ss_xy / ss_xx
            self.intercept = mean_y - self.slope * mean_x

        # R²
        ss_res = sum((y - (self.slope * x + self.intercept)) ** 2 for x, y in zip(xs, ys))
        ss_tot = sum((y - mean_y) ** 2 for y in ys)

        self.r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        # MAE
        self.mae = sum(abs(y - (self.slope * x + self.intercept)) for x, y in zip(xs, ys)) / n

        return self.summary()

    def predict(self, esm_score: float) -> dict:
        """对新的ESM分数进行校准预测"""
        if self.slope is None or self.intercept is None:
            return {"error": "模型未拟合，请先调用 fit()"}

        calibrated = self.slope * esm_score + self.intercept

        # 修正的预测区间 (prediction interval, 不是简单的置信区间)
        if self.n_samples >= 4:
            xs = [e["esm_score"] for e in self.entries]
            ys = [e["exp_ddg"] for e in self.entries]
            residuals = [y - (self.slope * x + self.intercept) for x, y in zip(xs, ys)]
            residual_std = math.sqrt(sum(r**2 for r in residuals) / (self.n_samples - 2))
            # t分布分位数 (90% CI, n-2自由度)
            # 小样本用保守近似: n<10时取3.0, n<30时取2.0, 否则1.645
            mean_x = sum(xs) / self.n_samples
            ss_xx = sum((x - mean_x)**2 for x in xs)
            if ss_xx > 0:
                se_pred = residual_std * math.sqrt(1 + 1.0/self.n_samples + (esm_score - mean_x)**2 / ss_xx)
            else:
                se_pred = residual_std
            t_val = 3.0 if self.n_samples < 10 else (2.0 if self.n_samples < 30 else 1.645)
            ci_half = t_val * se_pred
        else:
            ci_half = None

        return {
            "calibrated_ddg": round(calibrated, 3),
            "raw_esm_score": round(esm_score, 4),
            "ci_95": (
                [round(calibrated - ci_half, 3), round(calibrated + ci_half, 3)]
                if ci_half else None
            ),
            "r_squared": round(self.r_squared, 3) if self.r_squared else None,
            "n_samples": self.n_samples,
        }

    def summary(self) -> dict:
        return {
            "name": self.name,
            "n_samples": self.n_samples,
            "slope": round(self.slope, 4) if self.slope else None,
            "intercept": round(self.intercept, 4) if self.intercept else None,
            "r_squared": round(self.r_squared, 3) if self.r_squared else None,
            "mae": round(self.mae, 3) if self.mae else None,
            "reliability": (
                "good" if self.r_squared and self.r_squared > 0.6
                else "moderate" if self.r_squared and self.r_squared > 0.3
                else "low"
            ) if self.r_squared else "uncalibrated",
        }

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "entries": self.entries,
            "n_samples": self.n_samples,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "CalibrationDataset":
        ds = cls(data["name"])
        ds.entries = data.get("entries", [])
        ds.n_samples = len(ds.entries)
        if ds.n_samples >= 2:
            ds.fit()
        return ds


class CalibrationManager:
    """管理多个校准数据集"""

    def __init__(self, storage_dir: Path | None = None):
        self.storage_dir = storage_dir or (DATA_DIR / "calibrations")
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.datasets: dict[str, CalibrationDataset] = {}

    def _path(self, name: str) -> Path:
        # 防止路径遍历: 只允许字母数字下划线连字符
        import re
        if not re.match(r'^[a-zA-Z0-9_-]+$', name):
            raise ValueError(f"校准集名称只能包含字母数字下划线连字符: {name}")
        return self.storage_dir / f"{name}.json"

    def create(self, name: str) -> CalibrationDataset:
        ds = CalibrationDataset(name)
        self.datasets[name] = ds
        return ds

    def get(self, name: str) -> Optional[CalibrationDataset]:
        if name in self.datasets:
            return self.datasets[name]
        path = self._path(name)
        if path.exists():
            with open(path) as f:
                data = json.load(f)
            ds = CalibrationDataset.from_dict(data)
            self.datasets[name] = ds
            return ds
        return None

    def save(self, name: str) -> bool:
        ds = self.datasets.get(name)
        if not ds:
            return False
        with open(self._path(name), "w") as f:
            json.dump(ds.to_dict(), f, indent=2)
        return True

    def delete(self, name: str) -> bool:
        self.datasets.pop(name, None)
        path = self._path(name)
        if path.exists():
            path.unlink()
            return True
        return False

    def list_all(self) -> list[str]:
        names = set(self.datasets.keys())
        for p in self.storage_dir.glob("*.json"):
            names.add(p.stem)
        return sorted(names)


# 全局实例
calibration_mgr = CalibrationManager()
