"""Protein Design Studio — 配置"""
import os
from pathlib import Path

# 路径
BASE_DIR = Path(__file__).resolve().parent.parent  # protein-studio/
DATA_DIR = BASE_DIR / "data"
PDB_CACHE = DATA_DIR / "pdb_cache"
ESM_MODELS = Path(os.environ.get("ESM_MODEL_DIR", str(DATA_DIR / "esm_models")))
MODEL_NAME = "esm2_t33_650M_UR50D"

# 服务
HOST = "0.0.0.0"
PORT = int(os.environ.get("PORT", 8899))

# GPU
GPU_REQUIRED_VRAM_GB = 3.0  # ESM-2 650M 安全阈值
ESM_MAX_LEN = 1024  # ESM tokenizer 硬上限
ESM_FULL_SCAN_MAX_LEN = 500  # 全残基批量扫描上限

# DeepSeek API (直连)
DS_API_KEY = os.environ.get("DEEPSEEK_API_KEY", "")
DS_API_BASE = "https://api.deepseek.com/v1"
DS_MODEL = "deepseek-v4-pro"

# 结构分析
DEFAULT_SASA_THRESHOLD = 0.25  # relSASA > 25% = 表面残基
HBOND_DIST_CUTOFF = 3.5  # 氢键距离阈值 (Å)
HBOND_ANGLE_CUTOFF = 120  # 氢键角度阈值 (°) — 可选，Phase 1暂用距离

# 缓存
SCORE_CACHE_TTL = 3600  # 打分结果缓存1小时
