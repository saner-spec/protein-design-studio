#!/usr/bin/env python3
"""ESM-2 650M 模型预下载脚本

在中国网络环境下，首次自动下载可能失败或极慢。
此脚本提前下载模型权重到本地，避免启动时等待。

用法:
    python download_esm.py              # 下载到默认路径 (./data/esm_models/)
    python download_esm.py --dir /path  # 下载到指定路径

下载量: ~2.5GB
"""
import sys
import os
from pathlib import Path

# 设置模型缓存目录
script_dir = Path(__file__).resolve().parent
default_dir = script_dir / "data" / "esm_models"

target_dir = default_dir
for arg in sys.argv:
    if arg.startswith("--dir="):
        target_dir = Path(arg.split("=", 1)[1])
    elif arg == "--dir" and len(sys.argv) > sys.argv.index(arg) + 1:
        target_dir = Path(sys.argv[sys.argv.index(arg) + 1])

target_dir.mkdir(parents=True, exist_ok=True)
os.environ["TORCH_HOME"] = str(target_dir)

print(f"ESM模型下载目录: {target_dir}")
print("开始下载 ESM-2 650M (~2.5GB)...")
print("首次下载可能需要5-30分钟，取决于网络速度。")
print()

import torch
import esm

try:
    model, alphabet = esm.pretrained.load_model_and_alphabet("esm2_t33_650M_UR50D")
    print(f"\n✅ 下载完成！")
    print(f"   模型路径: {target_dir}")
    print(f"   启动时设置环境变量: export ESM_MODEL_DIR={target_dir}")
except Exception as e:
    print(f"\n❌ 下载失败: {e}")
    print("  尝试手动下载:")
    print("  1. 访问 https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt")
    print(f"  2. 下载后放入 {target_dir}/checkpoints/")
    sys.exit(1)
