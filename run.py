#!/usr/bin/env python3
"""Protein Design Studio — 一键启动

用法:
    python run.py            # 默认端口 8899
    python run.py --port 9000
    python run.py --offline  # 跳过ESM模型下载检查

环境变量:
    DEEPSEEK_API_KEY     DeepSeek API密钥 (OpenRouter)
    PORT                 服务端口 (默认 8899)
    ESM_MODEL_DIR        ESM模型缓存目录 (默认 ./data/esm_models)
"""
import sys
import os
from pathlib import Path

# 确保项目根在Python路径中
sys.path.insert(0, str(Path(__file__).resolve().parent))

from backend.diagnostics import full_diagnostics, format_diagnostics
from backend.config import PORT as DEFAULT_PORT


def main():
    offline = "--offline" in sys.argv

    print()
    print("  ╔══════════════════════════════════════════╗")
    print("  ║     Protein Design Studio v0.4.0         ║")
    print("  ╚══════════════════════════════════════════╝")
    print()

    # 环境变量覆盖配置
    port = int(os.environ.get("PORT", DEFAULT_PORT))
    for arg in sys.argv:
        if arg.startswith("--port="):
            port = int(arg.split("=")[1])

    # 诊断
    diag = full_diagnostics()
    print(format_diagnostics(diag))

    if not diag.get("all_ok"):
        print("\n  ❌ 核心依赖未就绪，无法启动。请修复上述问题后重试。")
        print(f"     提示: pip install -r requirements.txt")
        sys.exit(1)

    # ESM模型检查
    if not offline and diag.get("esm_model", {}).get("warning"):
        print("\n  ⏳ ESM模型未下载，首次运行将自动下载 (~2.5GB)...")
        print("     如需跳过: python run.py --offline")
        print()

    # 启动服务
    import os as _os
    if offline:
        _os.environ["ESM_OFFLINE"] = "1"
    
    from backend.main import app
    import uvicorn

    print(f"  ─────────────────────────────────────────")
    print(f"  服务启动: http://localhost:{port}")
    print(f"  按 Ctrl+C 停止")
    print(f"  ─────────────────────────────────────────")
    print()

    uvicorn.run(app, host="0.0.0.0", port=port, log_level="warning")


if __name__ == "__main__":
    main()
