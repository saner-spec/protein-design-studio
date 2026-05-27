"""Protein Design Studio — 系统诊断器

在服务启动前运行，逐一检查环境就绪状态。
CC审核要求: 诊断器必须在任何实质性操作之前运行。
"""
import os
import sys
import json
from pathlib import Path


def full_diagnostics() -> dict:
    """运行全部诊断，返回结果字典"""
    results = {}

    # 1. Python 版本
    py_ver = sys.version_info
    results["python"] = {
        "version": f"{py_ver.major}.{py_ver.minor}.{py_ver.micro}",
        "ok": py_ver >= (3, 9),
    }

    # 2. PyTorch + CUDA
    try:
        import torch
        results["pytorch"] = {
            "version": torch.__version__,
            "ok": True,
        }
        results["cuda"] = {
            "available": torch.cuda.is_available(),
            "ok": torch.cuda.is_available(),
        }
        if torch.cuda.is_available():
            props = torch.cuda.get_device_properties(0)
            total_gb = props.total_memory / 1e9
            free_gb = (props.total_memory - torch.cuda.memory_allocated()) / 1e9
            results["gpu"] = {
                "name": props.name,
                "total_vram_gb": round(total_gb, 1),
                "free_vram_gb": round(free_gb, 1),
                "ok": free_gb >= 3.0,  # ESM-2 650M 需要 ~2.5GB
            }
        else:
            results["gpu"] = {"ok": False, "warning": "CUDA不可用，将使用CPU推理(较慢)"}
    except ImportError:
        results["pytorch"] = {"ok": False, "error": "PyTorch未安装"}
        results["cuda"] = {"ok": False}
        results["gpu"] = {"ok": False}

    # 3. Biopython
    try:
        import Bio
        results["biopython"] = {"version": Bio.__version__, "ok": True}
    except ImportError:
        results["biopython"] = {"ok": False, "error": "Biopython未安装"}

    # 4. 磁盘空间
    data_dir = Path(__file__).resolve().parent.parent / "data"
    try:
        import shutil
        check_path = data_dir if data_dir.exists() else data_dir.parent
        usage = shutil.disk_usage(check_path)
        free_gb = usage.free / 1e9
        results["disk"] = {
            "free_gb": round(free_gb, 1),
            "ok": free_gb >= 5.0,
            "path": str(check_path),
        }
    except Exception:
        results["disk"] = {"ok": True, "warning": "无法检测磁盘空间，跳过"}

    # 5. 端口
    import socket
    from .config import PORT
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        port_free = sock.connect_ex(("127.0.0.1", PORT)) != 0
    finally:
        sock.close()
    results["port"] = {
        "port": PORT,
        "free": port_free,
        "ok": port_free,
    }

    # 6. API Key
    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    results["api_key"] = {
        "configured": bool(api_key),
        "ok": True,  # API key缺失不算失败，只是AI功能降级
        "warning": None if api_key else "DeepSeek API Key未配置，AI解释功能不可用",
    }

    # 7. ESM模型 (快速检查路径存在，不全量加载)
    from .config import ESM_MODELS, MODEL_NAME
    esm_dir = ESM_MODELS
    if esm_dir.exists() and any(esm_dir.iterdir()):
        results["esm_model"] = {"ok": True, "path": str(esm_dir)}
    else:
        results["esm_model"] = {
            "ok": True,
            "warning": "ESM模型未下载，首次运行将自动下载(~2.5GB)",
            "path": str(esm_dir),
        }

    # 汇总
    critical_checks = ["python", "pytorch", "biopython", "disk", "port"]
    all_ok = all(results.get(c, {}).get("ok", False) for c in critical_checks)
    results["all_ok"] = all_ok

    return results


def format_diagnostics(results: dict) -> str:
    """格式化诊断结果为终端输出"""
    lines = []
    lines.append("")
    lines.append("  环境诊断")
    lines.append("  " + "─" * 20)

    def status(ok):
        return "✓" if ok else "✗"

    py = results.get("python", {})
    lines.append(f"  {status(py.get('ok'))} Python {py.get('version', '?')}")

    pt = results.get("pytorch", {})
    lines.append(f"  {status(pt.get('ok'))} PyTorch {pt.get('version', '未安装')}")

    cuda = results.get("cuda", {})
    if cuda.get("available"):
        gpu = results.get("gpu", {})
        lines.append(f"  {status(gpu.get('ok'))} GPU: {gpu.get('name', '?')} ({gpu.get('free_vram_gb', 0)}GB 可用)")
    else:
        lines.append(f"  {status(False)} GPU 不可用 — 将使用CPU推理")

    bio = results.get("biopython", {})
    lines.append(f"  {status(bio.get('ok'))} Biopython {bio.get('version', '未安装')}")

    disk = results.get("disk", {})
    lines.append(f"  {status(disk.get('ok'))} 磁盘空闲 {disk.get('free_gb', '?')}GB")

    port = results.get("port", {})
    lines.append(f"  {status(port.get('ok'))} 端口 {port.get('port', '?')} 空闲")

    esm = results.get("esm_model", {})
    if esm.get("warning"):
        lines.append(f"  ⚠ {esm['warning']}")

    api = results.get("api_key", {})
    if api.get("warning"):
        lines.append(f"  ⚠ {api['warning']}")

    # Warnings
    for key, result in results.items():
        if isinstance(result, dict) and result.get("warning") and not result.get("ok", True):
            lines.append(f"  ⚠ [{key}] {result['warning']}")

    lines.append("")
    if results.get("all_ok"):
        lines.append("  ✅ 所有核心检查通过")
    else:
        lines.append("  ❌ 部分检查未通过，请查看上方详情")

    return "\n".join(lines)
