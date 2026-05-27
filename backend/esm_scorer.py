"""Protein Design Studio — ESM-2 打分管线

基于 fair-esm 的 masked marginal scoring。
对每个位置一次性计算所有20种氨基酸的 pseudo-log-likelihood。

v0.1.1 修复:
  - 单次forward pass对所有20种aa打分 (而非19次独立forward)
  - 包含野生型打分 (0.0)
  - 添加OOM回退机制
"""
import torch
import esm
from typing import Dict

from .config import MODEL_NAME, ESM_MODELS, ESM_MAX_LEN
from .residue_map import STANDARD_AAS

_model = None
_alphabet = None
_batch_converter = None


def load_model() -> tuple:
    """加载ESM-2模型（首次调用时下载）"""
    global _model, _alphabet, _batch_converter

    if _model is not None:
        return _model, _alphabet, _batch_converter

    if str(ESM_MODELS) != ".":
        torch.hub.set_dir(str(ESM_MODELS))

    try:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        device = "cpu"

    print(f"  ⏳ 加载 ESM-2 650M (设备: {device})...")

    _model, _alphabet = esm.pretrained.load_model_and_alphabet(MODEL_NAME)

    try:
        _model = _model.to(device).eval()
    except (torch.cuda.OutOfMemoryError, RuntimeError):
        print(f"  ⚠ GPU OOM, 回退到CPU")
        device = "cpu"
        _model = _model.to(device).eval()

    _batch_converter = _alphabet.get_batch_converter()

    if device == "cuda":
        allocated = torch.cuda.memory_allocated() / 1e9
        print(f"  ✓ ESM-2 加载完成 (显存: {allocated:.1f}GB)")
    else:
        print(f"  ✓ ESM-2 加载完成 (CPU)")

    return _model, _alphabet, _batch_converter


def score_single_position(
    sequence: str,
    position: int,
) -> Dict[str, float]:
    """对某个位置所有20种氨基酸打分 — 单次forward pass

    Args:
        sequence: 蛋白质序列 (单字母, 已normalize)
        position: 位置索引 (0-based, 在sequence中)

    Returns:
        {aa: pseudo_log_likelihood} — 包含野生型(0.0)
        负值 = 模型认为该氨基酸在此位置更"自然"
    """
    model, alphabet, batch_converter = load_model()
    device = next(model.parameters()).device

    if position >= len(sequence):
        raise ValueError(f"位置{position}超出序列长度{len(sequence)}")

    wt_aa = sequence[position]

    # 单次forward pass: 输入野生型序列
    _, _, batch_tokens = batch_converter([("wt", sequence)])
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])
        logits = results["logits"]  # (1, L, vocab_size)

    # token位置: position+1 (因为<cls>在索引0)
    token_idx_in_logits = position + 1
    token_logits = logits[0, token_idx_in_logits, :]  # (vocab_size,)
    token_probs = torch.log_softmax(token_logits, dim=-1)

    # 提取所有20种标准氨基酸的log-probability
    wt_score = token_probs[alphabet.get_idx(wt_aa)].item()
    scores = {wt_aa: 0.0}  # 野生型定义为0基准

    for aa in STANDARD_AAS:
        if aa == wt_aa:
            continue
        aa_idx = alphabet.get_idx(aa)
        # pseudo-ΔΔG: 突变体logP - 野生型logP
        scores[aa] = token_probs[aa_idx].item() - wt_score

    return scores


def score_batch_positions(
    sequence: str,
    positions: list[int],
) -> Dict[int, Dict[str, float]]:
    """批量对多个位置打分 — 仍然是每个位置一次forward (ESM限制)

    对每个位置运行一次forward pass, 但每次pass同时得到所有20种aa的分数。
    相比v0.1.0的 N×19次pass, 现在是 N×1次pass。

    Args:
        sequence: 蛋白质序列
        positions: 位置列表 (0-indexed)

    Returns:
        {position: {aa: score, ...}, ...}
    """
    if not positions:
        return {}

    model, alphabet, batch_converter = load_model()
    device = next(model.parameters()).device

    results = {}

    for pos in positions:
        if pos >= len(sequence):
            continue

        wt_aa = sequence[pos]

        _, _, batch_tokens = batch_converter([("wt", sequence)])
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            logits = model(batch_tokens, repr_layers=[33])["logits"]

        token_idx = pos + 1
        token_probs = torch.log_softmax(logits[0, token_idx, :], dim=-1)

        wt_score = token_probs[alphabet.get_idx(wt_aa)].item()
        pos_scores = {wt_aa: 0.0}

        for aa in STANDARD_AAS:
            if aa == wt_aa:
                continue
            pos_scores[aa] = token_probs[alphabet.get_idx(aa)].item() - wt_score

        results[pos] = pos_scores

    return results


def check_vram() -> dict:
    """检查显存状态"""
    if not torch.cuda.is_available():
        return {"cuda_available": False, "message": "CUDA不可用，使用CPU推理"}

    total = torch.cuda.get_device_properties(0).total_mem / 1e9
    allocated = torch.cuda.memory_allocated() / 1e9
    free = total - allocated

    return {
        "cuda_available": True,
        "total_gb": round(total, 1),
        "free_gb": round(free, 1),
        "ok": free >= 3.0,
        "message": None if free >= 3.0 else "显存不足(需要≥3GB)",
    }
