"""Protein Design Studio — DeepSeek API 客户端

通过 OpenRouter 调用 DeepSeek API。
仅单Provider，失败时返回降级消息，不崩溃。
"""
import json
import httpx
from typing import AsyncIterator

from .config import DS_API_KEY, DS_API_BASE, DS_MODEL

SYSTEM_PROMPT = """你是一位蛋白质设计导师。你的学生正在使用一个3D蛋白质结构查看器学习蛋白质设计。

你的角色:
1. 解释 — 用通俗语言解释蛋白质结构和突变效应
2. 引导 — 用苏格拉底式提问帮助学生自己思考
3. 诚实 — 清晰说明预测的不确定性，不假装精确

规则:
- 中文为主，专业术语保留英文
- 回答简洁 (200-500字)
- 如果学生选中了一个残基，先解释它的角色，再提问引导
- 引用具体数据时标注来源和局限性
- 不要给出危险建议 (如让非专业人士操作危险化学品)

当前上下文会随每次请求注入，包括蛋白信息和选中残基的详细数据。
"""


def build_context(
    pdb_id: str = "",
    residue_info: dict | None = None,
    scores: dict | None = None,
) -> str:
    """构建注入提示词的上下文"""
    parts = []
    if pdb_id:
        parts.append(f"当前蛋白质: PDB {pdb_id}")
    if residue_info:
        r = residue_info
        parts.append(
            f"学生选中了残基: {r.get('resname','')}{r.get('resnum','')} "
            f"(链{r.get('chain','')}, {r.get('aa','')}), "
            f"SASA: {r.get('sasa','?')}Å², "
            f"二级结构: {r.get('secondary_structure','?')}, "
            f"表面暴露: {'是' if r.get('surface') else '否'}"
        )
    if scores:
        top3 = sorted(scores.items(), key=lambda x: x[1])[:3]
        score_str = ", ".join(f"{aa}: {s:.2f}" for aa, s in top3)
        parts.append(f"ESM打分 (负值=可能稳定): {score_str}")
    return "\n".join(parts)


async def chat_stream(
    message: str,
    history: list[dict] | None = None,
    context: str = "",
    model: str = DS_MODEL,
) -> AsyncIterator[str]:
    """流式对话，逐token返回"""
    if not DS_API_KEY:
        yield "⚠️ DeepSeek API Key 未配置。AI解释功能不可用。"
        return

    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    if context:
        messages.append({"role": "system", "content": f"当前上下文:\n{context}"})
    if history:
        messages.extend(history)
    messages.append({"role": "user", "content": message})

    async with httpx.AsyncClient(timeout=60.0) as client:
        try:
            async with client.stream(
                "POST",
                f"{DS_API_BASE}/chat/completions",
                headers={
                    "Authorization": f"Bearer {DS_API_KEY}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": model,
                    "messages": messages,
                    "stream": True,
                    "temperature": 0.7,
                    "max_tokens": 1024,
                },
            ) as resp:
                if resp.status_code != 200:
                    yield f"⚠️ API 错误 ({resp.status_code})"
                    return

                async for line in resp.aiter_lines():
                    if line.startswith("data: "):
                        data = line[6:]
                        if data == "[DONE]":
                            break
                        try:
                            chunk = json.loads(data)
                            delta = chunk["choices"][0].get("delta", {})
                            content = delta.get("content", "")
                            if content:
                                yield content
                        except (json.JSONDecodeError, KeyError, IndexError):
                            continue

        except httpx.TimeoutException:
            yield "\n⚠️ API 超时，请稍后重试。"
        except Exception as e:
            yield f"\n⚠️ 请求失败: {str(e)}"


async def explain_residue(residue_info: dict, pdb_id: str = "") -> str:
    """生成残基解释（非流式，用于点击后自动显示）"""
    if not DS_API_KEY:
        return ""

    context = build_context(pdb_id=pdb_id, residue_info=residue_info)
    prompt = f"请简要解释这个残基在蛋白质中的角色（2-3句话），然后提一个引导性问题让学生思考。"

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "system", "content": f"当前上下文:\n{context}"},
        {"role": "user", "content": prompt},
    ]

    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            resp = await client.post(
                f"{DS_API_BASE}/chat/completions",
                headers={
                    "Authorization": f"Bearer {DS_API_KEY}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": DS_MODEL,
                    "messages": messages,
                    "temperature": 0.7,
                    "max_tokens": 512,
                },
            )
            if resp.status_code == 200:
                data = resp.json()
                return data["choices"][0]["message"]["content"]
        except (httpx.TimeoutException, httpx.ConnectError, httpx.HTTPStatusError):
            pass

    return ""
