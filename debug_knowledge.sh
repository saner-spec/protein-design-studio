#!/bin/bash
cd /mnt/e/AI_Agents/protein_designer
# Set your API key via environment variable before running:
#   export DEEPSEEK_API_KEY="sk-..."
if [ -z "$DEEPSEEK_API_KEY" ]; then
    echo "ERROR: DEEPSEEK_API_KEY not set. Export it first."
    exit 1
fi
.venv/bin/python << 'PYEOF'
import asyncio
from backend.ai_client import chat_stream
from backend.knowledge_engine import build_pdb_info_for_prompt, KNOWLEDGE_PROMPT

# Simulate 1ema structure
from backend.pdb_parser import parse_pdb
import asyncio

async def test():
    # Load structure
    pdb_path = "/mnt/e/AI_Agents/protein_designer/data/pdb_cache/1ema.pdb"
    structure = parse_pdb(pdb_path)
    pdb_info = build_pdb_info_for_prompt(structure)
    prompt = KNOWLEDGE_PROMPT.format(pdb_info=pdb_info)

    print("=== Prompt length:", len(prompt), "chars ===")
    print(prompt[:200])
    print("...")
    print("=== Streaming response ===")

    full = ""
    try:
        async for token in chat_stream(prompt, context="", model="deepseek-v4-pro"):
            full += token
            print(token, end="", flush=True)
    except Exception as e:
        print(f"\nEXCEPTION: {e}")

    print(f"\n=== Total: {len(full)} chars ===")

asyncio.run(test())
PYEOF
