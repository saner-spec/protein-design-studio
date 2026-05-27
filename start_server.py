"""Quick start script — skip diagnostics, go straight to uvicorn"""
import os
import sys

# API key from environment variable (required)
if not os.environ.get("DEEPSEEK_API_KEY"):
    print("ERROR: DEEPSEEK_API_KEY environment variable not set.")
    print("  export DEEPSEEK_API_KEY='sk-...'")
    sys.exit(1)

from backend.main import app
import uvicorn

print("Starting Protein Design Studio on :8899 ...")
uvicorn.run(app, host="0.0.0.0", port=8899, log_level="info")
