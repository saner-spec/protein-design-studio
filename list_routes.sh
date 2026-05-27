#!/bin/bash
cd /mnt/e/AI_Agents/protein_designer
.venv/bin/python << 'PYEOF'
from backend.main import app
for r in sorted(app.routes, key=lambda r: r.path):
    if "knowledge" in r.path.lower():
        print(r.methods, r.path)
PYEOF
