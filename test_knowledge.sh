#!/bin/bash
echo "=== Brief ==="
curl -s "http://localhost:8899/api/knowledge/1ema/brief" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('brief','ERROR'))"

echo ""
echo "=== Knowledge (first 500 chars) ==="
curl -s "http://localhost:8899/api/knowledge/1ema" 2>&1 | head -c 500
echo ""
echo "..."
