#!/bin/bash
echo "=== Test 1: No cache (should generate) ==="
curl -s "http://localhost:8899/api/knowledge/1ema" | python3 -c "
import sys,json
d=json.load(sys.stdin)
print('cached:', d.get('cached'), '| chars:', len(d.get('content','')))
"

echo "=== Test 2: Should hit cache ==="
curl -s "http://localhost:8899/api/knowledge/1ema" | python3 -c "
import sys,json
d=json.load(sys.stdin)
print('cached:', d.get('cached'), '| chars:', len(d.get('content','')))
"

echo "=== Test 3: Force refresh ==="
curl -s "http://localhost:8899/api/knowledge/1ema?refresh=true" | python3 -c "
import sys,json
d=json.load(sys.stdin)
print('cached:', d.get('cached'), '| chars:', len(d.get('content','')))
"
