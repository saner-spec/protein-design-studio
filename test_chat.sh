#!/bin/bash
# Test chat directly to see if DeepSeek responds
curl -s -X POST "http://localhost:8899/api/chat" \
  -H "Content-Type: application/json" \
  -d '{"message":"用一句话介绍绿色荧光蛋白GFP","pdb_id":"1ema"}' 2>&1 | head -c 500
