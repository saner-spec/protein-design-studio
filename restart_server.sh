#!/bin/bash
pkill -9 -f "start_server\|run.py\|uvicorn" 2>/dev/null
sleep 2
export TMPDIR=/mnt/e/AI_Agents/tmp
export TORCH_HOME=/mnt/e/AI_Agents/torch_cache
# Set your API key via environment variable before running:
#   export DEEPSEEK_API_KEY="sk-..."
if [ -z "$DEEPSEEK_API_KEY" ]; then
    echo "ERROR: DEEPSEEK_API_KEY not set. Export it first."
    exit 1
fi
cd /mnt/e/AI_Agents/protein_designer
setsid .venv/bin/python run.py > server.log 2>&1 < /dev/null &
disown
echo "Starting..."
for i in $(seq 1 12); do
  sleep 2
  if curl -s http://localhost:8899/api/health > /dev/null 2>&1; then
    echo "UP"
    exit 0
  fi
done
echo "TIMEOUT"
tail -10 server.log
exit 1
