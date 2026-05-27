FROM python:3.12-slim

LABEL description="Protein Design Studio — AI-guided protein mutation design tool"

WORKDIR /app

# 系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Python 依赖
RUN pip install --no-cache-dir \
    torch==2.5.1+cu121 --index-url https://download.pytorch.org/whl/cu121 \
    && pip install --no-cache-dir \
    fair-esm \
    fastapi \
    "uvicorn[standard]" \
    biopython \
    httpx \
    pydantic

# 复制代码
COPY backend/ ./backend/
COPY frontend/ ./frontend/
COPY run.py download_esm.py requirements.txt ./

# 数据目录
RUN mkdir -p /app/data/pdb_cache /app/data/esm_models /app/data/calibrations /app/data/knowledge

# 环境变量
ENV ESM_MODEL_DIR=/app/data/esm_models
ENV PIP_CACHE_DIR=/app/data/pip_cache
ENV TORCH_HOME=/app/data/torch_cache

EXPOSE 8899

CMD ["python", "run.py"]
