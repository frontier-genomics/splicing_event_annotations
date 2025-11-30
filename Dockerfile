FROM python:3.13-slim

# 1. Runtime env tweaks
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONPATH=/app:$PYTHONPATH

WORKDIR /app

# 2. Install system packages (ps from procps)
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# 3. Install Python deps with good caching
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copy application code
COPY . .

# 5. Non-root user
RUN adduser --disabled-password --gecos '' appuser && chown -R appuser /app
USER appuser

# 6. Entry point (works in Docker & Singularity)
ENTRYPOINT ["python", "-m", "splicing_event_annotator.main"]
CMD []