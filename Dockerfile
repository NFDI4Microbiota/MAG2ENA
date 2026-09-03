FROM python:3.13-slim

# Install system dependencies: Redis for job queue, Java for webin-cli jar
RUN apt-get update && apt-get install -y --no-install-recommends \
    redis-server \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Pull uv binary from the official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Create non-root user
RUN groupadd --gid 1001 appuser && \
    useradd --uid 1001 --gid 1001 --create-home appuser

WORKDIR /app

# Copy lockfile first so dependency install is cached independently of source changes
COPY pyproject.toml uv.lock ./

# Install into the image's system Python — scripts land in /usr/local/bin with stable shebangs
RUN uv pip install --system --no-cache-dir \
    "cryptography>=46.0.3" \
    "ipykernel>=7.1.0" \
    "pandas>=2.3.2" \
    "polars>=1.33.1" \
    "pydeck>=0.9.1" \
    "redis>=7.4.0" \
    "requests>=2.32.5" \
    "rq>=2.9.0" \
    "streamlit>=1.49.1" \
    "streamlit-option-menu>=0.4.0" \
    "tqdm>=4.67.1"

# Copy project source and transfer ownership to non-root user
COPY . .

# Directories written to at runtime (sqlite db, submitted jobs) are made group-writable
# and group-owned by root: some platforms run the container under a UID other than
# appuser's, but that UID is still placed in group 0, so this keeps writes working either way.
RUN chmod +x start.sh && \
    mkdir -p /app/App/task-results /app/App/jobs && \
    chown -R appuser:appuser /app && \
    chgrp -R 0 /app/App/task-results /app/App/jobs && \
    chmod -R g=u /app/App/task-results /app/App/jobs

USER appuser

EXPOSE 8501

CMD ["/app/start.sh"]
