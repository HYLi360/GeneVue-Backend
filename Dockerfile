# Pull a docker image.
ARG PYTHON_VERSION=3.9
FROM python:${PYTHON_VERSION}-slim

# Use it if you want running system-level command.
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    make \
    gffread \
    muscle \
    mafft \
    # and more \
    && rm -rf /var/lib/apt/lists/*

# make install here. Such as...
# RUN cd /tmp && curl -O http://example.com/soft.tar.gz && tar ... && make install

# Using UV
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Setup workdir.
WORKDIR /app

# Copy all file(s)
COPY . .
# Copy uv.lock or requirements.txt if you have
# COPY uv.lock .

# Install dependences.
RUN uv pip install --system --upgrade pip
RUN uv pip install --system -r pyproject.toml
RUN uv pip install --system -e .

# These are only used in CI, but not included in pyproject.toml
RUN uv pip install --system ruff pytest pytest-cov pytest-mpl

# Command as placeholder
CMD ["pytest"]
