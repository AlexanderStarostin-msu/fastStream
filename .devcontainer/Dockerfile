FROM mcr.microsoft.com/devcontainers/cpp:ubuntu-22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install core build tools, MPI, and essentials safely
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gdb \
    wget \
    curl \
    ca-certificates \
    libopenmpi-dev \
    openmpi-bin \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
