FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install python dependencies for scanpy and extensions
RUN pip install --no-cache-dir scanpy matplotlib seaborn pandas numpy leidenalg scrublet 'harmonypy<=0.0.9' celltypist scvelo multiqc
