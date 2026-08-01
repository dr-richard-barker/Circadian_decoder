FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    r-base \
    r-cran-metafor \
    r-cran-ggplot2 \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-fonts-recommended \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY . .

# Create output directories
RUN mkdir -p /results/figures /results/tables /results/data

# Set environment variables
ENV TF_CPP_MIN_LOG_LEVEL=3
ENV PYTHONPATH=/app/src

# Default command
CMD ["python", "src/run_analysis.py"]
