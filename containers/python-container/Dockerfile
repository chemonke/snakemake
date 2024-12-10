# Use a base image with Python 3.12
FROM python:3.12-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH="/opt/conda/bin:$PATH"

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && bash ~/miniconda.sh -b -p /opt/conda \
    && rm ~/miniconda.sh \
    && /opt/conda/bin/conda clean -afy

# Create a working directory
WORKDIR /workspace

# Copy the Conda environment file
COPY requirements.yaml .

# Create the Conda environment
RUN conda env create -f requirements.yaml \
    && conda clean -afy

# Set the Conda environment as the default
RUN echo "conda activate base" >> ~/.bashrc

# Default command
CMD ["bash"]
