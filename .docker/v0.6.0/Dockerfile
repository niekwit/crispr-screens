FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="d247ac7675afe5181718b120ea4371d8b8a5bf4ee1899bf205fdfe900c080302"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/a414a21b822baab7aecd9c5dd350be27
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.8
RUN mkdir -p /conda-envs/a414a21b822baab7aecd9c5dd350be27
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/cutadapt/se/environment.yaml /conda-envs/a414a21b822baab7aecd9c5dd350be27/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/hisat2/index/environment.yaml
#   prefix: /conda-envs/fb35c3e95b95b670b1aebc28877f773d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - hisat2 =2.2.1
#     - samtools =1.20
RUN mkdir -p /conda-envs/fb35c3e95b95b670b1aebc28877f773d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/hisat2/index/environment.yaml /conda-envs/fb35c3e95b95b670b1aebc28877f773d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/multiqc/environment.yaml
#   prefix: /conda-envs/f7d2de9338600e9135cf71b3d6ea89b1
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.21
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/f7d2de9338600e9135cf71b3d6ea89b1
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.10.2/bio/multiqc/environment.yaml /conda-envs/f7d2de9338600e9135cf71b3d6ea89b1/environment.yaml

# Conda environment:
#   source: workflow/envs/count.yaml
#   prefix: /conda-envs/37d3d3b0fcad4571e972404fb3bf1a98
#   name: count
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - hisat2=2.2.1
#     - pandas=2.1.1
RUN mkdir -p /conda-envs/37d3d3b0fcad4571e972404fb3bf1a98
COPY workflow/envs/count.yaml /conda-envs/37d3d3b0fcad4571e972404fb3bf1a98/environment.yaml

# Conda environment:
#   source: workflow/envs/stats.yaml
#   prefix: /conda-envs/db578f4a775b2f82d2ccf2dc03ebd055
#   name: stats
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.1
#     - r-viridis=0.6.4
#     - r-ggrepel=0.9.4
#     - r-ineq=0.2_13
#     - r-reshape2=1.4.4
#     - r-enrichr=3.2
#     - python=3.10.12
#     - mageck=0.5.9
#     - numpy<1.24
#     - pandas=2.0.3
#     - scipy=1.11.1
#     - scikit-learn=1.3.0
#     - seaborn=0.12.2
#     - click=8.1.7
RUN mkdir -p /conda-envs/db578f4a775b2f82d2ccf2dc03ebd055
COPY workflow/envs/stats.yaml /conda-envs/db578f4a775b2f82d2ccf2dc03ebd055/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/a414a21b822baab7aecd9c5dd350be27 --file /conda-envs/a414a21b822baab7aecd9c5dd350be27/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/fb35c3e95b95b670b1aebc28877f773d --file /conda-envs/fb35c3e95b95b670b1aebc28877f773d/environment.yaml && \
    mamba env create --prefix /conda-envs/f7d2de9338600e9135cf71b3d6ea89b1 --file /conda-envs/f7d2de9338600e9135cf71b3d6ea89b1/environment.yaml && \
    mamba env create --prefix /conda-envs/37d3d3b0fcad4571e972404fb3bf1a98 --file /conda-envs/37d3d3b0fcad4571e972404fb3bf1a98/environment.yaml && \
    mamba env create --prefix /conda-envs/db578f4a775b2f82d2ccf2dc03ebd055 --file /conda-envs/db578f4a775b2f82d2ccf2dc03ebd055/environment.yaml && \
    mamba clean --all -y
