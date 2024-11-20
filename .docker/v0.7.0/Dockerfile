FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="1200d5b46733c66c53ce0cf6ebfac79ea61f1163f663407b7842bc8d37b41620"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/afc705fc7d66595f67fba2816309d591
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.9
RUN mkdir -p /conda-envs/afc705fc7d66595f67fba2816309d591
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/cutadapt/se/environment.yaml /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/hisat2/index/environment.yaml
#   prefix: /conda-envs/b78cd6b273a8e5eaa856cc7165c924a1
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - hisat2 =2.2.1
RUN mkdir -p /conda-envs/b78cd6b273a8e5eaa856cc7165c924a1
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/hisat2/index/environment.yaml /conda-envs/b78cd6b273a8e5eaa856cc7165c924a1/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/863b76ea04987efcef2c54ec3473f2c3
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.25.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/863b76ea04987efcef2c54ec3473f2c3
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.1.0/bio/multiqc/environment.yaml /conda-envs/863b76ea04987efcef2c54ec3473f2c3/environment.yaml

# Conda environment:
#   source: workflow/envs/stats.yaml
#   prefix: /conda-envs/566b0c3f53a2f3378275090dc63dfacf
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
#     - hisat2=2.2.1
RUN mkdir -p /conda-envs/566b0c3f53a2f3378275090dc63dfacf
COPY workflow/envs/stats.yaml /conda-envs/566b0c3f53a2f3378275090dc63dfacf/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/afc705fc7d66595f67fba2816309d591 --file /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/b78cd6b273a8e5eaa856cc7165c924a1 --file /conda-envs/b78cd6b273a8e5eaa856cc7165c924a1/environment.yaml && \
    mamba env create --prefix /conda-envs/863b76ea04987efcef2c54ec3473f2c3 --file /conda-envs/863b76ea04987efcef2c54ec3473f2c3/environment.yaml && \
    mamba env create --prefix /conda-envs/566b0c3f53a2f3378275090dc63dfacf --file /conda-envs/566b0c3f53a2f3378275090dc63dfacf/environment.yaml && \
    mamba clean --all -y