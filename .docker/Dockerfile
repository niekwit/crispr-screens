FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bf953c0e7bbfc033ff3c54e84eb51d7cc3a749a26e64263e42e2cde37f42e8f0"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.18.0/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =5.2
RUN mkdir -p /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.18.0/bio/cutadapt/se/environment.yaml /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.18.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/43340c98ee16cbaa141fbba9b9dd2847
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.35
#     - snakemake-wrapper-utils =0.9.0
RUN mkdir -p /conda-envs/43340c98ee16cbaa141fbba9b9dd2847
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.18.0/bio/multiqc/environment.yaml /conda-envs/43340c98ee16cbaa141fbba9b9dd2847/environment.yaml

# Conda environment:
#   source: workflow/envs/stats.yaml
#   prefix: /conda-envs/a3bd78ddd9e1f75c8117e1071683fcc2
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
#     - bioconda::r-crisprcleanr=3.0.0
#     - python=3.10.12
#     - mageck=0.5.9
#     - numpy<1.24
#     - pandas=2.0.3
#     - scipy=1.11.1
#     - scikit-learn=1.3.0
#     - seaborn=0.12.2
#     - click=8.1.7
#     - bowtie=1.3.1
#     - conda-forge::git=2.47.1
#     - samtools=1.21
#     - conda-forge::bash=5.2.37
#     - conda-forge::gzip=1.14
#     - conda-forge::requests
#     - conda-forge::plotly
#     - conda-forge::jq
RUN mkdir -p /conda-envs/a3bd78ddd9e1f75c8117e1071683fcc2
COPY workflow/envs/stats.yaml /conda-envs/a3bd78ddd9e1f75c8117e1071683fcc2/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32 --file /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml && \
    conda env create --prefix /conda-envs/43340c98ee16cbaa141fbba9b9dd2847 --file /conda-envs/43340c98ee16cbaa141fbba9b9dd2847/environment.yaml && \
    conda env create --prefix /conda-envs/a3bd78ddd9e1f75c8117e1071683fcc2 --file /conda-envs/a3bd78ddd9e1f75c8117e1071683fcc2/environment.yaml && \
    conda clean --all -y