FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="0aae8a0ca0e8d368d2006d43ddda13c5acc70fd76d77ddeeacff6704a1d384d5"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.9
RUN mkdir -p /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/cutadapt/se/environment.yaml /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7/environment.yaml

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
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/hisat2/index/environment.yaml
#   prefix: /conda-envs/6bde2d6fbc33e1c614cc5b12e3191198
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - hisat2 =2.2.1
RUN mkdir -p /conda-envs/6bde2d6fbc33e1c614cc5b12e3191198
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/hisat2/index/environment.yaml /conda-envs/6bde2d6fbc33e1c614cc5b12e3191198/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/multiqc/environment.yaml
#   prefix: /conda-envs/9ab3161e947cc968897b9bb4302e9e45
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.25.2
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/9ab3161e947cc968897b9bb4302e9e45
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/multiqc/environment.yaml /conda-envs/9ab3161e947cc968897b9bb4302e9e45/environment.yaml

# Conda environment:
#   source: workflow/envs/stats.yaml
#   prefix: /conda-envs/7ec7a87a2587641ac871b018ae1e51df
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
#     - conda-forge::r-gprofiler2=0.2.3
#     - python=3.10.12
#     - mageck=0.5.9
#     - numpy<1.24
#     - pandas=2.0.3
#     - scipy=1.11.1
#     - scikit-learn=1.3.0
#     - seaborn=0.12.2
#     - click=8.1.7
#     - hisat2=2.2.1
RUN mkdir -p /conda-envs/7ec7a87a2587641ac871b018ae1e51df
COPY workflow/envs/stats.yaml /conda-envs/7ec7a87a2587641ac871b018ae1e51df/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7 --file /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7/environment.yaml && \
    conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/6bde2d6fbc33e1c614cc5b12e3191198 --file /conda-envs/6bde2d6fbc33e1c614cc5b12e3191198/environment.yaml && \
    conda env create --prefix /conda-envs/9ab3161e947cc968897b9bb4302e9e45 --file /conda-envs/9ab3161e947cc968897b9bb4302e9e45/environment.yaml && \
    conda env create --prefix /conda-envs/7ec7a87a2587641ac871b018ae1e51df --file /conda-envs/7ec7a87a2587641ac871b018ae1e51df/environment.yaml && \
    conda clean --all -y