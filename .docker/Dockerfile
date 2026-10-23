FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="870b60c16c6a118786ea136eccb3ba70e01856bae5177a407816e2cdf95d630e"

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
#   prefix: /conda-envs/2e2055326281268955f2572ee865118b
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
#     - conda-forge::python-kaleido
#     - conda-forge::jq
RUN mkdir -p /conda-envs/2e2055326281268955f2572ee865118b
COPY workflow/envs/stats.yaml /conda-envs/2e2055326281268955f2572ee865118b/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32 --file /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml && \
    conda env create --prefix /conda-envs/43340c98ee16cbaa141fbba9b9dd2847 --file /conda-envs/43340c98ee16cbaa141fbba9b9dd2847/environment.yaml && \
    conda env create --prefix /conda-envs/2e2055326281268955f2572ee865118b --file /conda-envs/2e2055326281268955f2572ee865118b/environment.yaml && \
    conda clean --all -y