# ==============================================================================
# Dockerfile for AD-GBM re-evaluation analysis
# Base: rocker/tidyverse with Bioconductor
# ==============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL maintainer="Ahmet Karakus <akarakus@bartin.edu.tr>"
LABEL project="AD-GBM cell-composition-controlled re-evaluation"

# System libraries needed by some Bioconductor packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhdf5-dev \
    libglpk-dev \
    libgmp-dev \
    libbz2-dev \
    liblzma-dev \
    libpng-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages
RUN R -e "install.packages(c( \
    'here', 'yaml', 'dplyr', 'tidyr', 'readr', 'purrr', 'stringr', 'tibble', 'glue', \
    'R.utils', 'data.table', 'car', 'ggplot2', 'patchwork', 'RColorBrewer', \
    'Seurat', 'renv' \
    ), repos = 'https://cloud.r-project.org')"

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'GEOquery', 'affy', 'oligo', 'hgu133plus2.db', 'hugene10sttranscriptcluster.db', \
    'sva', 'limma', 'edgeR', 'SingleCellExperiment', 'scran', 'scater', \
    'scDblFinder', 'WGCNA', 'fgsea', 'msigdbr', 'ComplexHeatmap' \
    ), update = FALSE, ask = FALSE)"

# GitHub-only packages
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"
RUN R -e "remotes::install_github('mcarona/BRETIGEA')"
RUN R -e "remotes::install_github('randel/MIND')"
RUN R -e "remotes::install_github('RRHO2/RRHO2')"
RUN R -e "remotes::install_github('saezlab/decoupleR')"

# Working directory
WORKDIR /workspace

# Default command: R shell
CMD ["R"]
