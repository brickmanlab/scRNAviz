FROM rocker/shiny:4.1.0

RUN apt update && apt install -y \
    --no-install-recommends \
    libglpk40 \
    libxml2-dev && \
    apt-get clean && \ 
    rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
    BiocManager \
    Seurat \
    plotly \
    ggplot2
