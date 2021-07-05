FROM rocker/shiny-verse:4.1.0

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('mutoss', 'metap'))"
RUN R -e "install.packages(c('Seurat', 'plotly', 'ggplot2'))"
