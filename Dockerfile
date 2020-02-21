FROM rocker/shiny-verse:3.6.1

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('mutoss', 'metap'))"
RUN R -e "install.packages(c('Seurat', 'plotly', 'ggplot2'))"

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
