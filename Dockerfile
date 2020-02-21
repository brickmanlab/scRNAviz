FROM rocker/shiny-verse:3.6.0

RUN R -e "install.packages(c('Seurat', 'plotly', 'ggplot2'))"

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
