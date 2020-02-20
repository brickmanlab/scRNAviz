library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

init_empty <- list(
  placeholder = 'Please select an option below',
  onInitialize = I('function() { this.setValue(""); }'))

rna_datasets <- list(messmer_et_al_2019 = "./data/messmer_et_al_2019.rds",
                     deng_et_al_2014 = "./data/deng_et_al_2014.rds",
                     mohammed_et_al_2017 = "./data/mohammed_et_al_2017.rds",
                     nakamura_et_al_2017 = "./data/nakamura_et_al_2017.rds",
                     nowotschin_et_al_2019 = "./data/nowotschin_et_al_2019.rds",
                     posfai_et_al_2017 = "./data/posfai_et_al_2017.rds",
                     stirparo_et_al_2018 = "./data/stirparo_et_al_2018.rds",
                     sara = "./data/sara-processed.rds"
                     )

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Single cell RNA visualizer"),
  
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for choosing dataset ----
      selectizeInput(inputId = "dataset",
                     label = "Dataset:",
                     options = init_empty,
                     multiple = FALSE,
                     choices = c("Messmer et al. (2019)",
                                 "Deng et al. (2014)",
                                 "Mohammned et al. (2017)",
                                 "Nakamura et al. (2017)",
                                 "Nowotschin et al. (2019)",
                                 "Posfai et al. (2017)",
                                 "Stirparo et al. (2018)",
                                 "Sara")),
      
      hr(),
      conditionalPanel(
        condition = "input.dataset != ''",
        # Input: Visualization method ----
        
        
        fluidRow(
          column(5,
                 selectizeInput(inputId = "viz",
                                options = init_empty,
                                label = "Choose a visualization method:",
                                choices = c("PCA","t-SNE", "UMAP"))),
          column(7,
                 selectizeInput(inputId = "metadata1",
                                label = "Choose metadata 1:",
                                options = NULL,
                                choices = NULL),
                 
                 selectizeInput(inputId = "metadata2",
                                label = "Choose metadata 2:",
                                options = NULL,
                                choices = NULL)
                 )
        ),
        
        hr(),
        # Input: Gene ----
        # Note: Changes made to the caption in the textInput control
        # are updated in the output area immediately as you type
        # textInput(inputId = "gene",
        #          label = "Check for gene expression:",
        #          value = "DNMT3L")
        # 
        fluidRow(
          column(8,
                 selectizeInput(inputId = "gene",
                                label = "Check for gene expression:",
                                choices = NULL,
                                options = NULL,
                                multiple = FALSE)),
          
          column(4,
                 checkboxInput("sort","Sort cells"))
        ))
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Verbatim text for data summary ----
      #verbatimTextOutput("summary"),
      
      # Output: HTML table with requested number of observations ----
      #tableOutput("view")
      
      # Output: Histogram ----
      fluidRow(
        column(6,
               plotlyOutput(outputId = "metadata1Plot")),
        
        column(6,
               plotlyOutput(outputId = "metadata2Plot"))
        
      ),
      
      hr(),
      
      fluidRow(
        column(3),
        column(6,
               plotlyOutput(outputId = "genePlot")),
        column(3)
      ),
      
      hr(),
      
      fluidRow(
        column(6,
               plotlyOutput(outputId = "metadata1VlnPlot")),
        
        column(6,
               plotlyOutput(outputId = "metadata2VlnPlot"))
        
      ),
      
      
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  
  # Return the requested dataset ----
  # By declaring datasetInput as a reactive expression we ensure
  # that:
  #
  # 1. It is only called when the inputs it depends on changes
  # 2. The computation and result are shared by all the callers,
  #    i.e. it only executes a single time
  datasetInput <- reactive({
    switch(input$dataset,
           "Messmer et al. (2019)" = readRDS(rna_datasets$messmer_et_al_2019),
           "Deng et al. (2014)" = readRDS(rna_datasets$deng_et_al_2014),
           "Mohammned et al. (2017)" = readRDS(rna_datasets$mohammed_et_al_2017),
           "Nakamura et al. (2017)" = readRDS(rna_datasets$nakamura_et_al_2017),
           "Nowotschin et al. (2019)" = readRDS(rna_datasets$nowotschin_et_al_2019),
           "Posfai et al. (2017)" = readRDS(rna_datasets$posfai_et_al_2017),
           "Stirparo et al. (2018)" = readRDS(rna_datasets$stirparo_et_al_2018),
           "Sara" = readRDS(rna_datasets$sara)
    )
  })
  
  vizInput <- reactive({
    switch(input$viz,
           "PCA" = "pca",
           "t-SNE" = "tsne",
           "UMAP" = "umap")
  })
  
  
  observeEvent(input$dataset, {
    req(input$dataset)
    updateSelectizeInput(session, 'gene',
                         choices = c(sort(rownames(datasetInput()))), 
                         server = TRUE)
    
    updateSelectizeInput(session, 'metadata1',
                         choices = c(sort(colnames(datasetInput()@meta.data))),
                         selected = "seurat_clusters",
                         server = TRUE)
    
    updateSelectizeInput(session, 'metadata2',
                         choices = c(sort(colnames(datasetInput()@meta.data))), 
                         server = TRUE)
  })
  
  
  # This expression that generates a plot is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  
  output$metadata1Plot <- renderPlotly({
    req(input$dataset, input$viz, input$metadata1)
    
    seurat_object <- datasetInput()
    viz_method <- vizInput()
    if (is.factor(seurat_object@meta.data[[input$metadata1]])) {
      ggplotly(DimPlot(seurat_object, reduction = viz_method, group.by = input$metadata1) + ggtitle(input$metadata1))
      
    } else {
      ggplotly(FeaturePlot(seurat_object, reduction = viz_method, features = input$metadata1) + ggtitle(input$metadata1))
    }
    
    
  }) 
  
  output$metadata2Plot <- renderPlotly({
    req(input$dataset, input$viz, input$metadata2)
    
    seurat_object <- datasetInput()
    viz_method <- vizInput()
    if (is.factor(seurat_object@meta.data[[input$metadata2]])) {
      ggplotly(DimPlot(seurat_object, reduction = viz_method, group.by = input$metadata2) + ggtitle(input$metadata2))
      
    } else {
      ggplotly(FeaturePlot(seurat_object, reduction = viz_method, features = input$metadata2) + ggtitle(input$metadata2))
    }
    
  })
  
  output$genePlot <- renderPlotly({
    req(input$gene, input$viz)
    seurat_object <- datasetInput()
    viz_method <- vizInput()
    plot_title <- paste0(input$gene," expression")
    ggplotly(FeaturePlot(seurat_object, reduction = viz_method, features = c(input$gene), sort = input$sort) + ggtitle(plot_title))
    
  })
  
  output$metadata1VlnPlot <- renderPlotly({
    req(input$dataset,input$gene)
    
    factor_metadata <- is.factor(datasetInput()@meta.data[[input$metadata1]])
    req(factor_metadata)
    
    seurat_object <- datasetInput()
    plot_title <- paste0(input$gene," expression by ", input$metadata1)
    ggplotly(VlnPlot(seurat_object, features = c(input$gene), group.by = input$metadata1) + ggtitle(plot_title))
  })
  
  output$metadata2VlnPlot <- renderPlotly({
    req(input$dataset,input$gene)
    
    factor_metadata <- is.factor(datasetInput()@meta.data[[input$metadata2]])
    req(factor_metadata)
    
    seurat_object <- datasetInput()
    
    plot_title <- paste0(input$gene," expression by ",input$metadata2)
    ggplotly(VlnPlot(seurat_object, features = c(input$gene), group.by = input$metadata2) + ggtitle(plot_title))
  })
  
  # Create caption ----
  # The output$caption is computed based on a reactive expression 
  # that returns input$caption. When the user changes the
  # "caption" field:
  #
  # 1. This function is automatically called to recompute the output
  # 2. New caption is pushed back to the browser for re-display
  #
  # Note that because the data-oriented reactive expressions
  # below don't depend on input$caption, those expressions are
  # NOT called when input$caption changes
  # output$caption <- renderText({
  #   input$caption
  # })
  
  # Generate a summary of the dataset ----
  # The output$summary depends on the datasetInput reactive
  # expression, so will be re-executed whenever datasetInput is
  # invalidated, i.e. whenever the input$dataset changes
  # output$summary <- renderPrint({
  #   dataset <- datasetInput()
  #   summary(dataset)
  # })
  
  # Show the first "n" observations ----
  # The output$view depends on both the databaseInput reactive
  # expression and input$obs, so it will be re-executed whenever
  # input$dataset or input$obs is changed
  #output$view <- renderTable({
  #  head(datasetInput(), n = input$obs)
  #})
}


# Create Shiny app ----
shinyApp(ui = ui, server = server)