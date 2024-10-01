# Clinical Data Visualization app

options(shiny.maxRequestSize=75*1024^2)
# Created a package to source supporting functions
source("PFRpackage.R")

library(shiny)
library(lme4)
library(doBy)
library(multcomp)
library(stringr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(tidyr)
library(ggsignif)
library(DT)
library(htmltools)
library(shinydashboard)
library(shinythemes)
library(bslib)

ui <- 
  tabsetPanel(
    
    # Create tab page to display results of submitted dataset
    tabPanel("Results",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          # Input bar for used-submitted CSV file
          fileInput("file", "Choose CSV File", accept = ".csv"),
          # Input for mixed linear model
          textInput("formula1", "Model 1:", 
                    value = "zScore ~ -1 + (1|Treatment) + (1|BioRep) + (1|TechRep)"),
          actionButton("update1", "Update Model 1"),
          textOutput("updated_code1"),
          textOutput("message1"),
          
          textInput("formula2", "Model 2:", 
                    value = "zScore ~ Treatment -1 + (1|BioRep) + (1|BioRep:TechRep)"),
          actionButton("update2", "Update Model 2"),
          textOutput("updated_code2"),
          textOutput("message2"),
          
          textInput("formula3", "Model 3:", 
                    value = "Log2Y ~ Treatment -1 + (1|BioRep) + (1|BioRep:TechRep) + (1|BioRep:TechRep)"),
          actionButton("update3", "Update Model 3"),
          textOutput("updated_code3"),
          textOutput("message3"),
          br(),
          # Create button to trigger functions for data analysis
          actionButton("runbutton", "Process"),
          # Select different datapoints to view qualitative analyses/visualizations
          uiOutput("proteoform_selector"),
          uiOutput("boxtreatment_selector"),
          uiOutput("volctreatment_selector"),
          radioButtons('show_biotech', 'Boxplot: Show biological/technical replicates?', c('NO', 'YES'),
                       inline = TRUE),
          fileInput("samplefile", "Sample file for PCA (optional)", accept = ".csv"),
          # Designate output file type
          radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                       inline = TRUE),
          # Create button to download onto local machine for email-sendable report
          downloadButton("downloadReport")
        ),
        mainPanel(
          # Value boxes to give quantitative data about measured data points 
          fluidRow(
            # Value Box 1 for number of proteoforms analyzed
            valueBoxOutput(outputId = "box_1", width = 3),
            # Value Box 2 for number of significant proteoforms found
            valueBoxOutput(outputId = "box_2", width = 3),
            # Value Box 3 for Missing Val percentage
            valueBoxOutput(outputId = "box_3", width = 3)
          ),
          br(),# HTML implementation to separate between UI sections
          hr(),
          column(6, style = 'border: 1px solid lightgrey; border-radius: 25px',
                 # title and info button
                 div(HTML('<b>Boxplot</b> '), style = 'display: inline-block;'),
                 # map plot within box
                 plotOutput('boxplot')
          ),
          column(6, style = 'border: 1px solid lightgrey; border-radius: 25px',
                 # title and info button
                 div(HTML('<b>Volcano Plot</b> '), style = 'display: inline-block;'),
                 # map plot within box
                 plotOutput('volcanoplot')
          ),
          # Use DT to create an interactive table for user to sort through columns
          DTOutput("signiftable"),
        )
      )
    )
  ),
  tabPanel("Heatmap", 
           fluidPage(
             plotOutput("heatmap")
           )),
  # Create tab for quality control of data
  tabPanel("QC",
           fluidPage(
             mainPanel(
               fluidRow(
                 column(6, style = 'border: 1px solid lightgrey; border-radius: 25px',
                        # ntitle and info button
                        div(HTML('<b>Mass Distribution</b> '), style = 'display: inline-block;'),
                        # map plot
                        plotOutput('massdist')
                 ),
                 column(6, style = 'border: 1px solid lightgrey; border-radius: 25px',
                        # ntitle and info button
                        div(HTML('<b>PCA</b> '), style = 'display: inline-block;'),
                        # map plot
                        plotOutput('PCA')
               )
             )
           )
        )
  )
)
server <- function(input, output, session) {
  # Update each model based on user input (3 models)
  observeEvent(input$update1, {
    updated_code <- update_model(input$formula1, 1)
    output$updated_code1 <- renderText(updated_code)
    update_script(updated_code, "model1")
    output$message1 <- renderText("Model 1 code updated successfully.")
    print(input$formula1)
  })
  
  observeEvent(input$update2, {
    updated_code <- update_model(input$formula2, 2)
    output$updated_code2 <- renderText(updated_code)
    update_script(updated_code, "model2")
    output$message2 <- renderText("Model 2 code updated successfully.")
  })
  
  observeEvent(input$update3, {
    updated_code <- update_model(input$formula3, 3)
    output$updated_code3 <- renderText(updated_code)
    update_script(updated_code, "model3")
    output$message3 <- renderText("Model 3 code updated successfully.")
  })
  # wrap input data in reactive expression to update when user input changes
  processed_data <- reactive({
    req(input$runbutton)
    # Read in file
    inFile <- input$file
    # Require input CSV file to process
    req(inFile)
    df <- read.csv(inFile$datapath)
    # Apply statistical models to data
  
    data <- MLMapply(df)
    # create list of returned data from package functions for easy access
    list(
      # data frame to make boxplots
      data4box = data$data4box,
      # data frame for various qualitative analyses
      finalData = data$finalData,
      # Percentage of missing values from data
      missing_vals = data$missing,
      boxdfs = MLMBoxPlot(data$data4box, data$finalData)
    )
    
  })
  
  observe({
    data <- processed_data()
    PFRs <- data$boxdfs$PFRs
    signif_PFRs <- subset(data$boxdfs$sub, qValue < 0.05, select = c(PFR))
    signif_PFRs <- unique(signif_PFRs$PFR)
    missing_vals <- data$missing
    # Display general dataset metrics on dashboard
    output$box_1 <- shinydashboard::renderValueBox({
      valueBox(length(PFRs), "Number of PFRs analyzed")
    })
    output$box_2 <- shinydashboard::renderValueBox({
      valueBox(length(signif_PFRs), "Number of significant PFRs found (p.adj < 0.05)")
    })
    output$box_3 <- shinydashboard::renderValueBox({
      valueBox(paste0(round(missing_vals, 3), "%"), "Missing values")
    })
    output$proteoform_selector <- renderUI({
      selectInput(
        inputId = "prots",
        label = "Choose a proteoform:",
        choices = PFRs
      )
    })
  })
  
  # Render the boxplot based on the selected proteoform
  output$boxplot <- renderPlot({
    data <- processed_data()
    sub <- data$boxdfs$sub
    
    req(input$prots)  # Ensure that a proteoform is selected
    pfr <- input$prots
    
    filtered <- sub[sub$PFR == pfr,]
    
    ggplot(filtered, aes(x = Treatment, y = Log2Y, fill = Treatment)) +
      geom_boxplot() +
      geom_jitter(shape = 16, position = position_jitter(0.2)) +
      scale_color_brewer(palette = "Dark2") +
      labs(title = pfr, x = "Treatment", y = "Log Intensity") +
      theme_classic() +
      geom_signif(comparisons = list(c("HigheGFR", "LoweGFR")), map_signif_level = TRUE, 
                  annotations = paste0("p.adj:",round(filtered$qValue[1], 4)))
  })
  # Volcano plot with proteoform display based on user input
  output$volcanoplot <- renderPlot({
    data <- processed_data()
    finalData <- data$finalData
    req(input$prots)  # Ensure that a proteoform is selected
    pfr <- input$prots
  
    comps <- unique(finalData %>% dplyr:: select(c('Treatment1','Treatment2')))
    
    MLMvolcanoPlot(comps[,"Treatment2"],comps[,"Treatment1"], finalData, pfr = pfr) 
  })
  # Create interactive table with only significant data points
  output$signiftable <- renderDT({
    data <- processed_data()
    finalData <- data$finalData
    signif_table <- subset(finalData, qValue < 0.05, select = c(PFR,
                                                                qValue,
                                                                TargetMass,
                                                                Description)) 
    signif_table <- signif_table[!duplicated(signif_table),] 
    
    # Use datatable() from DT package
    datatable(signif_table, options = list(autoWidth = TRUE))
  })
  # work in progress, interactive heatmap for user-defined treatment groups
  output$heatmap <- renderPlot({
    data <- processed_data()
    plotly()
  })
  reactive_data <- reactive({
    req(input$prots) 
    data <- processed_data()
    finalData <- data$finalData
    pfr <- input$prots
    sub <- data$boxdfs$sub 
    filtered <- sub[sub$PFR == pfr, ]
    
    list(pfr = pfr, filtered = filtered)
  })
  output$downloadReport <- downloadHandler(
    
    filename = function() {
      paste('my-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },
    # Upload content to remote file, input each input as params
    content = function(file) {
      src <- normalizePath('report.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      data <- reactive_data()
      filtered <- data$filtered
      # Pass data or input to the report
      params_list <- list(
        input_formula1 = input$formula1,
        input_formula2 = input$formula2,
        input_formula3 = input$formula3,
        pfr = data$pfr,
        in_prot = filtered,
        finalData = data$finalData
      )
      # Give multiple options for file output type
      library(rmarkdown)
      out <- render('report.Rmd',
                    params = params_list,
                    output_format = switch(
                      input$format,
                      PDF = pdf_document(),
                      HTML = html_document(),
                      Word = word_document()
                    ))
      file.rename(out, file)
    }
  )

  # Plot data distribution of mass values (sample file)
  output$massdist <- renderPlot({
    data <- processed_data()
    finalData <- data$finalData
    massdata <- finalData
    QCdist(massdata)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
