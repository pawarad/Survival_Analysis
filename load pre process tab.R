library(shiny)

ui <- fluidPage(
  titlePanel("Load Data and Pre-process Tabs"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Load Data", fileInput("file", "Choose CSV File")),
        tabPanel("Pre-process", 
                 checkboxInput("normalize", "Normalize Data", value = FALSE),
                 checkboxInput("impute", "Impute Missing Values", value = FALSE),
                 numericInput("filter", "Filter Genes by Count", value = 0, min = 0)
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Summary", verbatimTextOutput("summary")),
        tabPanel("Pre-processed Data", verbatimTextOutput("processed_data"))
      )
    )
  )
)

server <- function(input, output) {
  
  uploaded_data <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })
  
  output$summary <- renderPrint({
    summary(uploaded_data())
  })
  
  pre_processed_data <- reactive({
    data <- uploaded_data()
    if (input$normalize) {
      
    }
    if (input$impute) {
   
    }
    if (input$filter > 0) {
      
    }
    data
  })
  
  output$processed_data <- renderPrint({
    summary(pre_processed_data())
  })
}

shinyApp(ui = ui, server = server)
