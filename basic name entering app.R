library(shiny)

ui <- fluidPage(
  titlePanel("An shiny app"),
  sidebarLayout(
    sidebarPanel(
      textInput("txtInput","Input the text to display")
    ),
    mainPanel(
      paste("You are entering"),
        textOutput("textOutput")
    )
  )  
)

server <- shinyServer(function(input,output){
  output$txtOutput <- renderText({
    paste(input$txtInput)
  })
})

shinyApp(ui=ui, server = server)

