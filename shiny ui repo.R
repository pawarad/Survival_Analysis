library(shiny)
library(survival)
library(survminer)

data(lung)

ui <- fluidPage(
  titlePanel("Survival Analysis App"),
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", "Select Variable:",
                  choices = c("status", "sex"),
                  selected = "status"),
      actionButton("plot_button", "Generate Plot")
    ),
    mainPanel(
      plotOutput("km_plot")
    )
  )
)

server <- function(input, output) {
  
  km_plot <- eventReactive(input$plot_button, {
    if (input$variable == "status") {
      fit <- survfit(Surv(time, status) ~ 1, data = lung)
      ggsurvplot(fit, data = lung, palette = "jco", censor = FALSE)
    } else {
      fit <- survfit(Surv(time, status) ~ sex, data = lung)
      ggsurvplot(fit, data = lung, palette = "jco", censor = FALSE)
    }
  })
  
  output$km_plot <- renderPlot({
    km_plot()
  })
}

shinyApp(ui = ui, server = server)
