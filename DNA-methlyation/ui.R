# Load R packages
library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(fontawesome)

options(shiny.maxRequestSize=100*1024^2)
# Define UI
ui <- fluidPage(theme = shinytheme("united"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "The Ohio State University",
                  tabPanel("DNA Methylation survival",
                           sidebarPanel(
                             tags$h3("Input:"),
                             selectInput("Genes", label = "Genes:", 
                                         choices = NULL, 
                                         selected = "Rainy"),
                             
                             selectInput("Region", label = "Group:", 
                                         choices = NULL, 
                                         selected = "Rainy"),
                             
                             # selectInput("Island", label = "island region:", 
                             #             choices = NULL, 
                             #             selected = "Rainy"),
                             
                             selectInput("outlook", label = "cpgs:", 
                                         choices = NULL, 
                                         selected = "Rainy"),
                             
                             downloadButton("downloadPlot", "Download Plot as PDF")
                             
                             
                             
                           ), # sidebarPanel
                           mainPanel(
                             h2("The data used is publicly available TCGA methylation data using HM450K array."),
                             
                            h3("1) Select the gene of interest"),
                            h3("2) Select the group of interest"),
                                
                            h4("a)TS200 - region from the transcription start site (TSS) to -200 nucleotides upstrem of TSS"),
                            h4("b) TSS1500 - Region covering -200 to -1500 nucleotides upstream of TSS"),
                            h4("c) first exon"),
                            h4("d) 5' UTR"),
                            h4("e) 3' UTR"),
                            h4("f) body"),
                            h3("One CpG site can be annotated to multiple regions of a gene due to altenative TSSs"),
                          
                             h1("Survival Analysis"),
                             
                             h4("Kaplan Meier Curve"),
                             verbatimTextOutput("txtout"),
                             plotOutput(outputId = "KMplot")
                             
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  
                  
                  tabPanel("RNA-seq DEG",
                              type = "tabs",
                           tabsetPanel(
                             id = "tabsetPanelID",
                             type = "tabs",
                             tabPanel("Features", tabsetPanel(
                               tabPanel("Vocano Plot",
                                        plotOutput(outputId = "volcanoplot",click= "plot click"),
                                        DTOutput("dataTable")), 
                               tabPanel("pathway analysis")
                             ))
                           )), # Navbar 2, tabPanel
                  tabPanel("Navbar 3", "This panel is intentionally left blank")
                  
                ) # navbarPage
) # fluidPage