# Define server function  
server <- function(input, output) {
  
  ## Seletinput in server
  # updateSelectizeInput(inputId = 'outlook', choices = as.list(cpgs), server = TRUE)
  # 
  updateSelectizeInput(inputId = 'Genes', choices = total_genes, server = TRUE)
  
  ## Reactive function to select cpgs
  observeEvent(input$Genes,{
    ## Subset for specific cpg regions
    PAAD_meth_group <- subset(PAAD_meth, genes == input$Genes)
    subset_group <- unique(PAAD_meth_group$group)
    
    # PAAD_meth_island <- subset(PAAD_meth_group, group == input$Region)
    # subset_island <- unique(PAAD_meth_island$island)
    
    PAAD_meth_cpg <- subset(PAAD_meth_group, group == input$Region)
    subset_cpg <- unique(PAAD_meth_cpg$cpg_reg)
    
    # PAAD_meth_island <- subset(PAAD_meth_group, island == "N_Shore")
    updateSelectizeInput(inputId = 'Region', choices = subset_group, server = TRUE)
    # updateSelectizeInput(inputId = 'Island', choices = subset_island, server = TRUE)
    updateSelectizeInput(inputId = 'outlook', choices = as.factor(subset_cpg), server = TRUE)
    
  })
  
  
  ## Reactive function to select cpgs
  clinical_PAAD_ <- reactive({
    TCGA_cpg <-TCGA_PAAD_trans_ %>% select(sample, input$outlook)
    clinical_PAAD <- fread("survival_data.txt")
    
    ## subset based on samples
    #clinical_PAAD_ <- subset(clinical_PAAD, sample==TCGA_cpg$sample)
    clinical_PAAD_1 <- merge(clinical_PAAD,TCGA_cpg, by="sample")
    
    
    clinical_PAAD_ <- clinical_PAAD_1 %>%
      dplyr::mutate(
        bin = case_when(
          clinical_PAAD_1[[input$outlook]] >= mean(as.numeric(clinical_PAAD_1[[input$outlook]])) ~ "high",
          clinical_PAAD_1[[input$outlook]] < mean(as.numeric(clinical_PAAD_1[[input$outlook]])) ~ "low"
          ))
    
    observeEvent(input$outlook, {
      print(paste0("You have chosen: ", input$outlook))
      print(head(clinical_PAAD_)) })
    
    return(clinical_PAAD_)

  
  # observeEvent(input$outlook, {
  #   print(paste0("You have chosen firse: ", input$outlook))
  # })
  
  
  })
  
  
  
  output$KMplot <- renderPlot({
    fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~ bin, data = clinical_PAAD_())
    # 
    observeEvent(input$outlook, {
      print(paste0("You have chosen: ", input$outlook))
      print(fit)
    })
    ggsurvplot(fit, data = clinical_PAAD_(),
               title = paste0(input$Genes,":",input$outlook),
               pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(),
               palette = c("#990000", "#000099"))
  })
  
  # observeEvent(input$outlook, {
  #     print(paste0("You have chosen firse: ", head(data)))
  #   })
  
  # Download the plot as PDF
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("survival_plot", Sys.Date(), input$Genes,":",input$outlook,".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print( ggsurvplot(fit, data = clinical_PAAD_(),
               title = paste0(input$Genes,":",input$outlook),
               pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(),
               palette = c("#990000", "#000099")))
      dev.off()
    }
  )
  
  output$volcanoplot <- renderPlot({
      plot <- data %>%
        ggplot(aes(x = logFC,
                   y = minusLog10Pvalue,
                   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")),
                   text = Genes,
                   key = Genes)) +
        geom_point() +
        xlab("log fold change") +
        ylab("-log10(P-value)")
      
      plot
   })
  
  # Identify the clicked point and filter the data table accordingly
  observeEvent(input$plot_click, {
    click <- input$plot_click
    if (is.null(click)) return()
    
    # Identify the nearest point to the click
    selected <- data %>%
      mutate(distance = (logFC - click$x)^2 + (minusLog10Pvalue - click$y)^2) %>%
      filter(distance == min(distance)) %>%
      select(-distance)
    
    # Update the data table to show only the selected point
    output$dataTable <- renderDT({
      datatable(selected, options = list(pageLength = 10))
    })
  })
  
  # Initial rendering of the full data table
  output$dataTable <- renderDT({
    datatable(data, options = list(pageLength = 5))
  })

  
} # server
