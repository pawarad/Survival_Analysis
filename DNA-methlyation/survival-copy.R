# Load necessary libraries
setwd("/home/ubuntu/DNA-meth")

library(survival)
library(survminer)
library(data.table)
library(dplyr)

TCGA_PAAD <- fread("TCGA.PAAD.sampleMap%2FHumanMethylation450")
TCGA_PAAD_trans <- transpose(TCGA_PAAD[,-1])

# get row and colnames in order
colnames(TCGA_PAAD_trans) <- TCGA_PAAD$sample
#rownames(TCGA_PAAD_trans) <- colnames(TCGA_PAAD)
TCGA_PAAD_trans_ <- cbind(sample = colnames(TCGA_PAAD)[-1], TCGA_PAAD_trans)

 
TCGA_cpg <- TCGA_PAAD_trans_ %>% select(sample, cg13332474)
clinical_PAAD <- fread("survival_data.txt")

## subset based on samples
#clinical_PAAD_ <- subset(clinical_PAAD, sample==TCGA_cpg$sample)
clinical_PAAD_ <- merge(clinical_PAAD,TCGA_cpg, by="sample")

clinical_PAAD_ <- clinical_PAAD_ %>%
  mutate(
    bin = case_when(
      cg13332474 >= mean(cg13332474) ~ 1,
      cg13332474 < mean(cg13332474) ~ 2
    ))
clinical_PAAD_[["cg13332474"]]
# Create a survival object
surv_object <- Surv(time = clinical_PAAD_$OS.time, event = clinical_PAAD_$OS)

# Fit a Cox proportional hazards model
cox_model <- coxph(surv_object ~ TCGA_cpg$cg13332474, data = TCGA_cpg)
summary(cox_model)

fit <- survfit(Surv(clinical_PAAD_$OS.time, clinical_PAAD_$OS) ~ clinical_PAAD_$bin, data = clinical_PAAD_)
ggsurvplot(fit, 
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(),
           palette = c("#990000", "#000099"))



#######---------Survival for Methylation -----------------

load("PAAD_meth.RData")

names(PAAD_meth)
# Rename specific columns
names(PAAD_meth)[names(PAAD_meth) == "Row.names"] <- "cpg_reg"
names(PAAD_meth)[names(PAAD_meth) == "UCSC_RefGene_Name"] <- "genes"
names(PAAD_meth)[names(PAAD_meth) == "UCSC_RefGene_Group"] <- "group"
names(PAAD_meth)[names(PAAD_meth) == "Relation_to_UCSC_CpG_Island"] <- "island"

## list of genes
total_genes <- unique(PAAD_meth$genes)

## Subset for specific cpg regions
PAAD_meth_gene <- subset(PAAD_meth, genes == "MSX2")
PAAD_meth_group <- subset(PAAD_meth_gene, group == "TSS1500")
PAAD_meth_island <- subset(PAAD_meth_group, island == "N_Shore")

## Trsnpose the specific region data
PAAD_meth_island_trans <- transpose(PAAD_meth_island[,-(1:4)])
PAAD_meth_island_trans_ <- cbind(sample = colnames(PAAD_meth_island)[-(1:4)], PAAD_meth_island_trans)

TCGA_PAAD_trans_ <- cbind(sample = colnames(TCGA_PAAD)[-1], TCGA_PAAD_trans)


PAAD_meth_group <- subset(PAAD_meth, genes == "SFRP1")
subset_group <- unique(PAAD_meth_group$group)

PAAD_meth_island <- subset(PAAD_meth_group, group == "TSS200")
subset_cpg <- unique(PAAD_meth_island$cpg_reg)


####------------------------ back up survival code ------------------------------
# Define server function  
server <- function(input, output) {
  
  ## Seletinput in server
  updateSelectizeInput(inputId = 'outlook', choices = as.list(cpgs), server = TRUE)
  
  updateSelectizeInput(inputId = 'Genes', choices = as.list(total_genes), server = TRUE)
  
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
          clinical_PAAD_1[[input$outlook]] >= mean(as.numeric(clinical_PAAD_1[[input$outlook]])) ~ 1,
          clinical_PAAD_1[[input$outlook]] < mean(as.numeric(clinical_PAAD_1[[input$outlook]])) ~ 2
        ))
    
    observeEvent(input$outlook, {
      print(paste0("You have chosen: ", input$outlook))
      print(head(clinical_PAAD_)) })
    
    return(clinical_PAAD_)
    
    
    # Create a survival object
    #surv_object <- Surv(time = clinical_PAAD_$OS.time, event = clinical_PAAD_$OS)
    
    # Fit a Cox proportional hazards model
    #cox_model <- coxph(surv_object ~ TCGA_cpg[,2], data = TCGA_cpg)
    #summary(cox_model)
    
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
