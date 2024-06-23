# Load necessary libraries
library(survival)
library(survminer)
library(data.table)
library(dplyr)

TCGA_PAAD <- fread("TCGA.PAAD.sampleMap%2FHumanMethylation450")
TCGA_PAAD_trans <- transpose(TCGA_PAAD[,-1])

# get row and colnames in order
colnames(TCGA_PAAD_trans) <- TCGA_PAAD$sample
cpgs <- colnames(TCGA_PAAD_trans)
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

# Create a survival object
surv_object <- Surv(time = clinical_PAAD_$OS.time, event = clinical_PAAD_$OS)

# Fit a Cox proportional hazards model
cox_model <- coxph(surv_object ~ TCGA_cpg$cg13332474, data = TCGA_cpg)
summary(cox_model)

fit <- survfit(Surv(clinical_PAAD_$OS.time, clinical_PAAD_$OS) ~ clinical_PAAD_$bin, data = clinical_PAAD_)

########--------------------- Survival dynamic DNA meth ------------------------

load("PAAD_meth.RData")

names(PAAD_meth)
# Rename specific columns
names(PAAD_meth)[names(PAAD_meth) == "Row.names"] <- "cpg_reg"
names(PAAD_meth)[names(PAAD_meth) == "UCSC_RefGene_Name"] <- "genes"
names(PAAD_meth)[names(PAAD_meth) == "UCSC_RefGene_Group"] <- "group"
names(PAAD_meth)[names(PAAD_meth) == "Relation_to_UCSC_CpG_Island"] <- "island"

## list of genes
total_genes <- unique(PAAD_meth$genes)

#############------------------------- DEG --------------------------------------

setwd("/home/ubuntu/RNA-seq")

## Read deg data
data <- read.csv("DEG_results_Primary_Tumor vs. Solid_Tissue_Normal.csv",
                 header = TRUE, sep = ",")

## Rename column
names(data)[1] <- "Genes"
data <- data %>%
  mutate(
    minusLog10Pvalue = -log10(P.Value)
  )