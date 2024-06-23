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


plot <- data %>%
  ggplot(aes(x = logFC,
             y = minusLog10Pvalue,
             # colour = Pattern_number,
             text = Genes,
             key = Genes)) +
  geom_point() +
  xlab("log fold change") +
  ylab("-log10(P-value)")

plot %>%
  ggplotly(tooltip = "tooltip") %>%
  layout(dragmode = "select")



