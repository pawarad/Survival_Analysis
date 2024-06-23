head(PAAD_meth)
# Load necessary libraries
library(randomForestSRC)
library(gRbase)

# Remove rows where the name is "Doe"
PAAD_meth <- PAAD_meth[PAAD_meth$genes != "Unknown", ]

# Combine the two columns with a space in between
PAAD_meth$cpg_genes <- paste0(PAAD_meth$cpg_reg,"_", PAAD_meth$genes)

genes_PAAD <- sort(c(
  'TFAP2E', 'ASXL2', 'GSTM1', 'ONECUT2', 'TET1', 'HSPA5', 'SNRPF', 'PTGES2', 'CDA', 
  'SPARC', 'CES2', 'TSPAN1', 'ABCG2', 'BRCA1', 'CADM1', 'ABCB11', 'ABCB4', 'ABCC', 
  'ABCC10', 'ABCC3', 'ABCC5', 'SLC29A1', 'BCL2L1', 'ABCC2', 'FOXM1', 'TNFSF10', 
  'TM4SF1', 'DPYD', 'HMGA1', 'ABCB1', 'LCN2', 'SFN', 'RRM2', 'JAG1', 'NOTCH2', 
  'MUC5AC', 'SLC22A3', 'SLC29A3', 'MAP2', 'DCK', 'EGFR', 'MLH1', 'RRM1', 
  'SLC28A1', 'CHFR', 'MUC4', 'NRP1', 'SMARCA2', 'ABCC11', 'ARID1A', 'BNIP3', 
  'DKK3', 'IGFBP3', 'ISG15', 'LONRF2', 'MCL1', 'MUTYH', 'S100A4', 'SLC22A7', 
  'SLC28A2', 'SST', 'TYMS', 'NT5C1A', 'SLC22A2', 'SOX8', 'GSTM2', 'PROKR2', 
  'CTNNB1', 'ZEB1', 'VASH2', 'ABCB5', 'SLC28A3', 'MAP3K7', 'CXCL8', 'SLFN11', 
  'HNF1A', 'PYCARD',   
  'ADAMTS1', 'ACIN1', 'ALX4','BMP3', 'BNC1', 'CDKN2A', 'GSTP1', 'HIC1', 'MSX2', 
  'MEST', 'MGMT', 'MLH1', 'NEUROG1', 'NPTX2', 'PCDH10', 'PENK', 'POU4F1', 'SEMA5A', 
  'SEPTIN9', 'SFRP1', 'SFRP2', 'SPSB4', 'TFPI2', 'TNFRSF10C', 'WNT5A', 'UCHL1'))

gene.id <- PAAD_meth$genes
gene.num <- which(gene.id %in% genes_PAAD)
gene.nomatch <- genes_PAAD[-which(genes_PAAD %in% gene.id)]
gene.match <- genes_PAAD[which(genes_PAAD %in% gene.id)]
gene.match

gene.nomatch

## SElect only matched genes from DNA meth data
# Filter the data frame
PAAD_meth_match <- PAAD_meth %>%
  filter(genes %in% gene.match)


## transpose data for selected genes
columns_to_drop <- c("cpg_genes","cpg_reg","genes","group","island")
PAAD_meth_match_col <- PAAD_meth_match[, !(names(PAAD_meth_match) %in% columns_to_drop)]
PAAD_meth_match_mat <- as.data.frame(transpose(PAAD_meth_match_col))
colnames(PAAD_meth_match_mat) <- PAAD_meth_match$cpg_genes
PAAD_meth_match_mat <- cbind(sample = colnames(PAAD_meth_match_col), PAAD_meth_match_mat)

## Clean Clinical data
clinical_PAAD$sample <- clinical_PAAD$`_PATIENT`

## Combine survival and meth data
PAAD_meth_surv <- merge(clinical_PAAD,PAAD_meth_match_mat, by="sample")
columns_to_drop <- c("_PATIENT","sample","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time","Redaction")
PAAD_meth_surv <- select(PAAD_meth_surv, -any_of(columns_to_drop))


## Fit Random survival forest model
set.seed(1234)
rsf_model_multi <- rfsrc(Multivar(OS.time, OS) ~ ., data = PAAD_meth_surv,
                   save.memory = TRUE, perf.type = "none",
                   splitrule = "mahalanobis")

# Print the model summary
print(rsf_model_multi$forest)
plot(rsf_model_multi)

# Extract variable importance
var_importance <- vimp(rsf_model_multi)

# Print variable importance
print(var_importance)

# Get top genes based on variable importance
top_genes <- head(var_importance$importance, 10)
print(top_genes)

# Plot variable importance
plot(var_importance)

# If you want a more detailed visualization
library(ggplot2)

# Convert to data frame for ggplot2
var_importance_df <- data.frame(
  Gene = names(var_importance$importance),
  Importance = var_importance$importance
)

# Select top 10 genes for plotting
top_var_importance_df <- var_importance_df[order(-var_importance_df$Importance), ][1:10, ]

# Plot
ggplot(top_var_importance_df, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Genes by Variable Importance",
       x = "Genes",
       y = "Variable Importance")
# Fit Cox proportional hazards model with selected genes
cox_model <- coxph(Surv(OS.time, OS) ~ ., data = PAAD_meth_surv)

# Print the model summary
summary(cox_model)

# Function to fit Cox model for a given set of variables
fit_cox_model <- function(variables) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(variables, collapse = " + ")))
  cox_model <- coxph(formula, data = PAAD_meth_surv)
  return(cox_model)
}



variable_names <- colnames(PAAD_meth_match_mat)[-1]  
all_combinations <- lapply(1:4, function(i) RcppAlgos::comboGeneral(variable_names[1:100], i, repetition = FALSE))

###------------------------ elastic net ---------------

# Define a sequence of alpha values to test
alpha_values <- seq(0, 1, by = 0.1)
# lambdas <- 10^seq(-3,3,length=100) 

# Perform cross-validation for each alpha value
cv_results <- lapply(alpha_values, function(alpha) {
  cv.glmnet(as.matrix(X.train), S.train, family = "cox", alpha = alpha,type.measure = "deviance")
})

# Extract the mean cross-validated error for each alpha
mean_cv_errors <- sapply(cv_results, function(cv) min(cv$cvm))

# Find the alpha with the minimum cross-validated error
optimal_alpha <- alpha_values[which.min(mean_cv_errors)]
optimal_alpha

# Extract the cv.glmnet object for the optimal alpha
optimal_cv <- cv_results[[which.min(mean_cv_errors)]]

optimal_lambda <- optimal_cv$lambda.min
optimal_lambda

p.logrank.all <- rep(2,100)
fit.all <- rep(2,100)

for (i in 1:100) {
  print(paste0("iteration:",i))
  set.seed(i)
  X <- PAAD_meth_match_mat[,-1]
  time <- PAAD_meth_surv$OS.time
  time[time == 0] <- 1
  status <- as.numeric(PAAD_meth_surv$OS)
  S <- Surv(time,status)
  K <- 3
  fold <- split(sample(1:nrow(X)), rep(1:K, length=nrow(X)))
  id.train <- c(fold[[1]],fold[[2]])
  id.test <- fold[[K]]
  X.train <- X[id.train,]
  S.train <- S[id.train,]
  time.train <- time[id.train]
  status.train <- status[id.train]
  X.test <- X[id.test,]
  S.test <- S[id.test,]
  time.test <- time[id.test]
  status.test <- status[id.test]
  
  # Define a sequence of alpha values to test
  alpha_values <- seq(0, 1, by = 0.1)
  # lambdas <- 10^seq(-3,3,length=100) 
  
  # Perform cross-validation for each alpha value
  cv_results <- lapply(alpha_values, function(alpha) {
    cv.glmnet(as.matrix(X.train), S.train, family = "cox", alpha = alpha,type.measure = "deviance")
  })
  
  # Extract the mean cross-validated error for each alpha
  mean_cv_errors <- sapply(cv_results, function(cv) min(cv$cvm))
  
  # Find the alpha with the minimum cross-validated error
  optimal_alpha <- alpha_values[which.min(mean_cv_errors)]
  print(paste0("optimal_alpha",optimal_alpha))
  
  pen.lik <- function(alpha)
  {
    cv.pen <- try(cv.glmnet(as.matrix(X.train), S.train, family = "cox", alpha=alpha))
    print(paste0("alpha=",alpha))
    #genes selected with optimal model fitting
    coef.min <-  coef(cv.pen, s = "lambda.min")
    active.min <- which(coef.min != 0)
    index.min <- coef.min[active.min]
    sel.pen <- rownames(coef.min)[active.min]
    #risk score
    rs.train <- predict(cv.pen,newx=as.matrix(X.train),s="lambda.min")
    rs.test <- predict(cv.pen,newx=as.matrix(X.test),s="lambda.min")
    #lasso fit
    coxfit <- glmnet(X.train, S.train, family = "cox", alpha=alpha)
    return(list(cv.pen=cv.pen,coxfit=coxfit,rs.train=rs.train,rs.test=rs.test,sel.pen=sel.pen))
  }
  
  #elastic net method with all genes
  pen.elanet <- pen.lik(optimal_alpha)
  if (length(pen.elanet$sel.pen)==0) next
  cat(i)
  
  #median stratification
  good.prog <- (pen.elanet$rs.test < median(pen.elanet$rs.train))
  fit <- survfit(S.test ~ good.prog)
  # print(paste0("fit",fit))
  fit.all[i] <- fit
  
  
  #logrank test for test set
  logrank <- survdiff(S.test ~ good.prog)
  p.logrank <- 1-pchisq(logrank$chisq, 1)
  
  p.logrank.all[i] <- p.logrank
  
  #cox regression on risk score
  pi <- coxph(S.test ~ pen.elanet$rs.test, method = "breslow")
  p.pi <- 1-pchisq(pi$score, 1)
  print(paste0("p.logrank",p.logrank))
  
  if (p.logrank < 0.05){
    ## Plot KM
    #kaplan-meier curves
    
    # print(ggsurvplot(fit,data = PAAD_meth_match_mat,
    #            title = ("Kaplan-Meier Curve"),
    #            pval = TRUE, conf.int = TRUE,
    #            risk.table = TRUE, # Add risk table
    #            risk.table.col = "strata", # Change risk table color by groups
    #            linetype = "strata", # Change line type by groups
    #            surv.median.line = "hv", # Specify median survival
    #            ggtheme = theme_bw(),
    #            palette = c("#990000", "#000099")))
    # dev.off()
    
    plot(fit, lwd = 2, lty = c(1,1), col = c("red","blue"), xlab = 'Time', ylab = 'Estimated Survival Function')
    legend("topright", legend=c('high risk', 'low risk'), lty = c(1,1), col = c("red", "blue"), lwd = 2)
    text(1,0.05,paste("log-rank p=",round(p.logrank,4)))
    title("Kaplan-Meier Curve")
    
    print("Median OS:") 
    print(paste("high risk patient group:", round(unname(summary(fit)$table[,'median']*12)[1],2), "months"))
    print(paste("low risk patient group:", round(unname(summary(fit)$table[,'median']*12)[2],2), "months"))
    print("")
    print("")
    print(pen.elanet$sel.pen)
  }
}
