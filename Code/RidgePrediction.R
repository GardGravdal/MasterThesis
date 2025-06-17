##############################################################################
# Script to run Ridge regression prediction on genomic data.
# The script runs predictions on simulated phenotypes and calculates
# the accuracy score in Pearson correlation between predicted phenotype and 
# simulated phenotype. 
##############################################################################

library(glmnet)
library(data.table)
library(arrow)
library(dplyr)

# Choose architecture type and define globally
arch <- 1
# Choose data set 180k or 70k
dset <- "180k"

# Load the CSV file where 10 folds for ringnr train and validation is
wd_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/Sim_",dset)
setwd(wd_path)
cv_folds_path <- paste0("cv_folds_sim_", dset,"_arch_", arch, ".csv")
cv_folds <- read.csv(cv_folds_path)

# Load SNP_matrix as a regular data frame
X <- read_feather("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_180k.feather") %>%
  as.data.frame()

pheno_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",
                     dset,"/Sim_pheno_",dset,"_arch_", arch, ".csv")
phenotype <- read.csv(pheno_path) %>%
  as.data.frame()

X <- X %>% rename(ringnr = IID)


# Select relevant columns for y
y <- phenotype %>% select(pheno, ringnr)


# Initialize storage for prediction correlations
correlations <- numeric(10)
# Standardize SNPs?
standard = FALSE


# Loop through each fold
for (fold in 1:10) {
  cat("Fold ", fold, "\n")
  # Identify training and test indices
  test_indices <- cv_folds[cv_folds$Fold == fold & cv_folds$Set == "test", "ringnr"]
  train_indices <- cv_folds[cv_folds$Fold == fold & cv_folds$Set == "train", "ringnr"]
  
  # Split data, and remove 'ringnr' column
  X_train <- X[X$ringnr %in% train_indices, -which(names(X) == "ringnr")]  
  X_test <- X[X$ringnr %in% test_indices, -which(names(X) == "ringnr")]    
  y_train <- y[y$ringnr %in% train_indices, "pheno"]  
  y_test <- y[y$ringnr %in% test_indices, "pheno"]
  
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)
  y_test <- as.numeric(y_test)
  
  # Fit Ridge Regression (use cross-validation to select the best lambda on the training set)
  # Standardize = TRUE is set as deufalt
  ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, standardize = standard)  # alpha = 0 for Ridge
  
  # Best lambda from CV
  best_lambda <- ridge_cv$lambda.min
  
  # Fit the final Ridge model with the best lambda
  ridge_model <- glmnet(X_train, y_train, alpha = 0, lambda = best_lambda, standardize = standard)
  
  # Make predictions on the test set
  pred <- predict(ridge_model, X_test)
  pred <- as.vector(pred[, 1])
  
  # Calculate and store correlation for each fold
  correlations[fold] <- cor(pred, y_test)
  
  cat("Correlation fold ", fold, ": ", correlations[fold], "\n")
}


cat("Pearson Correlations:", correlations, "\n")
correlations_df <- data.frame(corr = correlations)  # Create a data frame with one column 'corr'

if (standard){
  save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/", dset,
                    "/corr_Ridge_arch_",arch,"_",dset,"_std",".csv")
} else {
  save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/", dset,
                    "/corr_Ridge_arch_",arch,"_",dset,".csv")
}
write.csv(correlations_df, save_path, row.names = FALSE)
