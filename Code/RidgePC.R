################################################################################
# Script for running Ridge regression with PCs as input instead of SNPs on
# simulated data
################################################################################
library(glmnet)
library(data.table)
library(arrow)
library(dplyr)

#################################################################################
# Prep by loading PCs

# Choose data set 180k or 70k
dset <- "180k"

# Load SNP_matrix as a regular data frame
X <- read_feather("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_180k.feather") %>%
  as.data.frame()


X <- X %>% rename(ringnr = IID)


# Excluding 'ringnr' column
SNPs <- X[, !colnames(X) %in% "ringnr"]

# Running PCA while standardizing the PCs
PCA <- prcomp(SNPs, center = TRUE)


# Compute variance explained by each PC
var_explained <- PCA$sdev^2 / sum(PCA$sdev^2)

# Cumulative variance explained
cum_var_explained <- cumsum(var_explained)

# Plot showing around 2500 PCs needed to explain 95% variance (almost all PCs)
plot(cum_var_explained)

# Store principal components
PCs <- PCA$x

# Testing re-scaling
# PCs <- PCs * PCA$sdev

# Convert PCs matrix to a data frame
PCs_df <- as.data.frame(PCs)

# Add 'ringnr' as the first column
PCs_df <- cbind(ringnr = X$ringnr, PCs_df)


X <- PCs_df

################################################################################
# Running Ridge models using PCs
{
# Choose architecture type and define globally
arch <- 9

# Load the CSV file where 10 folds for ringnr train and validation is
wd_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/Sim_",dset)
cv_folds_path <- paste0(wd_path,"/cv_folds_sim_", dset,"_arch_", arch, ".csv")
cv_folds <- read.csv(cv_folds_path)

pheno_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",
                     dset,"/Sim_pheno_",dset,"_arch_", arch, ".csv")
phenotype <- read.csv(pheno_path) %>%
  as.data.frame()

# Select relevant columns for y
y <- phenotype %>% select(pheno, ringnr)

# Initialize storage for prediction correlations
correlations <- numeric(10)

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
  ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, standardize = TRUE)  # alpha = 0 for Ridge
  
  # Best lambda from CV
  best_lambda <- ridge_cv$lambda.min
  
  # Fit the final Ridge model with the best lambda
  ridge_model <- glmnet(X_train, y_train, alpha = 0, lambda = best_lambda, standardize = TRUE)
  
  # Make predictions on the test set
  pred <- predict(ridge_model, X_test)
  pred <- as.vector(pred[, 1])
  
  # Calculate and store correlation for each fold
  correlations[fold] <- cor(pred, y_test)
  
  cat("Correlation fold ", fold, ": ", correlations[fold], "\n")
}


#boxplot(correlations)
correlations_df <- data.frame(corr = correlations)  # Create a data frame with one column 'corr'
save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/", dset,
                    "/RidgePC/corr_RidgePC_arch_",arch,"_std.csv")

write.csv(correlations_df, save_path, row.names = FALSE)
}
