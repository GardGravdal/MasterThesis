################################################################################
# Script for running Ridge regression with PC's from PCA as input. This script
# runs analysis on phenotypic data from the house sparrow 180k dataset.
# Generalized for phenotypes mass and wing (set manually). Returns accuracy in form
# of Pearson correlation between test and prediction.
################################################################################

library(glmnet)
library(data.table)
library(arrow)
library(dplyr)


#################################################################################
# Load phenotypic data, SNPs and CV-data

pheno <- "mass"
# 180K:
pheno_path <- sprintf("C:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/%sBV.feather", pheno)
data <- as.data.table(read_feather(pheno_path))

# Check that there are no duplicated ringnr's
print(any(duplicated(data$ringnr)))

cv_path <- paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/cv_folds_",pheno,"_180k.csv")
cv_folds <- read.csv(cv_path)

# 70k setup: Remove IID, 180k setup: remove FID. Else, everything is the same. 
# Sim: FID, but not hatchisland, ID, and mean_pheno
# X_CV is ringnr's + SNPs


X_CV <- data[, !c("ID", "mean_pheno", "FID", "MAT", "PAT", "SEX", "PHENOTYPE", "hatchisland"), with = FALSE] %>%
  as.data.frame()

# Some of the SNPs have NA-values. Set to 0
# This loop performs the equivalent of "X_CV[is.na(X_CV)] <- 0", but runs much faster for large matrices
for (j in seq_along(X_CV)) {
  set(X_CV, which(is.na(X_CV[[j]])), j, 0)
}

# Change from float to integer for all columns except 'ringnr' (i.e., all SNPs)
cols_to_convert <- setdiff(names(X_CV), "ringnr")
X_CV[cols_to_convert] <- lapply(X_CV[cols_to_convert], as.integer)

# y is ringnrs + pseudo phenotype (mean_pheno)
# rename mean_pheno to simply 'pheno' to generalize code
y <- data %>%
  select(mean_pheno, ringnr, ID) %>%
  as.data.frame()

################################################################################
# Prep by loading PCs

# Excluding 'ringnr' column
SNPs <- X_CV[, !colnames(X_CV) %in% "ringnr"]


# Running PCA while standardizing the PCs
PCA <- prcomp(SNPs, center = TRUE)


# Compute variance explained by each PC
var_explained <- PCA$sdev^2 / sum(PCA$sdev^2)

# Cumulative variance explained
cum_var_explained <- cumsum(var_explained)

# Plot showing around 2500 PCs needed to explain 95% variance (almost all PCs)
#plot(cum_var_explained)

# Store principal components
PCs <- PCA$x

# Testing re-scaling
# PCs <- PCs * PCA$sdev

# Convert PCs matrix to a data frame
PCs_df <- as.data.frame(PCs)

# Add 'ringnr' as the first column
PCs_df <- cbind(ringnr = X_CV$ringnr, PCs_df)


X <- PCs_df

################################################################################

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
  y_train <- y[y$ringnr %in% train_indices, "ID"]  
  y_test <- y[y$ringnr %in% test_indices, "mean_pheno"]
  
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)
  y_test <- as.numeric(y_test)
  
  # Fit Ridge Regression (use cross-validation to select the best lambda on the training set)
  ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, standardize = FALSE)  # alpha = 0 for Ridge
  
  # Best lambda from CV
  best_lambda <- ridge_cv$lambda.min
  
  # Fit the final Ridge model with the best lambda
  ridge_model <- glmnet(X_train, y_train, alpha = 0, lambda = best_lambda, standardize = FALSE)
  
  # Make predictions on the test set
  pred <- predict(ridge_model, X_test)
  pred <- as.vector(pred[, 1])
  
  # Calculate and store correlation for each fold
  correlations[fold] <- cor(pred, y_test)
  
  cat("Correlation fold ", fold, ": ", correlations[fold], "\n")
}


#boxplot(correlations)
correlations_df <- data.frame(corr = correlations)  # Create a data frame with one column 'corr'
save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno/corr_",pheno,"_RidgePC.csv")
write.csv(correlations_df, save_path, row.names = FALSE)
