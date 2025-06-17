################################################################################
                                        # Ridge pheno

# This script performs Ridge regression on phenotypic data from the house sparrow
# 180k dataset. The script  is generalizes such that it can be used for all three
# traits wing, mass and tarsus which are set globally at the start of the
# script. Accuracy in form of Pearson correlation from the 10-fold CV is stored
# at the end of the script. This script uses data from dataloader.R, which
# performs the so-called two-step procedure where environmental and other 
# effects are removed.
################################################################################


library(data.table)
library(arrow)
library(dplyr)
library(glmnet)


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
y_CV <- data %>%
  select(mean_pheno, ringnr, ID) %>%
  as.data.frame()


# Initialize storage for prediction correlations
correlations <- numeric(10)
# Decide whether to use standardization for SNPs or not
standard = FALSE

# Loop through each fold
for (fold in 1:1) {
  cat("Fold ", fold, "\n")
  # Identify training and test indices
  test_indices <- cv_folds[cv_folds$Fold == fold & cv_folds$Set == "test", "ringnr"]
  train_indices <- cv_folds[cv_folds$Fold == fold & cv_folds$Set == "train", "ringnr"]
  
  # Split data, and remove 'ringnr' column
  X_train <- X_CV[X_CV$ringnr %in% train_indices, -which(names(X_CV) == "ringnr")]  
  X_test <- X_CV[X_CV$ringnr %in% test_indices, -which(names(X_CV) == "ringnr")]    
  y_train <- y_CV[y_CV$ringnr %in% train_indices, "ID"]  
  y_test <- y_CV[y_CV$ringnr %in% test_indices, "mean_pheno"]
  
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)
  y_test <- as.numeric(y_test)
  
  # Fit Ridge Regression (use cross-validation to select the best lambda on the training set)
  # Generalized for use with and without standardization of SNPs
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
plot(order(ridge_model$beta))

cat("Pearson Correlations:", correlations, "\n")
correlations_df <- data.frame(corr = correlations)  # Create a data frame with one column 'corr'

# Generate different name for the results generated with standardization of SNPs
if (standard){
    save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno/corr_",pheno,"_Ridge_std.csv")
} else{
    save_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno/corr_",pheno,"_Ridge.csv")
}
#write.csv(correlations_df, save_path, row.names = FALSE)

