################################################################################
# Lasso coeffs
# This script generates coefs of lasso regression model used on the whole dataset.
# Can be used later for inference to be compared to e.g. SHAP, GWAS or phenotypic
# variance component. Can also be used to generate CFI plots. This script is 
# generalized such that it works for both simulated traits and real traits.
# Script is made for 70k dataset, 180k analysis not performed.
################################################################################



library(data.table)
library(arrow)
library(dplyr)
library(glmnet)

# Choose whether analyze simulated trait (TRUE) or real trait (FALSE):
sim <- TRUE

# Real traits
if (sim == FALSE){
  # Set phenotype (mass, wing, tarsus) globally
  pheno <- "wing"
  # 180K:
  pheno_path <- sprintf("C:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/%sBV_70k.feather", pheno)
  data <- as.data.table(read_feather(pheno_path))
  
  # Check that there are no duplicated ringnr's
  print(any(duplicated(data$ringnr)))
  
  
  # 70k setup: Remove IID, 180k setup: remove FID. Else, everything is the same. 
  # Sim: FID, but not hatchisland, ID, and mean_pheno
  # X_CV is ringnr's + SNPs
  X_CV <- data[, !c("ID", "mean_pheno", "IID", "MAT", "PAT", "SEX", "PHENOTYPE", "hatchisland"), with = FALSE] %>%
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
  
  # Train with ID effect
  y_train <- y_CV[,"ID"]
  # Exclude column ringnr in X
  X_train <- X_CV[, -which(names(X_CV) == "ringnr")]
  # Test set is the same as training in this case (task is not prediction)
  X_test <- X_CV[, -which(names(X_CV) == "ringnr")]
  # Test with mean_pheno
  y_test <- y_CV[,"mean_pheno"]
}

# Simulation study
if (sim == TRUE){
  # Choose architecture type and define globally
  arch <- 2
  # Choose data set 180k or 70k
  dset <- "70k"
  
  # Load the CSV file where 10 folds for ringnr train and validation is
  wd_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/Sim_",dset)
  setwd(wd_path)
  
  # Load SNP_matrix as a regular data frame
  X <- read_feather("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_70k.feather") %>%
    as.data.frame()
  
  pheno_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",
                       dset,"/Sim_pheno_",dset,"_arch_", arch, ".csv")
  phenotype <- read.csv(pheno_path) %>%
    as.data.frame()
  
  # Select relevant columns for y
  y <- phenotype %>% dplyr::select(pheno, ringnr)
  
  # Train with ID effect
  y_train <- y[,"pheno"]
  # Exclude column ringnr in X
  X_train <- X[, -which(names(X) == "ringnr")]
  # Test set is the same as training in this case (task is not prediction)
  X_test <- X[, -which(names(X) == "ringnr")]
  # Test with mean_pheno
  y_test <- y[,"pheno"]
}

X_train <- as.matrix(X_train)
X_test <- as.matrix(X_test)
y_train <- as.numeric(y_train)
y_test <- as.numeric(y_test)

# We do not standardize here
standard <- FALSE

# Find optimal lambda_R
lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1, standardize = standard)  # alpha = 1 for lasso

# Best lambda from CV
best_lambda <- lasso_cv$lambda.min

# Now we train and extract the coefficients

# Fit the final lasso model with the best lambda
lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = best_lambda, standardize = standard)

coefs <- as.vector(coef(lasso_model))
coef_names <- rownames(coef(lasso_model))  # This includes "(Intercept)" as the first entry

# Create dataframe
coef_df <- data.frame(
  SNP = coef_names,
  beta_L = coefs,
  row.names = NULL
)

# Remove intercept, only interested in SNPs
coef_df <- coef_df[coef_df$SNP != "(Intercept)", ]

head(coef_df)

file_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/Lasso/"

if (sim == FALSE){
  write.csv(coef_df, paste0(file_path,"lasso_coefs_",pheno,".csv"), row.names = FALSE)
} else {
  write.csv(coef_df, paste0(file_path,"lasso_coefs_arch",arch,".csv"), row.names = FALSE)
}
