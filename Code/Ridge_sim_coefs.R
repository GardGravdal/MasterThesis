################################################################################
                              # Ridge sim coeffs
# This script generates coefs of ridge regression model used on the whole dataset.
# Can be used later for inference to be compared to e.g. SHAP, GWAS or phenotypic
# variance component. Can also be used to generate CFI plots. This script is 
# meant for simulated traits. See RidgePheno_coefs.R for coefs of real traits.
# Script is also meant for 70k dataset.
################################################################################


library(data.table)
library(arrow)
library(dplyr)
library(glmnet)

library(glmnet)
library(data.table)
library(arrow)
library(dplyr)

# Choose architecture type and define globally
arch <- 1
# Choose data set 180k or 70k
dset <- "70k"

# Load the CSV file where 10 folds for ringnr train and validation is
wd_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/Sim_",dset)
setwd(wd_path)

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

# Train with ID effect
y_train <- y[,"pheno"]
# Exclude column ringnr in X
X_train <- X[, -which(names(X_CV) == "ringnr")]
# Test set is the same as training in this case (task is not prediction)
X_test <- X[, -which(names(X_CV) == "ringnr")]
# Test with mean_pheno
y_test <- y[,"pheno"]

#----------------

X_train <- as.matrix(X_train)
X_test <- as.matrix(X_test)
y_train <- as.numeric(y_train)
y_test <- as.numeric(y_test)

# We do not standardize here
standard <- FALSE

# Find optimal lambda_R
ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, standardize = standard)  # alpha = 0 for Ridge

# Best lambda from CV
best_lambda <- ridge_cv$lambda.min

# Now we train and extract the coefficients

# Fit the final Ridge model with the best lambda
ridge_model <- glmnet(X_train, y_train, alpha = 0, lambda = best_lambda, standardize = standard)

coefs <- as.vector(coef(ridge_model))
coef_names <- rownames(coef(ridge_model))  # This includes "(Intercept)" as the first entry

# Create dataframe
coef_df <- data.frame(
  SNP = coef_names,
  beta_R = coefs,
  row.names = NULL
)

# Remove intercept, only interested in SNPs
coef_df <- coef_df[coef_df$SNP != "(Intercept)", ]

head(coef_df)

file_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/Ridge/"

#----

write.csv(coef_df, paste0(file_path,"ridge_coeffs_",pheno,".csv"), row.names = FALSE)