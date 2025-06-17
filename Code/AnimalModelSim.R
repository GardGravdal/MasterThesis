################################################################################
# Markov animal model
# This script is based on BayesAnimalSim.R and is an adaptation to be used for 
# parallel running on Markov. Generates predictions from Bayesian Animal model
# with INLA and saves correlations.
################################################################################
library(nadiv)
library(pedigree)
library(MASS)
library(MCMCpack)
library(MCMCglmm)
# This is a self-made package that I send you to install locally:
#library(SMisc) # Take contact if you do not have this
library(dplyr)
library(INLA)
library(arrow)


{
# Set global architecture
arch <- 10
# Load the CSV file where 10 folds for ringnr train and validation is
cv_folds_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds/Sim_180k/","cv_folds_sim_180k_arch_", arch, ".csv")
cv_folds <- read.csv(cv_folds_path)

# Load simulated phenotypes and effect sizes
pheno_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/180k/", "Sim_pheno_180k_arch_", arch, ".csv")
#effects_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/180k/", "Sim_effect_180k_arch_", arch, ".csv")
# (n_animals x 2) First col phenotypes, second col ringnr
phenotypes <- read.csv(pheno_path)


################################################################################
Nanimals <- dim(phenotypes)[1]

# Creating the equivalent for d.sim in AnimalModel_INLA.R
d.sim <- phenotypes

# Some modifications (changing order) such that the columns are ordered correctly
d.sim["pheno"] <- phenotypes["ringnr"]
d.sim["ringnr"] <- phenotypes["pheno"]
# rename pheno to "sim"
colnames(d.sim) <- c("ringnr", "sim")
# Now d.sim is correctly ordered:
head(d.sim)

### Prepare for use in INLA -
d.sim$IDC <- 1:Nanimals

d.sim$IDC4 <- d.sim$IDC3 <- d.sim$IDC2 <- d.sim$IDC
head(d.sim)

# Creating the formula to be used in INLA, only intercept (population average) 
# and the inverse G-matrix included
formula.sim <- sim ~ 1 +
  f(IDC, model = "iid", hyper = list(
    prec = list(initial = log(1), prior = "pc.prec", param = c(1, 0.05))
  )) +
  f(IDC2,
    values = 1:Nanimals, model = "generic0",
    Cmatrix = Cmatrix,
    constr = TRUE,
    hyper = list(
      # The priors are relevant, need to discuss
      prec = list(initial = log(0.5), prior = "pc.prec", param = c(sqrt(2), 0.05))
    )
  )
}
################################################################################
{
# Relatedness matrix from get_Gmatrix.R (vanRaden method 1) where +0.01 was already added to diagnoal!
# Load G_sparse
G_sparse <- read.csv("C:/Users/gard_/Documents/MasterThesis/Code/Data/G_sparse.csv")

# G is a sparse matrix object. We can also verify that it is symmetric 
# (sometimes some numerical problems lead to non-symmetry)
G <- sparseMatrix(i = G_sparse[, 1], j = G_sparse[, 2], x = G_sparse[, 3], symmetric = T)
G[, ] <- as.numeric(G[, ])
isSymmetric(G)



# To ensure that the matrix is positive definite, we do a computational trick 
# (proposed by vanRaden 2008, see https://doi.org/10.3168/jds.2007-0980 :)
AAA <- diag(dim(G)[1])
GGG <- G * 0.99 + 0.01 * AAA # replace by Identity matrix

# Need to derive the inverse to give to INLA
Cmatrix <- solve(GGG)
if (!isSymmetric(Cmatrix)) {
  Cmatrix <- forceSymmetric(Cmatrix)
}
}
################################################################################
# Run 10-fold cross-validation

corr_csv_G <- c()

{
for (i in 0:9) {
  # get CV indices
  ringnr_train <- cv_folds %>%
    filter(Fold == (i + 1) & Set == "train") %>%
    pull(ringnr)
  ringnr_test <- cv_folds %>%
    filter(Fold == (i + 1) & Set == "test") %>%
    pull(ringnr)
  
  # make test and train set
  d.sim_train <- filter(d.sim, !ringnr %in% ringnr_test)
  d.sim_test <- filter(d.sim, ringnr %in% ringnr_test)
  
  
  n_train <- dim(d.sim_train)[1]
  n_test <- dim(d.sim_test)[1]
  N <- n_train + n_test
  
  # Save the phenotypic value in the test set
  pheno_test <- as.data.frame(d.sim_test["sim"])
  
  # However, INLA has no predict function, so have to fill the test-values with NAs and then merge it back into the train-set
  d.sim_test[, "sim"] <- NA
  d.sim_train <- union_all(d.sim_train, d.sim_test)
  
  
  # get the indices corresponding to the individuals in the test set
  idxs_test <- which(phenotypes$ringnr %in% ringnr_test)
  
  ##################################################################
  ### Run INLA based on the GBLUP approach
  ###
  ### To this end, use the GRM (genetic relatedness matrix) in the animal model
  ###
  ### !!! This is very slow - account for at least 20-30min waiting time before inla terminates !!!
  ##################################################################
  
  
  ##
  ## INLA formula
  ##
  # Here we use body mass as the response, and some fixed and random effects:
  
  cat("Starting INLA\n")
  model1.sim <- inla(
    formula = formula.sim, family = "gaussian",
    data = d.sim_train,
    control.family = list(hyper = list(theta = list(initial = log(0.5), prior = "pc.prec", param = c(sqrt(2), 0.05)))),
    control.compute = list(dic = F, return.marginals = FALSE), verbose = FALSE
    # control.compute=list(config = TRUE)
  )
  cat("INLA DONE\n")
  
  # Get breeding values
  preds <- model1.sim$summary.random$IDC2$mode[idxs_test]
  
  # calculate and save metrics
  corr <- cor(preds, pheno_test, method = "pearson")
  cat("result of fold", i, "corr_G:", corr, "corr_EG", "\n")
  corr_csv_G <- c(corr_csv_G, corr)
}

# save results
BV_csv <- data.frame(corr = corr_csv_G)

# Write the vector to a CSV file

BV_path <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/AnimalModel/corr_BV_arch_",arch,".csv")

write.csv(BV_csv, file = BV_path, row.names = FALSE)
}

BV <- read.csv(BV_path)

#boxplot(xgb_corrs,ylim = c(0.15, 0.35), main = "xgb", col = "blue")
boxplot(BV, ylim = c(0.25, 0.5), main = "BV", col = "yellow")

# read the result df
#result_df <- arrow::read_feather(path_to_results)
# we store both results of the predicted breeding value (G) and the predicted phenotype (EG)
#INLA_result_df <- data.frame(name = mod_name_EG, corr = corr_cvs_EG, phenotype = phenotype)
#INLA_result_df <- rbind(INLA_result_df, data.frame(name = mod_name_G, corr = corr_csv_G,phenotype = phenotype))
#result_df <- rbind(result_df, INLA_result_df)


#arrow::write_feather(result_df, path_to_results)