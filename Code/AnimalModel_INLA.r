
library(nadiv)
library(pedigree)
library(MASS)
library(MCMCpack)
library(MCMCglmm)
# This is a self-made package that I send you to install locally:
library(SMisc) # Take contact if you do not have this
library(dplyr)
library(INLA)
library(feather)


# Load the CSV file where 10 folds for ringnr train and validation is
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Data/CVfolds")

cv_folds <- read.csv('cv_folds_tarsus_180k.csv')

# Data preparation helper script:
setwd("C:/Users/gard_/Documents/MasterThesis/Code")
source("h_dataPrep.r")


# Some data wranging to ensure that the IDs in the data correspond to the IDs in the A and G-matrices (nothing to worry about):
# indicates that some IDs are missing:
#d.map[3110:3125, ]
# from this we see the number of animals
Nanimals <- 3116


# remove missing values (mass)
d.morph <- filter(d.morph, !is.na(eval(tarsus)))

# names(d.morph)
# In the reduced pedigree only Nanimals out of the 3147 IDs are preset.
d.map$IDC <- 1:nrow(d.map)

d.morph$IDC <- d.map[match(d.morph$ringnr, d.map$ringnr), "IDC"]

### Prepare for use in INLA -
d.morph$IDC4 <- d.morph$IDC3 <- d.morph$IDC2 <- d.morph$IDC


formula.tarsus <- tarsus ~ sex + FGRM + month + age + outer + other +
f(hatchisland, model = "iid", hyper = list(
    prec = list(initial = log(1), prior = "pc.prec", param = c(1, 0.05))
)) +
f(hatchyear, model = "iid", hyper = list(
    prec = list(initial = log(1), prior = "pc.prec", param = c(1, 0.05))
)) +
f(IDC, model = "iid", hyper = list(
    prec = list(initial = log(1), prior = "pc.prec", param = c(1, 0.05))
)) +
f(IDC2,
    values = 1:3116, model = "generic0",
    Cmatrix = Cmatrix,
    constr = TRUE,
    hyper = list(
        # The priors are relevant, need to discuss
        prec = list(initial = log(0.5), prior = "pc.prec", param = c(sqrt(2), 0.05))
    )
)

corr_cvs_EG <- c()
corr_cvs_G <- c()

# Relatedness matrix from Henrik (vanRaden method 1) where +0.01 was already added to diagnoal!
d.Gmatrix <- read.table(paste("C:/Users/gard_/Documents/MasterThesis/Code/Data/", "gghatvr3.triangle.g", sep = ""), header = F, sep = " ")


# keep only relatednesses that are relevant for animals in d.morph
# d.Gmatrix <- d.Gmatrix[d.Gmatrix[,1] %in% d.morph$ID & d.Gmatrix[,2] %in% d.morph$ID, ]

# G is a sparse matrix object. We can also verify that it is symmetric (sometimes some numerical problems lead to non-symmetry)
G <- sparseMatrix(i = d.Gmatrix[, 1], j = d.Gmatrix[, 2], x = d.Gmatrix[, 3], symmetric = T)
# Indicing to keep structure
G[, ] <- as.numeric(G[, ])
isSymmetric(G)

# Again extract the rows and columns for the individuals in the data set that we analyse
dim(G)
GG <- G[d.map[1:3116, 3], d.map[1:3116, 3]]


# To ensure that the matrix is positive definite, we do a computational trick (proposed by vanRaden 2008, see https://doi.org/10.3168/jds.2007-0980 :)
# diag(n) returns (nxn) identity matrix
AAA <- diag(dim(GG)[1])
GGG <- GG * 0.99 + 0.01 * AAA # replace by Identity matrix

# Need to derive the inverse to give to INLA
Cmatrix <- solve(GGG)
if (!isSymmetric(Cmatrix)) {
    Cmatrix <- forceSymmetric(Cmatrix)
}

# Run 10-fold cross-validation
# Make temp train/test set, to test the function of GBLUP/INLA
# Need to fix test/training set later
#head(cv_folds)

for (i in 0:9) {
    # get CV indices
    #ringnr_train <- pull(arrow::read_feather(paste(data_path, "temp/ringnr_train_", i, ".feather", sep = "")), "ringnr")
    #ringnr_test <- pull(arrow::read_feather(paste(data_path, "temp/ringnr_test_", i, ".feather", sep = "")), "ringnr")
    ringnr_train <- cv_folds %>%
        filter(Fold == (i + 1) & Set == "train") %>%
        pull(ringnr)
    ringnr_test <- cv_folds %>%
        filter(Fold == (i + 1) & Set == "test") %>%
        pull(ringnr)

    # make test and train set
    d.morph_train <- filter(d.morph, !ringnr %in% ringnr_test)
    d.morph_test <- filter(d.morph, ringnr %in% ringnr_test)

    n_train <- dim(d.morph_train)[1]
    n_test <- dim(d.morph_test)[1]
    N <- n_train + n_test

    # Save the phenotypic value in the test set, if we only looking at genetic effects (two-step) we take the average
    pheno_test_EG <- d.morph_test[, "tarsus"]
    pheno_test <- as.data.frame(d.morph_test %>%
        group_by(ringnr) %>%
        summarize(
            mean_pheno = mean(eval(tarsus))
        ))[, "mean_pheno"]

    # However, INLA has no predict function, so have to fill the test-values with NAs and then merge it back into the train-set
    d.morph_test[, "tarsus"] <- NA
    d.morph_train <- union_all(d.morph_train, d.morph_test)

    names(d.morph_train)
    # All individuals
    idxs <- 1:Nanimals
    # get the indices corresponding to the individuals in the test set
    idxs_test <- which(d.map$ringnr %in% unique(ringnr_test))

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
    model1.tarsus <- inla(
        formula = formula.tarsus, family = "gaussian",
        data = d.morph_train,
        control.family = list(hyper = list(theta = list(initial = log(0.5), prior = "pc.prec", param = c(sqrt(2), 0.05)))),
        control.compute = list(dic = F, return.marginals = FALSE), verbose = TRUE
        # control.compute=list(config = TRUE)
    )
    cat("INLA DONE\n")



    # get predicted phenotype
    #preds_EG <- model1.wing$summary.fitted.values$mean[(n_train + 1):N]

    # Get breeding values
    preds <- model1.tarsus$summary.random$IDC2$mode[idxs_test]

    # calculate and save metrics
    #corr_EG <- cor(preds_EG, pheno_test_EG, method = "pearson")
    corr <- cor(preds, pheno_test, method = "pearson")
    cat("result of fold", i, "corr_G:", corr, "\n")
    corr_cvs_G <- c(corr_cvs_G, corr)
    #corr_cvs_EG <- c(corr_cvs_EG, corr_EG)
}

# save results


#EG_csv <- data.frame(Value = corr_cvs_EG)
BV_csv <- data.frame(corr = corr_cvs_G)

#BV_csv
# Write the vector to a CSV file
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno")
BV_file_name <- paste0("corr_BV_tarsus.csv")
#write.csv(EG_csv, file = "corr_EG_WING.csv", row.names = FALSE)
write.csv(BV_csv, file = BV_file_name, row.names = FALSE)

BV <- read.csv(BV_file_name)
#EG <- read.csv("corr_EG_WING.csv")
#xgb <- read.csv("xgBoostCVHyperparams.csv")
#mean_pheno_corrs <- read.csv("mean_pheno_corrs.csv")
#xgb_corrs <- xgb$corr_val
#print(xgb_corrs)

#boxplot(xgb_corrs,ylim = c(0.15, 0.35), main = "xgb", col = "blue")
boxplot(BV, ylim = c(0.25, 0.5), main = "BV", col = "yellow")
#boxplot(EG, ylim = c(0.2,0.5), main = "EG")
#boxplot(mean_pheno_corrs, ylim = c(0.15, 0.35), main = "mean_pheno", col = "blue")

# read the result df
#result_df <- arrow::read_feather(path_to_results)
# we store both results of the predicted breeding value (G) and the predicted phenotype (EG)
#INLA_result_df <- data.frame(name = mod_name_EG, corr = corr_cvs_EG, phenotype = phenotype)
#INLA_result_df <- rbind(INLA_result_df, data.frame(name = mod_name_G, corr = corr_cvs_G,phenotype = phenotype))
#result_df <- rbind(result_df, INLA_result_df)


#arrow::write_feather(result_df, path_to_results)