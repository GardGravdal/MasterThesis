#######################################################################
# Script is based on code from Stefanie Muff
## This script is used to extract the data needed for to create the pseduo-response from the phenotype
# The script should not be run directly, but through the Makefile
# If ran manually, beware that two arguments are needed the phenotype ("mass" or "tarsus") and whether
# the 70K should be used or not (boolean argument)
#######################################################################

# CHANGE THIS TO YOUR OWN PATH: (i.e where the data is stored)
store_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/"
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Data")

# two inputs are needed for the script a phenotype (first argument, )
#args <- commandArgs(trailingOnly = TRUE)
 #phenotype <- args[1]

phenotype <- "tarsus"

#if (length(args) == 2) {
#    use_big_dataset <- TRUE
#    print("using big dataset")
#} else {
#   use_big_dataset <- FALSE
#}

# Packages needed for the script to run:

if (!require(nadiv)) {
    # install.packages("nadiv", repos = "http://cran.us.r-project.org", dependencies = TRUE)
   library(remotes); install_github("matthewwolak/nadiv", ref = "devel")
}
if (!require(pedigree)) {
    install.packages("pedigree", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}

if (!require(MASS)) {
    install.packages("MASS", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
if (!require(MCMCpack)) {
    install.packages("MCMCpack", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
if (!require(data.table)) {
    install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}

if (!require(SMisc)) {
    install.packages("C:/Users/gard_/OneDrive/Documents/MasterThesis/ProjectThesis/GenomicPrediction/src/data/SMisc.tar.gz", repos = NULL, type = "source")
}

if (!require(dplyr)) {
    install.packages("dplyr", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}

if (!require(lme4)) {
    install.packages("lme4", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}


if (!require(MCMCglmm)) {
    install.packages("MCMCglmm", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
if (!require(feather)) {
    install.packages("feather", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}






library(nadiv)
library(pedigree)
library(MASS)
library(MCMCpack)
library(MCMCglmm)
# This is a self-made package that I send you to install locally:
library(SMisc)

library(feather)

library(dplyr)


## Old packages no longer needed (but kept as comments, in case we need them in the future...)
# library(MasterBayes)
# library(pedigreemm)
# library(bdsmatrix)
# library(irlba)
# library(RSpectra)
# library(dplyr)

# big dataset refers to the 70k dataset. For 180k use_big_dataset = FALSE
use_big_dataset = TRUE
if (use_big_dataset == T) {
    print("using qc.raw")
    SNP_data_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/qc.raw"
    # loccation of the ringnr in the SNP matrix (qc.raw)
    ringnr_loc <- 1
    # where to save the phenotype-SNP data
    save_name <- paste(store_path, phenotype, "BV_70k.feather", sep = "")
    # save morphological data after some processing
    save_dd_name <- paste(store_path, phenotype, "Morph_70k.feather", sep = "")
    d.morph <- read.table("C:/Users/gard_/Documents/MasterThesis/Code/Data/AdultMorphology_20241117.csv", header = T, , sep = ",") # sep="\t")
    # Rename to match the previous dataset
    d.morph <- d.morph %>%
        rename(
            mass = body_mass,
            tarsus = thr_tarsus,
            wing = thr_wing,
            island_current = locality,
            hatchisland = first_locality,
            sex = adult_sex
        )
    d.morph$age <- d.morph$max_year - d.morph$hatch_year

    # Drop the southern_islands
    d.morph <- d.morph %>%
        dplyr::filter(!island_current %in% c(60, 61, 63, 67, 68))
} else {
    data_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/"
    SNP_data_path <- paste(data_path, "Helgeland_01_2018_QC.raw", sep = "")
    ringnr_loc <- 2
    save_name <- paste(store_path, phenotype, "BV.feather", sep = "")
    save_dd_name <- paste(store_path, phenotype, "Morph.feather", sep = "")
    # Data preparation helper script:
    # this creates d.morph and extracts the FGRM and the outer inner variable
    source("C:/Users/gard_/Documents/MasterThesis/Code/h_dataPrep.r")
    # Some data wranging to ensure that the IDs in the data correspond to the IDs in the A and G-matrices (nothing to worry about):
    # indicates that some IDs are missing:
    d.map[3110:3125, ]
    # from this we see the number of anmals
    Nanimals <- 3116


    # In the reduced pedigree only Nanimals out of the 3147 IDs are preset.
    d.map$IDC <- 1:nrow(d.map)

    d.morph$IDC <- d.map[match(d.morph$ringnr, d.map$ringnr), "IDC"]

    ### Prepare for use in INLA -
    d.morph$IDC4 <- d.morph$IDC3 <- d.morph$IDC2 <- d.morph$IDC
}




##############################################
### Preparations: Start with body mass and extract the id-value for each individual
##############################################

# keep all mass records that are not NA
d.pheno <- d.morph[!is.na(d.morph[phenotype]), c("ringnr", phenotype)]
names(d.pheno) <- c("ringnr", phenotype)

# RUN LMM ON PHENOTYPE TO SEPARETE ID EFFECT FROM ENVIRONMENTAL EFFECTS
d.mean.pheno <- as.data.frame(d.pheno %>%
    group_by(ringnr) %>%
    summarize(mean_pheno = mean(eval(as.symbol(phenotype)))))
if (use_big_dataset == T) {
    formula.pheno.lmm <- eval(as.symbol(phenotype)) ~   sex + month + age +
        (1 | island_current) +
        (1 | hatch_year) +
        (1 | ringnr)
} else {
    formula.pheno.lmm <- eval(as.symbol(phenotype)) ~   sex + FGRM + month + age + outer + other +
        (1 | island_current) +
        (1 | hatchyear) +
        (1 | ringnr)
}

library(lme4)
dd <- d.morph[!is.na(d.morph[phenotype]), ]
# save the morph data before the adjusting happens
write_feather(dd, save_dd_name)

# run LMM
r.pheno.lmer <- lmer(formula.pheno.lmm, data = dd,  control = lmerControl(optimizer ="Nelder_Mead"))
# Residuals
d.pheno.res <- data.frame(ringnr = dd$ringnr, pheno_res = residuals(r.pheno.lmer))
# Mean over repeats for each individual
d.mean.pheno.res <- as.data.frame(d.pheno.res %>%
    group_by(ringnr) %>%
    summarize(mean_pheno_res = mean(pheno_res)))

# ID effect
d.ID.pheno <- data.frame(ringnr = d.mean.pheno[, 1], ID.mass = ranef(r.pheno.lmer)$ringnr)

# We take as the new phenotype the estimated ID effect:
d.ID.pheno <- data.frame(ringnr = d.mean.pheno[, 1], ID = d.ID.pheno[, 2], mean_pheno = d.mean.pheno$mean_pheno)
tmp_hatchisland_ringnr_df <- d.morph[, c("ringnr", "hatchisland")]
tmp_hatchisland_ringnr_df <- tmp_hatchisland_ringnr_df %>%
    distinct(ringnr, .keep_all = T)
d.ID.pheno <- d.ID.pheno %>%
    merge(tmp_hatchisland_ringnr_df, by = "ringnr")

#############################################################
### Now we also load the raw SNP data matrix
#############################################################
library(data.table)

# Using the quality-controlled SNP matrix from Kenneth:
SNP.matrix <- data.frame(fread(SNP_data_path))

names(SNP.matrix)[ringnr_loc] <- "ringnr"
dim(SNP.matrix)
set.seed(323422)
sum(d.ID.pheno$ringnr %in% SNP.matrix$ringnr)



# Generate a data frame where individuals with ring numbers from d.ID.res.mass are contained, as well as the phenotype (here the residuals from the lmer analysis with mass as response)
d.dat.full <- merge(d.ID.pheno[, c("ringnr", "ID", "mean_pheno", "hatchisland")], SNP.matrix, by = "ringnr")
# SAVE THE FULL DATA SET:
write_feather(d.dat.full, save_name)

