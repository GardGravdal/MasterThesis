################################################################################
#                         SHAP vs Effects comparison

# This script generates nice Manhattan plots of the variance components of the
# phenotypes from simulated effect sizes. The variance components used are 
# \beta_j^2 cov(x_j) for each SNP j. These are compared to GWAS and SHAP results
# respectively. This script generates the plots that compare effects and SHAP 
# results. See GWAS_Effects_manhattan.r for the SHAP and effects comparison.
################################################################################

library(qqman)
library(tidyr)
library(dplyr)
library(stringr)
library(grid)
library(gridGraphics)
library(ggplot2)
library(dplyr)
library(arrow)
library(readr)
library(data.table)
library(smplot2)
library(latex2exp)



#setwd("C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/GWAS")

# Set architecture globally [1,2,...,10], where 1 is polygenic and 10 oliogenic
arch <- 10
# Choose dataset 180k or 70k
dset <- "70k"

SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_",dset,".feather"))
shap_vals <- read_feather(paste0('c:/Users/gard_/Documents/MasterThesis/Code/Results/',dset,'/SHAP/shap_',dset,'_arch_',arch,'.feather'))
effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",dset,"/Sim_effect_",dset,"_arch_",arch,".csv"))



# Some prepping is needed to match the SHAP values with chromosome locations
# NED A MAP FILE, gives where each SNP is located
if (dset == "70k"){
  map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
} else {
  map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
}

chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")


snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

SNP_cols <- effects["SNP"]
tmpdf <- data.frame(origsnp = SNP_cols, snpID = NA)
names(tmpdf) <- c("origsnp","snpID")

# Need to split the SNP name (in the shap df) to get the basepair
get_bp <- function(snp) {
  snp <- str_split(snp, "[_]")[[1]][1]
  return(snp[[1]])
}


tmpdf$snpID <- lapply(tmpdf$origsnp, get_bp)

# make df instead of list
tmpdf2 <- as.data.frame(lapply(tmpdf, unlist))

# merge the map-file stuff with the SNPname from the SNP columns in shap/gwas
tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")


snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]


# we now have a mapping between SHAP SNPs and the chromosome location
sum(chro$snpID %in% tmpdf$snpID)

# Use the desired subset
effects <- effects[effects$SNP %in% snp_to_keep, ]
shap_vals <- shap_vals[, names(shap_vals) %in% snp_to_keep, drop = FALSE]
SNP_subset <- SNPs[,c(snp_to_keep)]

# Re-order the names in shap_vals to align with the order in effects
shap_vals <- shap_vals[, effects$SNP, drop = FALSE]

# Compute variances for each column (Should not be necessary after first time - see read line next)
#variances <- apply(SNP_subset, 2, var)
#saveRDS(variances, "C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Read variances for 70k SNP data
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Check that they follow the same ordering before plotting
head(variances)
head(effects$SNP)
names(SNP_subset)[1:6]
names(shap_vals)[1:6]

# Create the variance terms to be compared to SHAP/GWAS values
var_term <- effects$effect^2 * variances

# Create variance term data frame and shap data frame
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term, row.names = NULL)
shap_frame <- data.frame(origsnp = effects$SNP, abs_shap = as.numeric(shap_vals[1, ]), row.names = NULL)

# Create the data frames with bp and chr information
variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)
shap_frame_df <- merge(shap_frame, tmpdfchro, by = "origsnp", sort = FALSE)

# Check that order is still retained
head(variance_term_df)
head(shap_frame_df)


# NOW WE CAN PLOT
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures/SHAP_Arch")


################################################################################
# Plot which highlights SNP with highest effect size

# Find the SNP with the highest effect size
max_effect_snp <- variance_term_df[which.max(variance_term_df$snp_effects), ]
max_effect_snp_name <- max_effect_snp$snpID
# Set the y-axis of effects plot relative to max effect size
y_limit <- max_effect_snp$snp_effects + 0.005


# Create Manhattan plots in a 1x2 layout
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)  # Adjust margins to fit both plots

# First plot: Manhattan plot for variance_term_df
manhattan(variance_term_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE, 
          suggestiveline = FALSE, yaxt = "n", highlight = max_effect_snp_name)
axis(2, at = seq(0, y_limit, by = 0.005), labels = format(seq(0, y_limit, by = 0.005),
                                                       digits = 2), cex.axis = 1.2)


# Add a title for the y-axis
title(ylab = expression(beta[j]^2 ~ Var(X[j])), line = 2, cex.lab = 2.3)

# Second plot: Manhattan plot for SHAP
manhattan(shap_frame_df,chr="chr",bp="bp",p="abs_shap",snp="snpID", ylim = c(0, 0.023), col = c("red1", "red4"),
          cex.axis = 1.5, cex.lab=2.3, cex.main=2.5, ylab="", highlight = max_effect_snp_name,
          logp = FALSE, genomewideline = FALSE, suggestiveline = FALSE, yaxt = "n")
axis(2, at = seq(0, y_limit, by = 0.005), labels = format(seq(0, y_limit, by = 0.005),
                                                          digits = 2), cex.axis = 1.2)
title(ylab = "mean |SHAP|", line = 4, cex.lab = 2.3)