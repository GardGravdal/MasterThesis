################################################################################
#                         GWAS vs Effects comparison

# This script generates nice Manhattan plots of the variance components of the
# phenotypes from simulated effect sizes. The variance components used are 
# \beta_j^2 cov(x_j) for each SNP j. These are compared to GWAS and SHAP results
# respectively. This script generates the plots that compare effects and GWAS 
# results. See get_SHAP_Effects.r for the SHAP and effects comparison.
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



setwd("C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/GWAS")

# Set architecture globally
arch <- 2
dset <- "70k"

SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_",dset,".feather"))
effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",dset,"/Sim_effect_",dset,"_arch_",arch,".csv"))


# Loading and preparing GWAS results for plotting
# ------------------------------------------------------------------------------
gemma_path <- paste0("GWASresults",arch,".lmm.assoc.txt")
resultGemma <- read_table(gemma_path)

# Some reformatting on the stored data:
# ---------------------------------
# Gemma stores lines in GWASresults across two lines... Need to fix
txt <- readLines(gemma_path)
str(txt)

txt <- txt[-1]

# merge every 2 subsequent lines into one to form a row of final dataframe
idx <- seq(1, length(txt), by=2)
txt[idx] <- paste(txt[idx], txt[idx+1])
txt <- txt[-(idx+1)]

# Final data
resultGemma <- read.table(text=txt, col.names=colnames(resultGemma))

# Exclude certain chromosomes in plot
resultGemma <- resultGemma[!(resultGemma$chr %in% c(0, 16,30,32)),]
# ------------------------------------------------------------------------------

head(resultGemma)
head(effects)

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
SNP_subset <- SNPs[,c(snp_to_keep)]

# Compute variances for each column (Should not be necessary after first time - see read line next)
#variances <- apply(SNP_subset, 2, var)
#saveRDS(variances, "C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Read variances for 70k SNP data
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

head(variances)
head(effects)
dim(variances)
# Check that they follow the same ordering before plotting
head(variances)
head(effects$SNP)
names(SNP_subset)[1:6]

# Create the variance terms to be compared to SHAP/GWAS values
var_term <- effects$effect^2 * variances

# Create variance term data frame
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term, row.names = NULL)


variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)

# Check that order is still retained
head(variance_term_df)
head(resultGemma)

# Check that dimensions are the same
dim(variance_term_df)
dim(resultGemma)


# NOW WE CAN PLOT
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures/GWAS_Arch")

# Compute the Bonferroni threshold
bonferroni <- -log10(0.05/ nrow(resultGemma))

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
title(ylab = expression(beta^2 ~ Var(X)), line = 2, cex.lab = 2.3)

# Second plot: Manhattan plot for resultGemma
manhattan(resultGemma, chr = "chr", bp = "ps", p = "p_lrt", snp = "rs", ylim = c(0, 7), 
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, 
          ylab = "", genomewideline = FALSE, suggestiveline = FALSE, 
          highlight = max_effect_snp_name)

# Add the Bonferroni line and a title for the y-axis for resultGemma
abline(h = bonferroni, col = "red", lty = 2)
title(ylab = TeX(r"($-log_{10}(p)$)"), line = 2, cex.lab = 2.3)









################################################################################
# Plot without highlight
#setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures/GWAS_Arch")

#par(mar = c(5, 6, 4, 1) + .1)
# Set up a 1x2 layout for the subplot
#par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)  # Adjust margins to fit both plots

#manhattan(variance_term_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
#          logp = FALSE, ylab = "", ylim = c(0,0.03), col = c("gray10", "gray50"),
#          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE, 
#          suggestiveline = FALSE, yaxt = "n")
#title(ylab = expression(beta[j]^2 ~ cov(X[j])), line = 2, cex.lab = 2.3)
#axis(2, at = seq(0, 0.03, by = 0.005), labels = format(seq(0, 0.03, by = 0.005),
#                                                       digits = 2), cex.axis = 1.2)

#manhattan(resultGemma,chr="chr",bp="ps",p="p_lrt",snp="rs", ylim = c(0, 7), col = c("gray10", "gray50"),
#          cex.axis = 1.5, cex.lab=2.3, cex.main=2.5, ylab="", genomewideline = FALSE, suggestiveline = FALSE)
#abline(h = bonferroni, col = "red", lty = 2)
#title(ylab = TeX(r"($-log_{10}(p)$)"),line=2, cex.lab=2.3)

# Reset par
#par(mar = c(5, 6, 4, 1) + .1)