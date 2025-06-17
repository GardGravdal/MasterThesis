################################################################################
#                         GWAS vs SHAP real trait

# This script generates nice Manhattan plots of GWAS and SHAP results from real
# traits. The script is generalized so that it works for all three traits: mass,
# tarsus and wing. This script only works for 70k dataset since GWAS and SHAP is
# only done on this dataset for real traits.
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



# Set phenotype globally
pheno <- "tarsus"

shap_vals <- read_feather(paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/SHAP/shap_",pheno,"_70k.feather"))

setwd("C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/GWAS")
# Loading and preparing GWAS results for plotting
# ------------------------------------------------------------------------------
gemma_path <- paste0("GWASresults_",pheno,".lmm.assoc.txt")
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

# Some prepping is needed to match the SHAP values with chromosome locations
# NED A MAP FILE, gives where each SNP is located
map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"

chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")


snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

SNP_cols <- names(shap_vals)

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


# Find SNPs in tmpdf$snpID that are not in resultGemma$rs
snps_not_in_gwas <- setdiff(tmpdf$snpID, resultGemma$rs)

# Make sure that GWAS and SHAP contains the same SNPs for comparison
snp_to_discard <- unique(c(snp_to_discard, snps_not_in_gwas))

snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]

# we now have a mapping between SHAP SNPs and the chromosome location
sum(chro$snpID %in% tmpdf$snpID)


# remove some SNP on chromosomes not desired for both SHAP and GEMMA
shap_vals <- shap_vals[, c(snp_to_keep)]


# Create dataframe with SHAP values
shap_frame <- data.frame(origsnp = rownames(t(shap_vals)), abs_shap = as.numeric(shap_vals[1, ]), row.names = NULL)

# Create the data frames with bp and chr information
shap_frame_df <- merge(shap_frame, tmpdfchro, by = "origsnp", sort = FALSE)


# Check that order is still retained
head(shap_frame_df)
head(resultGemma)

# Check that dimensions are the same
dim(shap_frame_df)
dim(resultGemma)


# NOW WE CAN PLOT
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures/SHAP_GWAS")

# Compute the Bonferroni threshold
bonferroni <- -log10(0.05/ nrow(resultGemma))

################################################################################
# Plot

# Set the y-axis of effects plot relative to max effect size
shap_y_limit <- max(shap_frame_df$abs_shap) * 1.1

# Create Manhattan plots in a 1x2 layout
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)  # Adjust margins to fit both plots

# First plot: Manhattan plot for variance_term_df
manhattan(shap_frame_df, chr = "chr", bp = "bp", p = "abs_shap", snp = "snpID", ylim = c(0, shap_y_limit),
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, ylab = "",
          logp = FALSE, genomewideline = FALSE, suggestiveline = FALSE, yaxt = "n")
axis(2, at = seq(0, shap_y_limit, by = shap_y_limit/3), labels = format(seq(0, shap_y_limit, by = shap_y_limit/3), digits = 2), cex.axis = 1.2)
title(ylab = "mean |SHAP|", line = 2, cex.lab = 2.3)


# Second plot: Manhattan plot for resultGemma
manhattan(resultGemma, chr = "chr", bp = "ps", p = "p_lrt", snp = "rs", ylim = c(0, 7),
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5,
          ylab = "", genomewideline = FALSE, suggestiveline = FALSE)
abline(h = bonferroni, col = "red", lty = 2)
title(ylab = TeX(r"($-log_{10}(p)$)"), line = 2, cex.lab = 2.3)