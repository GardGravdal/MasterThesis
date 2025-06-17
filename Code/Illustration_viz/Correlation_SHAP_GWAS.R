################################################################################
#                       Corr SHAP vs. GWAS real traits
# Script for generating correlation plot betweem GWAS p-values and mean |SHAP|
# importance measure for real traits. Generalized to work for all traits (mass,
# tarsus and wing).
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
pheno <- "wing"
pheno_name <- "Wing" # Capitalize first letter to be used in plotting

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

# trait 'wing' contains some NA's in p_lrt. Remove them
if (pheno == "wing"){
  resultGemma <- resultGemma[!is.na(resultGemma$p_lrt), ]
}
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

# Fit a linear model 
lm_model <- lm(shap_frame_df$abs_shap ~ I(-log10(resultGemma$p_lrt)))

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, length(shap_frame_df$abs_shap), length.out = 10)

# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

slope <- coef(lm_model)[2]
cor_val <- cor(-log10(resultGemma$p_lrt), shap_frame_df$abs_shap)

# Make the scatter plot with square aspect ratio
par(mar = c(5, 6, 4, 2) + 0.1, pty = "s")  # <- Added pty = "s"
plot(
  -log10(resultGemma$p_lrt), shap_frame_df$abs_shap,  
  xlab = TeX(r"($-log_{10}(p)$)"),  
  ylab = "mean |SHAP|",
  main = bquote(.(pheno_name)),
  cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.0,  # Adjust text sizes
  pch = 20, col = "black",  # Scatter point style
  xlim = range(-log10(resultGemma$p_lrt))
)

# Add the regression line to the plot
lines(x_vals, y_vals, col = "red", lwd = 2)

legend("topleft", 
       legend = c(
         paste("Corr:", round(cor_val, 3))
       ),
       cex = 1.2, 
       bty = "n", 
       text.col = c("black"),  
       seg.len = 2)