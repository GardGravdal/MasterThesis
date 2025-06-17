################################################################################
#                         GWAS vs Ridge coefs comparison

# This script generates nice Manhattan plots of comparison of GWAS p-values and Ridge mean coefs.
# The variance components used are  \beta_j^2 var(X_j) for each SNP j. These are compared to GWAS and SHAP results
# respectively. This script generates the plots that compare effects and GWAS results for real traits.
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

# Set phenotype globally
pheno <- "mass"
dset <- "70k"

# 180k:
#SNP_data <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/",pheno,"BV.feather"))
# 70k:
SNP_data <- as.data.table(read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/",pheno,"BV_70k.feather")))

SNPs <- SNP_data[, !c("ringnr","ID", "mean_pheno", "IID", "MAT", "PAT", "SEX", "PHENOTYPE", "hatchisland"), with = FALSE] %>%
  as.data.frame()
coefs <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/",dset,"/Ridge/","ridge_coeffs_",pheno,".csv"))

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
if (dset == "70k"){
  map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
} else {
  map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
}

chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")


snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

SNP_cols <- coefs["SNP"]
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
coefs <- coefs[coefs$SNP %in% snp_to_keep, ]
SNP_subset <- SNPs[,c(snp_to_keep)]

# Read variances for 70k SNP data
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Order to correspond to order in variances (in terms of SNP)
coefs <- coefs[match(names(variances), coefs$SNP), ]
SNP_subset <- SNP_subset[match(names(variances), names(SNP_subset))]

# Check that they follow the same ordering before plotting
head(resultGemma)
head(coefs$SNP)
names(SNP_subset)[1:6]
head(variances)


# Create coefficient term to be compared with GWAS
var_term <- coefs$beta_R^2 * variances
coefs <- data.frame(origsnp = coefs$SNP, snp_coefs = var_term, row.names = NULL)

coefs_df <- merge(coefs, tmpdfchro, by = "origsnp", sort = FALSE)

# Check that order is still retained
head(coefs_df)
head(resultGemma)

# Sometimes we have NAs in p-values:
sum(is.na(resultGemma$p_lrt))
# Remove SNPs with NA p-value:
resultGemma <- resultGemma[!is.na(resultGemma$p_lrt), ]

# Remove excessive SNPs in coefs_df
coefs_df <- coefs_df %>%
  filter(snpID %in% resultGemma$rs)

# Check that dimensions are the same
dim(coefs_df)
dim(resultGemma)


# NOW WE CAN PLOT
#setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures/Manhattan")
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures")

# Compute the Bonferroni threshold
bonferroni <- -log10(0.05/ nrow(resultGemma))

################################################################################
# Plot which highlights SNP with highest coef size

# Find the SNP with the highest coef size
max_coef_snp <- coefs_df[which.max(coefs_df$snp_coefs), ]
max_coef_snp_name <- max_coef_snp$snpID
# Set the y-axis of coefs plot relative to max coef size
y_limit <- max_coef_snp$snp_coefs + 0.0000005

# Create Manhattan plots in a 1x2 layout
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)  # Adjust margins to fit both plots

# First plot: Manhattan plot for variance_term_df
manhattan(coefs_df, chr = "chr", bp = "bp", p = "snp_coefs", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE, 
          suggestiveline = FALSE, yaxt = "n", highlight = max_coef_snp_name)
axis(2, at = seq(0, y_limit, by = 0.005), labels = format(seq(0, y_limit, by = 0.005),
                                                       digits = 2), cex.axis = 1.2)
# Add a title for the y-axis
title(ylab = expression(beta^2 ~ Var(X)), line = 2, cex.lab = 2.3)

# Second plot: Manhattan plot for resultGemma
manhattan(resultGemma, chr = "chr", bp = "ps", p = "p_lrt", snp = "rs", ylim = c(0, 7), 
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, 
          ylab = "", genomewideline = FALSE, suggestiveline = FALSE, 
          highlight = max_coef_snp_name)

# Add the Bonferroni line and a title for the y-axis for resultGemma
abline(h = bonferroni, col = "red", lty = 2)
title(ylab = TeX(r"($-log_{10}(p)$)"), line = 2, cex.lab = 2.3)


################# Variance vs p-vals #########################

# Fit a linear model
lm_model <- lm(coefs_df$snp_coefs ~ I(-log(resultGemma$p_lrt)))

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, length(-log(resultGemma$p_lrt)), length.out = length(resultGemma))


# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

slope <- coef(lm_model)[2]
slope <- -log(slope)
lm_model


# Make the scatter plot
par(mfrow = c(1,1))
par(mar = c(5, 7.5, 4, 2) + 0.1) 
plot(
  -log(resultGemma$p_lrt), coefs_df$snp_coefs,  
  xlab = "",  # Styled x-axis label
  ylab = expression(hat(beta)[R]^2 ~ Var(X)),  # Y-axis label
  main = paste("Correlation ridge coefficient term and GWAS p-values for", pheno),
  cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.0,  # Adjust text sizes
  pch = 20, col = "black",  # Scatter point style
  xlim = range(-log(resultGemma$p_lrt))
)
# Add the regression line to the plot
lines(x_vals, y_vals, col = "red", lwd = 2)
mtext(TeX(r"($-log_{10}(p)$)"), side = 1, line = 4, cex = 2.3)

correlation <- cor(coefs_df$snp_coefs, -log(resultGemma$p_lrt))
# Add a legend with the regression slope
legend("topleft", 
       legend = paste("Corr: ", round(correlation, 3)),  # Display the slope with 3 decimal places
       col = "red", 
       lwd = 2, 
       cex = 1.2, 
       bty = "n")  # bty = "n" removes the box around the legend
#cor(coefs_df$snp_coefs, -log(resultGemma$p_lrt))
