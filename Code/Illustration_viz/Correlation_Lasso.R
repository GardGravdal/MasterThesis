################################################################################
#                        Corr LASSO vs. variance term
# Script for generating correlation plot between phenotypic variance term and
# LASSO coefficient (|β_L|) importance measure for simulated traits.
# Generalized to work for all architecture types (1–10).
################################################################################

library(qqman)
library(tidyr)
library(dplyr)
library(stringr)
library(grid)
library(gridGraphics)
library(ggplot2)
library(arrow)
library(readr)
library(data.table)
library(smplot2)
library(latex2exp)

# Set architecture globally [1,2,...,10], where 1 is polygenic and 10 oligogenic
arch <- 10
pi_0_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
pi_0 <- pi_0_vals[arch]

# Choose dataset: "70k" or "180k"
dset <- "70k"

# Load SNPs and effect data
SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_", dset, ".feather"))
effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/", dset, "/Sim_effect_", dset, "_arch_", arch, ".csv"))

# Map file for SNP locations
map_path <- if (dset == "70k") {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
} else {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
}
chro <- read.table(map_path, header = FALSE)
names(chro) <- c("chr", "snpID", "value", "bp")

# Filter out unwanted chromosomes
snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

# Extract SNP mapping
SNP_cols <- effects[["SNP"]]
tmpdf <- data.frame(origsnp = SNP_cols, snpID = NA)

# Function to get base SNP name
get_bp <- function(snp) {
  snp <- str_split(snp, "[_]")[[1]][1]
  return(snp[[1]])
}
tmpdf$snpID <- lapply(tmpdf$origsnp, get_bp)

# Clean to proper data frame
tmpdf2 <- as.data.frame(lapply(tmpdf, unlist))

# Merge with chromosome info
tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")

# Filter SNPs to keep
snp_to_keep <- tmpdf2[!(tmpdf2$snpID %in% snp_to_discard), "origsnp"]

# Subset data
effects <- effects[effects$SNP %in% snp_to_keep, ]
SNP_subset <- SNPs[, snp_to_keep]

# Ensure correct order
effects <- effects[match(snp_to_keep, effects$SNP), ]

# Read SNP variances
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Calculate variance term (β² × Var(X))
var_term <- effects$effect^2 * variances

# Create variance term data frame
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term, row.names = NULL)


# --- NEW: Load LASSO Coefficients ---
lasso_path <- paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/", dset, "/Lasso/lasso_coefs_arch", arch, ".csv")
lasso_df <- read.csv(lasso_path)
lasso_df$SNP <- as.character(lasso_df$SNP)
lasso_df <- lasso_df[lasso_df$SNP %in% snp_to_keep, ]

# Create LASSO frame and merge with chr/bp info
lasso_frame <- data.frame(origsnp = lasso_df$SNP, var_lasso = lasso_df$beta_L^2 * variances, row.names = NULL)
lasso_frame_df <- merge(lasso_frame, tmpdfchro, by = "origsnp", sort = FALSE)

# Merge with variance term data for modeling
variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)

# Fit linear model
lm_model <- lm(variance_term_df$snp_effects ~ lasso_frame_df$var_lasso)

# X and Y values for line
x_vals <- seq(0, max(lasso_frame_df$var_lasso), length.out = length(variance_term_df$snp_effects))
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

# Correlation and slope
slope <- coef(lm_model)[2]
cor_val <- cor(variance_term_df$snp_effects, lasso_frame_df$var_lasso)

# --- Plotting ---
par(mar = c(5, 7, 4, 2) + 0.1, pty = "s", mgp = c(3.5, 1.5, 0))

plot(
  lasso_frame_df$var_lasso, variance_term_df$snp_effects,  
  xlab = "",
  ylab = expression(u^2 ~ Var(M)),
  main = bquote(pi[0] == .(pi_0)),
  cex.axis = 1.3, cex.lab = 2.0, cex.main = 2.0,
  pch = 20, col = "black",
  xlim = range(lasso_frame_df$var_lasso)
)
mtext(expression(hat(u)*phantom()[L]^2 ~ Var(M)), side = 1, line = 4.2, cex = 2.0)  
lines(x_vals, y_vals, col = "red", lwd = 2)

legend("topleft", 
       legend = c(
         paste("Corr:", round(cor_val, 3))
       ),
       cex = 1.2,
       bty = "n",
       text.col = c("black"),
       seg.len = 2)