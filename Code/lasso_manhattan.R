################################################################################
# Variance vs Lasso/Ridge coefs correlation plot

# This script generates 3 types of plots:
# 2) Manhattan plot of beta^2 * Var(M) for true effects and estimated effects.
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

# ========== CONFIGURATION ========== #
# Set architecture type
arch <- 1  # architecture [1:10]
dset <- "70k"  # "70k" or "180k"
method <- "Lasso"  # or "Ridge"

pi_0_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
pi_0 <- pi_0_vals[arch]

# ========== LOAD DATA ========== #
SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_", dset, ".feather"))
coefs <- read_csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/", dset, "/", method, "/", tolower(method), "_coefs_arch", arch, ".csv"))
effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/", dset, "/Sim_effect_", dset, "_arch_", arch, ".csv"))

# ========== MAP DATA ========== #
map_path <- if (dset == "70k") {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
} else {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
}

chro <- read.table(map_path, header = FALSE)
names(chro) <- c("chr", "snpID", "value", "bp")

snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

SNP_cols <- effects["SNP"]
tmpdf <- data.frame(origsnp = SNP_cols, snpID = NA)
names(tmpdf) <- c("origsnp", "snpID")

get_bp <- function(snp) str_split(snp, "[_]")[[1]][1]
tmpdf$snpID <- lapply(tmpdf$origsnp, get_bp)
tmpdf2 <- as.data.frame(lapply(tmpdf, unlist))

tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")
snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]

effects <- effects[effects$SNP %in% snp_to_keep, ]
coefs <- coefs[coefs$SNP %in% snp_to_keep, ]
SNP_subset <- SNPs[, c(snp_to_keep)]

variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# ========== COMPUTE VARIANCE TERMS ========== #
var_term <- effects$effect^2 * variances
beta_col <- if (method == "Ridge") "beta_R" else if (method == "Lasso") "beta_L" else stop("Invalid method")
var_term2 <- coefs[[beta_col]]^2 * variances

variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term)
coefs_frame <- data.frame(origsnp = effects$SNP, snp_effects = var_term2)

variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)
coefs_frame_df <- merge(coefs_frame, tmpdfchro, by = "origsnp", sort = FALSE)


# ========== PLOT 2: MANHATTAN PLOTS ========== #
top3 <- variance_term_df %>% arrange(desc(snp_effects)) %>% slice(1:3)
top3_snp_names <- top3$snpID
top3_colors <- c("blue", "green", "purple")
names(top3_colors) <- top3_snp_names

add_thicker_points <- function(df, y_col, snp_col, top_snps, chr_col, bp_col, col_map, cex = 2.5) {
  df_subset <- df[df[[snp_col]] %in% top_snps, ]
  if (grepl("-log10", y_col)) {
    raw_col <- gsub("-log10\\((.*)\\)", "\\1", y_col)
    df_subset$y <- -log10(df_subset[[raw_col]])
  } else {
    df_subset$y <- df_subset[[y_col]]
  }
  chr_df <- df_subset[[chr_col]]
  bp_df <- df_subset[[bp_col]]
  chr_levels <- sort(unique(df[[chr_col]]))
  chr_offsets <- cumsum(c(0, sapply(chr_levels, function(ch) max(df[df[[chr_col]] == ch, bp_col], na.rm = TRUE))))
  names(chr_offsets) <- chr_levels
  x_positions <- bp_df + chr_offsets[as.character(chr_df)]
  df_subset$color <- col_map[df_subset[[snp_col]]]
  points(x_positions, df_subset$y, pch = 21, bg = df_subset$color, col = "black", cex = cex, lwd = 1.2)
}

y_limit1 <- max(variance_term_df$snp_effects) * 1.1
y_limit2 <- max(coefs_frame_df$snp_effects) * 1.1

par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)

manhattan(variance_term_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit1), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE,
          suggestiveline = FALSE, yaxt = "n", highlight = top3_snp_names, highlight.col = top3_colors)
axis(2, at = seq(0, y_limit1, by = y_limit1/3), labels = format(seq(0, y_limit1, by = y_limit1/3), digits = 2), cex.axis = 1.2)
title(ylab = expression(u^2 ~ Var(M)), line = 2, cex.lab = 2.3)
add_thicker_points(variance_term_df, "snp_effects", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)

manhattan(coefs_frame_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit2), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE,
          suggestiveline = FALSE, yaxt = "n", highlight = top3_snp_names, highlight.col = top3_colors)
axis(2, at = seq(0, y_limit2, by = y_limit2/3), labels = format(seq(0, y_limit2, by = y_limit2/3), digits = 2), cex.axis = 1.2)
title(
  ylab = if (method == "Ridge") expression(hat(u)[R]^2 ~ Var(M)) else expression(hat(u)[L]^2 ~ Var(M)),
  line = 2, cex.lab = 2.3
)
add_thicker_points(coefs_frame_df, "snp_effects", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)
