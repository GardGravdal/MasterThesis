################################################################################
#                         SHAP vs Effects vs GWAS comparison

# This script generates nice Manhattan plots of the variance components of the
# phenotypes from simulated effect sizes. The variance components used are 
# \beta_j^2 cov(m_j) for each SNP j. These are compared to GWAS and SHAP results
# respectively. This script generates the plots that compare all the three:
# effects, SHAP and GWAS together. 
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

# Set architecture and dataset
arch <- 10
dset <- "70k"

# Load data
SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_",dset,".feather"))
shap_vals <- read_feather(paste0('c:/Users/gard_/Documents//MasterThesis/Code/Results/',dset,'/SHAP/Sim/shap_',dset,'_arch_',arch,'.feather'))
effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",dset,"/Sim_effect_",dset,"_arch_",arch,".csv"))


# Load GWAS results
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Results/70k/GWAS")
gemma_path <- paste0("GWASresults",arch,".lmm.assoc.txt")
resultGemma <- read_table(gemma_path)

# Fix GEMMA multiline output
txt <- readLines(gemma_path)[-1]
idx <- seq(1, length(txt), by=2)
txt[idx] <- paste(txt[idx], txt[idx+1])
txt <- txt[-(idx+1)]
resultGemma <- read.table(text=txt, col.names=colnames(resultGemma))
resultGemma <- resultGemma[!(resultGemma$chr %in% c(0, 16, 30, 32)),]

# Load map file
map_path <- if (dset == "70k") {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
} else {
  "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
}
chro <- read.table(map_path, header = F)
names(chro) <- c("chr", "snpID", "value", "bp")
snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

# Map SNP names to basepair positions
get_bp <- function(snp) str_split(snp, "[_]")[[1]][1]
tmpdf <- data.frame(origsnp = effects$SNP, snpID = sapply(effects$SNP, get_bp), stringsAsFactors = FALSE)
tmpdfchro <- merge(tmpdf, chro[, c("chr", "snpID", "bp")], by = "snpID")
snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]

# Subset data
effects <- effects[effects$SNP %in% snp_to_keep, ]
shap_vals <- shap_vals[, names(shap_vals) %in% snp_to_keep, drop = FALSE]
SNP_subset <- SNPs[, c(snp_to_keep)]
shap_vals <- shap_vals[, effects$SNP, drop = FALSE]

# Load precomputed SNP variances
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Compute variance components
var_term <- effects$effect^2 * variances
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term)
shap_frame <- data.frame(origsnp = effects$SNP, abs_shap = as.numeric(shap_vals[1, ]))

# Merge with chromosome info
variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)
shap_frame_df <- merge(shap_frame, tmpdfchro, by = "origsnp", sort = FALSE)

# Get the top 3 SNPs with the highest effect size
top3 <- variance_term_df %>% arrange(desc(snp_effects)) %>% slice(1:3)
top3_snp_names <- top3$snpID

# Define consistent, contrasting colors for the top 3 SNPs
top3_colors <- c("blue", "green", "purple")
names(top3_colors) <- top3_snp_names

# Custom function to add thicker points for highlighted SNPs
add_thicker_points <- function(df, y_col, snp_col, top_snps, chr_col, bp_col, col_map, cex = 2.5) {
  df_subset <- df[df[[snp_col]] %in% top_snps, ]
  
  # Evaluate y_col (supports expressions like "-log10(p_lrt)")
  if (grepl("-log10", y_col)) {
    raw_col <- gsub("-log10\\((.*)\\)", "\\1", y_col)
    df_subset$y <- -log10(df_subset[[raw_col]])
  } else {
    df_subset$y <- df_subset[[y_col]]
  }
  
  chr_df <- df_subset[[chr_col]]
  bp_df <- df_subset[[bp_col]]
  
  # Compute Manhattan plot x-axis positions
  chr_levels <- sort(unique(df[[chr_col]]))
  chr_offsets <- cumsum(c(0, sapply(chr_levels, function(ch) max(df[df[[chr_col]] == ch, bp_col], na.rm = TRUE))))
  names(chr_offsets) <- chr_levels
  x_positions <- bp_df + chr_offsets[as.character(chr_df)]
  
  df_subset$color <- col_map[df_subset[[snp_col]]]
  points(x_positions, df_subset$y, pch = 21, bg = df_subset$color, col = "black", cex = cex, lwd = 1.2)
}

# Compute y-limits based on max values
shap_y_limit <- max(shap_frame_df$abs_shap) * 1.1
effect_y_limit <- max(variance_term_df$snp_effects) * 1.1
# Compute the Bonferroni threshold
bonferroni <- -log10(0.05/ nrow(resultGemma))

# Set layout for 3 plots in a row and adjust margins
par(mfrow = c(1, 3), mar = c(5, 6.5, 4, 1), oma = c(0, 5.5, 0, 0), xpd = NA)



# First plot: SHAP Manhattan plot
manhattan(shap_frame_df, chr = "chr", bp = "bp", p = "abs_shap", snp = "snpID", ylim = c(0, shap_y_limit),
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, ylab = "",
          highlight = top3_snp_names, highlight.col = top3_colors,
          logp = FALSE, genomewideline = FALSE, suggestiveline = FALSE, yaxt = "n")
axis(2, at = seq(0, shap_y_limit, by = shap_y_limit/3), labels = format(seq(0, shap_y_limit, by = shap_y_limit/3), digits = 2), cex.axis = 1.2)
title(ylab = "mean |SHAP|", line = 2, cex.lab = 2.3)
add_thicker_points(shap_frame_df, "abs_shap", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)

# Second plot: Simulated effect sizes
manhattan(variance_term_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, effect_y_limit), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE,
          suggestiveline = FALSE, yaxt = "n", highlight = top3_snp_names, highlight.col = top3_colors)
axis(2, at = seq(0, effect_y_limit, by = effect_y_limit/3), labels = format(seq(0, effect_y_limit, by = effect_y_limit/3), digits = 2), cex.axis = 1.2)
title(ylab = expression(u^2 ~ Var(M)), line = 2, cex.lab = 2.3)
add_thicker_points(variance_term_df, "snp_effects", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)

# Third plot: GWAS Manhattan plot
manhattan(resultGemma, chr = "chr", bp = "ps", p = "p_lrt", snp = "rs", ylim = c(0, 7),
          col = c("red1", "red4"), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5,
          ylab = "", genomewideline = FALSE, suggestiveline = FALSE,
          highlight = top3_snp_names, highlight.col = top3_colors)
abline(h = bonferroni, col = "red", lty = 2)
title(ylab = TeX(r"($-log_{10}(p)$)"), line = 2, cex.lab = 2.3)
add_thicker_points(resultGemma, "-log10(p_lrt)", "rs", top3_snp_names, chr_col = "chr", bp_col = "ps", col_map = top3_colors)