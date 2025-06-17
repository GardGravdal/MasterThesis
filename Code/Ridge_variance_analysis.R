################################################################################
                  # Variance vs Ridge coefs correlation plot

# This script generates 3 types of plot:
# 1) Correlation plot for true variance vs. ridge estimated variance.
# 2) Manhattan beta^2 var(M) for true effects and ridge effects vs. chromosome 
# position.
# 3) Histogram plot for true effects and ridge effects.
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



# 1):
# Set architecture globally [1,2,...,10], where 1 is polygenic and 10 oliogenic
arch <- 10
pi_0_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
pi_0 <- pi_0_vals[arch]
# Choose dataset 180k or 70k
dset <- "70k"

SNPs <- read_feather(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_",dset,".feather"))
coefs <- read_csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/70k/Ridge/ridge_coefs_arch",arch,".csv"))
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
coefs <- coefs[coefs$SNP %in% snp_to_keep, ]
SNP_subset <- SNPs[,c(snp_to_keep)]


# Compute variances for each column (Should not be necessary after first time - see read line next)
#variances <- apply(SNP_subset, 2, var)
#saveRDS(variances, "C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Read variances for 70k SNP data
variances <- readRDS("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP_70k_variances.rds")

# Check that they follow the same ordering before plotting
head(variances)
head(effects$SNP)
names(SNP_subset)[1:6]
head(coefs$SNP)

# Create the variance terms to be compared to SHAP/GWAS values
var_term <- effects$effect^2 * variances
var_term2 <- coefs$beta_R^2 * variances

range(var_term)
range(var_term2)

# Create variance term data frame and shap data frame
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term, row.names = NULL)
coefs_frame <- data.frame(origsnp = effects$SNP, snp_effects = var_term2, row.names = NULL)

# Create the data frames with bp and chr information
variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)
coefs_frame_df <- merge(coefs_frame, tmpdfchro, by = "origsnp", sort = FALSE)

# Check that order is still retained
head(variance_term_df)
head(coefs_frame_df)


# NOW WE CAN PLOT
setwd("C:/Users/gard_/Documents/MasterThesis/Code/Figures")



# 1):
################################################################################

# Fit a linear model 
lm_model <- lm(variance_term_df$snp_effects ~ coefs_frame_df$snp_effects)

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, length(variance_term_df), length.out = length(variance_term_df))

# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

slope <- coef(lm_model)[2]
cor_val <- cor(variance_term_df$snp_effects, coefs_frame_df$snp_effects)

# Make the scatter plot
par(mar = c(5, 7, 4, 2) + 0.1, pty = "s", mgp = c(3.5, 1.5, 0))

plot(
  coefs_frame_df$snp_effects, variance_term_df$snp_effects,  
  xlab = "",  # Suppress default x-axis label
  ylab = expression(u^2 ~ Var(M)),  # Y-axis label
  main = bquote(pi[0] == .(pi_0)),
  cex.axis = 1.3, cex.lab = 2.0, cex.main = 2.0,
  pch = 20, col = "black",
  xlim = range(coefs_frame_df$snp_effects)
)

# Add x-axis label manually, moved lower with line = 3.5 or 4
mtext(expression(hat(u)*phantom()[R]^2 ~ Var(M)), side = 1, line = 4.2, cex = 2.0)

# Add the regression line
lines(x_vals, y_vals, col = "red", lwd = 2)

# Add the legend
legend("topleft", 
       legend = c(
         paste("Corr:", round(cor_val, 3))
       ),
       cex = 1.2,
       bty = "n",
       text.col = c("black"),
       seg.len = 2)



################################################################################
# 2) Plot which highlights SNP with highest effect size

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

# Set y-axis limits for both plots
y_limit1 <- max(variance_term_df$snp_effects) * 1.1
y_limit2 <- max(coefs_frame_df$snp_effects) * 1.1

# Set layout and margins for 2 plots
par(mfrow = c(1, 2), mar = c(5, 6, 4, 1) + 0.1)

# First plot: Manhattan for variance_term_df
manhattan(variance_term_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit1), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE,
          suggestiveline = FALSE, yaxt = "n", highlight = top3_snp_names, highlight.col = top3_colors)
axis(2, at = seq(0, y_limit1, by = y_limit1/3), labels = format(seq(0, y_limit1, by = y_limit1/3), digits = 2), cex.axis = 1.2)
title(ylab = expression(u^2 ~ Var(M)), line = 2, cex.lab = 2.3)
add_thicker_points(variance_term_df, "snp_effects", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)

# Second plot: Manhattan for coefs_frame_df
manhattan(coefs_frame_df, chr = "chr", bp = "bp", p = "snp_effects", snp = "snpID",
          logp = FALSE, ylab = "", ylim = c(0, y_limit2), col = c("red1", "red4"),
          cex.axis = 1.2, cex.lab = 2.3, cex.main = 2.5, genomewideline = FALSE,
          suggestiveline = FALSE, yaxt = "n", highlight = top3_snp_names, highlight.col = top3_colors)
axis(2, at = seq(0, y_limit2, by = y_limit2/3), labels = format(seq(0, y_limit2, by = y_limit2/3), digits = 2), cex.axis = 1.2)
title(ylab = expression(hat(u)[R]^2 ~ Var(M)), line = 2, cex.lab = 2.3)
add_thicker_points(coefs_frame_df, "snp_effects", "snpID", top3_snp_names, chr_col = "chr", bp_col = "bp", col_map = top3_colors)


################################################################################
# 3) Plot histofram of effect sizes (Ridge and true):

pi_0_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
par(mar = c(5, 5, 4, 2))

# Plot histogram ridge effects all architectures
for (arch in 1:10){
  pi_0 <- pi_0_vals[arch]
  coefs <- read_csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/70k/Ridge/ridge_coefs_arch",arch,".csv"))
  hist(coefs$beta_R, main = bquote(hat(beta)[R] ~ "for " ~ pi[0] == .(pi_0)),
       xlab = "beta_R", col = "skyblue", border = "white",
       cex.main = 2,    # main title size
       cex.lab = 2,   # axis label size
       cex.axis = 1.3)
}

# Plot histogram true effects all architectures
for (arch in 1:10){
  pi_0 <- pi_0_vals[arch]
  effects <- read.csv(paste0("c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/",dset,"/Sim_effect_",dset,"_arch_",arch,".csv"))
  hist(effects$effect, main = bquote(beta ~ "for" ~ pi[0] == .(pi_0)),
       xlab = "beta", col = "skyblue", border = "white",
       cex.main = 2,    # main title size
       cex.lab = 2,   # axis label size
       cex.axis = 1.3)
}

