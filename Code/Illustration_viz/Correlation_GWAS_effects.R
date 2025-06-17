################################################################################
#                             Corr GWAS vs. variance term
# Script for generating correlation plot betweem phenotypic variance term and
# p-values from GWAS for simulated traits. Generalized to work for all 
# architecture types (1-10).
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

# Set architecture globally [1,2,...,10], where 1 is polygenic and 10 oliogenic
arch <- 5
pi_0_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
pi_0 <- pi_0_vals[arch]

# Choose dataset 180k or 70k
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

# Check that they follow the same ordering before plotting
head(variances)
head(effects$SNP)
names(SNP_subset)[1:6]
head(resultGemma)


# Create the variance terms to be compared to SHAP/GWAS values
var_term <- effects$effect^2 * variances

# Create variance term data frame
variance_term <- data.frame(origsnp = effects$SNP, snp_effects = var_term, row.names = NULL)

# Create the data frames with bp and chr information
variance_term_df <- merge(variance_term, tmpdfchro, by = "origsnp", sort = FALSE)

# Remove excessive SNPs in variance_term
variance_term_df <- variance_term_df %>%
  filter(snpID %in% resultGemma$rs)

# Check that they are not the same dimension
dim(resultGemma)
dim(variance_term_df)


# Fit a linear model 
lm_model <- lm(variance_term_df$snp_effects ~ I(-log10(resultGemma$p_lrt)))

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, max(-log10(resultGemma$p_lrt)) + 100, length.out = length(variance_term_df$snpID))

# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

slope <- coef(lm_model)[2]
cor_val <- cor(variance_term_df$snp_effects, -log10(resultGemma$p_lrt))


# Set plot parameters with square plotting region
par(mar = c(5, 6, 4, 2) + 0.1, pty = "s")

# Make the scatter plot with square aspect ratio
plot(
  -log10(resultGemma$p_lrt), variance_term_df$snp_effects,  
  xlab = TeX(r"($-log_{10}(p)$)"),  # Styled x-axis label
  ylab = expression(u^2 ~ Var(M)),  # Y-axis label
  main = bquote(pi[0] == .(pi_0)),  # Simplified title
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