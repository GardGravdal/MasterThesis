################################################################################
#                               -log(u) plotting

# This is a script for generating plot of unif(p) on -log() scale, to be compared
# with GWAS. The purpose is to have an illustration and basis of comparison for
# false outliers in multiple hypothesis testing in a setting where H_0 is not
# true.

# We create dataframe [chr, bp, unif(p), snpID] and illustrate with manhattan plot.

################################################################################
library(dplyr)
library(qqman)
library(latex2exp)

#set.seed(3)
# Map for bp and chr
map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
chro <- read.table(map_path, header = F)
names(chro) <- c("chr", "snpID", "value", "bp")

# Number of covariates
p <- dim(chro)[1]

# Create p uniform variables
#set.seed(123)
u <- runif(p)

# Define dataframe [chr, bp, unif(p), snpID]
u_df <- data.frame(chr = chro$chr, bp = chro$bp, unif = u, snpID = chro$snpID, row.names = NULL)

# Exclude some chromosomes (correspond with other manhattan plots)
snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]
snp_to_keep <- chro[!(chro$snpID %in% snp_to_discard), "snpID"]
u_df <- u_df[u_df$snpID %in% snp_to_keep, ]


# Plot
bonferroni <- -log10(0.05/ nrow(u_df))


manhattan(u_df,chr="chr",bp="bp",p="unif",snp = "snpID", ylim = c(0, 9), col = c("gray10", "gray50"),
          cex.axis = 1.5, cex.lab=2.3, cex.main=2.5, ylab="",
          logp = TRUE, genomewideline = bonferroni, suggestiveline = FALSE, yaxt = "n")
axis(2, at = seq(0, 9, by = 1), labels = format(seq(0, 9, by = 1),
                                                          digits = 2), cex.axis = 1.2)
title(ylab = TeX(r"($-log_{10}(unif[0,1])$)"), line = 2, cex.lab = 1.8)

################################################################################
# Three realizations in one plot

# Read the map file
map_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Data/Helgeland_01_2018.map"
chro <- read.table(map_path, header = FALSE)
names(chro) <- c("chr", "snpID", "value", "bp")

# Number of covariates
p <- nrow(chro)

# Define chromosome exclusions
snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]
snp_to_keep <- chro[!(chro$snpID %in% snp_to_discard), "snpID"]

# Set up 1 row, 3 columns for plots
par(mfrow = c(1, 3))

# Plot 3 realizations
for (i in 1:3) {
  set.seed(123 + i)  # Different seeds for different realizations
  
  # Create uniform variables
  u <- runif(p)
  
  # Create dataframe [chr, bp, unif(p), snpID]
  u_df <- data.frame(chr = chro$chr, bp = chro$bp, unif = u, snpID = chro$snpID)
  
  # Filter unwanted SNPs
  u_df <- u_df[u_df$snpID %in% snp_to_keep, ]
  
  # Bonferroni threshold
  bonferroni <- -log10(0.05 / nrow(u_df))
  
  # Plot Manhattan
  manhattan(u_df, chr = "chr", bp = "bp", p = "unif", snp = "snpID", 
            ylim = c(0, 9), col = c("gray10", "gray50"),
            cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, ylab = "",
            logp = TRUE, genomewideline = min(bonferroni, 9), 
            suggestiveline = FALSE, yaxt = "n")
  
  axis(2, at = seq(0, 9, by = 1), labels = format(seq(0, 9, by = 1), digits = 2), cex.axis = 1.2)
  
  title(ylab = TeX(r"($-log_{10}(unif[0,1])$)"), line = 2, cex.lab = 1.8)
}


