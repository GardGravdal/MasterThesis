############################################################################################################
                                        # SNP subset selection

# This script generates subsets of SNPs from either of two strategies:
# 1) 50k Random subset selection
# 2) 10k Most correlated SNPs
# The script is generalized for either of the 3 phenotypes mass, wing and tarsus (define global variable).
# Returns vector 'SNP_subset' with the names of the SNPs in the subset from chosen strategy.
# Script assumes 180k dataset is used (preparation for prediction of real phenotypes).
############################################################################################################


library(data.table)
library(arrow)
library(dplyr)

# Define trait to run globally
pheno <- "wing"
# 180K processed pheno:
pheno_path <- sprintf("C:/Users/gard_/Documents/MasterThesis/Code/Data/Processed/%sBV.feather", pheno)
data <- as.data.table(read_feather(pheno_path))

# Check that there are no duplicated ringnr's
print(any(duplicated(data$ringnr)))


# X is only SNPs
X <- data[, !c("ID", "ringnr", "mean_pheno", "FID", "MAT", "PAT", "SEX", "PHENOTYPE", "hatchisland"), with = FALSE] %>%
  as.data.frame()

# Some of the SNPs have NA-values. Set to 0
# This loop performs the equivalent of "X[is.na(X)] <- 0", but runs much faster for large matrices
for (j in seq_along(X)) {
  set(X, which(is.na(X[[j]])), j, 0)
}

# Change from float to integer for all columns except 'ringnr' (i.e., all SNPs)
cols_to_convert <- names(X)
X[cols_to_convert] <- lapply(X[cols_to_convert], as.integer)

# ID effect is processed phenotype y^*
y <- data %>%
  select(ID) %>%
  as.data.frame()

random_subset <- function(SNP_names){
  # Returns 50k random subset of SNP names in vector 'SNP_subset'
  SNP_subset <- sample(SNP_names, size = 50000, replace = FALSE)
  return(SNP_subset)
}

most_correlated <- function(y, X){
  #####################
  # y: ID effect of a phenotype (processed pheno).
  # X: Only SNP data.
  # Returns vector 'SNP_subset' containing 10k subset of SNPs that
  # are the most correlated with the response y.
  #####################

  # Compute (absolute value of) correlation between each SNP and y
  cor_values <- abs(cor(X, y$ID))

  # Rank SNPs by absolute correlation and select the top 10 000
  top_snp_indices <- order(cor_values, decreasing = TRUE)[1:10000]

  # Get the corresponding SNP names
  SNP_subset <- names(X)[top_snp_indices]
  return(SNP_subset)
}

test <- most_correlated(y,X)

length(test)

main <- function(S = "random"){
  #########################################
  # S: Strategy ['random', 'correlated'].
  # S = 'random': Returns a random subset of 50k SNPs.
  # S = 'correlated': Returns 10k most correlated SNPs.
  #########################################
  # Load initial data


  # Run subset selection based on strategy
  if (S == "random"){
    return(random_subset(names(X)))
  } else if (S == "correlated"){
    return(most_correlated(y, X))
  }
}

random_SNPs <- main("random")
corr_SNPs <- main("correlated")

# Save subsets as csv
random_file_name <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Subset/",pheno,"_50k_random.csv")
correlated_file_name <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Data/Subset/",pheno,"_10k_correlated.csv")
write.csv(random_SNPs, file = random_file_name, row.names = FALSE)
write.csv(corr_SNPs, file = correlated_file_name, row.names = FALSE)

