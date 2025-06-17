################################################################################
# Script for generating G-matrix according to van Raden method 1
################################################################################

library(data.table)
library(arrow)
library(dplyr)
library(Matrix)

# Load SNP_matrix as a regular data frame
X <- read_feather("C:/Users/gard_/Documents/MasterThesis/Code/Data/SNP1_180k.feather") %>%
  as.data.frame()


X <- X %>% rename(ringnr = IID)


# Excluding 'ringnr' column
SNPs <- X[, !colnames(X) %in% "ringnr"]

# M_star is SNP matrix encoded {-1,0,1} instead
M_star <- SNPs - 1

M_star <- as.matrix(M_star)

# Compute allele frequencies (column means divided by 2)
p <- colMeans(SNPs) / 2

# Compute P matrix: P_ij = 2(p_j - 0.5)
P <- matrix(2 * (p - 0.5), nrow = nrow(SNPs), ncol = ncol(SNPs), byrow = TRUE)

#P <- as.data.frame(P)

# C as defined by van Raden method 1

C = M_star - P

denominator <- 2 * sum(p * (1 - p))

# Compute G matrix
G <- (C %*% t(C)) / denominator

# Adding +0.01 on diagonal of G
diag(G) <- diag(G) + 0.01

# Convert sparse matrix G to a three-column data frame (as Henrik did in his script)
#G_sparse <- as.data.frame(summary(G))
G_sparse <- as(G, "TsparseMatrix")  # Convert to triplet format

# Create a data frame with row index (i), column index (j), and values (x)
G_sparse <- data.frame(V1 = G_sparse@i + 1,  # +1 to adjust from 0-based to 1-based indexing
                        V2 = G_sparse@j + 1,
                        V3 = G_sparse@x)


write.csv(G_sparse, "C:/Users/gard_/Documents/MasterThesis/Code/Data/G_sparse.csv", row.names = FALSE)
