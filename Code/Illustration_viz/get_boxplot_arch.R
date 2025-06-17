################################################################################
# Script for plotting simulations in different architectures for:
# 1) Ridge regression with SNPs as input
# 2) Ridge regression with PCs as input
# 3) XGBoost
# 4) Animal model

library(ggplot2)

# Architecture vector in terms of pi_0
pi_values <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
# Define correct π₀ labels
pi_labels <- sapply(pi_values, function(x) bquote(pi[0] == .(x)))

# Read Ridge regression data
ridge_files <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/Ridge/corr_Ridge_arch_", 1:10, "_std.csv")
ridge_data <- lapply(ridge_files, read.csv)
ridge_values <- unlist(lapply(ridge_data, function(df) df$corr))

# Read RidgePC regression data
ridgePC_files <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RidgePC/corr_RidgePC_arch_", 1:10, ".csv")
ridgePC_data <- lapply(ridgePC_files, read.csv)
ridgePC_values <- unlist(lapply(ridgePC_data, function(df) df$corr))

# Read XGBoost regression data
xgb_files <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/XGB/corr_xgb_arch_", 1:10, ".csv")
xgb_data <- lapply(xgb_files, read.csv)
xgb_values <- unlist(lapply(xgb_data, function(df) df$corr))

# Read Animal Model (AM) regression data
AM_files <- paste0("C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/AnimalModel/corr_BV_arch_", 1:10, ".csv")
AM_data <- lapply(AM_files, read.csv)
AM_values <- unlist(lapply(AM_data, function(df) df$corr))


# Create data frame for Ridge, RidgePC, and XGBoost (Architectures 1-10)
data <- data.frame(
  Architecture = factor(rep(paste0("Arch ", 1:10), each = length(ridge_values) / 10), 
                        levels = paste0("Arch ", 1:10)),
  Values = c(ridge_values, ridgePC_values, xgb_values),
  Method = c(rep("RR", length(ridge_values)), 
             rep("PCRR", length(ridgePC_values)), 
             rep("XGB", length(xgb_values)))
)

# Create data frame for Animal Model (Only architectures 1-3)
AM_data_frame <- data.frame(
  Architecture = factor(rep(paste0("Arch ", 1:10), each = length(AM_values) / 10), 
                        levels = paste0("Arch ", 1:10)),
  Values = AM_values,
  Method = "Animal Model"
)

# Combine all models into a single data frame
data <- rbind(data, AM_data_frame)

# Define color palette for methods
method_colors <- c("RR" = "steelblue", "PCRR" = "darkorange", "XGB" = "forestgreen", "Animal Model" = "purple")

# Create grouped box plot
ggplot(data, aes(x = Architecture, y = Values, fill = Method)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 19, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = method_colors) +
  scale_x_discrete(labels = pi_labels) +  # Ensure proper π₀ formatting
  labs(
    title = "Prediction accuracy across architectures",
    y = expression(Corr(bar(y), hat(y)^"*")),
    x = expression(pi[0])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 21, hjust = 0.5),
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 21)
  )

