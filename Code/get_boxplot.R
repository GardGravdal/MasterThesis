# Script for generating nice multiple box-plots
library(ggplot2)

# Define path
real_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno/"

# Read data from CSV files
massBV <- read.csv(paste0(real_path, "corr_BV_mass.csv"))
wingBV <- read.csv(paste0(real_path, "corr_BV_wing.csv"))
xgb_mass <- read.csv(paste0(real_path, "corr_XGB_mass.csv"))
xgb_wing <- read.csv(paste0(real_path, "corr_XGB_wing.csv"))
RidgePC_mass <- read.csv(paste0(real_path, "corr_mass_RidgePC.csv"))
RidgePC_wing <- read.csv(paste0(real_path, "corr_wing_RidgePC.csv"))
Ridge_mass <- read.csv(paste0(real_path, "corr_mass_Ridge.csv"))
Ridge_wing <- read.csv(paste0(real_path, "corr_wing_Ridge.csv"))

# Include phenotype 'Tarsus'
tarsusBV <- read.csv(paste0(real_path, "corr_BV_tarsus.csv"))
Ridge_tarsus <- read.csv(paste0(real_path, "corr_tarsus_Ridge.csv"))
RidgePC_tarsus <- read.csv(paste0(real_path, "corr_tarsus_RidgePC.csv"))
xgb_tarsus <- read.csv(paste0(real_path, "corr_XGB_tarsus.csv"))

# Create a data frame including Ridge and RidgePC
data <- data.frame(
  Category = factor(rep(c("GEBV", "XGB", "Ridge", "RidgePC", 
                          "GEBV", "XGB", "Ridge", "RidgePC", 
                          "GEBV", "XGB", "Ridge", "RidgePC"), each = 10),
                    levels = c("GEBV", "XGB", "Ridge", "RidgePC")),
  Trait = factor(rep(c("Mass", "Mass", "Mass", "Mass", 
                       "Wing", "Wing", "Wing", "Wing", 
                       "Tarsus", "Tarsus", "Tarsus", "Tarsus"), each = 10),
                 levels = c("Mass", "Wing", "Tarsus")),
  Values = c(massBV$corr,
             xgb_mass$corr,
             Ridge_mass$corr,
             RidgePC_mass$corr,
             wingBV$corr,
             xgb_wing$corr,
             Ridge_wing$corr,
             RidgePC_wing$corr,
             tarsusBV$corr,
             xgb_tarsus$corr,
             Ridge_tarsus$corr,
             RidgePC_tarsus$corr)
)

# Define color mapping for categories
category_colors <- c("GEBV" = "aquamarine", 
                     "XGB" = "darkcyan", 
                     "Ridge" = "steelblue",
                     "RidgePC" = "orange")

# Create the box plot
ggplot(data, aes(x = Category, y = Values, fill = Category)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 19, outlier.size = 2) +
  scale_fill_manual(values = category_colors) +
  labs(
    title = "Prediction of genetic contribution",
    y = expression(Corr(bar(y), hat(y)^"*")),
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15, hjust = 0.5),   # Title size
    strip.text = element_text(size = 15, face = "bold"), # Subtitle text
    axis.title.x = element_text(size = 15),             # X-axis label size
    axis.title.y = element_text(size = 15),             # Y-axis label size
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 15)               # Y-axis tick label size
  ) +
  geom_hline(yintercept = 0.33, col = "blue", linetype = "dashed") +  # Dashed reference line
  facet_grid(~Trait, scales = "free_x", space = "free_x")  # Group by

