library(ggplot2)

# Define path
real_path <- "C:/Users/gard_/Documents/MasterThesis/Code/Results/180k/RealPheno/"

# Read data from CSV files
massBV <- read.csv(paste0(real_path, "corr_BV_mass.csv"))
wingBV <- read.csv(paste0(real_path, "corr_BV_wing.csv"))
PCRR_mass <- read.csv(paste0(real_path, "corr_mass_RidgePC.csv"))
PCRR_wing <- read.csv(paste0(real_path, "corr_wing_RidgePC.csv"))
RR_mass <- read.csv(paste0(real_path, "corr_mass_Ridge.csv"))
RR_wing <- read.csv(paste0(real_path, "corr_wing_Ridge.csv"))
xgb_mass <- read.csv(paste0(real_path, "corr_XGB_mass.csv"))
xgb_wing <- read.csv(paste0(real_path, "corr_XGB_wing.csv"))

# Tarsus models
tarsusBV <- read.csv(paste0(real_path, "corr_BV_tarsus.csv"))
RR_tarsus <- read.csv(paste0(real_path, "corr_tarsus_Ridge.csv"))
PCRR_tarsus <- read.csv(paste0(real_path, "corr_tarsus_RidgePC.csv"))
xgb_tarsus <- read.csv(paste0(real_path, "corr_XGB_tarsus.csv"))
xgb_tarsus_10k <- read.csv(paste0(real_path, "corr_XGB_tarsus_10k.csv"))
xgb_tarsus_50k <- read.csv(paste0(real_path, "corr_XGB_tarsus_50k.csv"))

xgb_mass_10k <- read.csv(paste0(real_path, "corr_XGB_mass_10k.csv"))
xgb_mass_50k <- read.csv(paste0(real_path, "corr_XGB_mass_50k.csv"))
xgb_wing_10k <- read.csv(paste0(real_path, "corr_XGB_wing_10k.csv"))
xgb_wing_50k <- read.csv(paste0(real_path, "corr_XGB_wing_50k.csv"))

# Create a data frame with correct row count
data <- data.frame(
  Method = factor(rep(c("Animal model", "RR", "PCRR", "XGB", "XGB_10k", "XGB_50k",
                        "Animal model", "RR", "PCRR", "XGB", "XGB_10k", "XGB_50k",
                        "Animal model", "RR", "PCRR", "XGB", "XGB_10k", "XGB_50k"),
                      each = 10),
                  levels = c("Animal model", "RR", "PCRR", "XGB", "XGB_10k", "XGB_50k")),
  Trait = factor(rep(c("Mass", "Mass", "Mass", "Mass", "Mass", "Mass",
                       "Wing", "Wing", "Wing", "Wing", "Wing", "Wing",
                       "Tarsus", "Tarsus", "Tarsus", "Tarsus", "Tarsus", "Tarsus"),
                     each = 10),
                 levels = c("Mass", "Wing", "Tarsus")),
  Values = c(massBV$corr,
             RR_mass$corr,
             PCRR_mass$corr,
             xgb_mass$corr,
             xgb_mass_10k$corr,
             xgb_mass_50k$corr,
             wingBV$corr,
             RR_wing$corr,
             PCRR_wing$corr,
             xgb_wing$corr,
             xgb_wing_10k$corr,
             xgb_wing_50k$corr,
             tarsusBV$corr,
             RR_tarsus$corr,
             PCRR_tarsus$corr,
             xgb_tarsus$corr,
             xgb_tarsus_10k$corr,
             xgb_tarsus_50k$corr)
)

# Define color mapping
method_colors <- c("Animal model" = "aquamarine", 
                   "XGB" = "darkcyan", 
                   "XGB_10k" = "purple", 
                   "XGB_50k" = "magenta",  
                   "RR" = "steelblue",
                   "PCRR" = "orange")

# Create the box plot
ggplot(data, aes(x = Method, y = Values, fill = Method)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 19, outlier.size = 2) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Prediction of Genetic Contribution",
    y = expression(Corr(bar(y), hat(y)^"*")),
    x = ""
  ) +
  theme_minimal(base_size = 18) +  # Base font size increased
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    strip.text = element_text(size = 18, face = "bold"), 
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 18),  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(1.5, "lines")
  ) +
  facet_grid(~Trait, scales = "free_x", space = "free_x")