# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(scales)  # Ensure the scales package is loaded

# Set file path
dset <- "70k"  # Change this to the correct dataset name
base_path <- "c:/Users/gard_/Documents/MasterThesis/Code/Results/70k/Ridge"

# Load the data
effects <- list()
for (i in 1:10) {
  file_path <- paste0(base_path, "/ridge_coefs_arch", i, ".csv")
  df <- read_csv(file_path)
  effects[[i]] <- df$beta_R  # Extracting just the beta_R values
}

# Define π values
pi_vec <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)

# Create histogram plots for each architecture with automatic scaling
plot_list <- list()
for (i in 1:10) {
  beta_values <- effects[[i]]  # Use the beta_R values directly
  
  # Calculate min and max of the beta values
  x_range <- range(beta_values, na.rm = TRUE)
  
  # Generate exactly 3 breaks, split the range into three intervals
  breaks <- seq(x_range[1], x_range[2], length.out = 3)  # 4 points = 3 intervals
  
  p <- ggplot(data.frame(Beta = beta_values), aes(x = Beta)) +
    geom_histogram(binwidth = 0.00015, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(bquote(pi[0] == .(pi_vec[i]))) +
    theme_minimal(base_size = 16) +  # Increase base font size
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 10),  # Larger tick labels
      plot.margin = margin(10, 10, 20, 10),  # Extra bottom margin
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
    ) +
    scale_x_continuous(
      breaks = breaks,  # Use the explicitly defined breaks
      labels = label_scientific(accuracy = 0.001)  # Use scientific notation for labels
    )
  
  plot_list[[i]] <- p
}

# Create an empty page
grid.newpage()

# Arrange plots in a 2x5 grid
grid_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = 5, bottom = textGrob("", gp = gpar(fontsize = 10)))

# Add common axis labels with improved spacing and size
grid.text("Frequency", x = 0.0085, y = 0.5, rot = 90, gp = gpar(fontsize = 20))  # shifted left
grid.text(expression(hat(u)[R]), x = 0.525, y = 0.035, gp = gpar(fontsize = 20))
