# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(arrow)

# Set file path
dset <- "70k"  # Change this to the correct dataset name
base_path <- paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/",dset,"/SHAP/Sim")

# Load the data
effects <- list()
for (i in 1:10) {
  file_path <- paste0(base_path, '/shap_',dset,'_arch_',i,'.feather')
  df <- read_feather(file_path)
  
  # Convert all columns to a single vector containing the SHAP values
  effects[[i]] <- as.vector(unlist(df))  # Flatten the data to a single vector
}

# Define π values
pi_vec <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)

# Create histogram plots for each architecture with automatic x-axis scaling
plot_list <- list()
for (i in 1:10) {
  shap_values <- effects[[i]]  # Use the SHAP values directly
  x_range <- range(shap_values, na.rm = TRUE)
  
  # Generate 2 evenly spaced breaks + 1 at the max
  core_breaks <- seq(x_range[1], x_range[2], length.out = 3)
  # Ensure the last break is exactly at the maximum value (avoid duplication if already included)
  breaks <- unique(c(core_breaks, x_range[2]))
  
  p <- ggplot(data.frame(Effect = shap_values), aes(x = Effect)) +
    geom_histogram(binwidth = 0.0002, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(bquote(pi[0] == .(pi_vec[i]))) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 10),
      plot.margin = margin(10, 10, 20, 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_continuous(
      breaks = breaks,
      labels = label_scientific(accuracy = 0.001)
    )
  
  plot_list[[i]] <- p
}

# Create an empty page for grid layout
grid.newpage()

# Arrange plots in a 2x5 grid
grid_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = 5)

# Add common axis labels with improved spacing and size
grid.text("Frequency", x = 0.0085, y = 0.5, rot = 90, gp = gpar(fontsize = 20))  # shifted left
grid.text(expression("mean |SHAP|"), x = 0.525, y = 0.02, gp = gpar(fontsize = 20))