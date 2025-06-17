# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(scales)
library(stringr)

# Set dataset and base path
dset <- "70k"
base_path <- paste0("c:/Users/gard_/Documents/MasterThesis/Code/Results/", dset, "/Lasso")

# List available coefficient files (lasso only)
files <- list.files(base_path, pattern = "^lasso_coefs_arch\\d+\\.csv$", full.names = TRUE)

# Extract architecture numbers from file names
arch_nums <- as.integer(str_extract(basename(files), "\\d+"))

# Sort architecture numbers (optional, for consistent order)
arch_nums <- sort(arch_nums)

# Define π values for all potential architectures
pi_vec_all <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)

# Load effect sizes into list
effects <- list()
for (i in seq_along(arch_nums)) {
  arch <- arch_nums[i]
  file_path <- paste0(base_path, "/lasso_coefs_arch", arch, ".csv")
  df <- read_csv(file_path, show_col_types = FALSE)
  effects[[i]] <- df$beta_L
}

# Generate plots
plot_list <- list()
for (i in seq_along(arch_nums)) {
  beta_values <- effects[[i]]
  arch <- arch_nums[i]
  pi_val <- pi_vec_all[arch]  # Use architecture number to index π
  
  # Calculate histogram range and breaks
  x_range <- range(beta_values, na.rm = TRUE)
  breaks <- seq(x_range[1], x_range[2], length.out = 3)
  
  p <- ggplot(data.frame(Beta = beta_values), aes(x = Beta)) +
    geom_histogram(binwidth = 0.0015, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(bquote(pi[0] == .(pi_val))) +
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

# Create and draw plot grid
grid.newpage()
grid.arrange(grobs = plot_list, nrow = 2, ncol = ceiling(length(plot_list) / 2))

# Add common axis labels
grid.text("Frequency", x = 0.0085, y = 0.5, rot = 90, gp = gpar(fontsize = 20))
grid.text(expression(hat(u)[L]), x = 0.525, y = 0.035, gp = gpar(fontsize = 20))
