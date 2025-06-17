# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)

# Set file path
dset <- "70k"  # Change this to the correct dataset name
base_path <- "c:/Users/gard_/Documents/MasterThesis/Code/Data/Phenotypes/"

# Load the data
effects <- list()
for (i in 1:10) {
  file_path <- paste0(base_path, dset, "/Sim_effect_", dset, "_arch_", i, ".csv")
  df <- read_csv(file_path)
  effects[[i]] <- df
}
head(effects)

# Define π values
pi_vec <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)

# Create histogram plots for each architecture with individual axis scaling
plot_list <- list()
for (i in 1:10) {
  effects_values <- effects[[i]][[1]]  # Use all effect values directly
  
  p <- ggplot(data.frame(Effect = effects_values), aes(x = Effect)) +
    geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(bquote(pi[0] == .(pi_vec[i]))) +
    theme_minimal(base_size = 16) +  # Increase base font size
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 10),  # Larger tick labels
      plot.margin = margin(10, 10, 20, 10)  # Extra bottom margin
    )
  
  plot_list[[i]] <- p
}

# Create an empty page
grid.newpage()

# Arrange plots in a 2x5 grid
grid_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = 5, bottom = textGrob("", gp = gpar(fontsize = 10)))

# Add common axis labels with improved spacing and size
grid.text("Frequency", x = 0.0085, y = 0.5, rot = 90, gp = gpar(fontsize = 20))  # shifted left
grid.text(expression(u), x = 0.525, y = 0.035, gp = gpar(fontsize = 20))
# Adjusted y position

