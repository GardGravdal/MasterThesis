####################################################################
# Script for making histogram plots of (observed) trait distributions
####################################################################

# Data preparation helper script:
setwd("C:/Users/gard_/Documents/MasterThesis/Code")
source("h_dataPrep.r")


# Some data wranging to ensure that the IDs in the data correspond to the IDs in the A and G-matrices (nothing to worry about):
# indicates that some IDs are missing:
#d.map[3110:3125, ]
# from this we see the number of animals
Nanimals <- 3116


# remove missing values (mass)
d.morph <- filter(d.morph, !is.na(eval(tarsus)))
d.morph <- filter(d.morph, !is.na(eval(mass)))
d.morph <- filter(d.morph, !is.na(eval(wing)))

var(d.morph$wing)
var(d.morph$tarsus)
var(d.morph$mass)

hist(d.morph$tarsus)
hist(d.morph$mass)
hist(d.morph$wing)

library(ggplot2)

# Calculate variances
var_wing <- var(d.morph$wing, na.rm = TRUE)
var_tarsus <- var(d.morph$tarsus, na.rm = TRUE)
var_mass <- var(d.morph$mass, na.rm = TRUE)

# Data frames
df_wing <- data.frame(value = d.morph$wing)
df_tarsus <- data.frame(value = d.morph$tarsus)
df_mass <- data.frame(value = d.morph$mass)

# Custom theme with large text, no title
large_text_theme <- theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.margin = margin(10, 10, 10, 10)
  )

#------------------

# Plot function WITH title and larger text
plot_hist <- function(data, var_value, fill_color, x_label, plot_title) {
  ggplot(data, aes(x = value)) +
    geom_histogram(fill = fill_color, color = "black", bins = 30) +
    labs(title = plot_title, x = x_label, y = "Frequency") +
    annotate("text", x = Inf, y = Inf, label = paste("Variance:", round(var_value, 3)),
             hjust = 1.1, vjust = 1.5, size = 8, fontface = "italic") +
    large_text_theme +
    theme(
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5),  # Centered and larger
      axis.title = element_text(size = 32),
      axis.text = element_text(size = 28),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22)
    )
}

# Create the plots WITH titles
plot_hist(df_wing, var_wing, "skyblue", "Wing length (mm)", "Wing")
plot_hist(df_tarsus, var_tarsus, "lightgreen", "Tarsus length (mm)", "Tarsus")
plot_hist(df_mass, var_mass, "salmon", "Body mass (g)", "Mass")
