# Load necessary libraries
library(ggplot2)
library(dplyr)

# Assuming you have the 'aggregated_data' from previous steps
# Create a smooth line plot
smooth_plot <- ggplot(aggregated_data, aes(x = bin_asset, y = bin_density)) +
  geom_line(group = 1, color = "blue", size = 1) +  # Line plot
  geom_point(color = "red") +  # Add points for clarity
  geom_smooth(method = "loess", color = "darkblue", se = FALSE) +  # Smooth line
  theme_minimal() +
  labs(title = "Smooth Distribution of Wealth",
       x = "Asset Value",
       y = "Aggregated Density") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(smooth_plot)

# Check if the directory 'plots' exists, if not, create it
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Save the plot
ggsave(filename = "plots/smooth_distribution_wealth.png", plot = smooth_plot, width = 10, height = 6, dpi = 300)
