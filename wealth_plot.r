# Load necessary library
library(ggplot2)
library(dplyr)

# Read data from the CSV file
data <- read.csv("output/wealth_distribution.csv")

# Aggregate data into fewer bins
# For example, create 10 bins from the range of asset values
data$bin <- cut(data$asset, breaks = 200, labels = FALSE)  # Assign each asset to a bin

# Sum densities within each bin
aggregated_data <- data %>%
  group_by(bin) %>%
  summarise(
    bin_density = sum(density),  # Sum of densities in each bin
    bin_asset = mean(asset)      # Mean asset value for label
  )

# Create the plot using ggplot2 with aggregated data
plot_g <- ggplot(aggregated_data, aes(x = bin_asset, y = bin_density)) +
  geom_col() +  # Still using geom_col because we manually aggregate densities
  theme_minimal() +
  labs(title = "Aggregated Histogram of Wealth Distribution",
       x = "Asset Value",
       y = "Aggregated Density") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(plot_g)

# Check if the directory 'plots' exists, if not, create it
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Save the plot
ggsave(filename = "plots/aggregated_wealth_distribution_histogram.png", plot = plot_g, width = 10, height = 6, dpi = 300)
