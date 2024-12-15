# Load required libraries
library(ggplot2) # For data visualization
library(devtools) # For installing packages from GitHub
devtools::install_github("JenniNiku/gllvm") # Install the latest gllvm package
library(gllvm) # Load gllvm package
library(reshape2) # For data reshaping
library(scales)
# Load microbial data
data("microbialdata") # Load example microbial data
data <- microbialdata$Y # Extract the response matrix

# Order columns based on sparseness (number of zeros in each column)
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum)) # Count zeros in each column
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1] # Sort columns by sparseness
data <- data[, ind] # Reorder columns accordingly

# Draw mean-variance plot, (shown in right picture of Figure 1)
# Compute mean and variance for each column
means <- apply(data, 2, mean) # Column-wise means
variances <- apply(data, 2, var) # Column-wise variances

# Define a Negative Binomial (NB) variance function
nb_model <- function(mean, theta) {
  mean + mean^2 * theta # NB variance formula
}

# Fit the NB model to the data using nonlinear least squares
fit <- nls(variances ~ nb_model(means, theta),
  start = list(theta = 0.1)
) # Use an initial guess for theta

# Extract the optimal dispersion parameter (theta)
theta_opt <- coef(fit)["theta"] # Extract the estimated theta value
cat("Optimal dispersion parameter (theta):", theta_opt, "\n")


# Save the overdispersion plot to a PDF
pdf(file = "overdispersion.pdf", width = 7, height = 6, useDingbats = FALSE)

# Set plot margins
par(mar = c(4, 4, 0.5, 0) + 0.1)

regions <- microbialdata$X$Region
# Extract unique regions and symbols
regions <- unique(microbialdata$X$Region)  # Unique region names
pch_values <- as.numeric(factor(regions)) + 15  # Define symbols manually: circle (15), triangle (16), plus (17)
colors <- c("black", "blue", "grey")[as.numeric(factor(regions))]  # Set specific colors for each region

# Plot mean-variance relationship
plot(means, variances,
  xlab = "Mean", # X-axis label
  ylab = "Variance", # Y-axis label
  pch = pch_values, # Point shape based on region
  cex = 0.7, # Point size
  col = colors, # Point color based on region
  xlim = c(0, 30), # Limit for X-axis
  ylim = c(0, 1000), # Limit for Y-axis
  cex.lab = 1.5, # Label font size
  cex.axis = 1.4, # Axis font size
  cex.main = 1.6
) # Title font size

# Add a Poisson reference line (variance = mean)
abline(a = 0, b = 1, col = "red", lwd = 2) # Poisson variance line

# Add the Negative Binomial relationship line
theta <- theta_opt # Use the estimated dispersion parameter
nb_variance <- means + means^2 * theta # Compute NB variance
lines(sort(means), sort(nb_variance),
  col = "red", lwd = 2, lty = 2
) # Add the NB variance line (dashed)

# Define full name of region as new region names
new_region_names <- c("Mayrhofen","Kilpisjarvi","Ny-Alesund")
# Add a legend with specified colors and symbols
legend("topleft",                              # Legend position
       legend = new_region_names,              # Legend labels
       pch = pch_values,                       # Corresponding symbols
       col = colors,                           # Corresponding colors (black, blue, grey)
       ncol = 3,
       cex = 1.1) 
# Close the PDF device
dev.off()


# Draw heatmap of data, (shown in left picture of Figure 1)

# Convert data to long format for plotting
long_data <- melt(data, varnames = c("Site", "Species"), value.name = "Richness")

# Set up the plot size and output as a PDF
pdf(file = "sparsity.pdf", width = 8, height = 7, useDingbats = FALSE)
par(mar = c(0, 0, 4, 4) + 0.1)

# Generate the heatmap
ggplot(long_data, aes(x = Species, y = Site, fill = Richness)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(hcl.colors(n = 8, palette = "Spectral")), # Reverse gradient colors
    values = rescale(c(0, 2, 5, 10, 15, 50, 100, 200, 1024)), # Emphasize 0-15 range
    oob = scales::squish, # Handle out-of-bounds values
    limits = c(0, 1024) # Set gradient range
  ) +
  theme_minimal() + # Minimal theme
  labs(
    x = "Species",
    y = "Sample",
    fill = "Count" # Axis and legend labels
  ) +
  theme(
    axis.text.x = element_blank(), # Remove x-axis labels
    axis.text.y = element_blank(), # Remove y-axis labels
    axis.ticks = element_blank(), # Remove axis ticks
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 9),
    legend.position = "bottom", # Move legend to the bottom
    legend.direction = "horizontal" # Align legend horizontally
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 32,
      barheight = 2,
      breaks = seq(0, 250, by = 50)
    )
  )

# Close the PDF device to save the plot
dev.off()
