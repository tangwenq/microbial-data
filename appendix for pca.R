#The following will show the ordination plot of classical plot, shown in Figure 8
# Load necessary libraries
library(gllvm)       # Generalized Linear Latent Variable Models
library(vegan)       # Community Ecology Package (for ordination and analysis)
library(robCompositions) # Robust Compositional Data Analysis
data("microbialdata") # Load the microbial dataset from the 'microbialdata' package

# Extract the data matrix (Y) from the microbialdata
data <- microbialdata$Y

# Order columns according to sparseness (number of zeros)
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum)) # Compute number of zeros per column
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]  # Sort columns by sparseness (zero count)
data <- data[, ind]  # Reorder columns based on sparseness

# Distance-based methods: nMDS (Non-metric Multidimensional Scaling)
# Perform nMDS analysis using Bray-Curtis distance
fit_mds <- try(metaMDS(data, distance = "bray", trace = FALSE))

# Extract scaled coordinates for plotting
mds_ords <- scale(fit_mds$points)

# Perform Singular Value Decomposition (SVD) to rotate the coordinates
svd_result <- svd(mds_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)

# Update the coordinates with the rotated values
mds_ords[, 1] = rotated_coords[, 1]
mds_ords[, 2] = rotated_coords[, 2]

# Save the nMDS ordination plot as a PDF
pdf(file = "nmds_ord.pdf", width = 6, height = 6, useDingbats = FALSE) 

# Set plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)

# Define colors and symbols for regions
regions <- unique(microbialdata$X$Region) # Unique region names
colors <- c("black", "blue", "grey")[as.numeric(factor(regions))] # Custom colors for regions
pch_values <- as.numeric(factor(regions)) + 15 # Custom symbols (e.g., circle, triangle, plus)

# Plot nMDS ordination with customized colors and symbols
plot(mds_ords[, 2], mds_ords[, 1], 
     col = colors[as.numeric(factor(microbialdata$X$Region))],  # Color by region
     pch = pch_values[as.numeric(factor(microbialdata$X$Region))],  # Symbol by region
     main = "", 
     ann = FALSE, 
     xlim = range(mds_ords[, 2]) * 1.1, 
     ylim = c(min(mds_ords[, 1]), max(mds_ords[, 1]) * 1.3)) # Set axis limits

# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)

# Define full region names for the legend
new_region_names <- c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund")

# Add legend with symbols and colors corresponding to regions
legend("topleft", 
       legend = new_region_names, 
       pch = pch_values, 
       col = colors, 
       ncol = 3, 
       cex = 1.3)

# Close the graphics device (finalize the plot)
dev.off()

# PCA on CLR-transformed data (Centered Log-Ratio Transformation)
clr_data <- cenLR(data + 1)  # Apply CLR transformation (adding 1 to avoid log(0))

# Perform nMDS on CLR-transformed data using Euclidean distance
fit_clrmds <- try(metaMDS(clr_data$x.clr, distance = "euclidean", autotransform = FALSE, noshare = FALSE, wascores = FALSE, trace = FALSE))
clrmds_ords <- scale(fit_clrmds$points)  # Extract scaled coordinates

# Perform Singular Value Decomposition (SVD) on CLR-based nMDS ordination
svd_result <- svd(clrmds_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)

# Update the coordinates with the rotated values
clrmds_ords[, 1] = rotated_coords[, 1]
clrmds_ords[, 2] = rotated_coords[, 2]

# Save CLR-based nMDS ordination plot as PDF
pdf(file = "clrmds_ord.pdf", width = 6, height = 6, useDingbats = FALSE)

# Set plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)

# Plot CLR-based nMDS ordination
plot(clrmds_ords[, 2], clrmds_ords[, 1], 
     col = colors[as.numeric(factor(microbialdata$X$Region))], 
     pch = pch_values[as.numeric(factor(microbialdata$X$Region))], 
     main = "", 
     ann = FALSE, 
     xlim = range(clrmds_ords[, 2]) * 1.1, 
     ylim = c(min(clrmds_ords[, 1]), max(clrmds_ords[, 1]) * 1.4))

# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 1", cex.lab = 1.3)

# Add legend for regions
legend("topleft", 
       legend = new_region_names, 
       pch = pch_values, 
       col = colors, 
       ncol = 3, 
       cex = 1.3)

# Close the graphics device
dev.off()

# PCA on CLR-transformed data using prcomp (Principal Component Analysis)
fit_pca <- prcomp(clr_data$x.clr, scale = TRUE)  # Perform PCA
pca_ords <- scale(fit_pca$x[, 1:2])  # Extract first two principal components

# Perform Singular Value Decomposition (SVD) on PCA ordination
svd_result <- svd(pca_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)

# Update the coordinates with the rotated values
pca_ords[, 1] = rotated_coords[, 1]
pca_ords[, 2] = rotated_coords[, 2]

# Save PCA ordination plot as PDF
pdf(file = "pca_ord.pdf", width = 6, height = 6, useDingbats = FALSE)

# Set plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)

# Plot PCA ordination
plot(pca_ords[, 2], pca_ords[, 1], 
     col = colors[as.numeric(factor(microbialdata$X$Region))], 
     pch = pch_values[as.numeric(factor(microbialdata$X$Region))], 
     main = "", 
     ann = FALSE, 
     xlim = range(pca_ords[, 2]) * 1.1, 
     ylim = c(min(pca_ords[, 1]), max(pca_ords[, 1]) * 1.2))

# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)

# Add legend for regions
legend("topleft", 
       legend = new_region_names, 
       pch = pch_values, 
       col = colors, 
       ncol = 3, 
       cex = 1.3)

# Close the graphics device
dev.off()

