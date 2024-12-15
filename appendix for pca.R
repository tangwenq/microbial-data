library(gllvm)
library(vegan)
library(robCompositions)
data("microbialdata")
data <- microbialdata$Y
# order columns according to sparseness
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum))
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]
data <- data[, ind]
# distance based methods; nMDS
# Perform nMDS analysis
fit_mds <- try(metaMDS(data, distance = "bray", trace = FALSE))

# Extract scaled coordinates for plotting
mds_ords <- scale(fit_mds$points)

# Perform Singular Value Decomposition (SVD)
svd_result <- svd(mds_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)
mds_ords[, 1] <- rotated_coords[, 1]
mds_ords[, 2] <- rotated_coords[, 2]
# Save ordination plot as PDF
pdf(file = "nmds_ord.pdf", width = 6, height = 6, useDingbats = FALSE) # Correct syntax

# Set plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)

# Define colors and symbols for regions
regions <- unique(microbialdata$X$Region) # Unique region names
colors <- c("black", "blue", "grey")[as.numeric(factor(regions))] # Custom colors
pch_values <- as.numeric(factor(regions)) + 15 # Custom symbols
# Plot nMDS ordination
plot(mds_ords[, 1], mds_ords[, 2],
  col = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE,
  xlim = range(mds_ords[, 1]) * 1.1,
  ylim = c(min(clrmds_ords[, 2]), max(clrmds_ords[, 2]) * 1.7)
)

# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)

# Define full region names for legend
new_region_names <- c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund")

# Add legend
legend("topleft",
  legend = new_region_names,
  pch = pch_values,
  col = colors,
  ncol = 3,
  cex = 1.3
)
# Close the graphics device
dev.off()

# PCA on clr-transformed data
clr_data <- cenLR(data + 1)
#
fit_clrmds <- try(metaMDS(clr_data$x.clr, distance = "euclidean", autotransform = FALSE, noshare = FALSE, wascores = FALSE, trace = FALSE))
clrmds_ords <- scale(fit_clrmds$points)

# Perform Singular Value Decomposition (SVD)
svd_result <- svd(clrmds_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)
clrmds_ords[, 1] <- rotated_coords[, 1]
clrmds_ords[, 2] <- rotated_coords[, 2]
pdf(file = "clrmds_ord.pdf", width = 6, height = 6, useDingbats = FALSE)

# Adjust plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)


# Plot CLR-based nMDS ordination
plot(clrmds_ords[, 1], clrmds_ords[, 2],
  col = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE,
  xlim = range(clrmds_ords[, 1]) * 1.1,
  ylim = c(min(clrmds_ords[, 2]), max(clrmds_ords[, 2]) * 1.5)
)

# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 1", cex.lab = 1.3)


# Add legend
legend("topleft",
  legend = new_region_names,
  pch = pch_values,
  col = colors,
  ncol = 3,
  cex = 1.3
)

# Close the graphics device
dev.off()

#
fit_pca <- prcomp(clr_data$x.clr, scale = TRUE)
pca_ords <- scale(fit_pca$x[, 1:2]) # extract simulated ordinations from PCA
# Perform Singular Value Decomposition (SVD)
svd_result <- svd(pca_ords)

# Reconstruct the rotated coordinates
rotated_coords <- svd_result$u %*% diag(svd_result$d)
pca_ords[, 1] <- rotated_coords[, 1]
pca_ords[, 2] <- rotated_coords[, 2]
# Save ordination plot as PDF
pdf(file = "pca_ord.pdf", width = 6, height = 6, useDingbats = FALSE)

# Adjust plot margins
par(mar = c(4.5, 4.5, 0, 0) + 0.1)

# Plot PCA ordination
plot(pca_ords[, 1], pca_ords[, 2],
  col = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE,
  xlim = range(pca_ords[, 1]) * 1.1,
  ylim = c(min(pca_ords[, 2]), max(pca_ords[, 2]) * 1.5)
)
# Add axis titles
title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)


# Add legend
legend("topleft",
  legend = new_region_names,
  pch = pch_values,
  col = colors,
  ncol = 3,
  cex = 1.3
)

# Close the graphics device
dev.off()
