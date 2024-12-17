# Load necessary R packages
source("clvm.R")
library(devtools) # Load devtools for package installation
devtools::install_github("JenniNiku/gllvm") # Install the latest gllvm package from GitHub
library(gllvm) # Load the gllvm package

# Load and preprocess data
data("microbialdata")
data <- microbialdata$Y # Extract response matrix

# Order columns by sparseness (number of zeros)
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum)) # Count zeros per column
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1] # Sort columns by sparseness
data <- data[, ind] # Reorder columns
dim(data) # Display dimensions of data
apply(data == 0, 2, sum) # Show zeros per column

# Extract and scale environmental covariates
X <- scale(microbialdata$Xenv[, c("SOM", "pH", "Phosp")])


# Fit constrained and unconstrained generalized linear latent variable models with negative binomial distribution
lvm_NB <- gllvm(
  y = data, X = X, family = "negative.binomial", sd.errors = FALSE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
) # Constrained model
save(lvm_NB, file = "lvm_NB.RData") # Save constrained model

unlvm_NB <- gllvm(
  y = data, X = X, family = "negative.binomial", sd.errors = FALSE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
) # Unconstrained model
save(unlvm_NB, file = "unlvmNB_model.RData") # Save unconstrained model

# Visualize residuals (reproduces Figure 8)

# the following plot will show the the upper plot in Figure 6, the D-S of concurrent constrained GLLVM_NB model with full data set

# Save D-S residuals against linear predictors, shown in upper left of Figure 8
pdf(file = "residuals against linear predictors.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(5, 5, 0, 0) + 0.1) # Set plot margins
plot(lvm_NB, which = 1 , caption =" ",var.colors = 1)
dev.off()

# Save Q-Q plot of residuals, shown in upper right of Figure 8
pdf(file = "qqnorm for residuals.pdf", width = 6, height = 6, useDingbats = FALSE)
plot(lvm_NB, which = 2 , caption =" ",var.colors = 1)
dev.off()

# the following plot will show the the bottom plot in Figure 6, the D-S of unconstrained GLLVM_NB model with full data set
# Plot D-S residuals against linear predictors, shown in bottom left of Figure 8
pdf(file = "un_residuals against linear predictors.pdf", width = 9, height = 9, useDingbats = FALSE)
# Set plot margins
par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(unlvm_NB, 
     which = 1, 
     caption = " ",       # Remove the default caption
     var.colors = 1,      # Set color for variables
     cex.lab = 2.5,       # Increase axis label text size
     cex.axis = 1.5,      # Increase axis tick label text size
     cex.main = 2.5,      # Increase main title text size
     cex = 1)             # Adjust the size of points or lines in the plot

dev.off()


# Save Q-Q plot of residuals, shown in bottom right of Figure 8
pdf(file = "un_qqnorm_for residuals.pdf", width = 9, height = 9, useDingbats = FALSE)
# Set plot margins
par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(unlvm_NB, 
     which = 2, 
     caption = " ",       # Remove the default caption
     var.colors = 1,      # Set color for variables
     cex.lab = 2.5,       # Increase axis label text size
     cex.axis = 1.5,      # Increase axis tick label text size
     cex.main = 2.5,      # Increase main title text size
     cex = 1)             # Adjust the size of points or lines in the plot

dev.off()


# Ordination plots for constrained and unconstrained models (reproduces Figure 7)
# Save ordination plot for the unconstrained model
pdf(file = "unconstrained_ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(4, 4, 0, 0) + 0.1) # Set plot margins

# Draw ordination plot
ordiplot(unlvm_NB,
  which.lvs = 1:2,
  s.colors = as.numeric(factor(microbialdata$X$Region)) + 1,
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15, # Use solid shapes
  ylim = c(-2, 4.9)
) # Adjust ylim to prevent legend overlap

# Extract unique regions and their corresponding symbols and colors
regions <- unique(microbialdata$X$Region) # Unique region names
pch_values <- as.numeric(factor(regions)) + 15 # Solid symbols: circle(16), triangle(17), diamond(18)
colors <- as.numeric(factor(regions)) + 1 # Assign colors to regions

# Add legend
legend("topleft", # Position of legend
  legend = regions, # Labels for regions
  pch = pch_values, # Symbols for regions
  col = colors
) # Colors for regions

dev.off()

# Save ordination plot for the constrained model
pdf(file = "ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(4, 4, 4, 0) + 0.1) # Set plot margins

# Draw ordination plot
ordiplot(lvm_NB,
  which.lvs = 1:2,
  s.colors = as.numeric(factor(microbialdata$X$Region)) + 1,
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15
) # Use solid shapes

# Add legend
legend("topleft", # Position of legend
  legend = regions, # Labels for regions
  pch = pch_values, # Symbols for regions
  col = colors
) # Colors for regions

dev.off()

# Compute AIC and BIC for other models (Table 2)
# Fit concurrent gllvm model with zero-inflated negative binomial (ZINB)
lvmZINB_model <- gllvm(
  y = data, X = X, family = "ZINB", sd.errors = TRUE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
)
save(lvmZINB_model, file = "lvmZINB_model.RData") # Save ZINB model

# Fit unconstrained gllvm model with zero-inflated negative binomial (ZINB)

unlvmZINB_model<- gllvm(y = data , family = "ZINB", sd.errors = TRUE, row.eff = "fixed",seed = 123)
save(lvmZINB_model, file = "unlvmZINB_model.RData") # Save ZINB model