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


# Plot D-S residuals against linear predictors, shown in left of Figure 8
pdf(file = "residuals against linear predictors.pdf", width = 9, height = 9, useDingbats = FALSE)
# Set plot margins
par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(lvm_NB,
  which = 1,
  caption = " ", # Remove the default caption
  var.colors = 1, # Set color for variables
  cex.lab = 2.5, # Increase axis label text size
  cex.axis = 1.5, # Increase axis tick label text size
  cex.main = 2.5, # Increase main title text size
  cex = 1
) # Adjust the size of points or lines in the plot

dev.off()


# Save Q-Q plot of residuals
pdf(file = "qqnorm_for residuals.pdf", width = 9, height = 9, useDingbats = FALSE)
# Set plot margins
par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(lvm_NB,
  which = 2,
  caption = " ", # Remove the default caption
  var.colors = 1, # Set color for variables
  cex.lab = 2.5, # Increase axis label text size
  cex.axis = 1.5, # Increase axis tick label text size
  cex.main = 2.5, # Increase main title text size
  cex = 1
) # Adjust the size of points or lines in the plot

dev.off()

# Ordination plots for constrained and unconstrained models (reproduces Figure 7)
# Save ordination plot for the unconstrained model
pdf(file = "unconstrained_ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(5.3, 5.3, 0, 0) + 0.1) # Set plot margins
# Draw ordination plot
ordiplot(unlvm_NB,
  which.lvs = 1:2,
  s.colors = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE
) # Set symbols based on region

title(xlab = "ordination coordinate 1", ylab = "ordination coordinate 2", cex.lab = 2)


# Extract unique regions and symbols
regions <- unique(microbialdata$X$Region) # Unique region names
pch_values <- as.numeric(factor(regions)) + 15 # Define symbols manually: circle (15), triangle (16), plus (17)
colors <- c("black", "blue", "grey")[as.numeric(factor(regions))] # Set specific colors for each region

# Define full name of region as new region names
new_region_names <- c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund")
# Add a legend with specified colors and symbols
legend("topleft", # Legend position
  legend = new_region_names, # Legend labels
  pch = pch_values, # Corresponding symbols
  col = colors, # Corresponding colors (black, blue, grey)
  ncol = 3,
  cex = 1.2
) # Arrange legend in 3 columns
dev.off()

# Save ordination plot for the constrained model
pdf(file = "ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(5.3, 5.3, 0, 0) + 0.1) # Set plot margins

# Draw ordination plot
ordiplot(lvm_NB,
  which.lvs = 1:2,
  s.colors = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE
) # Set symbols based on region

title(xlab = "ordination coordinate 1", ylab = "ordination coordinate 2", cex.lab = 2)


# Extract unique regions and symbols
regions <- unique(microbialdata$X$Region) # Unique region names
pch_values <- as.numeric(factor(regions)) + 15 # Define symbols manually: circle (15), triangle (16), plus (17)
colors <- c("black", "blue", "grey")[as.numeric(factor(regions))] # Set specific colors for each region


# Add a legend with specified colors and symbols
legend("topleft", # Legend position
  legend = new_region_names, # Legend labels with full name of sample site
  pch = pch_values, # Corresponding symbols
  col = colors, # Corresponding colors (black, blue, grey)
  ncol = 3, # Arrange legend in 3 columns
  cex = 1.2
)
dev.off()

# Compute AIC and BIC for other models (Table 2)
# Fit gllvm model with zero-inflated negative binomial (ZINB)
lvmZINB_model <- gllvm(
  y = data, X = X, family = "ZINB", sd.errors = TRUE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
)
save(lvmZINB_model, file = "lvmZINB_model.RData") # Save ZINB model

# Fit copula models
cNB_model <- fit_copula(data,
  X = X, reff = "fixed", gllvm.fam = "negative.binomial",
  sd.errors = TRUE, seed = 123, lv.n = 0
) # Copula NB model
save(cNB_model, file = "cNB_model.RData") # Save copula NB model

cZINB_model <- fit_copula(data,
  X = X, reff = "fixed", gllvm.fam = "ZINB",
  sd.errors = TRUE, seed = 123, lv.n = 0
) # Copula ZINB model
save(cZINB_model, file = "cZINB_model.RData") # Save copula ZINB model
