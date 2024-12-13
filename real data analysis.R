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
# Extract deviance residuals from the constrained model
lvm_res <- residuals(lvm_NB)$residuals

# Save histogram of residuals
pdf(file = "hist for residuals.pdf", width = 6, height = 6, useDingbats = FALSE)
hist(lvm_res,
  main = "Histogram of Residuals", xlab = "Residuals",
  col = "lightblue", border = "black"
)
dev.off()

# Save Q-Q plot of residuals
pdf(file = "qqnorm for residuals.pdf", width = 6, height = 6, useDingbats = FALSE)
qqnorm(lvm_res)
dev.off()

# Ordination plots for constrained and unconstrained models (reproduces Figure 7)
# Save ordination plot for the unconstrained model
pdf(file = "unconstrained_ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(4, 4, 4, 0) + 0.1) # Set plot margins

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
