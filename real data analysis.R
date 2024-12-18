# The following code will reproduce the Figure 6 and Figure 5
# Load necessary R packages
source("GCLVM.R")
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


# Fit constrained and unconstrained generalized linear latent variable models with negative binomial distribution(NB-GLLVM) with full data

# Constrained concurrent model
gllvm_NB <- gllvm(
  y = data, X = X, family = "negative.binomial", sd.errors = TRUE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
)
save(gllvm_NB, file = "gllvm_NB.RData") # Save constrained model

# Unconstrained model
ungllvm_NB <- gllvm(
  y = data, family = "negative.binomial", sd.errors = TRUE,
  row.eff = "fixed", num.lv = 2, seed = 123
)
save(ungllvm_NB, file = "ungllvmNB_model.RData") # Save unconstrained model

# Visualize residuals for two NB-GLLVM model(reproduces Figure 6)


# the D-S of constrained NB-GLLVM model with full data set, shown in bottom row in Figure 6

# D-S residuals against linear predictors
pdf(file = "residuals against linear predictors.pdf", width = 9, height = 9, useDingbats = FALSE)
par(mar = c(6, 6, 0, 0) + 0.1) # Set plot margins
# Create the plot
plot(gllvm_NB,
  which = 1,
  caption = " ",
  var.colors = 1,
  cex.lab = 2.5,
  cex.axis = 1.5,
  cex.main = 2.5,
  cex = 1
)

dev.off()

# Q-Q plot of residuals
pdf(file = "qqnorm for residuals1.pdf", width = 9, height = 9, useDingbats = FALSE)

par(mar = c(6, 6, 0, 0) + 0.1)

plot(gllvm_NB,
  which = 2,
  caption = " ", # Remove the default caption
  var.colors = 1, # Set color for variables
  cex.lab = 2.5, # Increase axis label text size
  cex.axis = 1.5, # Increase axis tick label text size
  cex.main = 2.5, # Increase main title text size
  cex = 1
) # Adjust the size of points or lines in the plot

dev.off()


# the D-S of unconstrained GLLVM_NB model with full data set, shown in top row in Figure 6

# D-S residuals against linear predictors
pdf(file = "un_residuals against linear predictors.pdf", width = 9, height = 9, useDingbats = FALSE)

par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(ungllvm_NB,
  which = 1,
  caption = " ",
  var.colors = 1,
  cex.lab = 2.5,
  cex.axis = 1.5,
  cex.main = 2.5,
  cex = 1
)

dev.off()


# Save Q-Q plot of residuals, shown in bottom right of Figure 6
pdf(file = "un_qqnorm_for residuals.pdf", width = 9, height = 9, useDingbats = FALSE)
# Set plot margins
par(mar = c(6, 6, 0, 0) + 0.1)

# Create the plot
plot(ungllvm_NB,
  which = 2,
  caption = " ", # Remove the default caption
  var.colors = 1, # Set color for variables
  cex.lab = 2.5, # Increase axis label text size
  cex.axis = 1.5, # Increase axis tick label text size
  cex.main = 2.5, # Increase main title text size
  cex = 1
) # Adjust the size of points or lines in the plot

dev.off()


# Ordination plots for constrained and unconstrained models (reproduces Figure 5)

# unconstrained NB-GLLVM model

pdf(file = "unconstrained_ord.pdf", width = 6, height = 6, useDingbats = FALSE)
par(mar = c(4.5, 4.5, 0, 0) + 0.1) 
# Draw ordination plot
ordiplot(ungllvm_NB,
  which.lvs = 1:2,
  s.colors = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE
) # Set symbols based on region

title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)


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
  cex = 1.3
) # Arrange legend in 3 columns
dev.off()

# constrained NB-GLLVM
pdf(file = "ord.pdf", width = 6, height = 6, useDingbats = FALSE)

par(mar = c(4.5, 4.5, 0, 0) + 0.1) 


ordiplot(gllvm_NB,
  which.lvs = 1:2,
  s.colors = c("black", "blue", "grey")[as.numeric(factor(microbialdata$X$Region))],
  rotate = TRUE,
  symbols = TRUE,
  pch = as.numeric(factor(microbialdata$X$Region)) + 15,
  main = "",
  ann = FALSE
) # Set symbols based on region

title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)


# Add a legend with specified colors and symbols
legend("topleft", # Legend position
  legend = new_region_names, # Legend labels
  pch = pch_values, # Corresponding symbols
  col = colors, # Corresponding colors (black, blue, grey)
  ncol = 3,
  cex = 1.3
) # Arrange legend in 3 columns

dev.off()


# Compute AIC and BIC for ZINB-GLLVM models
# Fit concurrent gllvm model with zero-inflated negative binomial
gllvm_ZINB <- gllvm(
  y = data, X = X, family = "ZINB", sd.errors = TRUE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
)
save(gllvm_ZINB, file = "gllvm_ZINB.RData") # Save ZINB model

# Fit unconstrained gllvm model with zero-inflated negative binomial (ZINB)

ungllvm_ZINB <- gllvm(y = data, family = "ZINB", sd.errors = TRUE, row.eff = "fixed", seed = 123)
save(ungllvm_ZINB, file = "ungllvm_ZINB.RData") # Save ZINB model
