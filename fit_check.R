# compare the performance of different modeling techniques (GLLVM vs Copula), shown in Table 1
# Load required libraries
source("GCLVM")
library(devtools) # Load devtools for package installation
devtools::install_github("JenniNiku/gllvm") # Install the latest gllvm package from GitHub
library(gllvm) # Generalized Linear Latent Variable Models
library(mvtnorm) # Multivariate normal and t-distributions

# Check the goodness of fit with 3 measures
# 1. MARNE: Mean Absolute Range Normalized Error.
# 2. Species-wise Correlations: mCOR
# 3. Overall Correlation (whole data): cor
# Global correlation between pred_y and y across all data points.
# Load microbial data
data("microbialdata")
data <- microbialdata$Y

# Order columns according to sparseness (number of zeros in each column)
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum))
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]
data <- data[, ind]

# Check dimensions and sparseness
dim(data) # Dimensions of the dataset
apply(data == 0, 2, sum) # Number of zeros in each column

# Compute goodness-of-fit measures for 4 models under different dimensions with microbial data
# The result is shown in Table 1

# 1. GLLVM with ZINB distribution

GOF_gllvmZINB <- matrix(NA, 4, 3) # Initialize results matrix
colnames(GOF_gllvmZINB) <- c("marne", "mCOR", "cor") # Column names for the measures
rownames(GOF_gllvmZINB) <- c("50", "100", "200", "400") # Rows correspond to different sample sizes
d <- c(50, 100, 200, 400) # Define the sample sizes
# Loop through different sample sizes (50, 100, 200, 400)
for (i in 1:4) {
  m <- d[i] # Select sample size
  y <- data[, 1:m] # Subset data for the selected sample size

  # Fit GLLVM with ZINB distribution
  model <- gllvm(y, family = "ZINB", num.lv = 2, row.eff = "fixed", sd.errors = FALSE, seed = 123)
  pred_y <- predict.gllvm(model, type = "response") # Predict values from the fitted model

  # 1. MARNE : Mean Absolute Range Normalized Error
  GOF_gllvmZINB[i, 1] <-  goodnessOfFit(y = y, pred = pred_y,  measure = c("MARNE"), species = FALSE)$MARNE
  
  # 2. Species-wise Correlations: mCOR
  GOF_gllvmZINB[i, 2] <- mean(goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = TRUE)$cor)

  # 3. Overall Correlation (whole data): cor
 
  GOF_gllvmZINB[i, 3] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = FALSE)$cor
}

# Print results for GLLVM with ZINB distribution
print(GOF_gllvmZINB)

# 2. GLLVM with NB distribution

GOF_gllvmNB <- matrix(NA, 4, 3) # Initialize results matrix for NB model
colnames(GOF_gllvmNB) <- c("marne", "mCOR", "cor")
rownames(GOF_gllvmNB) <- c("50", "100", "200", "400")

# Loop through different sample sizes (50, 100, 200, 400)
for (i in 1:4) {
  m <- d[i]
  y <- data[, 1:m]

  # Fit GLLVM with negative binomial(NB)  distribution
  model <- gllvm(y, family = "negative.binomial", num.lv = 2, row.eff = "fixed", sd.errors = FALSE, seed = 123)
  pred_y <- predict.gllvm(model, type = "response") # Predict values from the fitted model

  # 1. MARNE : Mean Absolute Range Normalized Error
  GOF_gllvmNB[i, 1] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("MARNE"), species = FALSE)$MARNE

  # 2. Species-wise Correlations: mCOR

  GOF_gllvmNB[i, 2] <-mean(goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = TRUE)$cor)

  # 3. Overall Correlation (whole data): cor

  GOF_gllvmNB[i, 3] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = FALSE)$cor
}

# Print results for GLLVM with NB distribution
print(GOF_gllvmNB)

# 3. Copula-based method with ZINB distribution

GOF_cZINB <- matrix(NA, 4, 3) # Initialize results matrix for copula ZINB model
colnames(GOF_cZINB) <- c("marne", "mCOR", "cor")
rownames(GOF_cZINB) <- c("50", "100", "200", "400")

# Loop through different sample sizes (50, 100, 200, 400)
for (i in 1:4) {
  m <- d[i]
  y <- data[, 1:m]

  # Fit Copula model with ZINB distribution
  model <- fit_copula(y, reff = "fixed", gllvm.fam = "ZINB", sd.errors = FALSE, seed = 123, lv.n = 0)
  true.ords <- model$scores # Extract latent variables (scores)
  true.load <- model$loadings # Extract factor loadings
  eta <- t(replicate(nrow(y), model$obj$params$beta0)) + c(0, model$obj$params$row.params) # Linear predictor

  # Predict values using the copula model
  pred_y <- exp(eta + (true.ords %*% t(true.load))) # Exponentiate to get predictions

  # 1. MARNE : Mean Absolute Range Normalized Error
  GOF_cZINB[i, 1] <-  goodnessOfFit(y = y, pred = pred_y,  measure = c("MARNE"), species = FALSE)$MARNE

  # 2. Species-wise Correlations: mCOR
  GOF_cZINB[i, 2] <- mean(goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = TRUE)$cor)

  # 3. Overall Correlation (whole data): cor
  
  GOF_cZINB[i, 3] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = FALSE)$cor
}

# Print results for Copula with ZINB distribution
print(GOF_cZINB)

# 4. Copula-based method with NB distribution

GOF_cNB <- matrix(NA, 4, 3) # Initialize results matrix for copula NB model
colnames(GOF_cNB) <- c("marne", "mCOR", "cor")
rownames(GOF_cNB) <- c("50", "100", "200", "400")

# Loop through different sample sizes (50, 100, 200, 400)
for (i in 1:4) {
  m <- d[i]
  y <- data[, 1:m]

  # Fit Copula model with NB distribution
  model <- fit_copula(y, reff = "fixed", gllvm.fam = "negative.binomial", sd.errors = FALSE, seed = 123, lv.n = 0)
  true.ords <- model$scores # Extract latent variables (scores)
  true.load <- model$loadings # Extract factor loadings
  eta <- t(replicate(nrow(y), model$obj$params$beta0)) + c(0, model$obj$params$row.params) # Linear predictor

  # Predict values using the copula model
  pred_y <- exp(eta + (true.ords %*% t(true.load))) # Exponentiate to get predictions

  # 1. MARNE : Mean Absolute Range Normalized Error
  GOF_cNB[i, 1] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("MARNE"), species = FALSE)$MARNE

  # 2. Species-wise Correlations: mCOR
  GOF_cNB[i, 2] <- mean(goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = TRUE)$cor)

  # 3. Overall Correlation (whole data): cor
  
  GOF_cNB[i, 3] <- goodnessOfFit(y = y, pred = pred_y,  measure = c("cor"), species = FALSE)$cor
}

# Print results for Copula with NB distribution
print(GOF_cNB)
