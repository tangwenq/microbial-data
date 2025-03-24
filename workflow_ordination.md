Workflow illustration to compare ordination methods in the analysis of
high-dimensional compositional count data using a  microbiome data set
================
Wenqi Tang, Pekka Korhonen, Jenni Niku, Klaus Nordhausen and Sara Taskinen
2025-03-21

This document reproduces the main content of the introduction and
simulation section of the paper **Comparing model-based unconstrained
ordination methods in the analysis of high-dimensional compositional
count data** by Wenqi Tang, Pekka Korhonen, Jenni Niku, Klaus Nordhausen
and Sara Taskinen. It includes microbiome data visualization and example
about how to handle compositional data and apply ordination methods in R
within the simulation section in paper. We primarily introduce the four
distinct ordination methods mentioned in the paper, which are applicable
to high-dimensional compositional count data. Specifically, we explore
classical ordination methods, Generalized Linear Latent Variable
Models(GLLVMs)-based ordination methods, and copula-based(GCLVM)
ordination methods, applying them to microbiome data in R package
`gllvm` and detailed introduction of dataset is in Kumar et al. (2017).
The dataset used in this analysis consists of a 56 × 985 matrix,
representing 985 bacterial species sampled from 56 soil sites across
three distinct climatic regions: NyÅlesund (high Arctic), Kilpisjärvi
(low Arctic), and Mayrhofen (European Alps) Kumar et al. (2017). Due to
biological or technical factors, such datasets often contain a
significant proportion of zero observations (e.g., 63% zeros in the
microbiome data). Additionally, microbiome data exhibit a compositional
nature, as highlighted in Gloor et al. (2016), meaning the data
represent relative abundances rather than absolute counts.

## Data preprocess

We preprocess data by ordering the columns according to sparseness and
use a heatmap to visualize the sparsity of the data and a mean-variance
plot to assess the overdisperison of the data (Figure 1).

``` r
library(robCompositions)
library(vegan)
library(devtools) # Load devtools for package installation
#devtools::install_github("JenniNiku/gllvm") # Install the latest gllvm package from GitHub
library(gllvm) # Load the gllvm package
source_url("https://raw.githubusercontent.com/tangwenq/microbial-data/main/GCLVM.R")
library(reshape2) # For data reshaping
library(scales)
library(TMB)

# Base the simulations to arctic microbialdata from Nissinen et al.
# Data included in gllvm package

data("microbialdata")
data <- microbialdata$Y
# order columns according to sparseness
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum))
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]
data <- data[, ind]
dim(data)
```

    ## [1]  56 985

``` r
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
```

    ## Optimal dispersion parameter (theta): 1.114091

``` r
# Plot mean-variance relationship
plot(means, variances,
  xlab = "Mean", # X-axis label
  ylab = "Variance", # Y-axis label
  pch = 16,
  cex = 0.7, # Point size
  col = "black",
  xlim = c(0, 30), # Limit for X-axis
  ylim = c(0, 1000), # Limit for Y-axis
  cex.lab = 1.5, # Label font size
  cex.axis = 1.4, # Axis font size
  cex.main = 1.6
) # Title font size

# Add a Poisson reference line (variance = mean)
abline(a = 0, b = 1, col = "blue", lwd = 2) # Poisson variance line

# Add the Negative Binomial relationship line
theta <- theta_opt # Use the estimated dispersion parameter
nb_variance <- means + means^2 * theta # Compute NB variance
lines(sort(means), sort(nb_variance),
  col = "blue", lwd = 2, lty = 2
) # Add the NB variance line (dashed)
```

<img src="workflow_ordination_files/figure-gfm/data-1.png" style="display: block; margin: auto;" />

``` r
# Draw heatmap of data, (shown in left picture of Figure 1)

# Convert data to long format for plotting
long_data <- melt(data, varnames = c("Site", "Species"), value.name = "Richness")


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
    axis.text.x = element_text(size = 14, face = "bold", angle = 90, hjust = 1), # Bold x-axis text
    axis.text.y = element_blank(), # Remove y-axis labels
    axis.ticks = element_blank(), # Remove axis ticks
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  scale_x_discrete(
    breaks = c(50, 100, 200, 400), # Specify discrete breaks if Species are numbered
    labels = c("50", "100", "200", "400") # Custom labels
  ) +
  geom_vline(
    xintercept = c(50, 100, 200, 400),
    color = "black",
    linetype = "solid", # Use solid lines for emphasis m = 50, 100, 200, 400
    size = 1
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 32,
      barheight = 2,
      breaks = seq(0, 250, by = 50)
    )
  )
```

<img src="workflow_ordination_files/figure-gfm/data-2.png" style="display: block; margin: auto;" />

## Classical ordination

Classical methods for unconstrained ordination, such as **non-metric
multidimensional scaling (nMDS)** and **principal component analysis
(PCA)**, are widely used due to their computational efficiency and ease
of implementation. Compositional data are subject to sum constraints,
require **centered log-ratio (clr) transformation** in R package
`robCompositions` to map compositional data from the Aitchison space to
real space before applying these methods. For nMDS, the Euclidean
distance can be used to evaluate the clr-transformed data by specifying
`distance = "euclidean"`. Alternatively, for the original compositional
data, dissimilarity measures such as the Bray-Curtis distance
`distance = "bray"` can be applied using the R package `vegan`.

## Example

In the following analysis, we demonstrate the application of classical
ordination methods to compositional data.

``` r
y = data # full dataset as example
# Distance-based methods: nMDS

    fit_mds <- metaMDS(y, distance = "bray", trace = FALSE)
    mds_ords <- scale(fit_mds$ points) # extract ordinations from nMDS
    
    # PCA and nMDS on CLR-transformed data
    clr_y <- cenLR(y + 1) # CLR transformation of data with added 1

    fit_clrmds <- metaMDS(clr_y$x.clr, distance = "euclidean", autotransform = FALSE, noshare = FALSE, wascores = FALSE, trace = FALSE)
    clrmds_ords <- scale(fit_clrmds$points)
    
    fit_pca <- prcomp(clr_y$x.clr, scale = TRUE) # PCA on the CLR-transformed data clr_y
    pca_ords <- scale(fit_pca$x[, 1:2]) # extract ordinations from PCA
```

The microbiome data was collected from three regions: “Mayrhofen,”
“Kilpisjarvi,” and “Ny-Alesund.” To analyze the community structure of
these data, we can visualize the ordinations extracted using classical
methods by plotting ordination plots. The following code demonstrates
how to perform ordination analysis on CLR-transformed data using PCA as
an example and visualize the results. We can find that classical methods
are lack of a probabilistic framework, challenges in handling zero
values.

``` r
unique_pch <- c(16, 17, 18)       # Mayrhofen=16, Kilpisjarvi=17, Ny-Alesund=18
unique_colors <- c("black", "blue", "grey") 

plot_pch <- unique_pch[as.numeric(factor(microbialdata$X$Region))]
plot_colors <- unique_colors[as.numeric(factor(microbialdata$X$Region))]


# draw ordination plot for nMDS
plot(clrmds_ords[, 1], clrmds_ords[, 2],
     col = plot_colors,
     pch = plot_pch,
     cex = 1.2,
     xlab = "ordination score 1",
     ylab = "ordination score 2",
     main = "nMDS of CLR-transformed Data")


legend("topleft",                      
       legend = c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund"),       
       pch = unique_pch,               
       col = unique_colors,                    
       ncol = 3,                        
       cex = 1)



```

<img src="workflow_ordination_files/figure-gfm/classical_vis-1.png" style="display: block; margin: auto;" />

``` r
# draw ordination plot for PCA
plot(pca_ords[, 1], pca_ords[, 2],
     col = plot_colors,
     pch = plot_pch,
     cex = 1.2,
     xlab = "ordination score 1",
     ylab = "ordination score 2",
     main = "PCA of CLR-transformed Data")


legend("topleft",                      
       legend = c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund"),       
       pch = unique_pch,               
       col = unique_colors,                    
       ncol = 3,                        
       cex = 1)
```

<img src="workflow_ordination_files/figure-gfm/classical_vis-2.png" style="display: block; margin: auto;" />

Results from different regions were represented using distinct colors
and symbols. From the ordination plots of PCA and nMDS, it is evident
that the sample points from the three regions are highly intermixed,
showing no clear clustering or separation trends. This indicates that
these methods fail to effectively capture any structure related to the
three regions. Specifically, under conditions of high sparsity,
classical ordination methods such as PCA and nMDS are unable to
adequately reveal the underlying regional patterns.

## Ordination based on latent variable models (GLLVM)

Generalized linear latent variable models (GLLVMs) provide a flexible
framework for jointly modeling multivariate response data, such as
microbial count data. GLLVMs extend generalized linear models (GLMs) by
incorporating correlations among response variables using a
factor-analytic approach. This allows GLLVMs to be used for model-based
ordination for any response type. Assume
![\boldsymbol{u}\_i = (u\_{i1}, \dots, u\_{id})^\top](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bu%7D_i%20%3D%20%28u_%7Bi1%7D%2C%20%5Cdots%2C%20u_%7Bid%7D%29%5E%5Ctop "\boldsymbol{u}_i = (u_{i1}, \dots, u_{id})^\top")
is a
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")-dimensional
latent variable following a standard multivariate normal distribution.
In GLLVMs, it is assumed that, conditional on the latent variables
![\boldsymbol{u}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bu%7D_i "\boldsymbol{u}_i"),
the response variables
![y\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bij%7D "y_{ij}")
are independently distributed according to some distribution
![F(\mu\_{ij}, \boldsymbol{\phi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F%28%5Cmu_%7Bij%7D%2C%20%5Cboldsymbol%7B%5Cphi%7D%29 "F(\mu_{ij}, \boldsymbol{\phi})"),
where
![\mu\_{ij} = \mathbb{E}(y\_{ij} \| \boldsymbol{u}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7Bij%7D%20%3D%20%5Cmathbb%7BE%7D%28y_%7Bij%7D%20%7C%20%5Cboldsymbol%7Bu%7D_i%29 "\mu_{ij} = \mathbb{E}(y_{ij} | \boldsymbol{u}_i)")
is the conditional mean, and
![\boldsymbol{\phi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cphi%7D "\boldsymbol{\phi}")
includes possible response-specific parameters (e.g., dispersion or
zero-inflation parameters). When GLLVMs are used for model-based
ordination,
![\mu\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7Bij%7D "\mu_{ij}")
is linked to the linear predictor via:

![g(\mu\_{ij}) = \alpha_i + \beta\_{0j} + \boldsymbol{\lambda}\_j^\top \boldsymbol{u}\_i,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Cmu_%7Bij%7D%29%20%3D%20%5Calpha_i%20%2B%20%5Cbeta_%7B0j%7D%20%2B%20%5Cboldsymbol%7B%5Clambda%7D_j%5E%5Ctop%20%5Cboldsymbol%7Bu%7D_i%2C "g(\mu_{ij}) = \alpha_i + \beta_{0j} + \boldsymbol{\lambda}_j^\top \boldsymbol{u}_i,")

where:

- ![g(\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Ccdot%29 "g(\cdot)")
  is a known link function (typically the log-link for count data).
- ![\beta\_{0j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B0j%7D "\beta_{0j}")
  is a column-specific intercept to account for differences in column
  totals.
- ![\boldsymbol{\lambda}\_j = (\lambda\_{j1}, \dots, \lambda\_{jd})^\top](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Clambda%7D_j%20%3D%20%28%5Clambda_%7Bj1%7D%2C%20%5Cdots%2C%20%5Clambda_%7Bjd%7D%29%5E%5Ctop "\boldsymbol{\lambda}_j = (\lambda_{j1}, \dots, \lambda_{jd})^\top")
  is
  a![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
  -dimensional vector of factor loadings.
- ![\alpha_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_i "\alpha_i")
  is a row-specific intercept to account for differences in row totals.

To ensure model identifiability, the upper triangular part of the factor
loading matrix
![\boldsymbol{\Lambda} = \[\boldsymbol{\lambda}\_1 \cdots \boldsymbol{\lambda}\_m\]^\top](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5CLambda%7D%20%3D%20%5B%5Cboldsymbol%7B%5Clambda%7D_1%20%5Ccdots%20%5Cboldsymbol%7B%5Clambda%7D_m%5D%5E%5Ctop "\boldsymbol{\Lambda} = [\boldsymbol{\lambda}_1 \cdots \boldsymbol{\lambda}_m]^\top")
is set to zero, and its diagonal elements are set to positive values.
The `gllvm` package implements the **zero-inflated negative binomial
(ZINB)** (`family = "ZINB"`) and **negative binomial (NB)**
(`family = "negative.binomial"`) GLLVMs to better handle datasets with a
large number of zeros. These models are particularly suited for
overdispersed count data, where excess zeros are common. For non-normal
response data, estimation is typically based on approximate marginal
likelihood using variational approximations (VA) or Laplace’s
approximation (LA).

## Example

It should be noted that the following code takes considerable time to
run, and we can speed up computations by enabling parallelization for
`gllvm` with the following command:
`TMB::openmp(n = parallel::detectCores() - 1, autopar = TRUE, DLL = "gllvm")`.
This command should be executed before fitting the model with
`gllvm()`.However, due to the algorithmic nature of the methods used,
enabling parallelization or using different versions that might
initialize parameters differently can lead to slight variations in the
final results. Additionally, as the `gllvm` package on GitHub is
continuously updated, installing the latest version via
`devtools::install_github("JenniNiku/gllvm")` may also introduce minor
differences. Despite these small variations, the overall results and
conclusions will remain consistent.

To fit the model:

``` r
# Fit GLLVMs with NB distribution with full data
# Extract and scale environmental covariates
X <- scale(microbialdata$Xenv[, c("SOM", "pH", "Phosp")])

# Uncomment to fit unconstrained model
#ungllvm_NB <- gllvm(
  #y = y, family = "negative.binomial", sd.errors = FALSE,
 # row.eff = "fixed", num.lv = 2, seed = 123
#)


# Constrained concurrent model
gllvm_NB <- gllvm(
  y = y, X=X,family = "negative.binomial", sd.errors = FALSE,
  row.eff = "fixed", num.lv.c = 2, seed = 123
)

scale(gllvm_NB$lvs) # extract ordinations
```

    ##             CLV1        CLV2
    ## AB2   0.26390421  0.19414480
    ## AB3   0.56364574 -0.45151891
    ## AB4   1.62391376 -0.18134564
    ## AB5   0.15627679 -1.36471469
    ## AT2   1.08875527 -0.02000472
    ## AT3   0.38455059 -0.83527610
    ## AT4   2.47477520  3.80100255
    ## AT5   1.62946930 -1.49061208
    ## CB2   1.74084310 -0.35200274
    ## CB3   1.66632776 -0.73071040
    ## CB4   1.42139411 -0.29505808
    ## CB5   1.10390801 -0.22121166
    ## JNB2 -0.77817423  1.28315681
    ## JNB3 -0.91167096  1.42826205
    ## JNB4 -1.03562918  1.39319531
    ## JNB5 -1.15865116  1.37865768
    ## JNT2 -0.05931051  0.45277657
    ## JNT3 -0.37874992  0.99715193
    ## JNT4 -0.63346301  0.84379640
    ## JNT5 -0.53881968  0.88120715
    ## JOB2 -2.03664195 -1.34072554
    ## JOB3 -1.31571823 -1.18588107
    ## JOB4 -1.53489148 -0.92974817
    ## JOB5 -0.71582063 -0.63951061
    ## JOT2  0.16619416 -1.00449510
    ## JOT3 -0.28065487 -1.01984326
    ## JOT4 -0.84716361 -0.20937951
    ## JOT5  0.08522503 -0.59217662
    ## KB2  -0.38519603  0.32517964
    ## KB3  -0.64533885 -0.20837786
    ## KB4  -1.07717807 -0.14760876
    ## KB5  -0.84431895 -0.04510399
    ## KT2   0.24097702  0.41458609
    ## KT3   0.60369817 -0.99300030
    ## KT4  -0.89002702 -1.36808813
    ## KT5   0.51291204 -1.47459452
    ## MB3  -0.38712704 -0.31189229
    ## MB4   0.66903928  0.10067609
    ## MB5   1.52157263  0.63241784
    ## MT3   0.37224930 -0.77222618
    ## MT4   1.70647726 -0.40534888
    ## MT5   0.60851315  1.29389120
    ## RRB2  0.11835907 -0.86170021
    ## RRB3 -0.62129574  0.38303645
    ## RRB4 -0.57434819  0.11786402
    ## RRB5 -1.31146222 -0.12168112
    ## RRT2 -0.19625430 -0.42890544
    ## RRT3 -0.24153547  0.36829809
    ## RRT4 -0.29067195 -0.13070128
    ## RRT5 -0.21096114 -0.89805273
    ## SB3  -1.54082736  1.14615595
    ## SB4  -1.12676959  1.50516406
    ## SB5  -0.42958645  1.72826575
    ## ST3   0.92996630 -0.47092264
    ## ST4   0.67158228  0.25981369
    ## ST5   0.67372825  0.57371910
    ## attr(,"scaled:center")
    ##         CLV1         CLV2 
    ##  0.002167298 -0.018585528 
    ## attr(,"scaled:scale")
    ##     CLV1     CLV2 
    ## 1.046706 1.009262

Through the `summary()`, `AIC()` or `BIC()` functions, we can observe
the AIC and BIC values of the model, which can also be found in **Table
2** of the paper.

``` r
AIC(gllvm_NB)
```

    ## [1] 126647.7

``` r
BIC(gllvm_NB)
```

    ## [1] 162310.8

For model-based methods, we can use the following code in `gllvm` to
generate residuals plots and ordination plots for analysis, assessing
the model’s fitting performance.

``` r
# Visualize residuals for  NB-GLLVM model(reproduces Figure 8)


# the D-S of concurrent NB-GLLVM model, shown in bottom row in Figure 8
par(mar = c(6, 6, 0, 0) + 0.1) # Set plot margins
plot(gllvm_NB,
  which = 1, #when which = 1, it produces a residuals plot
  caption = " ",
  var.colors = 1,
  cex.lab = 2.5,
  cex.axis = 1.5,
  cex.main = 2.5,
  cex = 1
)
```

<img src="workflow_ordination_files/figure-gfm/GLLVM_vis-1.png" style="display: block; margin: auto;" />

``` r
par(mar = c(6, 6, 0, 0) + 0.1) # Set plot margins
# Q-Q plot of residuals
plot(gllvm_NB,
  which = 2, #when which = 2, it generates a QQ norm plot
  caption = " ", 
  var.colors = 1, 
  cex.lab = 2.5, 
  cex.axis = 1.5, 
  cex.main = 2.5, 
  cex = 1
) 
```

<img src="workflow_ordination_files/figure-gfm/GLLVM_vis-2.png" style="display: block; margin: auto;" />

``` r
# Ordination plots models (reproduces Figure 5)
ordiplot(gllvm_NB,
  which.lvs = 1:2,
  s.colors = plot_colors,
  rotate = TRUE,
  symbols = TRUE,
  pch = plot_pch,
  main = "",
  ann = FALSE
) # Set symbols based on region

title(xlab = "ordination score 1", ylab = "ordination score 2", cex.lab = 1.3)
# Add a legend with specified colors and symbols
legend("topleft", # Legend position
  legend = c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund"), # Legend labels
  pch = unique_pch, # Corresponding symbols
  col = unique_colors, # Corresponding colors (black, blue, grey)
  ncol = 3,
  cex = 1
) # Arrange legend in 3 columns
```

<img src="workflow_ordination_files/figure-gfm/GLLVM_vis-3.png" style="display: block; margin: auto;" />

In ordination plot, longer arrows represent covariates with the largest
relative effects, and dark red arrows (here associated to pH) indicate
covariates with a significant effect to ordination.

## Ordination based on couplas (GCLVMs)

Copula models couple marginal models for the data with a multivariate
model that accounts for covariance across responses. Assume
![y\_{ij} \sim F_j(\mu\_{ij}, \boldsymbol{\phi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bij%7D%20%5Csim%20F_j%28%5Cmu_%7Bij%7D%2C%20%5Cboldsymbol%7B%5Cphi%7D%29 "y_{ij} \sim F_j(\mu_{ij}, \boldsymbol{\phi})"),
where
![\mu\_{ij} = \mathbb{E}(y\_{ij})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7Bij%7D%20%3D%20%5Cmathbb%7BE%7D%28y_%7Bij%7D%29 "\mu_{ij} = \mathbb{E}(y_{ij})"),
and
![\boldsymbol{\phi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cphi%7D "\boldsymbol{\phi}")
includes possible response-specific parameters. The mean
![\mu\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7Bij%7D "\mu_{ij}")
is linked to the linear predictor via generalized linear models (GLMs):

![g(\mu\_{ij}) = \alpha_i + \beta\_{0j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Cmu_%7Bij%7D%29%20%3D%20%5Calpha_i%20%2B%20%5Cbeta_%7B0j%7D "g(\mu_{ij}) = \alpha_i + \beta_{0j}")

where
![g(\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Ccdot%29 "g(\cdot)")
is a known link function,
![\alpha_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_i "\alpha_i")
and
![\beta\_{0j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B0j%7D "\beta_{0j}")
are row-specific and column-specific intercepts,respectively. In the
Gaussian copula model, count data
![y\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bij%7D "y_{ij}")
are mapped to copula values
![z\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bij%7D "z_{ij}")
that follow a multivariate normal distribution:

![F_j(y\_{ij} - 1) \leq \Phi(z\_{ij}) \< F_j(y\_{ij})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_j%28y_%7Bij%7D%20-%201%29%20%5Cleq%20%5CPhi%28z_%7Bij%7D%29%20%3C%20F_j%28y_%7Bij%7D%29 "F_j(y_{ij} - 1) \leq \Phi(z_{ij}) < F_j(y_{ij})")

where
![F_j(\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_j%28%5Ccdot%29 "F_j(\cdot)")
is the cumulative distribution function (cdf) assumed for the
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")-th
column in the data matrix under the marginal GLM, and
![\Phi(\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi%28%5Ccdot%29 "\Phi(\cdot)")
is the cdf of the standard normal distribution. For ordination analysis,
the copula values
![z\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bij%7D "z_{ij}")
are assumed to follow a factor-analytic model:

![z\_{ij} = \boldsymbol{\lambda}\_j^\top \boldsymbol{u}\_i + \epsilon\_{ij},](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bij%7D%20%3D%20%5Cboldsymbol%7B%5Clambda%7D_j%5E%5Ctop%20%5Cboldsymbol%7Bu%7D_i%20%2B%20%5Cepsilon_%7Bij%7D%2C "z_{ij} = \boldsymbol{\lambda}_j^\top \boldsymbol{u}_i + \epsilon_{ij},")

where:
![\boldsymbol{u}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bu%7D_i "\boldsymbol{u}_i")
is a
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")-dimensional
latent variable associated with the study unit.
![\boldsymbol{\lambda}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Clambda%7D_j "\boldsymbol{\lambda}_j")
is a
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")-dimensional
vector of factor loadings.
![\epsilon\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bij%7D "\epsilon_{ij}")
are independent Gaussian errors with variances
![\sigma_j^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_j%5E2 "\sigma_j^2").
The GCLVM data is generated based on modified code from Popovic (2021).

The parameters of the copula model are estimated using a two-step
procedure:

1.  Estimate the marginal distributions
    ![F_j(\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_j%28%5Ccdot%29 "F_j(\cdot)")
    using GLMs suitable for sparse, overdispersed count data.
2.  Use Monte Carlo expectation-maximization (MCEM) to estimate the
    covariance parameters in the copula model.

In place of marginal GLMs, due to the compositionality constraint, for
the base model, we use a multivariate GLM with row effects shared by the
different species in the data, estimated with the `gllvm` package using
the options `num.lv=0` and `row.eff="fixed`. Alternatively, one could
include the row effects by treating the data in the long format instead.


## Example

``` r
 # Uncomment for copula with ZINB distribution
    #c_ZINB <- fit_copula(y, gllvm.fam = "ZINB", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0)

  
  # Copula with NB distribtuion
    c_NB <- fit_copula(y, gllvm.fam = "negative.binomial", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0)
    
    c_NB_ords <- scale(c_NB$scores)# extract ordinations
c_NB_ords
```

    ##       
    ## X1        Factor1     Factor2
    ##   AB2   0.5086746 -0.40566805
    ##   AB3   0.3910097 -0.80810194
    ##   AB4   0.1120198 -1.68869799
    ##   AB5  -0.5901694 -0.85296243
    ##   AT2   1.1239948 -0.72255923
    ##   AT3   0.9562300 -0.99152217
    ##   AT4   1.3290799 -1.09588068
    ##   AT5   0.8605684 -1.52054317
    ##   CB2   1.8956152 -0.20279141
    ##   CB3   1.6794126 -0.04431647
    ##   CB4   1.8605660 -0.12784251
    ##   CB5   1.2186536  0.40141803
    ##   JNB2 -0.5104550  1.75144752
    ##   JNB3 -0.3644585  1.98641127
    ##   JNB4 -0.3182096  2.02576774
    ##   JNB5 -0.1808086  1.89219234
    ##   JNT2  0.1579177  1.32583888
    ##   JNT3  0.1548659  1.63914470
    ##   JNT4  0.2904853  1.62569651
    ##   JNT5  0.1907263  1.60763100
    ##   JOB2 -1.5145856  0.34481163
    ##   JOB3 -1.1408350  0.21751302
    ##   JOB4 -1.4985251  0.42707535
    ##   JOB5 -0.9589408  0.47005208
    ##   JOT2 -0.4184319 -0.40945170
    ##   JOT3 -0.4061251  0.06905622
    ##   JOT4 -0.4360063  0.52001080
    ##   JOT5  0.1381568  0.45971832
    ##   KB2  -1.1620381 -0.81350427
    ##   KB3  -0.7421644 -0.39881362
    ##   KB4  -1.0688430 -0.41773187
    ##   KB5  -1.0617087 -0.23077987
    ##   KT2  -0.6822235 -1.05316213
    ##   KT3   0.4274158 -0.54523074
    ##   KT4  -0.7820214 -0.56303634
    ##   KT5  -0.1393242 -1.05054080
    ##   MB3  -1.2531440 -0.70397574
    ##   MB4   0.6819726  0.11837939
    ##   MB5   1.5752083 -0.14361741
    ##   MT3  -0.5794767 -0.99746364
    ##   MT4   1.9070709 -0.71379382
    ##   MT5   1.7495239 -0.36437683
    ##   RRB2 -0.3918442 -0.74132585
    ##   RRB3 -1.1461971 -0.87157856
    ##   RRB4 -1.2239371 -0.68427180
    ##   RRB5 -1.3619212 -0.49149028
    ##   RRT2 -0.7148618 -0.05288195
    ##   RRT3 -0.8716803 -1.01064236
    ##   RRT4 -0.7736036 -1.19159798
    ##   RRT5 -0.7377571 -1.06169793
    ##   SB3  -0.2822288  1.27885145
    ##   SB4  -0.2679959  1.79841167
    ##   SB5   0.3327795  1.52382817
    ##   ST3   1.5561150  0.20241174
    ##   ST4   0.9442886  0.84033500
    ##   ST5   1.5381707  0.44584871
    ## attr(,"scaled:center")
    ##      Factor1      Factor2 
    ##  0.008505056 -0.017788075 
    ## attr(,"scaled:scale")
    ##  Factor1  Factor2 
    ## 1.006406 1.005247

we can use following code to draw ordination plot and extract model to
check AIC and BIC.

``` r
obj = c_NB$obj #extract model
obj
```

    ## Call: 
    ## gllvm(y = data, family = gllvm.fam, num.lv = lv.n, row.eff = reff, 
    ##     sd.errors = FALSE)
    ## family: 
    ## [1] "negative.binomial"
    ## method: 
    ## [1] "VA"
    ## 
    ## log-likelihood:  -70566.6 
    ## Residual degrees of freedom:  53135 
    ## AIC:  145183.2 
    ## AICc:  145337.6 
    ## BIC:  163242.1

``` r
# draw ordination plot for NB-GCLVM


# draw ordination plot for PCA
plot(c_NB_ords[, 1], c_NB_ords[, 2],
     col = plot_colors,
     pch = plot_pch,
     cex = 1.2,
     xlab = "ordination score 1",
     ylab = "ordination score 2",
     main = "NB_GCLVM ordination plot")


legend("topleft",                      
       legend = c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund"),       
       pch = unique_pch,               
       col = unique_colors,                    
       ncol = 3,                        
       cex = 1)


```

<img src="workflow_ordination_files/figure-gfm/GCLVM_ana-1.png" style="display: block; margin: auto;" />

From the ordination plot, it can be observed that while the NB-GCLVM
exhibits a stronger tendency for classification compared to PCA and
nMDS, it still performs less effectively when compared to the NB-GLLVM
visualization.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Glooretal2016" class="csl-entry">

Gloor, Gregory B., Jia Rong Wu, Vera Pawlowsky-Glahn, and Juan José
Egozcue. 2016. “It’s All Relative—Analyzing Microbiome Data as
Compositions” 26: 322–29.
<https://doi.org/10.1016/j.annepidem.2016.03.003>.

</div>

<div id="ref-Kumar2017" class="csl-entry">

Kumar, Manoj, Gerald Brader, Angela Sessitsch, Minna Mäki, Jan Dirk van
Elsas, and Riikka Nissinen. 2017. “Plants Assemble Species-Specific
Bacterial Communities from Common Core Taxa in Three Arcto-Alpine
Climate Zones” 8: 12. <https://doi.org/10.3389/fmicb.2017.00012>.

</div>

<div id="ref-Popovic2021" class="csl-entry">

Popovic, Gordana. 2021. “Fast Model-Based Ordination with
Copulas—Simulation Code (V1.0.0).”
<https://doi.org/10.5281/zenodo.5525716>.

</div>

</div>
