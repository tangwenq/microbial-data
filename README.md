# **Comparing model-based unconstrained ordination methods in the analysis of high-dimensional compositional count data**

This project contains code to reproduce the results of the paper with the same name by Tang, W., Korhonen, P., Niku, J., Nordhausen, K. and Taskinen, S.

The goal of the project was to compare different unconstraint ordination methods for high-dimensional compositional count data.

**Main functions**

All code here is written in R and requires the packages vegan, gllvm, ggplot2, mvtnorm, devtools, mvabund, ecoCoupla, distribution3, robCompositions, reshape2, scales and
Source of copula code.

The main files are:

- `GCLVM.R`: main functions in the file are `fit_copula` which performs model-based ordination along the lines of Popovic at al. (2022) and `simulate.copula` which simulates responses from a copula model fitted with `fit_copula`. The functions are modified versions of corresponding functions taken from https://github.com/gordy2x/ecoCopula and from the appendix of Popovic at al. (2022) available at https://doi.org/10.5281/zenodo.5525716.
- `simulation_gllvm.R`: The main function in this file is `sim_gllvm`, which simulates new response based on a fitted ZINB-GLLVM model with `gllvm` function. It applies seven ordination methods to the simulated data, extracts the estimated ordination scores, and calculates the Procrustes distance between the estimated and true ordinations extracted from `fit_gllvm`. The file conducts simulations using microbial data with 50, 100, 200, and 400 dimensions from the `gllvm` package and visualizes the results using boxplots(Figure 3).
- `simulation_copula.R` :The main function in this file is `sim_copula`, which simulates new response based on a fitted ZINB-GCLVM model with `fit_copula` function from `GCLVM` file. It applies seven ordination methods to the simulated data, extracts the estimated ordination scores, and calculates the Procrustes distance between the estimated and true ordinations extracted from `fit_ZINBcopula`. The file conducts simulations using microbial data with 50, 100, 200, and 400 dimensions from the `gllvm` package and visualizes the results using boxplots(Figure 4).
- `fit-check`: This file evaluates the goodness-of-fit of ZINB-GLLVM, NB-GLLVM, ZINB-GCLVM, and NB-GCLVM models across dimensions 50, 100, 200, and 400 of microbial data using the three measures: MAREN, mCOR, and gCOR (Table 1).
- `real data analysis`：This file uses the `gllvm` function to build unconstrained and concurrent GLLVM models incorporating three environmental variables (SOM, pH, and Phosphorus) with both NB and ZINB distributions (Table 2) using the full microbial dataset. Based on unconstrained and concurrent NB-GLLVM，generate D-S residual plots (Figure 6) are generated ordination plots (Figure 5) with `ordiplot` function from the `gllvm` package.
- `visulization of microbiome data`: The primary purpose of this file is to visualize the sparsity and overdispersion of microbial data using heatmaps and mean-variance plots (Figure 1).
- `appendix`: This file performs simulations when NB-GLLVM and NB-GCLVM are true models and visualizes the Procrustes distances using boxplots (Figures 7 and 8).

**Authors**

Tang, W., Korhonen, P., Niku, J., Nordhausen, K. and Taskinen, S.

**License**

LGPL >= 2.1

**References**

Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2014). Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6, 399–411.

Niku, J., Hui, F. K. C., Taskinen, S., and Warton, D. I. (2019). `gllvm` - Fast analysis of multivariate abundance data with generalized linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173–218. 

Popovic, G. C., Hui, F. K. C., and Warton, D. I. (2022). Fast model-based ordination with copulas. Methods in Ecology and Evolution, 13, 194–202.


