# **Comparing model-based unconstrained ordination methods in the analysis of high-dimensional compositional count data**

This project contains code to reproduce the results of the paper with the same name by Tang, W., Korhonen, P., Niku, J., Nordhausen, K. and Taskinen, S.

The goal of the project was to compare different unconstraint ordination methods for high-dimensional compositional count data.

**Main functions**

All code here is written in R and requires the packages vegan, gllvm, ggplot2, mvtnorm, ...
Source of copula code

The main files are:

- `clvm.R`: main functions in the file are `clvm` which performs model-based ordination along the lines of Popovic at al. (2022) and `simulate.copula` which simulates responses from a copula model fitted with `clvm`. The functions are modified versions of corresponding functions taken from https://github.com/gordy2x/ecoCopula and from the appendix of Popovic at al. (2022) available at https://doi.org/10.5281/zenodo.5525716.
- `simulation_copula.R`
- `simulation_gllvm.R`
- `visualization.R`

**Authors**

Tang, W., Korhonen, P., Niku, J., Nordhausen, K. and Taskinen, S.

**License**

LGPL >= 2.1

**References**

Niku, J., Hui, F. K. C., Taskinen, S., and Warton, D. I. (2019). `gllvm` - Fast analysis of multivariate abundance data with generalized linear latent variable models in R. Methods in Ecology and Evolution, 10, 2173–218. 

Popovic, G. C., Hui, F. K. C., and Warton, D. I. (2022). Fast model-based ordination with copulas. Methods in Ecology and Evolution, 13, 194–202.


