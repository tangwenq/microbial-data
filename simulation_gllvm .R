# packages used:
library(vegan)
library(mvabund)
library(MASS)
library(ecoCopula)
library(mvtnorm)
library(distributions3)
source("clvm.R")
library(ggplot2)
library(robCompositions)
library(gllvm)


# fit the ZINB_gllvm model with subset of data to obtain true model and true ordinations.
fit_gllvm <- function(data, m_species, seed_NO) {
  # Select the first m microbiome data
  y <- data[, 1:m_species]


  true_mod <- gllvm(y, family = "ZINB", sd.errors = FALSE, row.eff = "fixed", seed = seed_NO)

  true_ords <- scale(true_mod$lvs)


  return(list(true_mod = true_mod, true_ords = true_ords))
}

sim_gllvm <- function(gllvm_m, rep_k) {
  true_mod <- gllvm_m$true_mod
  true_ords <- gllvm_m$true_ords

  # Create an empty matrix to store Procrustes errors
  pro.errors <- matrix(NA, rep_k, 7)
  colnames(pro.errors) <- c("C-ZINB", "C-NB", "LVM-ZINB", "LVM-NB", "CLR+PCA", "CLR+nMDS", "nMDS")

  for (k in 1:rep_k) {
    print(k)
    set.seed(k)

    # Simulate the microbiome data using the true gllvm model
    sim_y <- simulate.gllvm(true_mod, conditional = TRUE)


    # fit the 6 methods with simulated y and calculate the Procrustes distance

    # Distance-based methods: nMDS

    fit_mds <- try(metaMDS(sim_y, distance = "bray", trace = FALSE))
    mds_ords <- scale(fit_mds$ points) # extract simluated ordinations from nMDS
    p_mds <- try(procrustes(true_ords, mds_ords)$ss) #
    if (class(p_mds)[1] == "try-error") {
      p_mds <- NA
    }

    # PCA and nMDS on CLR-transformed data
    clr_y <- cenLR(sim_y + 1) # CLR transformation of sim_y with added 1

    fit_clrmds <- try(metaMDS(clr_y$x.clr, distance = "euclidean", autotransform = FALSE, noshare = FALSE, wascores = FALSE, trace = FALSE))
    clrmds_ords <- scale(fit_clrmds$points)
    p_clrmds <- try(procrustes(true_ords, clrmds_ords)$ss)
    if (class(p_clrmds)[1] == "try-error") {
      p_clrmds <- NA
    }

    fit_pca <- prcomp(clr_y$x.clr, scale = TRUE) # PCA on the CLR-transformed data clr_y
    pca_ords <- scale(fit_pca$x[, 1:2]) # extract simluated ordinations from PCA

    p_pca <- try(procrustes(true_ords, pca_ords)$ss)
    if (class(p_pca)[1] == "try-error") {
      p_pca <- NA
    }

    # Ordination base on Copula with ZINB distribution
    c_ZINB <- try(fit_copula(sim_y, gllvm.fam = "ZINB", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0))
    c_ZINB_ords <- scale(c_ZINB$scores)
    p_cord1 <- try(procrustes(true_ords, c_ZINB_ords)$ss)
    if (class(p_cord1)[1] == "try-error") {
      p_cord1 <- NA
    }

    # Copula with NB distribtuion
    c_NB <- try(fit_copula(sim_y, gllvm.fam = "negative.binomial", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0))
    c_NB_ords <- scale(c_NB$scores)
    p_cord2 <- try(procrustes(true_ords, c_NB_ords)$ss)
    if (class(p_cord2)[1] == "try-error") {
      p_cord2 <- NA
    }

    # Ordination base on GLLVM model with ZINB distribution
    lvm_ZINB <- try(gllvm(sim_y, family = "ZINB", sd.errors = FALSE, row.eff = "fixed", num.lv = 2, seed = 123))
    lvm_ZINB_ords <- scale(lvm_ZINB$lvs)
    p_lvm1 <- try(procrustes(true_ords, lvm_ZINB_ords)$ss)
    if (class(p_lvm1)[1] == "try-error") {
      p_lvm1 <- NA
    }




    # GLLVM with NB distribution
    lvm_NB <- try(gllvm(sim_y, family = "negative.binomial", sd.errors = FALSE, row.eff = "fixed", num.lv = 2, seed = 123))
    lvm_NB_ords <- scale(lvm_NB$lvs)
    p_lvm2 <- try(procrustes(true_ords, lvm_NB_ords)$ss)
    if (class(p_lvm2)[1] == "try-error") {
      p_lvm2 <- NA
    }




    # Store the Procrustes errors for k simulation
    pro.errors[k, ] <- c(p_cord1, p_cord2, p_lvm1, p_lvm2, p_pca, p_clrmds, p_mds)
    print(pro.errors[k, ])
  }
  return(pro.errors)
}




# Base the simulations to arctic microbialdata from Nissinen et al.
# Data included in gllvm package

data("microbialdata")
data <- microbialdata$Y
# order columns according to sparseness
data.s <- cbind(1:dim(data)[2], apply(data == 0, 2, sum))
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]
data <- data[, ind]
dim(data)
# amount of zeros in each column
apply(data == 0, 2, sum)



# Do the simulation and save the results
# select first m microbes; need to cycle through 50, 100, 200, 400 with same random seed:123

# a) species/col = 50
gllvm_50 <- fit_gllvm(data, m_species = 50, seed_NO = 123)

sim50_gllvm <- sim_gllvm(gllvm_50, rep_k = 50)

write.table(sim50_gllvm, file = "sim50_gllvm.txt")

# b) species/col = 100
gllvm_100 <- fit_gllvm(data, m_species = 100, seed_NO = 123)

sim100_gllvm <- sim_gllvm(gllvm_100, rep_k = 50)


write.table(sim100_gllvm, file = "sim100_gllvm.txt")

# c) species/col = 200
gllvm_200 <- fit_gllvm(data, m_species = 200, seed_NO = 123)

sim200_gllvm <- sim_gllvm(gllvm_200, rep_k = 50)

write.table(sim200_gllvm, file = "sim200_gllvm.txt")

# d) species/col =400
gllvm_400 <- fit_gllvm(data, m_species = 400, seed_NO = 123)

sim400_gllvm <- sim_gllvm(gllvm_400, rep_k = 50)


write.table(sim400_gllvm, file = "sim400_gllvm.txt")
