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


# fit the ZINB_copula model with subset of data to obtain true model and true ordinations.
fit_copula <- function(data, m_species) {
  # Select the first m microbiome data
  y <- data[, 1:m_species]
  
  true_mod <- clvm(y, reff = "fixed", gllvm.fam = "ZINB", seed = 123, lv.n = 0)

 
  
  true_ords <- scale(true_mod$scores)
  
  
  return(list(true_mod = true_mod, true_ords = true_ords))
}

sim_copula <- function(copula_m, rep_k) {
  true_mod <- copula_m$true_mod
  true_ords <- copula_m$true_ords
  
  # Create an empty matrix to store Procrustes errors
  pro.errors <- matrix(NA, rep_k, 6)
  colnames(pro.errors) <- c("C-ZINB", "C-NB", "LVM-ZINB", "LVM-NB", "CLR+PCA", "nMDS")
  
  for (k in 1:rep_k) {
    print(k)
    set.seed(k)
    
    # Simulate the microbiome data using the true copula model
    sim_y <-  simulate.copula(true_mod)
    
   
    # fit the 6 methods with simulated y and calculate the Procrustes distance
    
    # Distance-based methods: nMDS
    
    fit_mds <- try(metaMDS(sim_y, distance = "bray", trace = FALSE))
    mds_ords <- scale(fit_mds$ points) # extract simulated ordinations from nMDS
    p_mds <- try(procrustes(true_ords, mds_ords)$ss) #
    if (class(p_mds)[1] == "try-error") {
      p_mds <- NA
    }
    
    # PCA on CLR-transformed data
    clr_y <- cenLR(sim_y + 1) # CLR transformation of sim_y with added 1
    fit_pca <- prcomp(clr_y$x.clr, scale = TRUE) # PCA on the CLR-transformed data clr_y
    pca_ords <- scale(fit_pca$x[, 1:2]) # extract simulated ordinations from PCA
    
    p_pca <- try(procrustes(true_ords, pca_ords)$ss)
    if (class(p_pca)[1] == "try-error") {
      p_pca <- NA
    }
    
    # Ordination base on Copula with ZINB distribution
    c_ZINB <- try(clvm(sim_y, gllvm.fam = "ZINB", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0))
    c_ZINB_ords <- scale(c_ZINB$scores)
    p_cord1 <- try(procrustes(true_ords, c_ZINB_ords)$ss)
    if (class(p_cord1)[1] == "try-error") {
      p_cord1 <- NA
    }
    
    # Copula with NB distribtuion
    c_NB <- try(clvm(sim_y, gllvm.fam = "negative.binomial", reff = "fixed", sd.errors = FALSE, seed = 123, lv.n = 0))
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
    pro.errors[k, ] <- c(p_cord1, p_cord2, p_lvm1, p_lvm2, p_pca, p_mds)
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



# Do the simulation rep_k times and save the results
# select first m microbes; need to cycle through 50, 100, 200, 400
# 

# a) species/col = 50
copula_50 <- fit_copula(data, m_species = 50)

sim50_copula <- sim_copula(copula_50, rep_k = 50)

write.table(sim50_copula, file = "sim50_copula.txt")

# b) species/col = 100
copula_100 <- fit_copula(data, m_species = 100)

sim100_copula <- sim_copula(copula_100, rep_k = 50)


write.table(sim100_copula, file = "sim100_copula.txt")

# c) species/col = 200
copula_200 <- fit_copula(data, m_species = 200)

sim200_copula <- sim_copula(copula_200, rep_k = 50)

write.table(sim200_copula, file = "sim200_copula.txt")

# d) species/col =400
copula_400 <- fit_copula(data, m_species = 400)
sim400_copula <- sim_copula(copula_400, rep_k = 50)


write.table(sim400_copula, file = "sim400_copula.txt")


style_file("//fileservices.ad.jyu.fi/homes/wentang/Downloads/simulation_copula.R")

