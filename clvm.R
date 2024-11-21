#' --- 
#' Function clvm performs model-based ordination along the lines of 
#' Popovic at al. (2022). Fast model‐based ordination with copulas.
#' Marginal model fitted using gllvm (not fast!) to allow ZINB model 
#' and fixed row effects for overdispersed compositional data.
#' 
#' Function simulate.copula simulates responses from copula model.  
#' ---
clvm <- function(data, gllvm.fam = "ZINB", reff = FALSE, lv.n = 0, sd.errors = FALSE) 
  {
  stacked_models <- gllvm(y = data, family = gllvm.fam, num.lv = lv.n, row.eff = reff, sd.errors = FALSE)  
  res_list <- simulate_res_S(obj = stacked_models, n.res = 500) # This should generally high e.g., like 500
  res_list <- plyr::alply(res_list, 3)
  S_list <- lapply(res_list, function(x) cov2cor(cov(x)))  #cov0 assumes 0 mean
  
  P <- dim(S_list[[1]])[1]
  
  A <- factor_opt(nlv = 2, 
                  S.list = S_list, 
                  full = TRUE, 
                  quick = FALSE, 
                  nobs = nrow(data))
  
  Th.out <- A$theta
  Sig.out <- A$sigma
  colnames(Sig.out) <- rownames(Sig.out) <- colnames(Th.out) <- rownames(Th.out) <- seq(1,dim(data)[2])
  
  res_list.mean <- plyr::aaply(plyr::laply(res_list,function(x) x), c(2,3), weighted.mean , weighs = A$weights)
  Scores <- t(as.matrix(A$loadings)) %*% A$theta %*% t(res_list.mean)
  
  make_cordobject <- list(loadings = A$loadings, 
                          scores = t(Scores), 
                          sigma = Sig.out, theta = Th.out,
                          obj = stacked_models)
  class(make_cordobject) <- "coord"
  
  return(make_cordobject)
} 

simulate_res_S <- function(obj, n.res = 500, seed = NULL) {
  
  single_res_fn <- function() {
    res <- residuals(obj)$residuals # Dunn-Smyth residuals
    out <- fix_inf(res)
    return(out)
  }
  
  out <- replicate(n.res, single_res_fn()) 
  
  set.seed(NULL)
  return(out)
}

fix_inf<-function(mat, lim=5){
  mat[mat >  lim] = lim
  mat[mat < (-lim)] = -lim
  mat
}

L.icov.prop = function(S, Theta) {
  # The trace of a product can be rewritten as the sum of entry-wise products of elements
  exp(-1/2 * sum(S * (Theta - diag(dim(Theta)[1]))))
}

factor_opt <- function(nlv, S.list, full = FALSE, quick = FALSE, nobs) {
  P <- dim(S.list[[1]])[1]
  J <- length(S.list)
  eps <- 1e-10
  maxit <- 10
  array.S <- array(unlist(S.list), c(P, P, J))
  
  #initialise weighted covariance and weights
  S.init <- cov2cor(apply(array.S, c(1, 2), mean))
  weights <- rep(1, J)/J
  
  #fit factor analysis
  A <- factanal(NA, nlv, covmat = S.init, nobs = nobs, nstart = 3)
  L <- A$loadings
  Fac <- A$factors
  
  # calculate covariance and precision matrices
  Psi <- diag(diag(S.init - L %*% t(L)))
  PsiInv <- diag(1/diag(Psi))
  Sest <- L %*% t(L) + Psi
  Test <- solve(Sest)
  
  Sigma.gl <- Theta.gl <- list()
  Sigma.gl[[1]] <- Sest
  Theta.gl[[1]] <- Test
  A$sigma <- Sest
  A$theta <- Test
  
  # if not quick, iterate till convergence
  if (!(quick)) {
    count <- 1
    diff <- eps + 1
    while ((diff > eps & count < maxit) & any(!is.na(Theta.gl[[count]]))) {
      
      #recalculate weights and weighted covariance
      weights <- plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
      weights <- weights/sum(weights)
      count <- count + 1
      Sigma.gl[[count]] <- cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))
      
      #factor analysis
      A <- factanal(NA, nlv, covmat = Sigma.gl[[count]], nobs = nobs, nstart = 3)
      L <- A$loadings
      Fac <- A$factors
      
      # calculate covariance and precision matrices and save
      Psi <- diag(diag(S.init - L %*% t(L)))
      PsiInv <- diag(1/diag(Psi))
      Sest <- L %*% t(L) + Psi
      Test <- solve(Sest)
      
      Sigma.gl <- Theta.gl <- list()
      Sigma.gl[[count]] <- Sest
      Theta.gl[[count]] <- Test
      A$sigma <- Sest
      A$theta <- Test
      A$weights <- weights
      
      if (any(!is.na(Theta.gl[[count]]))) {
        diff <- sum(((Theta.gl[[count]] - Theta.gl[[count - 1]])^2)/(P^2))
      } else {
        diff <- sum(((Sigma.gl[[count]] - Sigma.gl[[count - 1]])^2)/(P^2))
      }
    }
    
  }
  
  if (full) {
    return(A)
  } 
  else {
    return(A$theta)
  }
}

simulate.copula <- function(make_cordobject) {
  
  true.mod <- make_cordobject
  true.ords <- true.mod$scores

  # Simulate from NB model; for row.eff = "fixed" add to eta corresponding estimates
  if(true.mod$obj$family == "negative.binomial"){
    sig <- true.mod$sigma[[1]]
    eta <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$beta0)) + true.mod$obj$params$row.params
    phi.inv <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$inv.phi))
    true.load <- as.matrix(true.mod$loadings)
    Psi <- diag(diag(sig - true.load %*% t(true.load)))
    cx.z <- scale(true.ords%*%t(true.load) + rmvnorm(nrow(as.data.frame(true.mod$obj$data)),rep(0,ncol(as.data.frame(true.mod$obj$data))),Psi))
    #cw.y <- qnbinom(pnorm(cx.z),mu=exp(eta),size=phi.inv)
    cw.y <- matrix(NA, nrow  = nrow(eta), ncol = ncol(eta))
    for (i in 1:ncol(eta)) {
      cw.y[,i] <- qnbinom(pnorm(cx.z[,i]), mu = exp(eta[,i]), size = phi.inv[,i])
    }
  }
  
  # Simulate from ZINB model; for row.eff = "fixed" add to eta corresponding estimates
  if(true.mod$obj$family == "ZINB"){
    sig <- true.mod$sigma[[1]]
    eta <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$beta0)) + true.mod$obj$params$row.params
    phi.inv <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$ZINB.inv.phi))
    probs <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$phi))
    true.load <- as.matrix(true.mod$loadings)
    Psi <- diag(diag(sig - true.load %*% t(true.load)))
    cx.z <- scale(true.ords%*%t(true.load) + rmvnorm(nrow(as.data.frame(true.mod$obj$data)),rep(0,ncol(as.data.frame(true.mod$obj$data))),Psi))
    #cw.y <-qzinegbin(pnorm(cx.z), size = phi.inv, munb = exp(eta), pstr0 = probs)
    cw.y <- matrix(NA, nrow  = nrow(eta), ncol = ncol(eta))
    for (i in 1:ncol(eta)) {
      cw.y[,i] <- qzinbinom(pnorm(cx.z[,i]), mu = exp(eta[,i]), size = phi.inv[,i], pi = probs[,i])
    }
  }
  return(cw.y)
}



