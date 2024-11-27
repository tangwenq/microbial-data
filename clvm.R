# Define main function for model-based ordination using copulas with Zero-Inflated Negative Binomial(ZINB) models via gllvm

fit_ZINBcopula <- function(data, gllvm.fam = "ZINB", reff = FALSE, lv.n = 0, sd.errors = FALSE, seed = NULL)
# Set Random Seed to ensures reproducibility.
{
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  # Use ZINB distribution in "gllvm" to fit the marginal model.
  stacked_models <- gllvm(y = data, family = gllvm.fam, num.lv = lv.n, row.eff = reff, sd.errors = FALSE)

  # Simulate Dunn-Smyth residuals
  res_list <- simulate_res_S(obj = stacked_models, n.res = 500)
  res_list <- plyr::alply(res_list, 3)
  S_list <- lapply(res_list, function(x) cov2cor(cov(x))) # compute a list of correlation matrices (S_list), each with dimension P (number of variables)

  P <- dim(S_list[[1]])[1]
  # Factor Analysis
  A <- factor_opt(
    nlv = 2,
    S.list = S_list,
    full = TRUE,
    quick = FALSE,
    nobs = nrow(data)
  )

  Th.out <- A$theta
  Sig.out <- A$sigma
  colnames(Sig.out) <- rownames(Sig.out) <- colnames(Th.out) <- rownames(Th.out) <- seq(1, dim(data)[2])

  res_list.mean <- plyr::aaply(plyr::laply(res_list, function(x) x), c(2, 3), weighted.mean, weighs = A$weights)
  Scores <- t(as.matrix(A$loadings)) %*% A$theta %*% t(res_list.mean)

  make_cordobject <- list(
    loadings = A$loadings,
    scores = t(Scores),
    sigma = Sig.out, theta = Th.out,
    obj = stacked_models
  ) # calculate scores and construct a 'coord' object for model-based ordination
  class(make_cordobject) <- "coord"

  return(make_cordobject)
}


# Simulates residuals for a fitted gllvm model.
simulate_res_S <- function(obj, n.res = 500, seed = NULL) {
  single_res_fn <- function() {
    res <- residuals(obj)$residuals # Extract Dunn-Smyth residuals
    out <- fix_inf(res) # Cap Extreme residual.
    return(out)
  }
  # Simulate residuals n.res times
  out <- replicate(n.res, single_res_fn())

  set.seed(NULL)
  return(out)
}

# Caps extreme values in a matrix to ensure numerical stability.
fix_inf <- function(mat, lim = 5) {
  mat[mat > lim] <- lim
  mat[mat < (-lim)] <- -lim
  mat
}

# Calculates a proportionality measure based on the trace of the product of two matrices.
L.icov.prop <- function(S, Theta) {
  # The trace of a product can be rewritten as the sum of entry-wise products of elements
  exp(-1 / 2 * sum(S * (Theta - diag(dim(Theta)[1]))))
}


factor_opt <- function(nlv, S.list, full = FALSE, quick = FALSE, nobs) {
  P <- dim(S.list[[1]])[1]
  J <- length(S.list)
  eps <- 1e-10
  maxit <- 10
  array.S <- array(unlist(S.list), c(P, P, J))

  # initialise weighted covariance and weights
  S.init <- cov2cor(apply(array.S, c(1, 2), mean))
  weights <- rep(1, J) / J

  # fit factor analysis
  A <- factanal(NA, nlv, covmat = S.init, nobs = nobs, nstart = 3)
  L <- A$loadings
  Fac <- A$factors

  # calculate covariance and precision matrices
  Psi <- diag(diag(S.init - L %*% t(L)))
  PsiInv <- diag(1 / diag(Psi))
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
      # recalculate weights and weighted covariance
      weights <- plyr::laply(S.list, L.icov.prop, Theta = Theta.gl[[count]])
      weights <- weights / sum(weights)
      count <- count + 1
      Sigma.gl[[count]] <- cov2cor(apply(array.S, c(1, 2), weighted.mean, w = weights))

      # factor analysis
      A <- factanal(NA, nlv, covmat = Sigma.gl[[count]], nobs = nobs, nstart = 3)
      L <- A$loadings
      Fac <- A$factors

      # calculate covariance and precision matrices and save
      Psi <- diag(diag(S.init - L %*% t(L)))
      PsiInv <- diag(1 / diag(Psi))
      Sest <- L %*% t(L) + Psi
      Test <- solve(Sest)

      Sigma.gl <- Theta.gl <- list()
      Sigma.gl[[count]] <- Sest
      Theta.gl[[count]] <- Test
      A$sigma <- Sest
      A$theta <- Test
      A$weights <- weights

      if (any(!is.na(Theta.gl[[count]]))) {
        diff <- sum(((Theta.gl[[count]] - Theta.gl[[count - 1]])^2) / (P^2))
      } else {
        diff <- sum(((Sigma.gl[[count]] - Sigma.gl[[count - 1]])^2) / (P^2))
      }
    }
  }

  if (full) {
    return(A)
  } else {
    return(A$theta)
  }
}




# Simulates data from a copula model based on the fitted model object.

simulate_copula <- function(make_cordobject) {
  true.mod <- make_cordobject # Extract the model object
  true.ords <- true.mod$scores # Latent scores from the model
  n <- nrow(true.mod$obj$y) # Number of observations
  m <- ncol(true.mod$obj$y) # Number of variables
  sim_y <- matrix(0, n, m) # Initialize matrix for simulated responses

  # Simulate data from a Negative Binomial (NB) model
  if (true.mod$obj$family == "negative.binomial") {
    sig <- true.mod$sigma[[1]] # Covariance matrix
    eta <- t(replicate(n, true.mod$obj$params$beta0)) + c(0, true.mod$obj$params$row.params) # Linear predictor
    phi.inv <- t(replicate(n, true.mod$obj$params$inv.phi)) # Inverse dispersion parameter
    true.load <- as.matrix(true.mod$loadings) # Factor loadings matrix
    Psi <- diag(diag(sig - true.load %*% t(true.load))) # Covariance of residuals
    cx.z <- scale(true.ords %*% t(true.load) + rmvnorm(n, rep(0, m), Psi)) # Simulate latent variables

    # Generate responses using the NB distribution
    for (i in 1:m) {
      sim_y[, i] <- qnbinom(pnorm(cx.z[, i]), mu = exp(eta[, i]), size = phi.inv[, i]) # Simulate response
    }
  }

  # Simulate data from a Zero-Inflated Negative Binomial (ZINB) model
  if (true.mod$obj$family == "ZINB") {
    sig <- true.mod$sigma[[1]] # Covariance matrix
    eta <- t(replicate(n, true.mod$obj$params$beta0)) + c(0, true.mod$obj$params$row.params) # Linear predictor
    phi.inv <- t(replicate(nrow(as.data.frame(true.mod$obj$data)), true.mod$obj$params$ZINB.inv.phi)) # Inverse dispersion
    probs <- t(replicate(n, true.mod$obj$params$phi)) # Zero-inflation probabilities
    true.load <- as.matrix(true.mod$loadings) # Factor loadings matrix
    Psi <- diag(diag(sig - true.load %*% t(true.load))) # Covariance of residuals
    cx.z <- scale(true.ords %*% t(true.load) + rmvnorm(n, rep(0, m), Psi)) # Simulate latent variables

    # Generate responses using the Zero-Inflated Negative Binomial distribution
    for (i in 1:m) {
      sim_y[, i] <- qzinbinom(pnorm(cx.z[, i]), mu = exp(eta[, i]), size = phi.inv[, i], pi = probs[, i]) # Simulate response
    }
  }

  return(sim_y) # Return the simulated responses
}
