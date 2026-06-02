###############################################################################
## Description: Functions for generating contrastive multivariate functional data  
##              pairs with model components of Scenario 1 described in  
##              'Contrastive Latent Models for Structured Functional Data'.
###############################################################################

data_gen = function(n_x, n_y, s2_x, s2_y, seed = NULL){
  
  #############################################################################
  ## Description: This function generates the multivariate data sets described in  
  ##              the first simulation scenario of the paper.
  ## Args:        n_x: number of subjects within the first data set (integer).
  ##              n_y: number of subjects within the second data set (integer).
  ##              s2_x: common error variance for the first data set, sigma^2_{X,j}. (scalar)
  ##              s2_y: common error variance for the second data set, sigma^2_{Y,j}. (scalar)
  ## Returns:     A list containing generated data pairs (X.ls, Y.ls) (list).
  #############################################################################
  
  N = 50 # total number of time points
  tobs = seq(0, 1, length.out = N) # vector of total time points
  J = 12 # number of outcomes
  if(!is.null(seed)) set.seed(seed)
  s2_x = rep(s2_x, J)
  s2_y = rep(s2_y, J)
  beta.X <- c(1, 0.88, 0.84, 0.75, 0.71, 0.93, 0.81, 1.11, 0.80, 0.89, 0.90, 0.73)
  beta.Y <- c(1, 0.59, 0.63, 0.54, 0.35, 0.83, 0.73, 0.98, 0.78, 0.92, 0.69, 0.71)

  # Dimension of latent component spaces
  L = 1 # dimension of the subject-level shared space
  K = 1 # dimension of the subject-level unique space of the first group of data X
  R = 1 # dimension of the subject-level unique space of the second group of data Y
  V_x = 2 # dimension of the outcome-level latent space of X
  V_y = 3 # dimension of the outcome-level latent space of Y
  
  ## Subject-level
  Psi = eval.basis(tobs, suppressWarnings(create.fourier.basis(
    nbasis = 3, dropind = c(1,3))))
  Phi = eval.basis(tobs, suppressWarnings(create.fourier.basis(
    nbasis = 3, dropind = c(1,2))))
  Gamma = eval.basis(tobs, suppressWarnings(create.fourier.basis(
    nbasis = 5, dropind = c(1:3,5))))
  
  ## Outcome-level
  tau_x = cbind(
    matrix(sapply(tobs, FUN = function(x) 1), ncol = 1),
    matrix(sapply(tobs, FUN = function(x) (sqrt(3)*(2*x-1))), ncol = 1)
  ) 
  tau_y = cbind(
    matrix(sapply(tobs, FUN = function(x) 1), ncol = 1),
    matrix(sapply(tobs, FUN = function(x) (sqrt(5)*(6*x^2-6*x+1))), ncol = 1),
    matrix(sapply(tobs, FUN = function(x) (sqrt(7)*(20*x^3-30*x^2+12*x-1))), ncol = 1)
  ) 
  
  # Score variances
  ## Subject-level
  Omega_x = matrix(1)
  Omega_y = matrix(1)
  Theta = matrix(.5)
  Lambda = matrix(.5)
  
  ## Outcome-level
  varsigma_x = diag(c(1, .5))
  varsigma_y = diag(c(1, .5, .25))
  
  # Mean functions
  mu_x <- matrix(0, nrow = N, ncol = J)
  mu_y <- matrix(0, nrow = N, ncol = J)
  
  # Generate Dataset X
  eta_x  <- matrix(mvrnorm(n_x, mu = rep(0, L), Sigma = Omega_x),
                   nrow = n_x, ncol = L)
  if(K > 0){
    xi_x   <- matrix(mvrnorm(n_x, mu = rep(0, K), Sigma = Theta),
                     nrow = n_x, ncol = K)
  }else{
    xi_x <- matrix(0, nrow = n_x, ncol = K)
  }
  ## Between-outcome latent process
  if(K > 0){
    U_x <- eta_x %*% t(Psi) + xi_x %*% t(Phi)
  }else{
    U_x <- eta_x %*% t(Psi)
  }
  ## Generate each outcome j and store in list
  X.ls <- vector("list", J)
  for (j in 1:J) {
    nu_x <- matrix(mvrnorm(n_x, mu = rep(0, V_x), Sigma = varsigma_x),
                   nrow = n_x, ncol = V_x) # within-outcome scores
    W_xj <- nu_x %*% t(tau_x) # within-outcome process
    eps_xj <- matrix(rnorm(n_x * N, mean = 0, sd = sqrt(s2_x[j])),
                     nrow = n_x, ncol = N) # noise
    Xj <- sweep(beta.X[j] * (U_x + W_xj) + eps_xj, 2, mu_x[, j], "+") # full observed data
    
    # Convert to long-format data frame: (subj, argvals, y)
    X.ls[[j]] <- data.frame(
      subj    = rep(1:n_x, each = N),
      argvals = rep(tobs, times = n_x),
      y       = c(t(Xj))   # row-by-row vectorization
    )
  }
  
  # Generate Dataset Y
  eta_y  <- matrix(mvrnorm(n_y, mu = rep(0, L), Sigma = Omega_y),
                   nrow = n_y, ncol = L)
  if(R > 0){
    zeta_y <- matrix(mvrnorm(n_y, mu = rep(0, R), Sigma = Lambda),  
                     nrow = n_y, ncol = R)
  }else{
    zeta_y <- matrix(0, nrow = n_y, ncol = R)
  }
  ## Between-outcome latent process
  if(R > 0){
    U_y <- eta_y %*% t(Psi) + zeta_y %*% t(Gamma)
  }else{
    U_y <- eta_y %*% t(Psi)
  }

  ## Generate each outcome j and store in list
  Y.ls <- vector("list", J)
  for (j in 1:J) {
    nu_y <- matrix(mvrnorm(n_y, mu = rep(0, V_y), Sigma = varsigma_y), 
                   nrow = n_y, ncol = V_y) # within-outcome scores
    W_yj <- nu_y %*% t(tau_y) # within-outcome process
    eps_yj <- matrix(rnorm(n_y * N, mean = 0, sd = sqrt(s2_y[j])),
                     nrow = n_y, ncol = N) # noise
    Yj <- sweep(beta.Y[j] * (U_y + W_yj) + eps_yj, 2, mu_y[, j], "+") # full observed data
    
    # Convert to long-format data frame
    Y.ls[[j]] <- data.frame(
      subj    = rep(1:n_y, each = N),
      argvals = rep(tobs, times = n_y),
      y       = c(t(Yj))
    )
  }
  return(list(X.ls = X.ls, Y.ls = Y.ls))
}
