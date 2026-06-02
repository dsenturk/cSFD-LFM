###############################################################################
## Description: Functions for generating contrastive multilevel functional data  
##              pairs with model components of Scenario 1 described in  
##              'Contrastive Latent Models for Structured Functional Data'.
###############################################################################

data_gen = function(n_x, n_y, s2_x, s2_y, seed = NULL){
  
  #############################################################################
  ## Description: This function generates the multilevel data sets described in  
  ##              the first simulation scenario of the paper.
  ## Args:        n_x: number of subjects within the first data set (integer).
  ##              n_y: number of subjects within the second data set (integer).
  ##              s2_x: error variance for the first data set, sigma^2_X. (scalar)
  ##              s2_y: error variance for the second data set, sigma^2_Y. (scalar)
  ## Returns:     A list containing generated data pairs (X.df, Y.df) (dataframe: id, visit, data).
  #############################################################################
  
  N = 50 # total number of time points
  tobs = seq(0, 1, length.out = N) # vector of total time points
  J = 2 # number of repetitions
  if(!is.null(seed)) set.seed(seed)
  
  L_1 = 2 # dimension of the level-1 shared space
  K_1 = 1 # dimension of the level-1 unique space of the first group of data X
  R_1 = 1 # dimension of the level-1 unique space of the second group of data Y
  L_2 = 1 # dimension of the level-2 shared space
  K_2 = 2 # dimension of the level-2 unique space of the first group of data X
  R_2 = 1 # dimension of the level-2 unique space of the second group of data Y
  ## Level-1
  Psi_1 = eval.basis(tobs, suppressWarnings(
    create.fourier.basis(nbasis = 5, dropind = c(1,3,5)))) 
  Phi_1 = eval.basis(tobs, suppressWarnings(
    create.fourier.basis(nbasis = 3, dropind = c(1,2)))) 
  Gamma_1 = eval.basis(tobs, suppressWarnings(
    create.fourier.basis(nbasis = 5, dropind = c(1:4)))) 

  ## Level-2
  Psi_2 = matrix(sapply(tobs, FUN = function(x) 1), ncol = 1)
  Phi_2 = cbind(matrix(sapply(tobs, FUN = function(x) (sqrt(3)*(2*x-1))), ncol = 1),
                matrix(sapply(tobs, FUN = function(x) (sqrt(5)*(6*x^2-6*x+1))), ncol = 1)
                )
  Gamma_2 = matrix(sapply(tobs, FUN = function(x) (sqrt(7)*(20*x^3-30*x^2+12*x-1))), ncol = 1)

  ## Level-1
  Omega_x_1 = diag(c(1,.5))
  Omega_y_1 = diag(c(1,.5))
  Theta_1 = matrix(.25)
  Lambda_1 = matrix(.25)

  ## Level-2
  Omega_x_2 = matrix(1)
  Omega_y_2 = matrix(1)
  Theta_2 = diag(c(.5,.25))
  Lambda_2 = matrix(.5)
  

  # Generate X
  mu_x = matrix(data = 0, nrow = N, ncol = 1) 
  X_1 = matrix(data = NA, nrow = (n_x*J), ncol = N) # create an empty matrix to store generated data
  for (i in 1 : n_x) {
    ## level-1 scores
    if(L_1 > 0){
      eta_1 = matrix(data = rnorm(L_1, mean = rep(0, L_1), 
                                sd = sqrt(diag(Omega_x_1))), ncol = 1) 
    }
    if(K_1 > 0){
      xi_1 = matrix(data = rnorm(K_1, mean = rep(0, K_1), 
                               sd = sqrt(diag(Theta_1))), ncol = 1)
    }
    
    for(j in 1 : J){
      ## level-2 scores
      if(L_2 > 0){
        eta_2 = matrix(data = rnorm(L_2, mean = rep(0, L_2), 
                                    sd = sqrt(diag(Omega_x_2))), ncol = 1) 
      }
      if(K_2 > 0){
        xi_2 = matrix(data = rnorm(K_2, mean = rep(0, K_2), 
                                   sd = sqrt(diag(Theta_2))), ncol = 1)
      }
      
      ## generate X_ij
      epsilon = matrix(data = rnorm(N, mean = 0, sd = sqrt(s2_x)), ncol = 1)
      X_ij = mu_x + epsilon
      if(L_1 > 0){X_ij = X_ij + Psi_1 %*% eta_1}
      if(K_1 > 0){X_ij = X_ij + Phi_1 %*% xi_1}
      if(L_2 > 0){X_ij = X_ij + Psi_2 %*% eta_2}
      if(K_2 > 0){X_ij = X_ij + Phi_2 %*% xi_2}
      X_1[((i-1)*J+j),] = t(X_ij)
    }
  }
  X.df = data.frame(id = rep(1:n_x, each = J), visit = rep(1:J, n_x)) 
  X.df$data = X_1
  
  # Generate Y
  mu_y = matrix(data = 0, nrow = N, ncol = 1) 
  Y_1 = matrix(data = NA, nrow = (n_y*J), ncol = N)
  for (i in 1 : n_y) {
    ## level-1 scores
    if(L_1 > 0){
      eta_1 = matrix(data = rnorm(L_1, mean = rep(0, L_1), 
                                  sd = sqrt(diag(Omega_y_1))), ncol = 1) 
    }
    if(R_1 > 0){
      zeta_1 = matrix(data = rnorm(R_1, mean = rep(0, R_1), 
                                 sd = sqrt(diag(Lambda_1))), ncol = 1)
    }
    
    for(j in 1 : J){
      ## level-2 scores
      if(L_2 > 0){
        eta_2 = matrix(data = rnorm(L_2, mean = rep(0, L_2), 
                                    sd = sqrt(diag(Omega_y_2))), ncol = 1) 
      }
      if(R_2 > 0){
        zeta_2 = matrix(data = rnorm(R_2, mean = rep(0, R_2), 
                                   sd = sqrt(diag(Lambda_2))), ncol = 1)
      }
      
      ## generate Y_ij
      epsilon = matrix(data = rnorm(N, mean = 0, sd = sqrt(s2_y)), ncol = 1)
      Y_ij = mu_y + epsilon
      if(L_1 > 0){Y_ij = Y_ij + Psi_1 %*% eta_1}
      if(R_1 > 0){Y_ij = Y_ij + Gamma_1 %*% zeta_1}
      if(L_2 > 0){Y_ij = Y_ij + Psi_2 %*% eta_2}
      if(R_2 > 0){Y_ij = Y_ij + Gamma_2 %*% zeta_2}
      Y_1[((i-1)*J+j),] = t(Y_ij)
    }
  }
  Y.df = data.frame(id = rep(1:n_y, each = J), visit = rep(1:J, n_y)) 
  Y.df$data = Y_1
  
  return(list(X.df, Y.df))
}
