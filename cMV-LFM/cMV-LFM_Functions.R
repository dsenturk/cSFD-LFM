###############################################################################
## Description: Functions for conducting the proposed cMV-LFM estimation algorithm   
##              described in Algorithm 2 of 'Contrastive Latent Models for Structured Functional Data'.
###############################################################################
## Functions included:
## Main function:
##    cMV_LFM: Function for fitting a cMV-LFM model and executing the proposed estimation
##             algorithm described in Algorithm 2.
## Supporting functions used by the main function:
##    1. weight_fun: Supporting function that returns the trapezoidal weights based on the given time grid.
##    2. extract_MFMM: Supporting function that estimates scaling parameters beta_j and subject- and 
##                     outcome-level covariances.
###############################################################################


cMV_LFM = function(X.ls, Y.ls, tobs = NULL, 
                   pve_X = c(.90, .90), pve_Y = c(.90, .90)){
  
  #############################################################################
  ## Description: (main function) This function executes the proposed cMV-LFM estimation algorithm described in Algorithm 2.
  ##              The two input lists (X.ls, Y.ls) are constructed with J dataframes (outcomes) each containing three variables: 
  ##              subj (id), argvals (time), y (data), representing `argvals` measurements of `subj` subjects for the `y`th outcome.
  ## Args:        X.ls: list of the first group of data (J dataframes with three arguments: (1) subj (id): subject indices; 
  ##                    (2) argvals (time): observation times; (3) y (data): values of observations for each outcome).
  ##              Y.ls: list of the second group of data (J dataframes: subj (id), argvals (time), y (data)).
  ##              tobs: total time grid T (vector, N*1).
  ##              pve_X: proportion of variance explained for retaining leading eigenfunctions of X in Stage I
  ##              pve_Y: proportion of variance explained for retaining leading eigenfunctions of Y in Stage I
  ## Returns:     A list containing:
  ##                mu_x: estimate of the outcome-specific mean function mu^x_j(t) (matrix, N*J)
  ##                mu_y: estimate of the outcome-specific mean function mu^y_j(t) (matrix, N*J)
  ##                s2_x: estimate of the outcome-specific error variance sigma^2_{X,j} (vector, J*1)
  ##                s2_y: estimate of the outcome-specific error variance sigma^2_{Y,j} (vector, J*1)
  ##                beta.X: estimate of the outcome-specific scaling parameters beta^x_j (vector, J*1)
  ##                beta.Y: estimate of the outcome-specific scaling parameters beta^y_j (vector, J*1)
  ##                Psi: estimate of the subject-level shared latent components {psi_l(t), l = 1,...,L} (matrix, N*L).
  ##                Phi: estimate of the subject-level unique latent components for group X {phi_k(t), k = 1,...,K} (matrix, N*K).
  ##                Gamma: estimate of the subject-level unique latent components for group Y {gamma_r(t), r = 1,...,R} (matrix, N*R).
  ##                tau_x: estimate of the outcome-level unique latent components for group X {tau^x_v(t), v = 1,...,V^x} (matrix, N*V^x).
  ##                tau_y: estimate of the outcome-level unique latent components for group Y {tau^y_v(t), v = 1,...,V^y} (matrix, N*V^y).
  ##                Omega_x: diagonal matrix of estimates of the score variances {omega^{x}_l, l = 1,...,L} (matrix, L*L).
  ##                Omega_y: diagonal matrix of estimates of the score variances {omega^{y}_l, l = 1,...,L} (matrix, L*L).
  ##                Theta: diagonal matrix of estimates of the score variances {theta_k, k = 1,...,K} (matrix, K*K).
  ##                Lambda: diagonal matrix of estimates of the score variance {lambda_r, r = 1,...,R} (matrix, R*R).
  ##                varsigma_x: diagonal matrix of estimates of the score variances {varsigma^x_v, v = 1,...,V^x} (matrix, V^x*V^x).
  ##                varsigma_y: diagonal matrix of estimates of the score variance {varsigma^y_v, v = 1,...,V^y} (matrix, V^y*V^y).
  ##                X_score: BLUP estimate of the subject-specific latent scores for group X (dataframe: id, outcome, subject-level score, outcome-level score).
  ##                Y_score: BLUP estimate of the subject-specific latent scores for group Y (dataframe: id, outcome, subject-level score, outcome-level score).
  #############################################################################
  
  stopifnot((!is.null(X.ls) && !is.null(Y.ls)))
  
  # Convert the input data lists into dataframes for contrastive analysis
  X.df <- do.call(rbind, lapply(seq_along(X.ls), function(j) {
    nm <- if (!is.null(names(X.ls))) names(X.ls)[j] else j
    cbind(X.ls[[j]], outcome = nm)
  }))
  X.df$outcome <- factor(X.df$outcome)
  colnames(X.df) <- c("id", "tobs", "data", "outcome")
  
  Y.df <- do.call(rbind, lapply(seq_along(Y.ls), function(j) {
    nm <- if (!is.null(names(Y.ls))) names(Y.ls)[j] else j
    cbind(Y.ls[[j]], outcome = nm)
  }))
  Y.df$outcome <- factor(Y.df$outcome)
  colnames(Y.df) <- c("id", "tobs", "data", "outcome")
  
  if(setequal(X.df$tobs, Y.df$tobs)){ # check the consistency of time grids from two groups
    if(is.null(tobs)){tobs = sort(unique(X.df$tobs))}  # total time grid T (vector, N*1).
    N = length(tobs) # total number of time points
  }else{
    stop(paste("The total distinct time points are inconsistent between two",
               "groups of data."))
  }
  if(setequal(X.df$outcome, Y.df$outcome)){ # check the consistency of outcomes from two groups
    J = length(unique(X.df$outcome))  # total number of outcomes
  }else{
    stop(paste("The number of outcomes are inconsistent between two",
               "groups of data."))
  }
  n_x = length(unique(X.df$id)) ## number of subjects in X
  n_y = length(unique(Y.df$id)) ## number of subjects in Y
  
  X.wide <- X.df %>%
    pivot_wider(id_cols = c(id, outcome), names_from  = tobs, values_from = data)
  X.df <- data.frame(id = X.wide$id, outcome = X.wide$outcome)
  X.df$data <- as.matrix(X.wide[, -(1:2)]) # attach matrix as a variable
  nOutcomes_x <- data.frame(table(X.df$id)) # number of outcomes for each subject X_i
  colnames(nOutcomes_x) <- c("id", "numOutcomes")
  rm(X.wide)
  Y.wide <- Y.df %>%
    pivot_wider(id_cols = c(id, outcome), names_from = tobs, values_from = data)
  Y.df <- data.frame(id = Y.wide$id, outcome = Y.wide$outcome)
  Y.df$data <- as.matrix(Y.wide[, -(1:2)])
  nOutcomes_y <- data.frame(table(Y.df$id)) ## number of outcomes for each subject Y_i
  colnames(nOutcomes_y) <- c("id", "numOutcomes")
  rm(Y.wide)
  
  
  #################################
  # Step 1. Estimate the mean functions mu^x_j(t) and mu^y_j(t) (by penalized smoothing spline), 
  #         error variances sigma^2_{X,j} and sigma^2_{Y,}, and auto- and cross-covariances 
  #         K_jj(s,t), K_jj'(s,t) using fast covariance method (Li et al., 2020).
  #################################
  
  # X
  ## Estimate mu^x(t) and sigma^2_{X,j}
  mu_x <- matrix(nrow = N, ncol = J)
  for (j in 1:J) {
    fit_df <- data.frame(
      y       = as.vector(X.ls[[j]]$y),
      argvals = as.vector(X.ls[[j]]$argvals)
    )
    gam0 <- gam(y ~ s(argvals, k = 10), data = fit_df)
    mu_x[, j] <- predict(gam0, newdata = data.frame(argvals = tobs))
    tobs_idx <- match(X.ls[[j]]$argvals, tobs)
    X.ls[[j]]$y <- X.ls[[j]]$y - mu_x[tobs_idx, j] # centralize the data
  }
  for(ind in 1:NROW(X.df)) {
    X.df$data[ind,] = X.df$data[ind,] - mu_x[,X.df$outcome[ind]]
  }
  
  ## Fast cov estimation
  fit_1 <- mface.sparse(X.ls, center = FALSE, argvals.new = tobs, knots = 7) 
  s2_x = fit_1$var.error.new[1,]
  if(!all(s2_x > 0)){stop("Estimate of error variance for X is non-positive.")}
  
  
  # Y
  ## mu^y(t) and sigma^2_Y
  mu_y <- matrix(nrow = N, ncol = J)
  for (j in 1:J) {
    fit_df <- data.frame(
      y       = as.vector(Y.ls[[j]]$y),
      argvals = as.vector(Y.ls[[j]]$argvals)
    )
    gam0 <- gam(y ~ s(argvals, k = 10), data = fit_df)
    mu_y[, j] <- predict(gam0, newdata = data.frame(argvals = tobs))
    tobs_idx <- match(Y.ls[[j]]$argvals, tobs)
    Y.ls[[j]]$y <- Y.ls[[j]]$y - mu_y[tobs_idx, j]
  }
  for(ind in 1:NROW(Y.df)) {
    Y.df$data[ind,] = Y.df$data[ind,] - mu_y[,Y.df$outcome[ind]]
  }
  
  ## Fast cov estimation
  fit_2 <- mface.sparse(Y.ls, center = FALSE, argvals.new = tobs, knots = 7) 
  s2_y = fit_2$var.error.new[1,]
  if(!all(s2_y > 0)){stop("Estimate of error variance for Y is non-positive.")}

  
  #################################
  # Step 2. Estimate outcome-specific scaling parameters \beta_j and subject- and outcome-level 
  #         covariance functions G^x_1, G^y_1, G^x_\tau, G^y_\tau
  #################################
  
  ## X
  res_1 <- extract_MFMM(
    Chat_new = as.matrix(fit_1$Chat.new),
    N    = N,
    J    = J,
    tobs = tobs
  )
  
  beta.X = res_1$beta
  G_x_0 = res_1$C0
  G_x_1 = res_1$C1
  
  ## Y
  res_2 <- extract_MFMM(
    Chat_new = as.matrix(fit_2$Chat.new),
    N    = N,
    J    = J,
    tobs = tobs
  )
  
  beta.Y = res_2$beta
  G_y_0 = res_2$C0
  G_y_1 = res_2$C1
  
  
  #################################
  # Step 3. Perform FPCA on subject- and outcome-level covariance functions and retain
  #         leading eigenfunctions based on pre-specified pve thresholds.
  #################################
  
  # X
  W_diag = weight_fun(tobs) # construct the trapezoidal weights
  W = diag(x = W_diag, nrow = N)
  Wsqrt_diag = sqrt(W_diag)
  Wsqrt = diag(x = Wsqrt_diag, nrow = N)
  Winvsqrt_diag = 1 / Wsqrt_diag
  Winvsqrt <- diag(x = Winvsqrt_diag, nrow = N)
  
  S_tilde <- sweep(sweep(G_x_0, 1, Wsqrt_diag, "*"), 2, Wsqrt_diag, "*") # S_tilde[i,j] = Wsqrt[i] * G_x_0[i,j] * Wsqrt[j]
  eig <- eigen(S_tilde, symmetric = TRUE) # eigen-decomposition of the symmetric transformed matrix
  pos_idx <- eig$values > 0 # Keep only positive eigenvalues
  if (!any(pos_idx)) {stop("No positive eigenvalues;")}
  vals <- eig$values[pos_idx]
  psi  <- eig$vectors[, pos_idx, drop = FALSE]
  
  # Select number of components by pre-specified PVE
  cum_pve <- cumsum(vals) / sum(vals)
  npc <- which.max(cum_pve >= pve_X[1])
  vals <- vals[1:npc]
  psi  <- psi[, 1:npc, drop = FALSE]
  
  Varphi_x <- sweep(psi, 1, Winvsqrt_diag, "*") # recover subject-level eigenfunctions
  M_x = NCOL(Varphi_x)
  
  
  # Y
  S_tilde <- sweep(sweep(G_y_0, 1, Wsqrt_diag, "*"), 2, Wsqrt_diag, "*")
  eig <- eigen(S_tilde, symmetric = TRUE)
  pos_idx <- eig$values > 0
  if (!any(pos_idx)) {stop("No positive eigenvalues;")}
  vals <- eig$values[pos_idx]
  psi  <- eig$vectors[, pos_idx, drop = FALSE]
  
  # Select number of components by pre-specified PVE
  cum_pve <- cumsum(vals) / sum(vals)
  npc <- which.max(cum_pve >= pve_Y[1])
  vals <- vals[1:npc]
  psi  <- psi[, 1:npc, drop = FALSE]
  
  Varphi_y <- sweep(psi, 1, Winvsqrt_diag, "*") 
  M_y = NCOL(Varphi_y)
  
  
  #################################
  # Step 4. Construct the group-specific projection operators and their weighted average projection operators
  #################################
  
  # X
  P_x = Wsqrt %*% tcrossprod(Varphi_x) %*% Wsqrt 
  # Y
  P_y = Wsqrt %*% tcrossprod(Varphi_y) %*% Wsqrt 
  # Weighted average 
  P_w = (P_x * n_x + P_y * n_y) / (n_x + n_y)
  
  #################################
  # Step 5. Identify the rank of shared space and construct the shared projection operator P_s(s,t) 
  #################################
  
  eigen_res = eigen(P_w, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  L = sum(eigen_res$values > .9) # identify the rank of the subject-level shared space
  P_s = tcrossprod(as.matrix(eigen_res$vectors[, 1:L])) # projector of shared subspace
  
  
  #################################
  # Step 6. Target the subject-level shared latent components {\psi_l(t), l = 1, ..., L}
  #################################
  
  G_x_s = Winvsqrt %*% P_s %*% Wsqrt %*% G_x_0 %*% Wsqrt %*% P_s %*% Winvsqrt
  G_y_s = Winvsqrt %*% P_s %*% Wsqrt %*% G_y_0 %*% Wsqrt %*% P_s %*% Winvsqrt
  eigen_res = eigen(Wsqrt %*% (G_x_s + G_y_s) %*% Wsqrt, symmetric = TRUE)
  Psi = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:L])
 
  
  #################################
  # Step 7. Construct group-specific complementary projection operators and identify the ranks of unique spaces
  #################################
  
  ## X
  P_x.c = P_x %*% (diag(N) - P_s) %*% P_x
  eigen_res = eigen(P_x.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  K = sum(eigen_res$values > .9) # identify the rank of the subject-level unique subspace for X
  if(K > 0){P_x_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:K]))}
  
  ## Y
  P_y.c = P_y %*% (diag(N) - P_s) %*% P_y
  eigen_res = eigen(P_y.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  R = sum(eigen_res$values > .9) # identify the rank of the subject-level unique subspace for Y
  if(R > 0){P_y_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:R]))}
  
  
  #################################
  # Step 8. Target group-specific unique latent components 
  #################################
  
  if(K > 0){
    G_x_u = Winvsqrt %*% P_x_u %*% Wsqrt %*% G_x_0 %*% Wsqrt %*% P_x_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_x_u %*% Wsqrt, symmetric = TRUE)
    Phi = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:K])
  }else{Phi = NULL}
  
  if(R > 0){
    G_y_u = Winvsqrt %*% P_y_u %*% Wsqrt %*% G_y_0 %*% Wsqrt %*% P_y_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_y_u %*% Wsqrt, symmetric = TRUE)
    Gamma = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:R])
  }else{Gamma = NULL}
  
  #################################
  # Step 9. Estimate outcome-level latent components
  #################################
  
  # X
  S_tilde <- sweep(sweep(G_x_1, 1, Wsqrt_diag, "*"), 2, Wsqrt_diag, "*")
  eig <- eigen(S_tilde, symmetric = TRUE) # eigen-decomposition of the symmetric transformed matrix
  pos_idx <- eig$values > 0 # Keep only positive eigenvalues
  if (!any(pos_idx)) {stop("No positive eigenvalues;")}
  vals <- eig$values[pos_idx]
  psi  <- eig$vectors[, pos_idx, drop = FALSE]
  
  ## Select number of components by pre-specified PVE
  cum_pve <- cumsum(vals) / sum(vals)
  npc <- which.max(cum_pve >= pve_X[2])
  vals <- vals[1:npc]
  psi  <- psi[, 1:npc, drop = FALSE]
  
  tau_x <- sweep(psi, 1, Winvsqrt_diag, "*") # recover outcome-level eigenfunctions
  varsigma_x = vals
  V_x = NCOL(tau_x)
  
  # Y
  S_tilde <- sweep(sweep(G_y_1, 1, Wsqrt_diag, "*"), 2, Wsqrt_diag, "*")
  eig <- eigen(S_tilde, symmetric = TRUE)
  pos_idx <- eig$values > 0
  if (!any(pos_idx)) {stop("No positive eigenvalues;")}
  vals <- eig$values[pos_idx]
  psi  <- eig$vectors[, pos_idx, drop = FALSE]
  
  ## Select number of components by pre-specified PVE
  cum_pve <- cumsum(vals) / sum(vals)
  npc <- which.max(cum_pve >= pve_Y[2])
  vals <- vals[1:npc]
  psi  <- psi[, 1:npc, drop = FALSE]
  
  tau_y <- sweep(psi, 1, Winvsqrt_diag, "*")
  varsigma_y = vals
  V_y = NCOL(tau_y)
  
  
  #################################
  # Step 10. Target variance components
  #################################

  # For X
  if(L > 0){
    Omega_x = diag(x = diag(as.matrix(t(Psi) %*% W %*% G_x_0 %*% W %*% Psi)), nrow = L)
  }else{Omega_x = NULL}
  if(K > 0){
    Theta = diag(x = diag(as.matrix(t(Phi) %*% W %*% G_x_0 %*% W %*% Phi)), nrow = K)
  }else{
    Theta = NULL
  }
  
  # Inversion of data covariance matrix Sigma
  F_1 = cbind(Psi, Phi)
  F_2 = tau_x
  D_1 = diag(x = c(diag(Omega_x), diag(Theta)), nrow = (L + K))
  D_2 = diag(x = varsigma_x, nrow = length(varsigma_x))

  ## Calculate inverse of covariance matrix Sigma for each subject
  FDF_1 = F_1 %*% D_1 %*% t(F_1)
  FDF_2 = F_2 %*% D_2 %*% t(F_2)
  Sigma_inv = solve(Matrix::kronecker(matrix(1, J, J), FDF_1) +
                      Matrix::kronecker(diag(J), FDF_2) + 
                      Matrix::kronecker(diag((s2_x/(beta.X^2)), nrow = J), diag(N)))

  ## Calculate score estimates for each subject
  X_score = tibble(X.df[,1:2])
  X_score$subj_score = vector("list", NROW(X_score))
  X_score$oc_score = vector("list", NROW(X_score))
  for (ind in unique(X.df$id)) {
    X_i = matrix(t(X.df$data[X.df$id == ind,]), ncol = 1) # flatten to a 1-d vector
    J_i = nOutcomes_x$numOutcomes[nOutcomes_x$id == ind]
    outcome_i = X.df$outcome[X.df$id == ind]
    if(is.null(Sigma_inv)){
      Sigma_inv = solve(Matrix::kronecker(matrix(1, J_i, J_i), FDF_1) + 
                          Matrix::kronecker(diag(J_i), FDF_2) + 
                          Matrix::kronecker(diag((s2_x/(beta.X^2)), nrow = J)[outcome_i, outcome_i], diag(N)))
    }
    beta_i = rep(beta.X[outcome_i], each = N)
    score_mat = rbind(Matrix::kronecker(matrix(1, 1, J_i), (D_1 %*% t(F_1))),
                      Matrix::kronecker(diag(J_i), (D_2 %*% t(F_2)))) %*%
      Sigma_inv %*% (X_i / beta_i)
    for(j in 1:J_i){
      X_score$subj_score[X_score$id == ind & X_score$outcome == outcome_i[j]][[1]] = score_mat[1:(L+K)]
      X_score$oc_score[X_score$id == ind & X_score$outcome == outcome_i[j]][[1]] = score_mat[(L+K+(j-1)*(V_x)+1):(L+K+j*(V_x))]
    }
  }


  # For Y
  if(L > 0){Omega_y = diag(x = diag(as.matrix(t(Psi) %*% W %*% G_y_0 %*% W %*% Psi)), nrow = L)}else{Omega_y = NULL}
  if(R > 0){Lambda = diag(x = diag(as.matrix(t(Gamma) %*% W %*% G_y_0 %*% W %*% Gamma)), nrow = R)}else{Lambda = NULL}
  

  # Inversion of data covariance matrix Sigma
  F_1 = cbind(Psi, Gamma)
  F_2 = cbind(tau_y)
  D_1 = diag(x = c(diag(Omega_y), diag(Lambda)), nrow = (L + R))
  D_2 = diag(x = varsigma_y, nrow = length(varsigma_y))
  
  ## Calculate inverse of covariance matrix Sigma for each subject
  FDF_1 = F_1 %*% D_1 %*% t(F_1)
  FDF_2 = F_2 %*% D_2 %*% t(F_2)
  Sigma_inv = solve(Matrix::kronecker(matrix(1, J, J), FDF_1) +
                      Matrix::kronecker(diag(J), FDF_2) + 
                      Matrix::kronecker(diag((s2_y/(beta.Y^2)), nrow = J), diag(N)))

  
  ## Calculate score estimates for each subject
  Y_score = tibble(Y.df[,1:2])
  Y_score$subj_score = vector("list", NROW(Y_score))
  Y_score$oc_score = vector("list", NROW(Y_score))
  for (ind in unique(Y.df$id)) {
    Y_i = matrix(t(Y.df$data[Y.df$id == ind,]), ncol = 1)
    J_i = nOutcomes_y$numOutcomes[nOutcomes_y$id == ind]
    outcome_i = Y.df$outcome[Y.df$id == ind]
    if(is.null(Sigma_inv)){
      Sigma_inv = solve(Matrix::kronecker(matrix(1, J_i, J_i), FDF_1) + 
                          Matrix::kronecker(diag(J_i), FDF_2) + 
                          Matrix::kronecker(diag((s2_y/(beta.Y^2)), nrow = J)[outcome_i, outcome_i], diag(N)))
    }
    beta_i = rep(beta.Y[outcome_i], each = N)
    score_mat = rbind(Matrix::kronecker(matrix(1, 1, J_i), (D_1 %*% t(F_1))),
                      Matrix::kronecker(diag(J_i), (D_2 %*% t(F_2)))) %*%
      Sigma_inv %*% (Y_i / beta_i)
    for(j in 1:J_i){
      Y_score$subj_score[Y_score$id == ind & Y_score$outcome == outcome_i[j]][[1]] = score_mat[1:(L+R)]
      Y_score$oc_score[Y_score$id == ind & Y_score$outcome == outcome_i[j]][[1]] = score_mat[(L+R+(j-1)*(V_y)+1):(L+R+j*(V_y))]
    }
  }
  
  
  # Return the estimates
  print("The estimated latent dimensions (L, K, R, V^x, V^y) are:")
  print(c(L, K, R, V_x, V_y))
  
  res_ls = list(mu_x, mu_y, s2_x, s2_y, beta.X, beta.Y,
                Psi, Phi, Gamma, tau_x, tau_y,
                Omega_x, Omega_y, Theta, Lambda, 
                varsigma_x, varsigma_y,
                X_score, Y_score)
  names(res_ls) = c("mu_x", "mu_y", "s2_x", "s2_y", "beta.X", "beta.Y",
                    "Psi", "Phi", "Gamma", "tau_x", "tau_y",
                    "Omega_x", "Omega_y", "Theta", "Lambda",
                    "varsigma_x", "varsigma_y",
                    "X_score", "Y_score")
  return(res_ls)
}




weight_fun = function(tobs, regular = TRUE){
  
  #############################################################################
  ## Description: (supporting function) This function returns the trapezoidal weights based on 
  ##              the input time grid.
  ## Args:        tobs: total time grid T (vector, N*1).
  ##              regular: indicator whether the input time grid is equidistant or not. (logical).
  ## Returns:     an array of trapezoidal weights.
  #############################################################################
  
  n <- length(tobs)
  if(regular){
    # grid step
    dt <- (tail(tobs, n = 1) - tobs[1]) / (n - 1)
    # trapezoidal weights
    w <- rep(1, n); w[c(1, n)] <- 0.5; w <- w * dt
  }else{
    stopifnot(n >= 2)
    w <- numeric(n)
    w[1] <- (tobs[2] - tobs[1]) / 2
    w[n] <- (tobs[n] - tobs[n-1]) / 2
    if (n > 2) {
      w[2:(n-1)] <- (tobs[3:n] - tobs[1:(n-2)]) / 2
    }
  }
  return(w)
}




extract_MFMM <- function(Chat_new, N, J, tobs) {
  
  #############################################################################
  ## Description: (supporting function) This function refers to the MFMM in "Joint model for survival and multivariate  
  ##              sparse functional data with application to a study of Alzheimer’s Disease" (Li et al. 2020) and
  ##              is adapted to only estimate the scaling parameters beta_j and the subject- and outcome-level covariance functions.
  ## Args:        Chat_new: cross-covariance functions obtained from the fast covariance estimation method (matrix, NJ*NJ).
  ##              N: total number of distinct times across all subjects and outcomes (scalar)
  ##              J: total number of outcomes (scalar)
  ##              tobs: total time grid T (vector, N*1).
  ## Returns:     a list containing:
  ##                beta: estimate of the outcome-specific scaling parameters (vector, J*1)
  ##                C0: estimate of subject-level covariance surface (matrix, N*N)
  ##                C1: estimate of outcome-level covariance surface. (matrix, N*N)
  #############################################################################

  stopifnot(nrow(Chat_new) == N * J, ncol(Chat_new) == N * J)
  stopifnot(length(tobs) == N)
  wt <- weight_fun(tobs)
  
  # Extract all auto- or cross-covariance blocks 
  extract_block <- function(C_full, j, jp, N) {
    rows <- ((j  - 1) * N + 1):(j  * N)
    cols <- ((jp - 1) * N + 1):(jp * N)
    C_full[rows, cols]
  }
  Cblk <- vector("list", J)
  for (j in 1:J) {
    Cblk[[j]] <- vector("list", J)
    for (jp in 1:J) {
      Cblk[[j]][[jp]] <- extract_block(Chat_new, j, jp, N)
    }
  }
  
  
  # Estimate beta_j
  beta <- rep(NA_real_, J)
  beta[1] <- 1 # set beta_1 to be 1 for identifiability
  
  weighted_L2_inner <- function(A, B, wt) {
    W2 <- wt %o% wt # W2 = w %o% w is the N x N outer product of weights
    sum(W2 * A * B) # sum(W2 * A * B) = sum over all (i,j) of w_i * w_j * A_ij * B_ij
  }
  
  if(J == 2){
    # J = 2: ratio approach
    C11 <- Cblk[[1]][[1]]
    C22 <- Cblk[[2]][[2]]
    C12 <- Cblk[[1]][[2]]
    
    beta2_sq   <- weighted_L2_inner(C11, C22, wt) / weighted_L2_inner(C11, C11, wt)
    beta2_sign <- sign(sum(diag(C12)))
    beta[2]    <- beta2_sign * sqrt(abs(beta2_sq))
  }else{
    # J >= 3: pool over j' != j, j' != 1
    for (j in 2:J) {
      numer <- 0
      denom <- 0
      for (jp in 2:J) {
        if (jp == j) next
        Cjjp <- Cblk[[j]][[jp]]
        Cjp1 <- Cblk[[jp]][[1]]
        numer <- numer + weighted_L2_inner(Cjjp, Cjp1, wt)
        denom <- denom + weighted_L2_inner(Cjp1, Cjp1, wt)
      }
      beta[j] <- numer / denom
    }
  }
  

  # Estimate C_0
  C0_numer <- matrix(0, N, N)
  C0_denom <- 0
  for (j in 1:J) {
    for (jp in 1:J) {
      if (j == jp) next
      w <- beta[j] * beta[jp]
      C0_numer <- C0_numer + w * Cblk[[j]][[jp]]
      C0_denom <- C0_denom + w^2
    }
  }
  C0_raw <- C0_numer / C0_denom
  C0_raw <- as.matrix(forceSymmetric(C0_raw))
  
  
  functional_eigen_recon <- function(S, wt) {
    N <- nrow(S)
    sqrt_w <- sqrt(wt)
    S_tilde <- sweep(sweep(S, 1, sqrt_w, "*"), 2, sqrt_w, "*")
    eig <- eigen(S_tilde, symmetric = TRUE)
    pos_idx <- eig$values > 0 # keep only positive eigenvalues (ensure PSD)
    if (!any(pos_idx)) {stop("No positive eigenvalues.")}
    
    # Reconstruct covariance
    vals <- eig$values[pos_idx]
    psi  <- eig$vectors[, pos_idx, drop = FALSE]
    npc <- which.max((cumsum(vals) / sum(vals)) >= .99)
    vals <- vals[1:npc]
    psi  <- psi[, 1:npc, drop = FALSE]
    phi <- sweep(psi, 1, (1 / sqrt_w), "*")
    S_recon <- phi %*% diag(vals, nrow = npc, ncol = npc) %*% t(phi)
    S_recon
  }
  C0_fit <- functional_eigen_recon(C0_raw, wt) # Reconstruct C0 to ensure PSD
  

  # Estimate C_1
  C1_raw <- matrix(0, N, N)
  for (j in 1:J) {
    C1_raw <- C1_raw + (Cblk[[j]][[j]] / beta[j]^2 - C0_fit)
  }
  C1_raw <- C1_raw / J
  C1_raw <- as.matrix(forceSymmetric(C1_raw))
  C1_fit <- functional_eigen_recon(C1_raw, wt) # Reconstruct C1 to ensure PSD
  
  list(beta = beta, C0 = C0_fit, C1 = C1_fit)
}










