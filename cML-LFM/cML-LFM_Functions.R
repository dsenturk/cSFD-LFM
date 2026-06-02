###############################################################################
## Description: Functions for conducting the proposed cML-LFM estimation algorithm   
##              described in Algorithm 1 of 'Contrastive Latent Models for Structured Functional Data'.
###############################################################################
## Functions included:
## Main function:
##    cML_LFM: Function for fitting a cML-LFM model and executing the proposed estimation
##             algorithm described in Algorithm 1.
## Supporting functions used by the main function:
##    1. weight_fun: Supporting function that returns the trapezoidal weights based on the given time grid.
##    2. mfpca: Supporting function that conducts marginal multilevel FPCA.
###############################################################################

cML_LFM = function(X.df, Y.df, tobs = NULL, 
                   pve_X = c(.90, .90), pve_Y = c(.90, .90)){
  
  #############################################################################
  ## Description: (main function) This function executes the proposed cML-LFM estimation algorithm described in Algorithm 1.
  ##              The two input dataframes (X.df, Y.df) are constructed with three variables, id, visit, and data,
  ##              representing N measurements of id-th subject at the visit-th repetition.
  ## Args:        X.df: dataframe of the first group of data (three variables: id, visit, data).
  ##              Y.df: dataframe of the second group of data (three variables: id, visit, data).
  ##              tobs: total time grid T (vector, N*1).
  ##              pve_X: proportion of variance explained for retaining leading eigenfunctions of X in Stage I
  ##              pve_Y: proportion of variance explained for retaining leading eigenfunctions of Y in Stage I
  ## Returns:     A list containing:
  ##                mu_x: estimate of the repetition-specific mean function mu^x_j(t) (matrix, N*J)
  ##                mu_y: estimate of the repetition-specific mean function mu^y_j(t) (matrix, N*J)
  ##                s2_x: estimate of the error variance sigma^2_X (scalar)
  ##                s2_y: estimate of the error variance sigma^2_Y (scalar)
  ##                Psi_1: estimate of the level-1 shared latent components {psi^{(1)}_l(t), l = 1,...,L_1} (matrix, N*L_1).
  ##                Phi_1: estimate of the level-1 unique latent components for group X {phi^{(1)}_k(t), k = 1,...,K_1} (matrix, N*K_1).
  ##                Gamma_1: estimate of the level-1 unique latent components for group Y {gamma^{(1)}_r(t), r = 1,...,R_1} (matrix, N*R_1).
  ##                Psi_2: estimate of the level-2 shared latent components {psi^{(2)}_l(t), l = 1,...,L_2} (matrix, N*L_2).
  ##                Phi_2: estimate of the level-2 unique latent components for group X {phi^{(2)}_k(t), k = 1,...,K_2} (matrix, N*K_2).
  ##                Gamma_2: estimate of the level-2 unique latent components for group Y {gamma^{(2)}_r(t), r = 1,...,R_2} (matrix, N*R_2).
  ##                Omega_x_1: diagonal matrix of estimates of the score variances {omega^{x(1)}_l, l = 1,...,L_1} (matrix, L_1*L_1).
  ##                Omega_y_1: diagonal matrix of estimates of the score variances {omega^{y(1)}_l, l = 1,...,L_1} (matrix, L_1*L_1).
  ##                Theta_1: diagonal matrix of estimates of the score variances {theta^{(1)}_k, k = 1,...,K_1} (matrix, K_1*K_1).
  ##                Lambda_1: diagonal matrix of estimates of the score variance {lambda^{(1)}_r, r = 1,...,R_1} (matrix, R_1*R_1).
  ##                Omega_x_2: diagonal matrix of estimates of the score variances {omega^{x(2)}_l, l = 1,...,L_2} (matrix, L_2*L_2).
  ##                Omega_y_2: diagonal matrix of estimates of the score variances {omega^{y(2)}_l, l = 1,...,L_2} (matrix, L_2*L_2).
  ##                Theta_2: diagonal matrix of estimates of the score variances {theta^{(2)}_k, k = 1,...,K_2} (matrix, K_2*K_2).
  ##                Lambda_2: diagonal matrix of estimates of the score variance {lambda^{(2)}_r, r = 1,...,R_2} (matrix, R_2*R_2).
  ##                X_score: BLUP estimate of the subject-specific latent scores for group X (dataframe: id, visit, level-1 score, level-2 score).
  ##                Y_score: BLUP estimate of the subject-specific latent scores for group Y (dataframe: id, visit, level-1 score, level-2 score).
  #############################################################################
  
  stopifnot((!is.null(X.df$data) && !is.null(X.df$id) && !is.null(Y.df$data) && !is.null(Y.df$id)))
  
  if(NCOL(X.df$data) != NCOL(Y.df$data)){
    stop(paste("The total distinct time points are inconsistent between two",
               "groups of data."))
  }else{
    N = NCOL(X.df$data) # total number of time points
    if(is.null(tobs)){tobs = seq(0, 1, length.out = N)}  # total time grid T (vector, N*1).
  }
  if (!is.null(X.df$visit)){
    X.df$visit = as.integer(factor(X.df$visit))
  }else{X.df$visit = ave(X.df$id, X.df$id, FUN=seq_along)}
  if (!is.null(Y.df$visit)){
    Y.df$visit = as.integer(factor(Y.df$visit))
  }else{Y.df$visit = ave(Y.df$id, Y.df$id, FUN=seq_along)}
  J = length(unique(X.df$visit))  ## total number of repetitions/visits
  n_x = length(unique(X.df$id)) ## number of subjects in X
  n_y = length(unique(Y.df$id)) ## number of subjects in Y
  nVisits_x = data.frame(table(X.df$id))  ## number of repetitions/visits for each subject X_i
  colnames(nVisits_x) = c("id", "numVisits")
  nVisits_y = data.frame(table(Y.df$id))  ## number of repetitions/visits for each subject Y_i
  colnames(nVisits_y) = c("id", "numVisits")
  
  
  
  #################################
  # Step 1. Estimate the mean functions mu^x_j(t) and mu^y_j(t) (by penalized smoothing spline), 
  #         error variances sigma^2_X and sigma^2_Y, and between- and within-subject covariances.
  #################################
  
  # X
  ## Estimate mu^x(t) and sigma^2_X
  mfpca_res = mfpca(X.df, tobs)
  
  mu_x = mfpca_res$mu.j  ## matrix to store repetition-specific means
  s2_x = mfpca_res$sigma2
  if(s2_x <= 0){stop("Estimate of error variance for X is non-positive.")}
  ## Centralize the first group of data, (X(t) - mu^x(t))
  for(ind in 1:NROW(X.df)) {
    X.df$data[ind,] = X.df$data[ind,] - mu_x[,X.df$visit[ind]]
  }
  ## Estimate level-1 and level-2 covariance functions
  G_x_1 = mfpca_res$G_B
  G_x_2 = mfpca_res$G_W
  

  # Y
  ## Estimate mu^y(t) and sigma^2_Y
  mfpca_res = mfpca(Y.df, tobs)
  
  mu_y = mfpca_res$mu.j  
  s2_y = mfpca_res$sigma2
  if(s2_y <= 0){stop("Estimate of error variance for Y is non-positive.")}
  ## Centralize the second group of data, (Y(t) - mu^y(t))
  for(ind in 1:NROW(Y.df)) {
    Y.df$data[ind,] = Y.df$data[ind,] - mu_y[,Y.df$visit[ind]]
  }
  ## Estimate level-1 and level-2 covariance functions
  G_y_1 = mfpca_res$G_B
  G_y_2 = mfpca_res$G_W
  
  
  
  #################################
  # Step 2. Perform FPCA on between- and within-subject covariances to acquire leading eigenfunctions
  #################################
  
  W = diag(x = weight_fun(tobs), nrow = N) # construct the trapezoidal weight functions
  Wsqrt = diag(x = sqrt(diag(W)), nrow = N)
  Winvsqrt = diag(x = 1/sqrt(diag(W)), nrow = N)
  
  # X
  G_x = list(level1 = G_x_1, level2 = G_x_2) # put G_B and G_W together.
  V = lapply(G_x, function(x) Wsqrt %*% x %*% Wsqrt)
  evalues = lapply(V, function(x) eigen(x, symmetric = TRUE, 
                                        only.values = TRUE)$values)
  evalues = lapply(evalues, function(x) replace(x, which(x <= 0), 0)) 
  ## Determine the number of retained eigenfunctions
  npc = list()
  npc = append(npc, min(which(cumsum(evalues[[1]])/sum(evalues[[1]]) > pve_X[1])))
  npc = append(npc, min(which(cumsum(evalues[[2]])/sum(evalues[[2]]) > pve_X[2])))
  names(npc) = names(evalues)
  
  ## Calculate the retained eigenfunctions
  efunctions = lapply(names(V), function(x) 
    matrix(Winvsqrt %*% eigen(V[[x]], symmetric = TRUE)$vectors[, seq(len = npc[[x]])], 
           nrow = N, ncol = npc[[x]]))
  names(efunctions) = c("level1", "level2")
  
  Varphi_x_1 = efunctions$level1
  Varphi_x_2 = efunctions$level2
  M_x_1 = NCOL(Varphi_x_1)
  M_x_2 = NCOL(Varphi_x_2)
  
  
  # Y
  G_y = list(level1 = G_y_1, level2 = G_y_2) 
  V = lapply(G_y, function(x) Wsqrt %*% x %*% Wsqrt)
  evalues = lapply(V, function(x) eigen(x, symmetric = TRUE, 
                                        only.values = TRUE)$values)
  evalues = lapply(evalues, function(x) replace(x, which(x <= 0), 0))
  ## Determine the number of retained eigenfunctions
  npc = list()
  npc = append(npc, min(which(cumsum(evalues[[1]])/sum(evalues[[1]]) > pve_Y[1])))
  npc = append(npc, min(which(cumsum(evalues[[2]])/sum(evalues[[2]]) > pve_Y[2])))
  names(npc) = names(evalues)
  
  ## Calculate the retained eigenfunctions
  efunctions = lapply(names(V), function(x) 
    matrix(Winvsqrt %*% eigen(V[[x]], symmetric = TRUE)$vectors[, seq(len = npc[[x]])], 
           nrow = N, ncol = npc[[x]]))
  names(efunctions) = c("level1", "level2")
  
  Varphi_y_1 = efunctions$level1
  Varphi_y_2 = efunctions$level2
  M_y_1 = NCOL(Varphi_y_1)
  M_y_2 = NCOL(Varphi_y_2)
  
  
  
  #################################
  # Step 3. Construct the group- and level-specific projection operators and their weighted average projection operators
  #################################
  
  # X
  P_x_1 = Wsqrt %*% tcrossprod(Varphi_x_1) %*% Wsqrt 
  P_x_2 = Wsqrt %*% tcrossprod(Varphi_x_2) %*% Wsqrt 
  # Y
  P_y_1 = Wsqrt %*% tcrossprod(Varphi_y_1) %*% Wsqrt 
  P_y_2 = Wsqrt %*% tcrossprod(Varphi_y_2) %*% Wsqrt 
  # Weighted average 
  P_w_1 = (P_x_1 * n_x + P_y_1 * n_y) / (n_x + n_y)
  P_w_2 = (P_x_2 * n_x + P_y_2 * n_y) / (n_x + n_y)
  
  
  
  #################################
  # Step 4. Identify the rank of shared spaces and construct the corresponding shared projection operators
  #################################
  
  # Level-1
  eigen_res = eigen(P_w_1, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  L_1 = sum(eigen_res$values > .9) # identify the rank of the level-1 shared space
  P_s_1 = tcrossprod(as.matrix(eigen_res$vectors[, 1:L_1])) # Construct the corresponding shared projection operators

  
  # Level-2
  eigen_res = eigen(P_w_2, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  L_2 = sum(eigen_res$values > .9) # identify the rank of the level-2 shared space
  P_s_2 = tcrossprod(as.matrix(eigen_res$vectors[, 1:L_2]))
  
  #################################
  # Step 5. Target shared latent components through projection of level-specific covariance functions onto shared subspaces
  #################################
  
  # Level-1
  G_x_1_s = Winvsqrt %*% P_s_1 %*% Wsqrt %*% G_x_1 %*% Wsqrt %*% P_s_1 %*% Winvsqrt
  G_y_1_s = Winvsqrt %*% P_s_1 %*% Wsqrt %*% G_y_1 %*% Wsqrt %*% P_s_1 %*% Winvsqrt
  eigen_res = eigen(Wsqrt %*% (G_x_1_s + G_y_1_s) %*% Wsqrt, symmetric = TRUE)
  Psi_1 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:L_1])
  
  # Level-2
  G_x_2_s = Winvsqrt %*% P_s_2 %*% Wsqrt %*% G_x_2 %*% Wsqrt %*% P_s_2 %*% Winvsqrt
  G_y_2_s = Winvsqrt %*% P_s_2 %*% Wsqrt %*% G_y_2 %*% Wsqrt %*% P_s_2 %*% Winvsqrt
  eigen_res = eigen(Wsqrt %*% (G_x_2_s + G_y_2_s) %*% Wsqrt, symmetric = TRUE)
  Psi_2 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:L_2])
  
  
  #################################
  # Step 6. Construct complementary projection operators and identify dimensions of unique latent components
  #################################
  
  # Level-1
  ## X
  P_x_1.c = P_x_1 %*% (diag(N) - P_s_1) %*% P_x_1
  eigen_res = eigen(P_x_1.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  K_1 = sum(eigen_res$values > .9) # identify the rank of the level-1 unique subspace for X
  if(K_1 > 0){P_x_1_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:K_1]))} # Construct the group-specific complementary projection operators
  
  ## Y
  P_y_1.c = P_y_1 %*% (diag(N) - P_s_1) %*% P_y_1
  eigen_res = eigen(P_y_1.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  R_1 = sum(eigen_res$values > .9) # identify the rank of the level-1 unique subspace for Y
  if(R_1 > 0){P_y_1_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:R_1]))}
  
  # Level-2
  ## X
  P_x_2.c = P_x_2 %*% (diag(N) - P_s_2) %*% P_x_2
  eigen_res = eigen(P_x_2.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  K_2 = sum(eigen_res$values > .9) # identify the rank of the level-2 unique subspace for X
  if(K_2 > 0){P_x_2_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:K_2]))}
  
  ## Y
  P_y_2.c = P_y_2 %*% (diag(N) - P_s_2) %*% P_y_2
  eigen_res = eigen(P_y_2.c, symmetric = TRUE)
  pos_n = sum(eigen_res$values > 0)
  eigen_res$values = eigen_res$values[1:pos_n]
  eigen_res$vectors = as.matrix(eigen_res$vectors[,1:pos_n])
  R_2 = sum(eigen_res$values > .9) # identify the rank of the level-2 unique subspace for Y
  if(R_2 > 0){P_y_2_u = tcrossprod(as.matrix(eigen_res$vectors[, 1:R_2]))}
  
  
  
  #################################
  # Step 7. Construct complementary projection operators and identify dimensions of unique latent components
  #################################
  
  # Level-1
  if(K_1 > 0){
    G_x_1_u = Winvsqrt %*% P_x_1_u %*% Wsqrt %*% G_x_1 %*% Wsqrt %*% P_x_1_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_x_1_u %*% Wsqrt, symmetric = TRUE)
    Phi_1 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:K_1])
  }else{
    Phi_1 = NULL
  }
  if(R_1 > 0){
    G_y_1_u = Winvsqrt %*% P_y_1_u %*% Wsqrt %*% G_y_1 %*% Wsqrt %*% P_y_1_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_y_1_u %*% Wsqrt, symmetric = TRUE)
    Gamma_1 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:R_1])
  }else{
    Gamma_1 = NULL
  }
  
  # Level-2
  if(K_2 > 0){
    G_x_2_u = Winvsqrt %*% P_x_2_u %*% Wsqrt %*% G_x_2 %*% Wsqrt %*% P_x_2_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_x_2_u %*% Wsqrt, symmetric = TRUE)
    Phi_2 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:K_2])
  }else{
    Phi_2 = NULL
  }
  if(R_2 > 0){
    G_y_2_u = Winvsqrt %*% P_y_2_u %*% Wsqrt %*% G_y_2 %*% Wsqrt %*% P_y_2_u %*% Winvsqrt
    eigen_res = eigen(Wsqrt %*% G_y_2_u %*% Wsqrt, symmetric = TRUE)
    Gamma_2 = Winvsqrt %*% as.matrix(eigen_res$vectors[, 1:R_2])
  }else{
    Gamma_2 = NULL
  }
  
  
  
  #################################
  # Step 8. Target variance components and subject-specific latent scores
  #################################
  
  # Variance components
  ## X
  Omega_x_1 = diag(x = diag(as.matrix(t(Psi_1) %*% W %*% G_x_1 %*% W %*% Psi_1)), nrow = L_1)
  Omega_x_2 = diag(x = diag(as.matrix(t(Psi_2) %*% W %*% G_x_2 %*% W %*% Psi_2)), nrow = L_2)
  if(K_1 > 0){Theta_1 = diag(x = diag(as.matrix(t(Phi_1) %*% W %*% G_x_1 %*% W %*% Phi_1)), nrow = K_1)}else{Theta_1 = NULL}
  if(K_2 > 0){Theta_2 = diag(x = diag(as.matrix(t(Phi_2) %*% W %*% G_x_2 %*% W %*% Phi_2)), nrow = K_2)}else{Theta_2 = NULL}
  
  ## Y
  Omega_y_1 = diag(x = diag(as.matrix(t(Psi_1) %*% W %*% G_y_1 %*% W %*% Psi_1)), nrow = L_1)
  Omega_y_2 = diag(x = diag(as.matrix(t(Psi_2) %*% W %*% G_y_2 %*% W %*% Psi_2)), nrow = L_2)
  if(R_1 > 0){Lambda_1 = diag(x = diag(as.matrix(t(Gamma_1) %*% W %*% G_y_1 %*% W %*% Gamma_1)), nrow = R_1)}else{Lambda_1 = NULL}
  if(R_2 > 0){Lambda_2 = diag(x = diag(as.matrix(t(Gamma_2) %*% W %*% G_y_2 %*% W %*% Gamma_2)), nrow = R_2)}else{Lambda_2 = NULL}
  
  # Score estimation
  ## X
  ### Inversion of data covariance matrix Sigma 
  id_mat = diag(N)
  F_1 = cbind(Psi_1, Phi_1)
  F_2 = cbind(Psi_2, Phi_2)
  D_1 = diag(x = c(diag(Omega_x_1), diag(Theta_1)), nrow = (L_1 + K_1))
  D_2 = diag(x = c(diag(Omega_x_2), diag(Theta_2)), nrow = (L_2 + K_2))
  
  ### First, inverse G_W
  G_w_inv = (1/s2_x) * (id_mat - F_2 %*% Matrix::solve(s2_x * Matrix::solve(D_2) + crossprod(F_2)) %*% t(F_2))  

  ### Second, calculate G_B and inverse (J_i * G_B + G_W)
  G_BW_inv = tibble(visit_num = unique(nVisits_x$numVisits),
                    inv_mat = lapply(visit_num, function(i) matrix(NA, nrow = N, ncol = N)))
  
  for (J_i in unique(nVisits_x$numVisits)) {
    G_BW_inv$inv_mat[G_BW_inv$visit_num == J_i][[1]] = G_w_inv %*% 
      (id_mat - F_1 %*% Matrix::solve(1/J_i * Matrix::solve(D_1) + 
                                        t(F_1) %*% G_w_inv %*% F_1) %*% 
         t(F_1) %*% G_w_inv)  
  }
  G_B = F_1 %*% D_1 %*% t(F_1)
  
  ### Third, calculate inverse of covariance matrix Sigma for each subject
  Sigma_inv = tibble(visit_num = unique(nVisits_x$numVisits),
                     inv_mat = vector("list", length(unique(nVisits_x$numVisits))))
  for (J_i in unique(nVisits_x$numVisits)) {
    Sigma_inv$inv_mat[Sigma_inv$visit_num == J_i][[1]] = Matrix::kronecker(diag(J_i), G_w_inv) - 
      Matrix::kronecker(matrix(1, J_i, J_i), (G_w_inv %*% G_B %*% G_BW_inv$inv_mat[G_BW_inv$visit_num == J_i][[1]]))
  }
  
  ### Calculate score estimates for each subject
  X_score = tibble(X.df[,1:2])
  X_score$level_1_score = vector("list", NROW(X_score))
  X_score$level_2_score = vector("list", NROW(X_score))
  for (ind in unique(X.df$id)) {
    X_i = matrix(t(X.df$data[X.df$id == ind,]), ncol = 1) # flatten to a 1-d vector
    J_i = nVisits_x$numVisits[nVisits_x$id == ind]
    visit_i = X.df$visit[X.df$id == ind]
    score_mat = rbind(Matrix::kronecker(matrix(1, 1, J_i), (D_1 %*% t(F_1))), 
                      Matrix::kronecker(diag(J_i), (D_2 %*% t(F_2)))) %*%
      Sigma_inv$inv_mat[Sigma_inv$visit_num == J_i][[1]] %*% X_i
    for(j in 1:J_i){
      X_score$level_1_score[X_score$id == ind & X_score$visit == visit_i[j]][[1]] = score_mat[1:(L_1+K_1)]
      X_score$level_2_score[X_score$id == ind & X_score$visit == visit_i[j]][[1]] = score_mat[(L_1+K_1+(j-1)*(L_2+K_2)+1):(L_1+K_1+j*(L_2+K_2))]
    }
  }
  
  
  ## Y
  ### Inversion of data covariance matrix Sigma 
  id_mat = diag(N)
  F_1 = cbind(Psi_1, Gamma_1)
  F_2 = cbind(Psi_2, Gamma_2)
  D_1 = diag(x = c(diag(Omega_y_1), diag(Lambda_1)), nrow = (L_1 + R_1))
  D_2 = diag(x = c(diag(Omega_y_2), diag(Lambda_2)), nrow = (L_2 + R_2))
  
  ### First, inverse G_W
  G_w_inv = (1/s2_y) * (id_mat - F_2 %*% Matrix::solve(s2_y * Matrix::solve(D_2) + crossprod(F_2)) %*% t(F_2))  
  
  ### Second, calculate G_B and inverse (J_i * G_B + G_W)
  G_BW_inv = tibble(visit_num = unique(nVisits_y$numVisits),
                    inv_mat = lapply(visit_num, function(i) matrix(NA, nrow = N, ncol = N)))
  
  for (J_i in unique(nVisits_y$numVisits)) {
    G_BW_inv$inv_mat[G_BW_inv$visit_num == J_i][[1]] = G_w_inv %*% 
      (id_mat - F_1 %*% Matrix::solve(1/J_i * Matrix::solve(D_1) + 
                                        t(F_1) %*% G_w_inv %*% F_1) %*% 
         t(F_1) %*% G_w_inv)  
  }
  G_B = F_1 %*% D_1 %*% t(F_1)
  
  ### Third, calculate inverse of covariance matrix Sigma for each subject
  Sigma_inv = tibble(visit_num = unique(nVisits_y$numVisits),
                     inv_mat = vector("list", length(unique(nVisits_y$numVisits))))
  for (J_i in unique(nVisits_y$numVisits)) {
    Sigma_inv$inv_mat[Sigma_inv$visit_num == J_i][[1]] = Matrix::kronecker(diag(J_i), G_w_inv) - 
      Matrix::kronecker(matrix(1, J_i, J_i), (G_w_inv %*% G_B %*% G_BW_inv$inv_mat[G_BW_inv$visit_num == J_i][[1]]))
  }
  
  ### Calculate score estimates for each subject
  Y_score = tibble(Y.df[,1:2])
  Y_score$level_1_score = vector("list", NROW(Y_score))
  Y_score$level_2_score = vector("list", NROW(Y_score))
  for (ind in unique(Y.df$id)) {
    Y_i = matrix(t(Y.df$data[Y.df$id == ind,]), ncol = 1)
    J_i = nVisits_y$numVisits[nVisits_y$id == ind]
    visit_i = Y.df$visit[Y.df$id == ind]
    score_mat = rbind(Matrix::kronecker(matrix(1, 1, J_i), (D_1 %*% t(F_1))), 
                      Matrix::kronecker(diag(J_i), (D_2 %*% t(F_2)))) %*%
      Sigma_inv$inv_mat[Sigma_inv$visit_num == J_i][[1]] %*% Y_i
    for(j in 1:J_i){
      Y_score$level_1_score[Y_score$id == ind & Y_score$visit == visit_i[j]][[1]] = score_mat[1:(L_1+R_1)]
      Y_score$level_2_score[Y_score$id == ind & Y_score$visit == visit_i[j]][[1]] = score_mat[(L_1+R_1+(j-1)*(L_2+R_2)+1):(L_1+R_1+j*(L_2+R_2))]
    }
  }
  
  print("The estimated latent dimensions (L_1, K_1, R_1, L_2, K_2, R_2) are:")
  print(c(L_1, K_1, R_1, L_2, K_2, R_2))
  
  res_ls = list(mu_x, mu_y, s2_x, s2_y, 
                Psi_1, Phi_1, Gamma_1, Psi_2, Phi_2, Gamma_2, 
                Omega_x_1, Omega_y_1, Theta_1, Lambda_1, 
                Omega_x_2, Omega_y_2, Theta_2, Lambda_2, 
                X_score, Y_score)
  names(res_ls) = c("mu_x", "mu_y", "s2_x", "s2_y", 
                    "Psi_1", "Phi_1", "Gamma_1", "Psi_2", "Phi_2", "Gamma_2",
                    "Omega_x_1", "Omega_y_1", "Theta_1", "Lambda_1",
                    "Omega_x_2", "Omega_y_2", "Theta_2", "Lambda_2",
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
  
  n = length(tobs)
  if(regular){
    # grid step
    dt = (tail(tobs, n = 1) - tobs[1]) / (n - 1)
    # trapezoidal weights
    w = rep(1, n); w[c(1, n)] = 0.5; w = w * dt
  }else{
    stopifnot(n >= 2)
    w = numeric(n)
    w[1] = (tobs[2] - tobs[1]) / 2
    w[n] = (tobs[n] - tobs[n-1]) / 2
    if (n > 2) {
      w[2:(n-1)] = (tobs[3:n] - tobs[1:(n-2)]) / 2
    }
  }
  return(w)
}



mfpca = function (Y.df, tobs){
  
  #############################################################################
  ## Description: (supporting function) This function refers to multilevel FPCA in "Multilevel functional principal
  ##              component analysis" (Di et al. 2009) and is adapted to only estimate the mean functions, error variances 
  ##              and between- and within-subject covariance functions.
  ## Args:        Y.df: dataframe being analysed (three variables: id, visit, data).
  ##              tobs: total time grid T (vector, N*1).
  ## Returns:     a list containing:
  ##                mu.j: estimate of the repetition-specific mean function (matrix, N*J)
  ##                sigma2: estimate of the error variance sigma^2 (scalar)
  ##                G_T: estimate of total covariance surface (matrix, N*N)
  ##                G_B: estimate of between-subject covariance surface. (matrix, N*N)
  ##                G_W: estimate of within-subject covariance surface. (matrix, N*N)
  #############################################################################
  
  J = length(unique(Y.df$visit))
  n = length(unique(Y.df$id))
  N = NCOL(Y.df$data)
  nVisits = data.frame(table(Y.df$id))
  colnames(nVisits) = c("id", "numVisits")
  d.vec = rep(tobs, each = NROW(Y.df$data))
  
  ######## Calculate repetition-specific mean deviations.
  mu.j = matrix(0, N, J)
  Y.split = as.list(rep(NA, J))
  for (j in 1:J) {
    Y.split[[j]] = subset(Y.df, visit == j)
    d.vec.split = rep(tobs, each = NROW(Y.split[[j]]))
    fit.eta = gam(as.vector(Y.split[[j]]$data) ~ s(d.vec.split, k = 4))
    mu.j[, j] = predict(fit.eta, newdata = data.frame(d.vec.split = tobs))
    Y.split[[j]]$Y.tilde = Y.split[[j]]$data - 
      matrix(mu.j[, j], NROW(Y.split[[j]]), N, byrow = TRUE)
  }
  Y.df.new = Reduce(function(...) merge(..., by = c("id", "visit", "data", 
                                                    "Y.tilde"), 
                                        all = TRUE, sort = FALSE), Y.split)
  Y.df.new = Y.df.new[order(Y.df.new$id, Y.df.new$visit), ]
  Y.tilde = Y.df.new$Y.tilde
  
  
  ######## Estimation of the total covariance function
  cov.sum = cov.count = cov.mean = matrix(0, N, N) # create three containers
  row.ind = 0
  for (m in 1:n) {
    for (j in 1:nVisits[m, 2]) {
      row.ind.temp = row.ind + j
      obs.points = which(!is.na(Y.df$data[row.ind.temp, ])) # filter out those non-missing data points
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1 # count the number of non-missing data points for each entry ij.
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + 
        tcrossprod(Y.tilde[row.ind.temp, obs.points]) # calculate the entry ij of the cov matrix by summing all non-missing values
    }
    row.ind = row.ind + nVisits[m, 2]
  }
  G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count) # run over all cov entries and divide them with their counts of non-missing obs
  diag.G0 = diag(G.0) # save the diagonal of the cov
  diag(G.0) = NA # remove the diagonal
  
  # 2D smoothing of this cov
  row.vec = rep(tobs, each = N)
  col.vec = rep(tobs, N)
  s.npc.0 = predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = 7), 
                        weights = as.vector(cov.count)), 
                    newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
  # Construct the total covariance function npc.0 (i.e. G_T)
  npc.0 = matrix(s.npc.0, N, N)
  npc.0 = (npc.0 + t(npc.0))/2
  
  
  ######## Estimation of the between-subject covariance function
  cov.sum = cov.count = cov.mean = matrix(0, N, N)
  row.ind = 0
  ids.KB = nVisits[nVisits$numVisits > 1, c("id")]
  for (m in 1:n) {
    if (Y.df$id[m] %in% ids.KB) {
      for (j in 1:(nVisits[m, 2] - 1)) {
        row.ind1 = row.ind + j
        obs.points1 = which(!is.na(Y.df$data[row.ind1, ]))
        for (k in (j + 1):nVisits[m, 2]) {
          row.ind2 = row.ind + k
          obs.points2 = which(!is.na(Y.df$data[row.ind2, ]))
          cov.count[obs.points1, obs.points2] = 
            cov.count[obs.points1, obs.points2] + 1
          cov.sum[obs.points1, obs.points2] = 
            cov.sum[obs.points1, obs.points2] + tcrossprod(Y.tilde[row.ind1, obs.points1], 
                                                           Y.tilde[row.ind2, obs.points2])
          cov.count[obs.points2, obs.points1] = 
            cov.count[obs.points2, obs.points1] + 1
          cov.sum[obs.points2, obs.points1] = 
            cov.sum[obs.points2, obs.points1] + tcrossprod(Y.tilde[row.ind2, obs.points2], 
                                                           Y.tilde[row.ind1, obs.points1])
        }
      }
    }
    row.ind = row.ind + nVisits[m, 2]
  }
  G.0b = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  row.vec = rep(tobs, each = N)
  col.vec = rep(tobs, N)
  s.npc.0b = predict(gam(as.vector(G.0b) ~ te(row.vec, col.vec, k = 10), 
                         weights = as.vector(cov.count)), 
                     newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
  # npc.0b is the estimated between-subject covariance function (i.e. G_B)
  npc.0b = matrix(s.npc.0b, N, N)
  npc.0b = (npc.0b + t(npc.0b))/2
  
  # npc.0w is the estimated within-subject covariance function (i.e. G_W)
  s.npc.0w = s.npc.0 - s.npc.0b
  npc.0w = matrix(s.npc.0w, N, N)
  npc.0w = (npc.0w + t(npc.0w))/2
  
  ######## Estimate the error variance.
  T.len = tobs[N] - tobs[1]
  T1.min = min(which(tobs >= tobs[1] + 0.25 * T.len))
  T1.max = max(which(tobs <= tobs[N] - 0.25 * T.len))
  DIAG = (diag.G0 - diag(npc.0b) - diag(npc.0w))[T1.min:T1.max]
  w2 = weight_fun(tobs[T1.min:T1.max])
  sigma2 = max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
  
  
  ret.objects = c("mu.j", "sigma2", "npc.0", "npc.0b", "npc.0w")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = c("mu.j", "sigma2", "G_T", "G_B", "G_W")
  class(ret) = "mfpca"
  return(ret)
}






