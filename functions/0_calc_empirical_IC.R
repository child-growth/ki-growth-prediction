
cal_IC_for_beta <- function(X, Y, Y_hat, beta_n){
  n <- dim(X)[1]
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  # 1. calculate score: X(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(X, 1, res, `*`)
  # 2. calculate the derivative of phi:
  d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
  d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  # 4. calculate influence curves
  IC <- t(tmat %*% t(score))
  return(IC)
}

cal_IC_for_EY <- function(X_new, beta_n, IC_beta){
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
  d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  IC = diag(d_phi_new %*% t(IC_beta))
  return(IC)
}

cal_IC_for_ATE <- function(X_new_a, X_new_0, beta_n, IC_beta){
  if (!is.matrix(X_new_a)) X_new_a <- as.matrix(X_new_a)
  if (!is.matrix(X_new_0)) X_new_0 <- as.matrix(X_new_0)
  d_phi_scaler_new_a <- as.vector(exp(- beta_n %*% t(X_new_a)) / ((1 + exp(- beta_n %*% t(X_new_a)))^2))
  d_phi_new_a <- sweep(X_new_a, 1, d_phi_scaler_new_a, `*`)
  d_phi_scaler_new_0 <- as.vector(exp(- beta_n %*% t(X_new_0)) / ((1 + exp(- beta_n %*% t(X_new_0)))^2))
  d_phi_new_0 <- sweep(X_new_0, 1, d_phi_scaler_new_0, `*`)
  d_phi_new <- d_phi_new_a - d_phi_new_0
  IC = diag(d_phi_new %*% t(IC_beta))
  return(IC)
}

cal_IC_for_phi <- function(X_new, beta_n, IC_beta){
  
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
  d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  
  IC = diag(d_phi_new %*% t(IC_beta))
  
  return(IC)
}


get_pred_empirical_ci <- function(fit, X, Y, n){
  
  if(!is.null(fit)){
    Y_hat_init <- predict(fit, new_data = X)
    init_coef <- fit$coefs[-1]
    nonzero_col <- which(init_coef != 0)
    init_coef_nonzero <- init_coef[nonzero_col]
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    
    if(length(unique(Y))>2){
      IC_beta <- get_ic_beta_cont(X = basis_mat, Y = Y, Y_hat = Y_hat_init,
                                 beta_n = init_coef_nonzero)
    }else{
      IC_beta <- cal_IC_for_beta(X = basis_mat, Y = Y, Y_hat = Y_hat_init,
                                 beta_n = init_coef_nonzero)      
    }

    IC_phi <- NULL
    try(IC_phi <- cal_IC_for_phi(X_new = basis_mat, beta_n = init_coef_nonzero, IC_beta))
    
    if(!is.null(IC_phi)){
      se_IC <- sqrt(var(IC_phi)/n)
      
      res <- data.frame(predY = mean(Y_hat_init))
      res$se <- se_IC
      res$lb <- res$predY - 1.96 * res$se
      res$ub <- res$predY + 1.96 * res$se

      return(res)
    }else{
      return(NULL)
    }
  }  
  
}

#ic beta function
get_ic_beta <- function(fit, X, Y, n){
  
  if(!is.null(fit)){
    Y_hat_init <- predict(fit, new_data = X)
    init_coef <- fit$coefs[-1]
    nonzero_col <- which(init_coef != 0)
    init_coef_nonzero <- init_coef[nonzero_col]
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    IC_beta <- cal_IC_for_beta(X = basis_mat, Y = Y, Y_hat = Y_hat_init,
                               beta_n = init_coef_nonzero)

    return(IC_beta)
  }else{
    return(NULL)
  }  
}


get_counterfactual_empirical_ci <- function(fit, X, Y, n, IC_beta){
  
  if(!is.null(fit)){
    Y_hat_init <- predict(fit, new_data = X)
    
    init_coef <- fit$coefs[-1]
    nonzero_col <- which(init_coef != 0)
    init_coef_nonzero <- init_coef[nonzero_col]
    
    x_basis <- make_design_matrix(as.matrix(X), fit$basis_list, p_reserve = 0.75)
    x_basis <- as.matrix(x_basis[, nonzero_col])

    IC_phi <- NULL
    try(IC_phi <- cal_IC_for_phi(X_new = x_basis, beta_n = init_coef_nonzero, IC_beta))
    
    if(!is.null(IC_phi)){
      se_IC <- sqrt(var(IC_phi)/n)
      
      res <- data.frame(predY = mean(Y_hat_init))
      res$se <- se_IC
      res$lb <- res$predY - 1.96 * res$se
      res$ub <- res$predY + 1.96 * res$se
      
      return(res)
    }else{
      return(NULL)
    }
  }  
  
}




# init_coef <- fit_init$coefs[-1]
# nonzero_col <- which(init_coef != 0)
# init_coef_nonzero <- init_coef[nonzero_col]
# basis_mat <- as.matrix(fit_init$x_basis)
# basis_mat <- as.matrix(basis_mat[, nonzero_col])
# 
# IC_phi <- cal_IC_for_phi(X_new = basis_mat, 
#                          beta_n = init_coef_nonzero, IC_beta)
# se_IC <- sqrt(var(IC_phi)/n)
# 
# 
# X_a1 <- X
# X_a1$A = 1



########################
# calculating efficient influence curves for continuous outcome
########################


get_ic_beta_cont <- function(fit, X, Y, n){
  
  if(!is.null(fit)){
    Y_hat_init <- predict(fit, new_data = X)
    if(length(unique(Y_hat_init))>1){
    init_coef <- fit$coefs[-1]
    nonzero_col <- which(init_coef != 0)
    init_coef_nonzero <- init_coef[nonzero_col]
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    IC_beta <- NULL
    try(IC_beta <- cal_IC_for_beta_cont(X = basis_mat, Y = Y, Y_hat = Y_hat_init, beta_n = init_coef_nonzero))
    
      return(IC_beta)
    }else{
      return(NA)
    }
    
  }else{
    return(NA)
  }  
}


cal_IC_for_beta_cont <- function(X, Y, Y_hat, beta_n){
  
  if (!is.matrix(X)) X <- as.matrix(X)
  # 1. calculate score: X'(Y-Y_hat)
  score <- sweep(X, 1, (Y - Y_hat), `*`)
  # 2. calculate E_{P_n}(X'X)^(-1)
  d_scaler = solve(t(X) %*% X)
  # 3. calculate influence curves
  IC <- d_scaler %*% t(score)
  return(IC)
}

cal_IC_for_ATE_cont <- function(Xa_X0_diff, IC_beta){
  if (!is.matrix(Xa_X0_diff)) Xa_X0_diff <- as.matrix(Xa_X0_diff)
  IC <- Xa_X0_diff %*% IC_beta
  return(IC)
}

