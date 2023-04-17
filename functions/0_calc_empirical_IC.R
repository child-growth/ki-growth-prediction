
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