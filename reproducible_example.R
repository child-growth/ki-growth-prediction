
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/functions/0_calc_empirical_IC.R"))
library(hal9001)
library(washb)


#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS")) %>% filter(studyid=="MAL-ED", country=="INDIA")

#Subset to control arms
#d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")



y_name = "whz_24" 
x_names = c("sex", "W_mage", "W_meducyrs", "W_nrooms", "W_mwtkg","W_mhtcm","W_mbmi",
            "miss_W_mage","miss_W_meducyrs", "miss_W_nrooms",   "miss_W_mwtkg",    "miss_W_mhtcm",    "miss_W_mbmi" , "arm", "waz_birth")
dfull <- d %>% ungroup() %>% select(all_of(c("studyid","country",y_name, x_names)))

#temp: complete cases
dim(dfull)
table(!is.na(dfull$waz_birth))
dfull <- dfull[complete.cases(dfull),]
dfull <- droplevels(dfull)
dim(dfull)


n=nrow(dfull)
  Y <- as.numeric(as.matrix(dfull %>% select(all_of(y_name))))
  X <- dfull %>% ungroup() %>%
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.character), as.factor)
  
  NZV <- nearZeroVar(X, saveMetrics = TRUE)
  X <- X[,!NZV$zeroVar]
  
  X <-design_matrix(as.data.frame(X))
  


# hal_fit_list <- readRDS(paste0(model_dir,"res_wazBirth_whz24_hal.RDS"))
# 
# hal_fit<-hal_fit_list$`MAL-ED-INDIA`

fit <- fit_hal(X = X,
                    Y = Y,
                    smoothness_orders = 0,
                    return_x_basis = TRUE,
                    family = "gaussian",
                    num_knots = hal9001:::num_knots_generator(
                      max_degree = ifelse(ncol(X) >= 20, 2, 3),
                      smoothness_orders = 0,
                      base_num_knots_0 = max(100, ceiling(sqrt(n)))
                    ))

fit$coefs[fit$coefs!=0]




#t-statistic based mean outcome
library(distributions3)
n <- length(Y)
T_stat <- StudentsT(df = n-1)
mean(Y) + quantile(T_stat, 0.12 / 2) * sd(Y) / sqrt(n)
mean(Y) + quantile(T_stat, 1 - 0.12 / 2) * sd(Y) / sqrt(n)

#washb mean + robust CI
res_mean <- data.frame(washb::washb_mean(Y=Y, id=1:length(Y), print=FALSE))
res_mean


#Get HAL predictions - just with the data used to fit the model
#res_hal <- get_pred_empirical_ci(fit=hal_fit_list[[i]]$fit_init, X=hal_fit_list[[i]]$X, Y=hal_fit_list[[i]]$Y, n=hal_fit_list[[i]]$n))

Y_hat <- predict(fit, new_data = X)
init_coef <- fit$coefs[-1]
nonzero_col <- which(init_coef != 0)
init_coef_nonzero <- init_coef[nonzero_col]
basis_mat <- as.matrix(fit$x_basis)
basis_mat <- as.matrix(basis_mat[, nonzero_col])


#cal_IC_for_beta_cont <- function(X, Y, Y_hat, beta_n){
  # 1. calculate score: X'(Y-Y_hat)
  score <- sweep(basis_mat, 1, (Y - Y_hat), `*`)
  # 2. calculate E_{P_n}(X'X)^(-1)
  d_scaler = solve(t(basis_mat) %*% basis_mat)
  # 3. calculate influence curves
  IC_beta <- t(d_scaler %*% t(score))
#   return(IC)
# }

#IC_phi <- cal_IC_for_phi(X_new = basis_mat, beta_n = init_coef_nonzero, IC_beta)
beta_n = init_coef_nonzero
X_new=basis_mat
  d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
  d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  
  IC_phi = diag(d_phi_new %*% t(IC_beta))
  

  
  se_IC <- sqrt(var(IC_phi)/n)
  
  res <- data.frame(predY = mean(Y_hat))
  res$se <- se_IC
  res$lb <- res$predY - 1.96 * res$se
  res$ub <- res$predY + 1.96 * res$se
  res





  # get_counterfactual_empirical_ci
  # function(fit, X, Y, n, IC_beta){
  #   
  #   if(!is.null(fit)){
  #     Y_hat_init <- predict(fit, new_data = X)
  #     
  #     init_coef <- fit$coefs[-1]
  #     nonzero_col <- which(init_coef != 0)
  #     init_coef_nonzero <- init_coef[nonzero_col]
  #     
  #     x_basis <- make_design_matrix(as.matrix(X), fit$basis_list, p_reserve = 0.75)
  #     x_basis <- as.matrix(x_basis[, nonzero_col])
  #     
  #     IC_phi <- NULL
  #     try(IC_phi <- cal_IC_for_phi(X_new = x_basis, beta_n = init_coef_nonzero, IC_beta))
  #     
  #     if(!is.null(IC_phi)){
  #       se_IC <- sqrt(var(IC_phi)/n)
  #       
  #       res <- data.frame(predY = mean(Y_hat_init))
  #       res$se <- se_IC
  #       res$lb <- res$predY - 1.96 * res$se
  #       res$ub <- res$predY + 1.96 * res$se
  #       
  #       return(res)
  #     }else{
  #       return(NULL)
  #     }
  #   }  
  #   
  # }
  
  X1 <- X %>% mutate(waz_birth=waz_birth+0.5)

  
  # function(fit, X, Y, n, IC_beta){
    
      Y_hat_init <- predict(fit, new_data = X1)
      
      init_coef <- fit$coefs[-1]
      nonzero_col <- which(init_coef != 0)
      init_coef_nonzero <- init_coef[nonzero_col]
      
      x_basis1 <- make_design_matrix(as.matrix(X1), fit$basis_list, p_reserve = 0.75)
      x_basis1 <- as.matrix(x_basis1[, nonzero_col])
      
      IC_phi <- cal_IC_for_phi(X_new = x_basis1, beta_n = init_coef_nonzero, IC_beta)
      
        se_IC <- sqrt(var(IC_phi)/n)
        
        res_counterfactual <- data.frame(predY = mean(Y_hat_init))
        res_counterfactual$se <- se_IC
        res_counterfactual$lb <- res_counterfactual$predY - 1.96 * res_counterfactual$se
        res_counterfactual$ub <- res_counterfactual$predY + 1.96 * res_counterfactual$se
        
        res
        res_counterfactual

  