##########################################
# ki longitudinal manuscripts
# growth failure prediction analysis

# Clean primary analysis datasets from 3 epi
# manuscripts into single dataset for prediction
# analysis by subsetting to included studies
# and merge outcome-specific datasets together
##########################################

# Mode function
Mode <- function(x, na.rm = TRUE){
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}



impute_missing_values <- function (data, type = "standard", add_indicators = TRUE, prefix = "miss_", 
          skip_vars = NULL, all_vars = FALSE, remove_constant = TRUE, 
          remove_collinear = TRUE, values = NULL, h2o_glrm = NULL, 
          glrm_k = 10L, verbose = FALSE){
  missing_indicators = NULL
  new_data = data
  non_skipped_vars = !colnames(data) %in% skip_vars
  results = list(type = type, add_indicators = add_indicators, 
                 skip_vars = skip_vars, prefix = prefix)
  any_nas = which(sapply(colnames(data), function(col) !col %in% 
                           skip_vars && anyNA(data[[col]])))
  if (verbose) {
    cat("Found", length(any_nas), "variables with NAs.\n")
  }
  if (type == "standard") {
    if (verbose) {
      cat("Running standard imputation.\n")
    }
    impute_values = vector("list", sum(non_skipped_vars))
    names(impute_values) = colnames(data[non_skipped_vars])
    if (all_vars) {
      loop_over = which(non_skipped_vars)
      names(loop_over) = colnames(data)[non_skipped_vars]
    }
    else {
      loop_over = any_nas
    }
    sum_nas = sapply(loop_over, function(col_i) sum(is.na(data[[col_i]])))
    col_classes = sapply(loop_over, function(col_i) class(data[[col_i]]))
    for (i in loop_over) {
      colname = names(loop_over)[loop_over == i]
      nas = sum_nas[colname]
      col_class = col_classes[colname]
      if (verbose) {
        cat("Imputing", colname, paste0("(", i, " ", 
                                        col_class, ")"), "with", prettyNum(nas, big.mark = ","), 
            "NAs.")
      }
      if (colname %in% names(values)) {
        impute_value = values[[colname]]
        if (verbose) {
          cat(" Pre-filled.")
        }
      }
      else if (col_class %in% c("factor")) {
        impute_value = Mode(data[[i]], exclude_na = TRUE)[1]
      }
      else if (col_class %in% c("integer", "numeric", "logical", 
                                "labelled", "integer64")) {
        impute_value = median(data[[i]], na.rm = TRUE)
      }
      else {
        warning(paste(colname, "should be numeric or factor type. But its class is", 
                      col_class))
      }
      if (verbose) {
        cat(" Impute value:", impute_value, "\n")
      }
      impute_values[[colname]] = impute_value
      if (nas == nrow(data)) {
        if (verbose) {
          cat("Note: cannot impute", colname, "because all values are NA.\n")
        }
        next
      }
      else if (nas == 0) {
        next
      }
      else {
        new_data[is.na(data[[i]]), i] = impute_value
      }
    }
    if (!all_vars) {
      impute_values = impute_values[names(any_nas)]
    }
    results$impute_values = impute_values
  }
  else if (type == "knn") {
    if (verbose) {
      cat("Running knn imputation. NOTE: this will standardize your data!\n")
    }
    if (!"RANN" %in% installed.packages()) {
      stop("knn imputation requires the RANN package. Please run install.packages(\"RANN\")")
    }
    impute_info = caret::preProcess(new_data, method = c("knnImpute"))
    new_data = predict(impute_info, new_data)
    results$impute_info = impute_info
  }
  else if (type == "glrm") {
    if (verbose) {
      cat("Running glrm imputation via h2o.\n")
    }
    capture.output({
      h2o::h2o.init(nthreads = -1)
    }, split = verbose)
    capture.output({
      df_h2o = h2o::as.h2o(new_data[, !names(new_data) %in% 
                                      skip_vars])
    }, split = verbose)
    if (is.null(h2o_glrm)) {
      capture.output({
        model_glrm = h2o::h2o.glrm(training_frame = df_h2o, 
                                   k = min(ncol(df_h2o), glrm_k), loss = "Quadratic", 
                                   init = "SVD", svd_method = "GramSVD", regularization_x = "None", 
                                   regularization_y = "None", min_step_size = 1e-06, 
                                   max_iterations = 1000)
      }, split = verbose)
    }
    else {
      model_glrm = h2o_glrm
    }
    capture.output({
      imp_h2o = predict(model_glrm, df_h2o)
    }, split = verbose)
    results$h2o_glrm = model_glrm
    capture.output({
      glrm_data = as.data.frame(imp_h2o)
    }, split = verbose)
    names(glrm_data) = setdiff(names(data), skip_vars)
    for (colname_i in names(data)[any_nas]) {
      missing_val = is.na(new_data[, colname_i])
      new_data[missing_val, colname_i] = glrm_data[missing_val, 
                                                   colname_i]
    }
  }
  if (add_indicators) {
    if (length(any_nas) > 0L) {
      if (verbose) {
        cat("Generating missingness indicators.\n")
      }
      missing_indicators = missingness_indicators(data[, 
                                                       any_nas, drop = FALSE], prefix = prefix, remove_constant = remove_constant, 
                                                  remove_collinear = remove_collinear, verbose = verbose)
      if (verbose) {
        cat(paste0("Indicators added (", ncol(missing_indicators), 
                   "):"), paste(colnames(missing_indicators), 
                                collapse = ", "), "\n")
      }
      results$indicators_added = colnames(missing_indicators)
      new_data = cbind(new_data, missing_indicators)
    }
  }
  results$data = new_data
  results
}

missingness_indicators<-function (data, prefix = "miss_", remove_constant = TRUE, 
          remove_collinear = TRUE, skip_vars = c(), verbose = FALSE) 
{
  any_nas = which(sapply(data[, !colnames(data) %in% skip_vars, 
                              drop = FALSE], function(col) anyNA(col)))
  if (verbose) {
    cat("Generating", length(any_nas), "missingness indicators.\n")
  }
  indicators = 1L * is.na(data[, names(any_nas), drop = FALSE])
  if (length(any_nas) > 0) {
    colnames(indicators) = paste0(prefix, names(any_nas))
  }
  if (remove_constant) {
    col_means = colMeans(indicators)
    if (verbose) {
      num_removed = sum(col_means %in% c(0, 1))
      if (num_removed > 0) {
        cat("Removing", num_removed, "indicators that are constant.\n")
      }
    }
    indicators = indicators[, !col_means %in% c(0, 1), drop = FALSE]
  }
  if (remove_collinear) {
    if (verbose) {
      cat("Checking for collinearity of indicators.\n")
    }
    linear_combos = caret::findLinearCombos(indicators)
    remove_columns = linear_combos$remove
    if (length(linear_combos$remove) > 0L) {
      if (verbose) {
        cat("Removing", length(linear_combos$remove), 
            "indicators due to collinearity:\n")
        cat(paste0(colnames(indicators)[linear_combos$remove], 
                   collapse = ", "), "\n")
      }
      indicators = indicators[, -linear_combos$remove, 
                              drop = FALSE]
    }
  }
  if (verbose) {
    cat("Final number of indicators:", ncol(indicators), 
        "\n")
  }
  return(indicators)
}


#-------------------------------------------
# set up SL components
#-------------------------------------------



#Set up SL components
# lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
# lrnr_glmnet <- make_learner(Lrnr_glmnet)
# lrnr_ranger100 <- make_learner(Lrnr_ranger, num.trees = 100)
# lrnr_hal_simple <- make_learner(Lrnr_hal9001, degrees = 1, n_folds = 2)
# lrnr_gam <- Lrnr_pkg_SuperLearner$new("SL.gam")
# lrnr_bayesglm <- Lrnr_pkg_SuperLearner$new("SL.bayesglm")
# 
# 
# stack <- make_learner(
#   Stack,
#   lrnr_glm, lrnr_mean, lrnr_ranger100, lrnr_glmnet,
#   lrnr_gam, lrnr_bayesglm
# )


lrnr_glmnet <- Lrnr_glmnet$new()
random_forest <- Lrnr_randomForest$new()
glm_fast <- Lrnr_glm_fast$new()
nnls_lrnr <- Lrnr_nnls$new()

xgboost_lrnr <- Lrnr_xgboost$new()
ranger_lrnr <- Lrnr_ranger$new()
#gbm_lrnr <- Lrnr_gbm$new()
earth_lrnr <- Lrnr_earth$new()
#dbarts_lrnr <- Lrnr_dbarts$new()
#hal_lrnr <- Lrnr_hal9001$new()
gam_lrnr <- Lrnr_gam$new()
polyspline<-Lrnr_polspline$new()






screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")
screen_glmnet <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
#screen_corRank <- Lrnr_screener_corRank$new(rank=3)


cor_nnls_lrnr <- make_learner(Pipeline, screen_cor, glm_fast)
cor_glm <- make_learner(Pipeline, screen_cor, glm_fast)
cor_spline <- make_learner(Pipeline, screen_cor, polyspline)
#screened_hal <- make_learner(Pipeline, screen_glmnet, hal_lrnr)



stack <- make_learner(
  Stack,
  lrnr_mean,
  glm_fast,
  cor_glm,
  lrnr_glmnet,
  ranger_lrnr,
  xgboost_lrnr,
  gam_lrnr,
  cor_spline#,
  #screened_hal
)

stack <- make_learner(
  Stack,
  lrnr_mean,
  glm_fast,
  # cor_glm,
   lrnr_glmnet,
  # ranger_lrnr,
  xgboost_lrnr,
  gam_lrnr#,
  # cor_spline#,
  #screened_hal
)


weighted_stack <- make_learner(
  Stack,
  glm_fast,
  cor_glm,
  lrnr_glmnet,
  xgboost_lrnr
)


cor_pipeline <- try(make_learner(Pipeline, screen_cor, stack))


fancy_stack <- make_learner(Stack, cor_pipeline, stack)


metalearner <- make_learner(Lrnr_nnls)



sl <- make_learner(Lrnr_sl,
                   learners = stack,
                   metalearner = metalearner
)

glm_sl <- make_learner(Lrnr_sl,
                   learners = make_learner(Stack, glm_fast),
                   metalearner = metalearner
)


weighted_sl <- make_learner(Lrnr_sl,
                            learners = weighted_stack,
                            metalearner = make_learner(Lrnr_solnp) 
)





#function to drop obs with missing earlier growth measures
clean_prediction_dataset <- function(d, outcome = "haz_12", covars=covars, covars2=NULL){
  if(is.null(covars2)){
    d <- d[!is.na(d[,outcome]),]
    covars <- covars[grepl("haz", covars)|grepl("waz", covars)|grepl("whz", covars)|
                       grepl("stunt", covars)|grepl("wast", covars)|grepl("underwt", covars)]
    d$n_missing_anthro <- rowSums(is.na(d[,covars]))
    d <- d %>% filter(n_missing_anthro < length(covars))
  }else{
    d <- d[!is.na(d[,outcome]),]
    covars <- covars[grepl("haz", covars)|grepl("waz", covars)|grepl("whz", covars)|
                       grepl("stunt", covars)|grepl("wast", covars)|grepl("underwt", covars)]
    covars2 <- covars2[grepl("haz", covars2)|grepl("waz", covars2)|grepl("whz", covars2)|
                       grepl("stunt", covars2)|grepl("wast", covars2)|grepl("underwt", covars2)]
    d$n_missing_anthro <- rowSums(is.na(d[,covars]))
    d$n_missing_anthro2 <- rowSums(is.na(d[,covars2]))
    d <- d %>% filter(n_missing_anthro < length(covars) & n_missing_anthro2 < length(covars2))    
  }
  return(d)
}


## ------------------------------------------------------
## Function to fit SL model
## ------------------------------------------------------

#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
#HAODONG CHECK THIS FUNCTION
#Notes: 1)Check if you can get the SL cross-validation working
# 2) Check the calculations on lines 412-480:
# Am I correctoly calculating MSE, R2, and AUC?
#XXXXXXXXXXXXXXXXXXXXXXXXXXXX


fit_SL_fun <- function(dat,  
                       outcome = "haz",
                       id="subjid",
                       family="gaussian",
                       covars,
                       slmod=sl,
                       CV=TRUE,
                       folds=5){
  
  Y_miss <- sum(is.na(dat[,outcome]))
  if(Y_miss>0){
    dat <- dat[!is.na(dat[,outcome]),]
    cat("Dropping ",Y_miss," missing outcome observations\n")
  }
  
  start_time <- Sys.time()
  full_covars <- covars
  #Drop missing variables or all-NA variables 
  missing_covars <- covars[!(covars %in% colnames(dat))]
  
  
  #Drop near-zero variance predictors
  not_cov <- dat[,!(colnames(dat) %in% covars)]
  cov <- dat[,covars]
  #impute missing
  cov <- cov %>% 
    do(impute_missing_values(., type = "standard", add_indicators = T, prefix = "missing_")$data) %>%
    as.data.frame()
  nzv_cols <- nearZeroVar(cov)
  dropped_covars <-  colnames(cov)[nzv_cols]
  cat("Dropping for low variance: ", dropped_covars)
  if(length(nzv_cols) > 0){ 
    covars <- covars[!(covars %in% colnames(cov)[nzv_cols])]
    cov <- cov[, -nzv_cols]
  }
  not_cov <- data.frame(not_cov)
  cov <- data.frame(cov)
  dat <-cbind(not_cov,cov)
  dat <- data.table(dat)
  
  
  
  n= nrow(dat)
  y <- as.matrix(dat[, outcome, with=FALSE])
  colnames(y) <- NULL
  
  
  
  
  #parallelize
  uname <- system("uname -a", intern = TRUE)
  os <- sub(" .*", "", uname)
  if(os=="Darwin"){
    cpus_logical <- as.numeric(system("sysctl -n hw.logicalcpu", intern = TRUE))
  } else if(os=="Linux"){
    cpus_logical <- system("lscpu | grep '^CPU(s)'", intern = TRUE)
    cpus_logical <- as.numeric(gsub("^.*:[[:blank:]]*","", cpus_logical))
  } else {
    stop("unsupported OS")
  }

  plan(multicore, workers=floor(cpus_logical/4))
  
  # create the sl3 task
  set.seed(12345)
  cv_folds = make_folds(dat, fold_fun = folds_vfold, V = folds)
  SL_task <- make_sl3_Task(
    data = dat,
    covariates = covars,
    outcome = outcome,
    folds = cv_folds
  )
  
  
  
  
  ##3. Fit the full model
  
  
  #get cross validated fit
  # if(CV){
  #   cv_sl <- make_learner(Lrnr_cv, slmod, full_fit = TRUE)
  #   cv_glm <- make_learner(Lrnr_cv, glm_sl, full_fit = TRUE)
  #   sl_fit <- cv_sl$train(SL_task)
  #   glm_fit <- cv_glm$train(SL_task)
  # }else{
  #   sl_fit <- slmod$train(SL_task)
  #   glm_fit <- glm_sl$train(SL_task)
  # }
  
  #0903 test
  sl_fit <- slmod$train(SL_task)
  glm_fit <- glm_sl$train(SL_task)
  
  if (CV){
    
    # start_time <- Sys.time()
    # CVsl <- CV_lrnr_sl(sl_fit, SL_task, loss_squared_error)
    # end_time <- Sys.time()
    # print(end_time - start_time)
    
    # start_time <- Sys.time()
    cv_sl <- make_learner(Lrnr_cv, sl_fit, full_fit = TRUE)
    fit_cv_sl <- cv_sl$train(SL_task)
    yhat_full <- fit_cv_sl$predict_fold(SL_task,"validation")
    # end_time <- Sys.time()
    # print(end_time - start_time)
    
    # start_time <- Sys.time()
    cv_glm <- make_learner(Lrnr_cv, glm_fit, full_fit = TRUE)
    fit_cv_glm <- cv_glm$train(SL_task)
    yhat_glm <- fit_cv_glm$predict_fold(SL_task,"validation")
    # end_time <- Sys.time()
    # print(end_time - start_time)
    
    lrnr_cv_null <- make_learner(Lrnr_cv, make_learner(Lrnr_mean))
    fit_null <- lrnr_cv_null$train(SL_task)
    yhat_null <- fit_null$predict_fold(SL_task,"validation")
    
  }else{
    yhat_full <- sl_fit$predict()
    yhat_glm <- glm_fit$predict()
    yhat_null <- mean(y)
  }
  
  #save residuals
  SL_residuals <- y-yhat_full
  glm_residuals <- y-yhat_glm
  
  
  
  
  if(family=="gaussian"){
    
    mse_full <- 1/n * sum((yhat_full-y)^2)
    mse_glm <- 1/n * sum((yhat_glm-y)^2)
    mse_null <- 1/n * sum((yhat_null-y)^2)
    
    
    
    ##5. Construct CI for mse_full/mse_null
    IC <- 1/mse_null * ((yhat_full-y)^2 - mse_full) - mse_full/(mse_null)^2 * ((yhat_null-y)^2- mse_null)
    IC_glm <- 1/mse_null * ((yhat_glm-y)^2 - mse_glm) - mse_glm/(mse_null)^2 * ((yhat_null-y)^2- mse_null)
    
    psi <- mse_full/mse_null
    se <- sqrt(var(IC)/n)
    CI_up <- psi + 1.96*se
    CI_low <- psi - 1.96*se
    R2 <- 1 - psi
    R2.ci1 <- 1 - CI_up
    R2.ci2 <- 1 - CI_low
    
    glm.psi <- mse_glm/mse_null
    glm.se <- sqrt(var(IC_glm)/n)
    glm.CI_up <- glm.psi + 1.96*glm.se
    glm.CI_low <- glm.psi - 1.96*glm.se
    glm.R2 <- 1 - glm.psi
    glm.R2.ci1 <- 1 - glm.CI_up
    glm.R2.ci2 <- 1 - glm.CI_low
    
    
    R2_full <- data.frame(R2=R2, R2.se=se, R2.ci1=R2.ci1, R2.ci2=R2.ci2)
    R2_glm <- data.frame(R2=glm.R2, R2.se=glm.se, R2.ci1=glm.R2.ci1, R2.ci2=glm.R2.ci2)
    
    ##6. add more results to pool R2
    
    # ic_mse_full <- (yhat_full-y)^2 - mse_full
    # ic_mse_null <- (yhat_null-y)^2 - mse_null
    # 
    # se_mse_full <- sqrt(var(ic_mse_full)/n)
    # se_mse_null <- sqrt(var(ic_mse_null)/n)
    # 
    # mse_ic <- list(ic_mse_full=ic_mse_full,
    #                ic_mse_null=ic_mse_null,
    #                se_mse_full=se_mse_full,
    #                se_mse_null=se_mse_null)
    
    perf_metrics <- list(mse_full = mse_full,
                         mse_glm = mse_glm,
                         mse_null = mse_null,
                         R2_full = R2_full,
                         R2_glm = R2_glm,
                         mse_ic=mse_ic,
                         mse = mse_full/mse_null,
                         mse_var = se^2
    )
    
    fold_index <- origami::folds2foldvec(SL_task$folds) 
    
  }else{
    fold_index <- origami::folds2foldvec(SL_task$folds) 
    
    outcome <- SL_task$Y
    
    auc.plotdf <- cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index)
    auc.ci<-ci.cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index, confidence=0.95)
    glm.auc.ci<-ci.cvAUC(predictions=yhat_glm, labels=outcome, folds=fold_index, confidence=0.95)
    
    perf_metrics <- list(auc.plotdf = auc.plotdf,
                         auc.ci = auc.ci,
                         glm.auc.ci = glm.auc.ci,
                         fold_index = fold_index)
  }
  
  
  end_time <- Sys.time()
  
  run_time <- end_time - start_time
  
  
  return(
    list(
      sl_fit = sl_fit,
      covars = list(specified_covars=full_covars,
                    missing_covars=missing_covars, 
                    dropped_covars=dropped_covars,
                    used_covars=covars),
      N_obs=n,
      yhat_full = yhat_full,
      yhat_glm = yhat_glm,
      Y= SL_task$Y,
      fold_index=fold_index,
      perf_metrics=perf_metrics,
      SL_residuals = SL_residuals,
      glm_residuals = glm_residuals,
      run_time=run_time,
      Y_miss=Y_miss,
      id=dat$subjid
    )
  )
  #}
}


#wrapper function to detect errors
fit_SL <- purrr::safely(fit_SL_fun)

clean_res <- function(res){

  clean_res <- res %>% 
    mutate(no_error = fit %>% purrr::map_lgl(.f = ~ is.null(.x$error)))
  return(clean_res)

}



## ------------------------------------------------------
## Function to fit SL model
## ------------------------------------------------------


fit_weighted_SL_fun <- function(dat,  
                       outcome = "haz",
                       id="subjid",
                       covars,
                       slmod=weighted_sl,
                       CV=TRUE,
                       folds=5){
  
  Y_miss <- sum(is.na(dat[,outcome]))
  if(Y_miss>0){
    dat <- dat[!is.na(dat[,outcome]),]
    cat("Dropping ",Y_miss," missing outcome observations\n")
  }
  
  start_time <- Sys.time()
  full_covars <- covars
  #Drop missing variables or all-NA variables 
  missing_covars <- covars[!(covars %in% colnames(dat))]
  
  
  #Drop near-zero variance predictors
  not_cov <- dat[,!(colnames(dat) %in% covars)]
  cov <- dat[,covars]
  #impute missing
  cov <- cov %>% 
    do(impute_missing_values(., type = "standard", add_indicators = T, prefix = "missing_")$data) %>%
    as.data.frame()
  nzv_cols <- nearZeroVar(cov)
  dropped_covars <-  colnames(cov)[nzv_cols]
  cat("Dropping for low variance: ", dropped_covars)
  if(length(nzv_cols) > 0){ 
    covars <- covars[!(covars %in% colnames(cov)[nzv_cols])]
    cov <- cov[, -nzv_cols]
  }
  not_cov <- data.frame(not_cov)
  cov <- data.frame(cov)
  dat <-cbind(not_cov,cov)
  dat <- data.table(dat)
  
  
  #Truncate weights
  dat$trunc_weights <- dat$weights
  dat$trunc_weights[dat$weights<0.1] <- 0.1 
  dat$trunc_weights[dat$weights>0.9] <- 0.9 
  summary(dat$trunc_weights)
  dat$inv_trunc_weights <- 1/(1-dat$trunc_weights)
  
  
  n= nrow(dat)
  y <- as.matrix(dat[, outcome, with=FALSE])
  colnames(y) <- NULL
 
  
  #parallelize
  uname <- system("uname -a", intern = TRUE)
  os <- sub(" .*", "", uname)
  if(os=="Darwin"){
    cpus_logical <- as.numeric(system("sysctl -n hw.logicalcpu", intern = TRUE))
  } else if(os=="Linux"){
    cpus_logical <- system("lscpu | grep '^CPU(s)'", intern = TRUE)
    cpus_logical <- as.numeric(gsub("^.*:[[:blank:]]*","", cpus_logical))
  } else {
    stop("unsupported OS")
  }
  
  plan(multicore, workers=floor(cpus_logical/4))
  
  # create the sl3 task
  set.seed(12345)
  SL_task <- make_sl3_Task(
    data = dat,
    covariates = covars,
    outcome = outcome,
    weights = "inv_trunc_weights",
    folds = make_folds(dat, fold_fun = folds_vfold, V = folds)
  )
  
  
  
  
  ##3. Fit the full model
  
  
  #get cross validated fit
  if(CV){
    cv_sl <- make_learner(Lrnr_cv, slmod, full_fit = TRUE)
    sl_fit <- cv_sl$train(SL_task)
  }else{
    sl_fit <- slmod$train(SL_task)
  }
  
  #get outcome predictions
  yhat_full <- sl_fit$predict_fold(SL_task,"validation")
  
  #save SL residuals
  SL_residuals <- y-yhat_full
  
  
  ##4. Fit the null model
  lrnr_cv_null <- make_learner(Lrnr_cv, make_learner(Lrnr_mean))
  fit_null <- lrnr_cv_null$train(SL_task)
  yhat_null <- fit_null$predict_fold(SL_task,"validation")
  
    
    mse_full <- 1/n * sum((yhat_full-y)^2)
    mse_null <- 1/n * sum((yhat_null-y)^2)
    
    
    
    ##5. Construct CI for mse_full/mse_null
    IC <- 1/mse_null * ((yhat_full-y)^2 - mse_full) - mse_full/(mse_null)^2 * ((yhat_null-y)^2- mse_null)
    
    
    psi <- mse_full/mse_null
    se <- sqrt(var(IC)/n)
    CI_up <- psi + 1.96*se
    CI_low <- psi - 1.96*se
    
    R2  <- 1 - psi
    R2.ci1 <- 1 - CI_up
    R2.ci2 <- 1 - CI_low
    
    
    R2_full <- data.frame(R2=R2, R2.se=se, R2.ci1=R2.ci1, R2.ci2=R2.ci2)
    
    ##6. add more results to pool R2
    
    ic_mse_full <- (yhat_full-y)^2 - mse_full
    ic_mse_null <- (yhat_null-y)^2 - mse_null
    
    se_mse_full <- sqrt(var(ic_mse_full)/n)
    se_mse_null <- sqrt(var(ic_mse_null)/n)
    
    mse_ic <- list(ic_mse_full=ic_mse_full,
                   ic_mse_null=ic_mse_null,
                   se_mse_full=se_mse_full,
                   se_mse_null=se_mse_null)
    
    perf_metrics <- list(mse_full = mse_full,
                         mse_null = mse_null,
                         R2_full = R2_full,
                         mse_ic=mse_ic,
                         SL_residuals = SL_residuals)
    
    fold_index <- origami::folds2foldvec(SL_task$folds) 
    
 
  
  
  end_time <- Sys.time()
  
  run_time <- end_time - start_time
  

    return(
      list(
        sl_fit = sl_fit,
        covars = list(specified_covars=full_covars,
                      missing_covars=missing_covars, 
                      dropped_covars=dropped_covars,
                      used_covars=covars),
        N_obs=n,
        yhat_full = yhat_full,
        Y= SL_task$Y,
        fold_index=fold_index,
        perf_metrics=perf_metrics,
        SL_residuals = SL_residuals,
        run_time=run_time,
        Y_miss=Y_miss,
        id==dat$subjid
      )
    )
}


#wrapper function to detect errors
fit_weighted_SL <- purrr::safely(fit_weighted_SL_fun)






fit.rma <- function(data, ni, xi = NULL, yi = NULL, vi = NULL, measure = "PLO", nlab = "", method = "REML", age=NULL) {
  if(!is.null(age)){
    data <- data[data$agecat==age,]
    data$age <- age
  }
  
  mode_continuous <- !is.null(yi) | !is.null(vi)
  mode_binary <- !is.null(xi)
  if (mode_binary & mode_continuous) stop("can only do binary or continuous")
  # check if measure=="PLO" - default for all parameters bounded between 0 and 1 (prevalence, cumulative incidence)
  # because then output needs back transformation
  if (measure == "PLO" & mode_binary) {
    # If only one row of the data, no need for pooling (single study, often seens when looping over agecats),
    # so just use the escalc() function that underlies RMA to calculate the parameter for the single study
    if (nrow(data) == 1) {
      fit <- NULL
      try(fit <- escalc(data = data, ni = data[[ni]], xi = data[[xi]], method = method, measure = measure, append = T))
      data <- fit
      data$se <- sqrt(data$vi)
      out <- data %>%
        ungroup() %>%
        mutate(nstudies = 1, nmeas = data[[ni]]) %>%
        mutate(
          est = plogis(yi),
          lb = plogis(yi - 1.96 * se),
          ub = plogis(yi + 1.96 * se),
          nmeas.f = paste0("N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE), " ", nlab),
          nstudy.f = paste0("N=", nstudies, " studies"),
          method.used=method,
          ptest.f = sprintf("%0.0f", est)
        ) %>%
        subset(., select =c(nstudies, nmeas, est, se, lb, ub, nmeas.f, nstudy.f, method.used, agecat, ptest.f)) %>%
        as.tibble()
      rownames(out) <- NULL
      # If input is more than 1 row (multiple studies), pool across studies with rma() function from metafor package
    } else {
      # Check if 0 cases of the outcome
      # Use FE model if all 0 counts because no heterogeneity and rma.glmm fails
      if(!is.null(xi)){
        if (sum(data[[xi]]) == 0) method <- "FE"
      }
      fit <- NULL
      method_fit <- method
      try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = method,
        xi = data[[xi]],
        measure = measure
      ))
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = "ML",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "ML"
      }
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ni = data[[ni]],
        method = "DL",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "DL"
      }
      # Create formatted dataset to return
      out <- data %>%
        ungroup() %>%
        summarise(
          nstudies = length(unique(paste0(studyid," ",country))),
          nmeas = sum(data[[ni]])
        ) %>%
        mutate(
          est = plogis(fit$beta),
          se = fit$se,
          lb = plogis(fit$beta - 1.96 * fit$se),
          ub = plogis(fit$beta + 1.96 * fit$se),
          nmeas.f = paste0("N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE), " ", nlab),
          nstudy.f = paste0("N=", nstudies, " studies"),
          method.used=method_fit
        ) %>%
        as.tibble()
      rownames(out) <- NULL
    }
  } else {
    
    if(measure == "IR"){
      # If measure if IR for incidence rate, use `to="if0all"` argument to add .5 to all cells of the 2x2 table if one is 0 so rma() works
      to <- "if0all"
      fit <- NULL
      method_fit <- method
      try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = method,
        xi = data[[xi]],
        measure = measure,
        to=to
      ))
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = "ML",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "ML"
      }
      if(is.null(fit)){try(fit <- rma(
        data = data,
        ti = data[[ni]],
        method = "DL",
        xi = data[[xi]],
        measure = measure
      ))
        method_fit <- "DL"
      }
      
    }else{
      
      
      # If measure other than PLO or IR is chosen:
      to <- "only0"
      
      fit <- NULL
      method_fit <- method
      
      if (mode_binary) {
        try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = method,
          xi = data[[xi]],
          measure = measure,
          to=to
        ))
        if(is.null(fit)){try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = "ML",
          xi = data[[xi]],
          measure = measure
        ))
          method_fit <- "ML"
        }
        if(is.null(fit)){try(fit <- rma(
          data = data,
          ni = data[[ni]],
          method = "DL",
          xi = data[[xi]],
          measure = measure
        ))
          method_fit <- "DL"
        }
      }
      if (mode_continuous) {
        try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = method,
          measure = "MN"
        ))
        if(is.null(fit)){try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = "ML",
          measure = "MN"))
          method_fit <- "ML"
        }
        if(is.null(fit)){try(fit <- rma(
          data = data,
          mi = data[[yi]],
          sdi = sqrt(data[[vi]]),
          ni = data[[ni]],
          method = "DL",
          measure = "MN"))
          method_fit <- "DL"
        }
      }
    }
    out <- data %>%
      ungroup() %>%
      summarise(
        nstudies = length(unique(studyid)),
        nmeas = sum(data[[ni]])
      ) %>%
      mutate(
        est = fit$beta,
        se = fit$se,
        lb = fit$ci.lb,
        ub = fit$ci.ub,
        nmeas.f = paste0(
          "N=", format(sum(data[[ni]]), big.mark = ",", scientific = FALSE),
          " ", nlab
        ),
        nstudy.f = paste0("N=", nstudies, " studies"),
        method.used=method_fit
      )
    
  }
  return(out)
}





#Function to identify high and low risk children and estimate treatment effect within subgroup
risk_quartile_tmle <- function(d, model, study, Y, contrast, Avar="intervention"){
  d <- d[d[[Avar]] %in% contrast,]
  d <- d %>% filter(studyid==study)
  dim(d)
  
  #extract model
  model <- model %>% filter(studyid==study) %>% select(fit)
  model <- model[[1]]
  model <- model[[1]]$result$sl_fit
  
  d <- d[!is.na(d[,Y]),]
  task <- make_sl3_Task(data = d, covariates = covars12,
                        outcome = Y, outcome_type="continuous",
                        drop_missing_outcome=TRUE)
  d$preds<-model$predict(task)
  summary(d$preds)
  
  #categorize risk
  d$pred_quart <- ntile(d$preds, 3)
  d <- d %>% mutate(pred_quart = case_when(pred_quart==1 ~ "High",
                                           pred_quart==2 ~ "Med.",
                                           pred_quart==3 ~ "Low"))
  d$pred_quart <- factor(d$pred_quart, levels = c("High","Med.","Low"))
  d <- d %>% subset(., select = c(Y,"pred_quart",Avar,covars))
  
  #specify TMLE
  node_list <- list(
    W=c(covars,"pred_quart"),
    A=Avar,
    Y=Y
  )
  processed <- process_missing(data=d, node_list,  max_p_missing = 0.5)
  df_processed <- processed$data
  df_processed <- droplevels(df_processed)
  node_list <- processed$node_list
  
  ate_spec <- tmle_ATE(
    treatment_level = contrast[2],
    control_level = contrast[1]
  )
  learner_list <- list(A = tmle_mod, Y = tmle_mod)
  
  
  
  # tmle_fit1 <- tmle3(ate_spec, df_processed[df_processed$pred_quart==1,], node_list, learner_list)
  # tmle_fit2 <- tmle3(ate_spec, df_processed[df_processed$pred_quart==2,], node_list, learner_list)
  # tmle_fit3 <- tmle3(ate_spec, df_processed[df_processed$pred_quart==3,], node_list, learner_list)
  # tmle_fit4 <- tmle3(ate_spec, df_processed[df_processed$pred_quart==4,], node_list, learner_list)
  # 
  tmle_fit1 <- tmle3(ate_spec, df_processed[df_processed$pred_quart=="Low",], node_list, learner_list)
  tmle_fit2 <- tmle3(ate_spec, df_processed[df_processed$pred_quart=="Med.",], node_list, learner_list)
  tmle_fit3 <- tmle3(ate_spec, df_processed[df_processed$pred_quart=="High",], node_list, learner_list)
  
  res <- bind_rows(tmle_fit1$summary, tmle_fit2$summary, tmle_fit3$summary) #, tmle_fit4$summary)
  res$pred_quart <- c("Low", "Med.", "High")
  res$studyid <- study
  res$intervention <- contrast[2]
  res$Y <- Y
  
  return(res)
}
