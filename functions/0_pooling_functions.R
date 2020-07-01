
#-------------------------------------------
# List extraction function
#-------------------------------------------

extract_results = function(x) x[["result"]]
extract_perf_metrics = function(x) x[["perf_metrics"]]
extract_R2_full = function(x) x[["R2_full"]]
extract_R2_glm = function(x) x[["R2_glm"]]
extract_N_obs = function(x) x[["N_obs"]]

extract_mse_full = function(x) x[["mse_full"]]
extract_mse_null = function(x) x[["mse_null"]]
extract_mse_ic = function(x) x[["mse_ic"]]
extract_se_mse_full = function(x) x[["se_mse_full"]]
extract_se_mse_null = function(x) x[["se_mse_null"]]

extract_yhat_full = function(x) x[["yhat_full"]]
extract_Y = function(x) x[["Y"]]
extract_fold_index = function(x) x[["fold_index"]]
extract_id = function(x) x[["id"]]

extract_AUC_full = function(x) x[["auc.ci"]]
extract_cvAUC = function(x) x[["cvAUC"]]
extract_ci.AUC = function(x) x[["ci"]]


#-------------------------------------------
# Extract AUC
#-------------------------------------------
extract_AUC <- function(res){
  results <- Map(extract_results, res$fit)
  N_obs <- Map(extract_N_obs, results)
  perf_metrics <- Map(extract_perf_metrics, results)
  AUC_full <- Map(extract_AUC_full, perf_metrics)
  cvAUC <- Map(extract_cvAUC, AUC_full)
  ci.AUC <- Map(extract_ci.AUC, AUC_full)
  names(ci.AUC) <- letters[1:length(ci.AUC)]
  ci.AUC <- t(bind_rows(ci.AUC))
  AUC_full <- data.frame(studyid=res$studyid, N_obs=unlist(N_obs, use.names=FALSE), cvAUC=unlist(cvAUC), ci.lb=ci.AUC[,1],  ci.ub=ci.AUC[,2])
  rownames(AUC_full) <- NULL
  return(AUC_full)
}


extract_preds <- function(res){
  res <- res %>% filter(no_error)
  results <- Map(extract_results, res$fit)
  yhat_full <- Map(extract_yhat_full, results)
  N_obs <- Map(extract_N_obs, results)
  Y <- Map(extract_Y, results)
  fold_index <- Map(extract_fold_index, results)
  
  studyid=rep(res$studyid, lengths(Y))
  country=rep(res$country, lengths(Y))
  N_obs= rep(unlist(N_obs, use.names=FALSE), lengths(Y))
  df <- data.frame(studyid=studyid,
                   country=country,
                   N_obs=N_obs,
                   Y=unlist(Y), 
                   yhat_full=unlist(yhat_full), 
                   fold_index=unlist(fold_index))
  df$Ybin <- 1*(df$Y < -2)
  df$prob <- pnorm(-(df$yhat_full - (-2)))
  return(df)
}

#-------------------------------------------
# Pool AUC
#-------------------------------------------
calc.ci.cvAUC <- function(d){
  cat(as.character(d$studyid)[1], " ", as.character(d$country)[1],"\n")
  res <- NULL
  if(length(unique(d$Ybin))==2){
    suppressMessages({try(res <- ci.cvAUC(predictions=d$prob, labels=d$Ybin, folds=d$fold_index, confidence=0.95))})
    if(is.null(res)){
      try(res <- ci.cvAUC(predictions=d$prob, labels=d$Ybin,confidence=0.95))
    }
  resdf <- data.frame(N_obs=d$N_obs[1], AUC=res$cvAUC, se=res$se, AUC.ci1=res$ci[1], AUC.ci2=res$ci[2])
  }else{
    resdf <- data.frame(N_obs=d$N_obs[1], AUC=NA, se=NA, AUC.ci1=NA, AUC.ci2=NA)
  }
  return(resdf)
}

calc_predY = function(d){
 #cat(d$studyid[1], " ",d$country[1],"\n")
 #  analysis <- roc(response=d$Ybin, predictor=d$prob)
 #  
 #  
 #  #Find t that minimizes error
 #  e <- cbind(analysis$thresholds,analysis$sensitivities+analysis$specificities)
 #  e <- matrix(e[is.finite(rowSums(e)),])
 #  if(ncol(e)==1) e<-t(e)
 #  
 #  opt_t <- subset(e,e[,2]==max(e[,2]))[,1]
 
  if(length(unique(d$Ybin))==2){
  #From https://ethen8181.github.io/machine-learning/unbalanced/unbalanced.html
  # user-defined different cost for false negative and false positive
  cost_fp <- 100
  cost_fn <- 200
  roc_info <- ROCInfo(data = d, predict = "prob", 
                       actual = "Ybin", cost.fp = cost_fp, cost.fn = cost_fn )
  
  #Predict outcome classes
  # d$pred_Y <- ifelse(d$prob >= max(opt_t), 1, 0)
  # d$opt_t <- max(opt_t)
  
  d$pred_Y <- ifelse(d$prob >= roc_info$cutoff, 1, 0)
  d$opt_t <- roc_info$cutoff
  }
  return(d)
}


calc_tp = function(d){
  #Get TP, FP, TN, FP
  tp = sum(d$pred_Y == 1 & d$Ybin == 1)
  fp = sum(d$pred_Y == 1 & d$Ybin == 0)
  tn = sum(d$pred_Y == 0 & d$Ybin == 0)
  fn = sum(d$pred_Y == 0 & d$Ybin == 1)
  
  return(data.frame(TP=tp, FP=fp, TN=tn, FN=fn))
}


get_pooled_AUC <- function(d, study_AUC){
  
  fit.phm.homo <- phm(d, hetero = FALSE)
  sum.phm.homo <- summary(fit.phm.homo)
  
  #details on AUC vs. pAUC
  #https://rdrr.io/rforge/mada/man/phm-class.html
  #pAUC is partial AUC
  #https://www.ncbi.nlm.nih.gov/pubmed/15900606
  res.phm.homo <- data.frame(
    studyid="Pooled-FE",
    country=NA,
    AUC=sum.phm.homo$AUC$AUC,
    AUC.ci1=sum.phm.homo$AUC$ci[2],
    AUC.ci2=sum.phm.homo$AUC$ci[1])
  
  fit.phm.hetero <- phm(d, hetero = TRUE)
  sum.phm.hetero <- summary(fit.phm.hetero)
  res.phm.hetero <- data.frame(
    studyid="Pooled-RE",
    country=NA,
    AUC=sum.phm.hetero$AUC$AUC,
    AUC.ci1=sum.phm.hetero$AUC$ci[2],
    AUC.ci2=sum.phm.hetero$AUC$ci[1])
  
  res <- bind_rows(
    data.frame(study_AUC, pooled=0),
    data.frame(res.phm.homo, pooled=1), 
    data.frame(res.phm.hetero, pooled=1))
  
  res$pooled <- factor(res$pooled, levels = c("0","1"))
  return(res)
}


pool_AUC <- function(res, Y, age, covars){
  objectname=deparse(substitute(res))
  
  res <- res %>% filter(no_error)
  preds <- extract_preds(res)
  preds <- preds %>% filter(N_obs>=50)
  study_AUC <- preds %>% group_by(studyid, country) %>%
    do(calc.ci.cvAUC(.))
  
  preds <- preds %>% group_by(studyid, country) %>%
    do(calc_predY(.))
  preds <- preds %>% filter(!is.na(pred_Y))
  study_AUC <- study_AUC%>%filter(!is.na(AUC))
  
  study_confusion_matrix <- preds %>% group_by(studyid, country) %>%
    do(calc_tp(.))
  
  study_confusion_matrix <- study_confusion_matrix %>% filter(!is.na(FP))
  AUC_res <- get_pooled_AUC(study_confusion_matrix, study_AUC)
  AUC_res$Y = Y
  AUC_res$age = age
  AUC_res$covars = covars
  
  saveRDS(preds, here(paste0("prediction-models/binary_preds/",objectname,"_preds.RDS")))
  return(AUC_res)
}


pool_bin_AUC <- function(res, Y, age, covars, save=T){
  res <- res %>% filter(no_error)
  preds <- extract_preds(res)
  preds$prob <- preds$yhat_full
  preds$Ybin <- preds$Y

  study_AUC <- preds %>% group_by(studyid, country) %>%
    do(calc.ci.cvAUC(.))
  
  preds <- preds %>% group_by(studyid, country) %>%
    do(calc_predY(.))
  study_confusion_matrix <- preds %>% group_by(studyid, country) %>%
    do(calc_tp(.))
  
  AUC_res <- get_pooled_AUC(study_confusion_matrix, study_AUC)
  AUC_res$Y = Y
  AUC_res$age = age
  AUC_res$covars = covars
  if(save){
    saveRDS(preds, here(paste0("prediction-models/binary_preds/",deparse(substitute(res)),"_preds.RDS")))
  }
  return(AUC_res)
}


#-------------------------------------------
# Extract MSE
#-------------------------------------------
extract_mse <- function(res){
  res <- res %>% filter(no_error)
  results <- Map(extract_results, res$fit)
  N_obs <- Map(extract_N_obs, results)
  perf_metrics <- Map(extract_perf_metrics, results)
  mse_ic <- Map(extract_mse_ic, perf_metrics)
  mse_full <- Map(extract_mse_full, perf_metrics)
  mse_null <- Map(extract_mse_null, perf_metrics)
  se_mse_full <- Map(extract_se_mse_full, mse_ic)
  se_mse_null <- Map(extract_se_mse_null, mse_ic)
  mse_full <- data.frame(studyid=res$studyid, 
                         country=res$country, 
                         N_obs=unlist(N_obs, use.names=FALSE),
                         mse_full=unlist(mse_full),
                         mse_null=unlist(mse_null),
                         se_mse_full=unlist(se_mse_full),
                         se_mse_null=unlist(se_mse_null))
  return(mse_full)
}


#-------------------------------------------
# Pool R2
#-------------------------------------------

extract_R2 <- function(res){
  res <- res %>% filter(no_error)
  results <- Map(extract_results, res$fit)
  N_obs <- Map(extract_N_obs, results)
  perf_metrics <- Map(extract_perf_metrics, results)
  R2_full <- Map(extract_R2_full, perf_metrics)
  R2_full <- data.frame(studyid=res$studyid, country=res$country, N_obs=unlist(N_obs, use.names=FALSE), bind_rows(R2_full))
  return(R2_full)
}

pool_R2 <- function(res, Y, age, covars){
  
  R2 <- extract_R2(res)
  R2 <- R2 %>% filter(N_obs>=50)
  
  #(temp) calculate SE
  R2$se <- ((1-R2$R2.ci1) - (1-R2$R2))/1.96
  
  
  mse <- extract_mse(res) %>% filter(N_obs>=50)
  #calc variance
  mse <- mse %>%
    mutate(var_mse_full=se_mse_full^2,var_mse_null=se_mse_null^2, N=n()) %>%  
    dplyr::select (-c(se_mse_full, se_mse_null))
  
  pooled_mse_full <-fit.rma(data=mse, ni="N",  yi="mse_full", vi= "var_mse_full",
                                age=NULL, method = "FE", measure="MN")
  pooled_mse_null <-fit.rma(data=mse, ni="N", yi="mse_null", vi= "var_mse_null",
                                age=NULL, method = "FE", measure="MN")
  pooled_mse <- as.matrix(pooled_mse_full$est)/as.matrix(pooled_mse_null$est)
  pooled_lb <- as.matrix(pooled_mse_full$lb)/as.matrix(pooled_mse_null$lb)
  pooled_ub <- as.matrix(pooled_mse_full$ub)/as.matrix(pooled_mse_null$ub)
  pooled_r2  <- 1-pooled_mse 
  pooled_r2_ub  <- 1-pooled_lb 
  pooled_r2_lb  <- 1-pooled_ub 
  cat("Pooled R2: ", pooled_r2, "\n")
  
  #delta method
  #https://migariane.github.io/DeltaMethodEpiTutorial.nb.html
  
  pooled_mse_full <-fit.rma(data=mse, ni="N",  yi="mse_full", vi= "var_mse_full",
                            age=NULL, method = "REML", measure="MN")
  pooled_mse_null <-fit.rma(data=mse, ni="N", yi="mse_null", vi= "var_mse_null",
                            age=NULL, method = "REML", measure="MN")
  pooled_mse <- as.matrix(pooled_mse_full$est)/as.matrix(pooled_mse_null$est)
  pooled_lb <- as.matrix(pooled_mse_full$lb)/as.matrix(pooled_mse_null$lb)
  pooled_ub <- as.matrix(pooled_mse_full$ub)/as.matrix(pooled_mse_null$ub)
  pooled_r2_re  <- 1-pooled_mse 
  pooled_r2_ub_re  <- 1-pooled_lb 
  pooled_r2_lb_re  <- 1-pooled_ub 
  cat("Pooled R2: ", pooled_r2_re, "\n")
  
  
  # # Calculate pooled R2 + CI
  # R2$N<-1
  # R2$var<-R2$se^2
  # pooled_r2_fe <- fit.rma(data=R2, ni="N", yi="R2", vi= "var",
  #                      age=NULL, method = "FE", measure="MN")
  # pooled_R2_re <- fit.rma(data=R2, ni="N", yi="R2", vi= "var",
  #                         age=NULL, method = "REML", measure="MN")
  # 
  # pooled_r2_fe <- data.frame(studyid="Pooled - FE", R2=pooled_r2_fe$est, R2.ci1=pooled_r2_fe$lb, R2.ci2=pooled_r2_fe$ub)
  # pooled_r2_re <- data.frame(studyid="Pooled - RE", R2=pooled_R2_re$est, R2.ci1=pooled_R2_re$lb, R2.ci2=pooled_R2_re$ub)
  
  
  pooled_r2_fe <- data.frame(studyid="Pooled - FE", R2=pooled_r2, R2.ci1=pooled_r2_lb, R2.ci2=pooled_r2_ub)
  pooled_r2_re <- data.frame(studyid="Pooled - RE", R2=pooled_r2_re, R2.ci1=pooled_r2_lb_re, R2.ci2=pooled_r2_ub_re)
  
  
  
  R2.df <- bind_rows(R2, pooled_r2_fe, pooled_r2_re)
  R2.df$Y <- Y
  R2.df$age <- age
  R2.df$covars <- covars
  
  return(R2.df)
}


#-------------------------------------------
# Pool GLM R2
#-------------------------------------------

extract_glm_R2 <- function(res){
  res <- res %>% filter(no_error)
  results <- Map(extract_results, res$fit)
  N_obs <- Map(extract_N_obs, results)
  perf_metrics <- Map(extract_perf_metrics, results)
  R2_glm <- Map(extract_R2_glm, perf_metrics)
  R2_glm <- data.frame(studyid=res$studyid, country=res$country, N_obs=unlist(N_obs, use.names=FALSE), bind_rows(R2_glm))
  return(R2_glm)
}

pool_glm_R2 <- function(res, Y, age, covars){
  
  R2 <- extract_glm_R2(res)
  R2 <- R2 %>% filter(N_obs>=50)
  
  #(temp) calculate SE
  R2$se <- ((1-R2$R2.ci1) - (1-R2$R2))/1.96
  
  
  mse <- extract_mse(res)
  #calc variance
  mse <- mse %>%
    mutate(var_mse_full=se_mse_full^2,var_mse_null=se_mse_null^2, N=n()) %>%  
    dplyr::select (-c(se_mse_full, se_mse_null))
  
  pooled_mse_full <-fit.rma(data=mse, ni="N",  yi="mse_full", vi= "var_mse_full",
                            age=NULL, method = "FE", measure="MN")
  pooled_mse_null <-fit.rma(data=mse, ni="N", yi="mse_null", vi= "var_mse_null",
                            age=NULL, method = "FE", measure="MN")
  pooled_mse <- as.matrix(pooled_mse_full$est)/as.matrix(pooled_mse_null$est)
  pooled_r2  <- 1-pooled_mse 
  cat("Pooled R2: ", pooled_r2, "\n")
  
  # Calculate pooled R2 + CI
  R2$N<-1
  R2$var<-R2$se^2
  pooled_r2 <- fit.rma(data=R2, ni="N", yi="R2", vi= "var",
                       age=NULL, method = "FE", measure="MN")
  pooled_R2_re <- fit.rma(data=R2, ni="N", yi="R2", vi= "var",
                          age=NULL, method = "REML", measure="MN")
  
  pooled_r2 <- data.frame(studyid="Pooled - FE", R2=pooled_r2$est, R2.ci1=pooled_r2$lb, R2.ci2=pooled_r2$ub)
  pooled_r2_re <- data.frame(studyid="Pooled - RE", R2=pooled_R2_re$est, R2.ci1=pooled_R2_re$lb, R2.ci2=pooled_R2_re$ub)
  
  R2.df <- bind_rows(R2, pooled_r2, pooled_r2_re)
  R2.df$Y <- Y
  R2.df$age <- age
  R2.df$covars <- covars
  
  return(R2.df)
}


#-------------------------------------------
# Extract binary weights
#-------------------------------------------
extract_weights <- function(res){
  res <- res %>% filter(no_error)
  results <- Map(extract_results, res$fit)
  yhat <- Map(extract_yhat_full, results)
  id <- Map(extract_id, results)
  
  studyid=rep(res$studyid, lengths(yhat))
  country=rep(res$country, lengths(yhat))
  
  df <- data.frame(studyid=studyid,
                   country=country,
                   subjid=unlist(id),
                   weights=unlist(yhat))
  return(df)
}









# ------------------------------------------------------------------------------------------
# [ROCInfo] : 
# Pass in the data that already consists the predicted score and actual outcome.
# to obtain the ROC curve 
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cost.fp : associated cost for a false positive 
# @cost.fn : associated cost for a false negative 
# return   : a list containing  
#			 1. plot        : a side by side roc and cost plot, title showing optimal cutoff value
# 				 	   		  title showing optimal cutoff, total cost, and area under the curve (auc)
# 		     2. cutoff      : optimal cutoff value according to the specified fp/fn cost 
#		     3. totalcost   : total cost according to the specified fp/fn cost
#			 4. auc 		: area under the curve
#		     5. sensitivity : TP / (TP + FN)
#		     6. specificity : TN / (FP + TN)

ROCInfo <- function( data, predict, actual, cost.fp, cost.fn )
{
  # calculate the values using the ROCR library
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost 
  # false postive rate * number of negative instances * false positive cost + 
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 0 ) + 
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 1 )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)   
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) + 
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )				
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )	
  
  # the main title for the two arranged plot
  sub_title <- sprintf( "Cutoff at %.2f - Total Cost = %f, AUC = %.3f", 
                        best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2, 
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  return( list( plot 		  = plot, 
                cutoff 	  = best_cutoff, 
                totalcost   = best_cost, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr ) )
}

