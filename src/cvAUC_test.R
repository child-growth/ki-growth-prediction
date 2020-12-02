

source(paste0(here::here(), "/0-config.R"))
library(cvAUC)

#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))
#Subset to control arms
d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")
d <- droplevels(d)

df <- d %>% filter(studyid=="MAL-ED", country=="BANGLADESH")

# test <- fit_SL(df, outcome = "stunt_18", family="gaussian", 
#                covars=covars12, slmod=sl, CV=T)



dat = df
outcome = "stunt_18"
id="subjid"
family="binomial"
covars=covars12
slmod=sl
CV=TRUE
folds=5


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
#0903 test
sl_fit <- slmod$train(SL_task)
glm_fit <- glm_sl$train(SL_task)

# start_time <- Sys.time()
cv_sl <- make_learner(Lrnr_cv, sl_fit, full_fit = TRUE)
fit_cv_sl <- cv_sl$train(SL_task)
yhat_full <- fit_cv_sl$predict_fold(SL_task,"validation")
# end_time <- Sys.time()
# print(end_time - start_time)
  
  
  
  
  


fold_index <- origami::folds2foldvec(SL_task$folds) 
outcome <- SL_task$Y

auc.plotdf <- cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index)
auc.ci <-ci.cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index, confidence=0.95)
naiveAUC <- AUC(predictions=yhat_full, labels=outcome)

foldAUC = rep(NA, 5)
for (i in 1:5){
  preds <- yhat_full[cv_folds[[i]]$validation_set]
  obs <- outcome[cv_folds[[i]]$validation_set]
  foldAUC[i] <- AUC(predictions=preds, labels=obs)
}

naiveAUC
auc.plotdf$cvAUC
auc.ci$cvAUC
mean(foldAUC)

# fold_index2 <- list(cv_folds[[1]]$validation_set,
#                     cv_folds[[2]]$validation_set,
#                     cv_folds[[3]]$validation_set,
#                     cv_folds[[4]]$validation_set,
#                     cv_folds[[5]]$validation_set)


# .cvFolds <- function(Y, V){
#   # Create CV folds (stratify by outcome)	
#   Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
#   Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
#   folds <- vector("list", length=V)
#   for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}  	
#   return(folds)
# }
# folds <- .cvFolds(Y = outcome, V = 5)













