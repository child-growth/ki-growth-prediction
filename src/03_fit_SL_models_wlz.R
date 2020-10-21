




source(paste0(here::here(), "/0-config.R"))




#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))
#Subset to control arms
d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")
d <- droplevels(d)

#Set parameter
CV_setting = TRUE


#wrapper function to fit and clean SL
fit_SuperLearner <- function(d, outcome, covars, slmod=sl, CV_setting=CV_setting, include_anthro=T, covars2=NULL){
  if(include_anthro){
    res <- d %>% group_by(studyid, country) %>% 
      do(clean_prediction_dataset(., outcome = outcome, covars=covars, covars2=covars2)) %>%
      do(fit=try(fit_SL(., outcome = outcome, covars=unique(c(covars,covars2)), slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }else{
    res <- d %>% group_by(studyid, country) %>% 
      do(fit=try(fit_SL(., outcome = outcome, covars=covars, slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }
  return(res)
}

# source(paste0(here::here(), "/0-config.R"))
# d<- d %>% filter(studyid=="iLiNS-DOSE")
# res <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars12, slmod=sl, CV=CV_setting)
# 
# res$fit[[1]]$result$run_time
# 
# res$fit[[1]]$result$perf_metrics$R2_full

#d <- d %>% filter(country=="TANZANIA, UNITED REPUBLIC OF")

start_time <- Sys.time()


#12 month whz
res_whz12_covars <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_whz12_covars_birth <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars_birth, slmod=sl, CV=CV_setting)
res_whz12_covars3 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars3, slmod=sl, CV=CV_setting)
res_whz12_covars6 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars6, slmod=sl, CV=CV_setting)
res_whz12_covars9 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars9, slmod=sl, CV=CV_setting)

saveRDS(res_whz12_covars, paste0(model_dir,"res_whz12_covars.RDS"))
saveRDS(res_whz12_covars_birth, paste0(model_dir,"res_whz12_covars_birth.RDS"))
saveRDS(res_whz12_covars3, paste0(model_dir,"res_whz12_covars3.RDS"))
saveRDS(res_whz12_covars6, paste0(model_dir,"res_whz12_covars6.RDS"))
saveRDS(res_whz12_covars9, paste0(model_dir,"res_whz12_covars9.RDS"))


#18 month whz
res_whz18_covars <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_whz18_covars_birth <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars_birth, slmod=sl, CV=CV_setting)
res_whz18_covars3 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars3, slmod=sl, CV=CV_setting)
res_whz18_covars6 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars6, slmod=sl, CV=CV_setting)
res_whz18_covars9 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars9, slmod=sl, CV=CV_setting)
res_whz18_covars12 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_whz18_covars, paste0(model_dir,"res_whz18_covars.RDS"))
saveRDS(res_whz18_covars_birth, paste0(model_dir,"res_whz18_covars_birth.RDS"))
saveRDS(res_whz18_covars3, paste0(model_dir,"res_whz18_covars3.RDS"))
saveRDS(res_whz18_covars6, paste0(model_dir,"res_whz18_covars6.RDS"))
saveRDS(res_whz18_covars9, paste0(model_dir,"res_whz18_covars9.RDS"))
saveRDS(res_whz18_covars12, paste0(model_dir,"res_whz18_covars12.RDS"))

#24 month whz
res_whz24_covars <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_whz24_covars_birth <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars_birth, slmod=sl, CV=CV_setting)
res_whz24_covars3 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars3, slmod=sl, CV=CV_setting)
res_whz24_covars6 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars6, slmod=sl, CV=CV_setting)
res_whz24_covars9 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars9, slmod=sl, CV=CV_setting)
res_whz24_covars12 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_whz24_covars, paste0(model_dir,"res_whz24_covars.RDS"))
saveRDS(res_whz24_covars_birth, paste0(model_dir,"res_whz24_covars_birth.RDS"))
saveRDS(res_whz24_covars3, paste0(model_dir,"res_whz24_covars3.RDS"))
saveRDS(res_whz24_covars6, paste0(model_dir,"res_whz24_covars6.RDS"))
saveRDS(res_whz24_covars9, paste0(model_dir,"res_whz24_covars9.RDS"))
saveRDS(res_whz24_covars12, paste0(model_dir,"res_whz24_covars12.RDS"))



#Predictions with growth velocity measures included
res_whz12_covars_birth_3 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_whz12_covars_3_6 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_whz12_covars_6_9 <- fit_SuperLearner(d,  outcome = "whz_12", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)

res_whz18_covars_birth_3 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_whz18_covars_3_6 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_whz18_covars_6_9 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)
res_whz18_covars_9_12 <- fit_SuperLearner(d,  outcome = "whz_18", covars=covars9, covars2=covars12, slmod=sl, CV=CV_setting)

res_whz24_covars_birth_3 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_whz24_covars_3_6 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_whz24_covars_6_9 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)
res_whz24_covars_9_12 <- fit_SuperLearner(d,  outcome = "whz_24", covars=covars9, covars2=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_whz12_covars_birth_3, paste0(model_dir,"res_whz12_covars_birth_3.RDS"))
saveRDS(res_whz12_covars_3_6, paste0(model_dir,"res_whz12_covars3_6.RDS"))
saveRDS(res_whz12_covars_6_9, paste0(model_dir,"res_whz12_covars6_9.RDS"))

saveRDS(res_whz18_covars_birth_3, paste0(model_dir,"res_whz18_covars_birth_3.RDS"))
saveRDS(res_whz18_covars_3_6, paste0(model_dir,"res_whz18_covars3_6.RDS"))
saveRDS(res_whz18_covars_6_9, paste0(model_dir,"res_whz18_covars6_9.RDS"))
saveRDS(res_whz18_covars_9_12, paste0(model_dir,"res_whz18_covars9_12.RDS"))

saveRDS(res_whz24_covars_birth_3, paste0(model_dir,"res_whz24_covars_birth_3.RDS"))
saveRDS(res_whz24_covars_3_6, paste0(model_dir,"res_whz24_covars3_6.RDS"))
saveRDS(res_whz24_covars_6_9, paste0(model_dir,"res_whz24_covars6_9.RDS"))
saveRDS(res_whz24_covars_9_12, paste0(model_dir,"res_whz24_covars9_12.RDS"))


end_time <- Sys.time()

timerun <- end_time - start_time
timerun


