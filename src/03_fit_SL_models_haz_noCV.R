




source(paste0(here::here(), "/0-config.R"))




#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))
#Subset to control arms
d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")
d <- droplevels(d)

#Set parameter
CV_setting = FALSE


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
# res <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars12, slmod=sl, CV=CV_setting)
# 
# res$fit[[1]]$result$run_time
# 
# res$fit[[1]]$result$perf_metrics$R2_full

#d <- d %>% filter(country=="TANZANIA, UNITED REPUBLIC OF")

start_time <- Sys.time()


#12 month haz
res_haz12_covars <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_haz12_covars_birth <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars_birth, slmod=sl, CV=CV_setting)
res_haz12_covars3 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars3, slmod=sl, CV=CV_setting)
res_haz12_covars6 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars6, slmod=sl, CV=CV_setting)
res_haz12_covars9 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars9, slmod=sl, CV=CV_setting)

saveRDS(res_haz12_covars, paste0(model_dir,"res_haz12_covars_noCV.RDS"))
saveRDS(res_haz12_covars_birth, paste0(model_dir,"res_haz12_covars_birth_noCV.RDS"))
saveRDS(res_haz12_covars3, paste0(model_dir,"res_haz12_covars3_noCV.RDS"))
saveRDS(res_haz12_covars6, paste0(model_dir,"res_haz12_covars6_noCV.RDS"))
saveRDS(res_haz12_covars9, paste0(model_dir,"res_haz12_covars9_noCV.RDS"))


#18 month haz
res_haz18_covars <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_haz18_covars_birth <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars_birth, slmod=sl, CV=CV_setting)
res_haz18_covars3 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars3, slmod=sl, CV=CV_setting)
res_haz18_covars6 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars6, slmod=sl, CV=CV_setting)
res_haz18_covars9 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars9, slmod=sl, CV=CV_setting)
res_haz18_covars12 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_haz18_covars, paste0(model_dir,"res_haz18_covars_noCV.RDS"))
saveRDS(res_haz18_covars_birth, paste0(model_dir,"res_haz18_covars_birth_noCV.RDS"))
saveRDS(res_haz18_covars3, paste0(model_dir,"res_haz18_covars3_noCV.RDS"))
saveRDS(res_haz18_covars6, paste0(model_dir,"res_haz18_covars6_noCV.RDS"))
saveRDS(res_haz18_covars9, paste0(model_dir,"res_haz18_covars9_noCV.RDS"))
saveRDS(res_haz18_covars12, paste0(model_dir,"res_haz18_covars12_noCV.RDS"))

#24 month haz
res_haz24_covars <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_haz24_covars_birth <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars_birth, slmod=sl, CV=CV_setting)
res_haz24_covars3 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars3, slmod=sl, CV=CV_setting)
res_haz24_covars6 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars6, slmod=sl, CV=CV_setting)
res_haz24_covars9 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars9, slmod=sl, CV=CV_setting)
res_haz24_covars12 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_haz24_covars, paste0(model_dir,"res_haz24_covars_noCV.RDS"))
saveRDS(res_haz24_covars_birth, paste0(model_dir,"res_haz24_covars_birth_noCV.RDS"))
saveRDS(res_haz24_covars3, paste0(model_dir,"res_haz24_covars3_noCV.RDS"))
saveRDS(res_haz24_covars6, paste0(model_dir,"res_haz24_covars6_noCV.RDS"))
saveRDS(res_haz24_covars9, paste0(model_dir,"res_haz24_covars9_noCV.RDS"))
saveRDS(res_haz24_covars12, paste0(model_dir,"res_haz24_covars12_noCV.RDS"))



#Predictions with growth velocity measures included
res_haz12_covars_birth_3 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_haz12_covars_3_6 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_haz12_covars_6_9 <- fit_SuperLearner(d,  outcome = "haz_12", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)

res_haz18_covars_birth_3 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_haz18_covars_3_6 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_haz18_covars_6_9 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)
res_haz18_covars_9_12 <- fit_SuperLearner(d,  outcome = "haz_18", covars=covars9, covars2=covars12, slmod=sl, CV=CV_setting)

res_haz24_covars_birth_3 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars_birth, covars2=covars3, slmod=sl, CV=CV_setting)
res_haz24_covars_3_6 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars3, covars2=covars6, slmod=sl, CV=CV_setting)
res_haz24_covars_6_9 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars6, covars2=covars9, slmod=sl, CV=CV_setting)
res_haz24_covars_9_12 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars9, covars2=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_haz12_covars_birth_3, paste0(model_dir,"res_haz12_covars_birth_3_noCV.RDS"))
saveRDS(res_haz12_covars_3_6, paste0(model_dir,"res_haz12_covars3_6_noCV.RDS"))
saveRDS(res_haz12_covars_6_9, paste0(model_dir,"res_haz12_covars6_9_noCV.RDS"))

saveRDS(res_haz18_covars_birth_3, paste0(model_dir,"res_haz18_covars_birth_3_noCV.RDS"))
saveRDS(res_haz18_covars_3_6, paste0(model_dir,"res_haz18_covars3_6_noCV.RDS"))
saveRDS(res_haz18_covars_6_9, paste0(model_dir,"res_haz18_covars6_9_noCV.RDS"))
saveRDS(res_haz18_covars_9_12, paste0(model_dir,"res_haz18_covars9_12_noCV.RDS"))

saveRDS(res_haz24_covars_birth_3, paste0(model_dir,"res_haz24_covars_birth_3_noCV.RDS"))
saveRDS(res_haz24_covars_3_6, paste0(model_dir,"res_haz24_covars3_6_noCV.RDS"))
saveRDS(res_haz24_covars_6_9, paste0(model_dir,"res_haz24_covars6_9_noCV.RDS"))
saveRDS(res_haz24_covars_9_12, paste0(model_dir,"res_haz24_covars9_12_noCV.RDS"))


end_time <- Sys.time()

timerun <- end_time - start_time

saveRDS(timerun, paste0(model_dir,"time_benchmark_noCV_noCV.RDS"))


