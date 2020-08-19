



source(paste0(here::here(), "/0-config.R"))




#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))
#Subset to control arms
d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")
d <- droplevels(d)

#Set cross-valisation parameter
CV_setting = FALSE

#-------------------------------------
#Fit Superlearner model for a single cohort
#-------------------------------------


df <- d %>% filter(studyid=="MAL-ED", country=="BANGLADESH")

#fit model
start_time <- Sys.time()
test <- fit_SL(df, outcome = "haz_12", family="gaussian", covars=covars12, slmod=sl, CV=FALSE)
end_time <- Sys.time()

end_time - start_time


start_time2 <- Sys.time()
test2 <- fit_SL(df, outcome = "haz_12", family="gaussian", covars=covars12, slmod=sl, CV=TRUE)
end_time2 <- Sys.time()

end_time2 - start_time2

#look at performance metrics:
test$result$perf_metrics



#-------------------------------------
# Loop over cohorts
#-------------------------------------

#wrapper function to fit and clean SL
fit_SuperLearner <- function(d, outcome, covars, slmod=sl, CV=CV_setting, family="gaussian", include_anthro=T, covars2=NULL){
  if(include_anthro){
    res <- d %>% group_by(studyid, country) %>% 
      do(clean_prediction_dataset(., outcome = outcome, covars=covars, covars2=covars2)) %>%
      do(fit=try(fit_SL(., outcome = outcome, family=family, covars=unique(c(covars,covars2)), slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }else{
    res <- d %>% group_by(studyid, country) %>% 
      do(fit=try(fit_SL(., outcome = outcome, family=family, covars=covars, slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }
  return(res)
}

#12 month haz, birth covariates
res_haz12_covars_birth <- fit_SuperLearner(d,  outcome = "haz_9", covars=covars12, slmod=sl, CV=CV_setting)

#12 month stunting, birth covariates
res_stunt12_covars_birth <- fit_SuperLearner(d,  outcome = "stunt_12", covars=covars12, family="binomial", slmod=sl, CV=CV_setting)



#Test for pooled R2
res_haz12_covars_birth <- fit_SuperLearner(d,  outcome = "haz_9", covars=covars12, slmod=sl, CV=CV_setting)
pR2 = pool_R2(res=res_haz12_covars_birth, Y="haz_9", age= NULL, covars =covars12)



