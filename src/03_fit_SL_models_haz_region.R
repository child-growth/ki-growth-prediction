





source(paste0(here::here(), "/0-config.R"))




#Load data
#d <- readRDS(here("data/clean_data.RDS"))
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))
#Subset to control arms
d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")
d <- droplevels(d)

d <- d %>% mutate(
  region =  case_when(
    country %in% c("BANGLADESH","INDIA","NEPAL") ~ "South Asia",
    country %in% c("BRAZIL","PERU") ~ "Latin America",
    country %in% c("MALAWI","SOUTH AFRICA","TANZANIA, UNITED REPUBLIC OF") ~ "Africa"
  ))

#Set parameter
CV_setting = TRUE


#wrapper function to fit and clean SL
fit_SuperLearner <- function(d, outcome, covars, slmod=sl, CV=CV_setting, include_anthro=T, covars2=NULL){
  if(include_anthro){
    res <- d %>% group_by(region) %>% 
      do(clean_prediction_dataset(., outcome = outcome, covars=covars, covars2=covars2)) %>%
      do(fit=try(fit_SL(., outcome = outcome, covars=unique(c(covars,covars2)), slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }else{
    res <- d %>% group_by(region) %>% 
      do(fit=try(fit_SL(., outcome = outcome, covars=covars, slmod=slmod, CV=CV_setting)))%>%
      ungroup() %>% do(clean_res(.))
  }
  return(res)
}


#24 month haz
res_haz24_covars <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars, slmod=sl, CV=CV_setting, include_anthro=F)
res_haz24_covars_birth <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars_birth, slmod=sl, CV=CV_setting)
res_haz24_covars3 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars3, slmod=sl, CV=CV_setting)
res_haz24_covars6 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars6, slmod=sl, CV=CV_setting)
res_haz24_covars9 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars9, slmod=sl, CV=CV_setting)
res_haz24_covars12 <- fit_SuperLearner(d,  outcome = "haz_24", covars=covars12, slmod=sl, CV=CV_setting)

saveRDS(res_haz24_covars, here("prediction-models/region_res_haz24_covars.RDS"))
saveRDS(res_haz24_covars_birth, here("prediction-models/region_res_haz24_covars_birth.RDS"))
saveRDS(res_haz24_covars3, here("prediction-models/region_res_haz24_covars3.RDS"))
saveRDS(res_haz24_covars6, here("prediction-models/region_res_haz24_covars6.RDS"))
saveRDS(res_haz24_covars9, here("prediction-models/region_res_haz24_covars9.RDS"))
saveRDS(res_haz24_covars12, here("prediction-models/region_res_haz24_covars12.RDS"))

