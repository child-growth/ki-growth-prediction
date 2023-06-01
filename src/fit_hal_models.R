

rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/functions/0_calc_empirical_IC.R"))
library(hal9001)
library(washb)
library(caret)

#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))

#Subset to control arms
#d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")



y_name = "whz_24" 
x_names = c(covars, "arm", "waz_birth")
dfull <- d %>% ungroup() %>% select(all_of(c("studyid","country",y_name, x_names)))

#temp: complete cases
dim(dfull)
dfull <- dfull[complete.cases(dfull),]
dfull <- droplevels(dfull)
dim(dfull)

dfull$cohort <- paste0(dfull$studyid,"-",dfull$country)

hal_fit_list <- list()

cohorts <- unique(dfull$cohort)
i=cohorts[1]

for(i in cohorts){
  df <- dfull %>% filter(cohort==i)
  
  n=nrow(df)
  
  if(n>=50){
  Y <- as.numeric(as.matrix(df %>% select(all_of(y_name))))
  X <- df %>% ungroup() %>%
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.character), as.factor)
  
  NZV <- nearZeroVar(X, saveMetrics = TRUE)
  X <- X[,!NZV$zeroVar]

  X <-design_matrix(as.data.frame(X))
  
  fit_init <- fit_hal(X = X,
                      Y = Y,
                      smoothness_orders = 0,
                      return_x_basis = TRUE,
                      family = "gaussian",
                      num_knots = hal9001:::num_knots_generator(
                        max_degree = ifelse(ncol(X) >= 20, 2, 3),
                        smoothness_orders = 0,
                        base_num_knots_0 = max(100, ceiling(sqrt(n)))
                      ))
  
  res_list = list(Y=Y, X=X, n=n, fit_init=fit_init)
  hal_fit_list[[i]] <- res_list
  }
}


saveRDS(hal_fit_list, paste0(model_dir,"res_wazBirth_whz24_hal.RDS"))

