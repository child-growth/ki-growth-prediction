
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/functions/0_calc_empirical_IC.R"))
library(hal9001)
library(washb)


#Load data
d <- readRDS(paste0(data_dir,"analysis_data.RDS"))

#Subset to control arms
#d <- d %>% filter(arm=="" | arm=="Control" | arm=="Iron and folic acid supplementation" | arm=="Iron Folic Acid")



y_name = "wast_24" 
x_names = c(covars, "arm", "waz_birth")
dfull <- d %>% ungroup() %>% select(all_of(c("studyid","country",y_name, x_names)))

#temp: complete cases
dim(dfull)
dfull <- dfull[complete.cases(dfull),]
dfull <- droplevels(dfull)
dim(dfull)

dfull$cohort <- paste0(dfull$studyid,"-",dfull$country)
cohorts <- unique(dfull$cohort)


hal_fit_list <- readRDS(paste0(model_dir,"res_wazBirth_wast24_hal.RDS"))



observed_SE_list <- list()

for(i in cohorts){
  if(!is.null(hal_fit_list[[i]]$Y)){
    res <- data.frame(washb::washb_mean(Y=hal_fit_list[[i]]$Y, id=1:length(hal_fit_list[[i]]$Y), print=FALSE))
    res_hal=NULL
    try(res_hal <- get_pred_empirical_ci(fit=hal_fit_list[[i]]$fit_init, X=hal_fit_list[[i]]$X, Y=hal_fit_list[[i]]$Y, n=hal_fit_list[[i]]$n))
    res <- bind_cols(res, res_hal)
    observed_SE_list[[i]] <- res
  }
}

predicted_means = rbindlist(observed_SE_list, fill=TRUE)
predicted_means = data.frame(study=names(observed_SE_list), predicted_means)
predicted_means



beta_ic_list <- list()

for(i in cohorts){
  if(!is.null(hal_fit_list[[i]]$Y)){
    res_ic_beta=NULL
    try(res_ic_beta <- get_ic_beta(fit=hal_fit_list[[i]]$fit_init, X=hal_fit_list[[i]]$X, Y=hal_fit_list[[i]]$Y, n=hal_fit_list[[i]]$n))
    beta_ic_list[[i]] <- res_ic_beta
  }
}



#counterfactual predictions
counterfactual_SE_list <- list()

#shift birth waz for all
shift_function <- function(Xnew, shift_var, shift=1){
  Xnew[[shift_var]] = Xnew[[shift_var]] + shift
  return(Xnew)
}
shift= 0


#no low birthweight
# #waz for male 2500 g at birth
# -1.89878
# #waz for female 2500 g at birth
# -1.726807

# shift_function <- function(Xnew, shift_var, shift=1){
#   
#   Xnew[[shift_var]] <- ifelse(Xnew$sexMale==1 & Xnew[[shift_var]] < (-1.89878),  -1.89878, Xnew[[shift_var]])
#   Xnew[[shift_var]] <- ifelse(Xnew$sexMale==0 & Xnew[[shift_var]] < (-1.726807),  -1.726807, Xnew[[shift_var]])
#   return(Xnew)
# }
# shift= 1


#LNS for moms:
#40.96 g birthweight increase



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# To do:
# Need to use the rma escalc function to get proportion estimates
# then check the CI calculations for the counterfactuals
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


shift_var="waz_birth"
i=cohorts[1]
IC_beta=beta_ic_list[[1]]

for(i in cohorts){
  if(!is.null(hal_fit_list[[i]]$Y)){
    res <- data.frame(washb::washb_mean(Y=hal_fit_list[[i]]$Y, id=1:length(hal_fit_list[[i]]$Y), print=FALSE))
    res_hal=NULL
    Xnew <- hal_fit_list[[i]]$X 
    Xnew <- shift_function(Xnew=Xnew, shift_var=shift_var, shift=shift)

    try(res_hal <- get_counterfactual_empirical_ci(
      fit=hal_fit_list[[i]]$fit_init, 
      X=Xnew, 
      Y=hal_fit_list[[i]]$Y, 
      n=hal_fit_list[[i]]$n,
      IC_beta=beta_ic_list[[i]])) 
    res <- bind_cols(res, res_hal)
    counterfactual_SE_list[[i]] <- res
  }
}

counterfactual_means = rbindlist(counterfactual_SE_list, fill=TRUE)
counterfactual_means = data.frame(study=names(counterfactual_SE_list), counterfactual_means)
counterfactual_means


df <- counterfactual_means %>% filter(!is.na(se) & se!=0, lb!=ub)

method_fit="REML"
fit <- NULL
fit <- rma(yi=df$predY, sei=df$se, measure="MN",method=method_fit)
if(is.null(fit)){fit <- rma(yi=df$predY, sei=df$se, measure="MN",method="DL")}

res = data.frame(
  est = (fit$beta),
  se = fit$se,
  lb = (fit$beta - 1.96 * fit$se),
  ub = (fit$beta + 1.96 * fit$se),
  nmeas.f = paste0("N=", format(sum(df$N), big.mark = ",", scientific = FALSE)),
  nstudy.f = paste0("N=", nrow(df), " studies"),
  method.used=method_fit
) 
res

method_fit="FE"
fit_FE <- rma(yi=df$predY, sei=df$se, measure="MN",method=method_fit)

res_FE = data.frame(
  est = (fit_FE$beta),
  se = fit_FE$se,
  lb = (fit_FE$beta - 1.96 * fit_FE$se),
  ub = (fit_FE$beta + 1.96 * fit_FE$se),
  nmeas.f = paste0("N=", format(sum(df$N), big.mark = ",", scientific = FALSE)),
  nstudy.f = paste0("N=", nrow(df), " studies"),
  method.used=method_fit
) 
res_FE


method_fit="REML"
fit_obs <- rma(yi=df$Mean, sei=df$Robust.SE, measure="MN",method=method_fit)

res_obs = data.frame(
  est = (fit_obs$beta),
  se = fit_obs$se,
  lb = (fit_obs$beta - 1.96 * fit_obs$se),
  ub = (fit_obs$beta + 1.96 * fit_obs$se),
  nmeas.f = paste0("N=", format(sum(df$N), big.mark = ",", scientific = FALSE)),
  nstudy.f = paste0("N=", nrow(df), " studies"),
  method.used=method_fit
) 
res_obs


fit_obs_FE <- rma(yi=df$Mean, sei=df$Robust.SE, measure="MN",method="FE")
res_obs_FE = data.frame(
  est = (fit_obs$beta),
  se = fit_obs$se,
  lb = (fit_obs$beta - 1.96 * fit_obs$se),
  ub = (fit_obs$beta + 1.96 * fit_obs$se),
  nmeas.f = paste0("N=", format(sum(df$N), big.mark = ",", scientific = FALSE)),
  nstudy.f = paste0("N=", nrow(df), " studies"),
  method.used="FE"
) 
res_obs_FE


#Get the pooled observed estimate



plotdf <- data.frame(
  bind_rows(
    df %>% filter(!is.na(predY)) %>%
      select(study, Mean,Lower.95.CI, Upper.95.CI) %>%
      rename(est=Mean, 
             lb=Lower.95.CI, 
             ub=Upper.95.CI) %>%
      mutate(estimate="observed"),
    df %>% select(study, predY, lb, ub) %>% 
      rename(est=predY) %>%
      mutate(estimate="counterfactual"),
    res %>% select(est, lb, ub) %>% mutate(estimate="counterfactual", study="Pooled"),
    res_obs %>% select(est, lb, ub) %>% mutate(estimate="observed", study="Pooled"),
    res_FE %>% select(est, lb, ub) %>% mutate(estimate="counterfactual", study="Pooled-FE"),
    res_obs %>% select(est, lb, ub) %>% mutate(estimate="observed", study="Pooled-FE"))
) %>% filter(lb!=ub)
plotdf$study = factor(plotdf$study, levels=c("Pooled-FE", "Pooled", unique(df$study)))


#----------------------------------------------------------------------------------------------
# plot
#----------------------------------------------------------------------------------------------

#plot parameters
tableau11 <- c("Black","#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")


# plot
p1 <- ggplot(plotdf, aes(x=study, y=est, color=estimate, group=estimate)) + 
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin=lb, ymax=ub), position = position_dodge(0.5)) +
  geom_vline(xintercept = 2.5, linetype="dashed", color="grey50") +
  scale_color_manual(values = tableau11[2:3]) +
  coord_flip() +
  theme(legend.position = "bottom")
p1

saveRDS(list(p=p1, plotdf=plotdf), paste0(here::here(),"/figures/waz_wast_shift_example.RDS"))
#saveRDS(list(p=p1, plotdf=plotdf), paste0(here::here(),"/figures/birthweight_shift_example.RDS"))

#Add the FE for both observed and counterfactuals






#can also get a mean difference between the observed and predicted
#also get the proportion wasted before and after the shift from 
# a new function that returns Y_hat_init

# #to do: add region to the output
# #convert waz to birthweight
# library(zscorer)
# https://nutriverse.io/zscorer/reference/getWGS.html
# getWGS(sexObserved=2, firstPart=2.5, secondPart=0, index="wfa")
# 
# library(growthstandards)
# 
# zscorer::run_zscorer()
