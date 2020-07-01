##########################################
# ki longitudinal manuscripts
# growth failure prediction analysis

# Clean primary analysis datasets from 3 epi
# manuscripts into single dataset for prediction
# analysis by subsetting to included studies
# and merge outcome-specific datasets together
##########################################

source(paste0(here::here(), "/0-config.R"))


#-------------------------------------------
# Load and pre-process data
#-------------------------------------------


#Load data
#d <- readRDS(here("data/mock_data.RDS"))

#load covariate dataset (one row per child)
cov <- readRDS(paste0(ghapdata_dir,"FINAL_clean_covariates.rds"))
cov <- cov %>% subset(., select = c(studyid, subjid, country, tr, arm,sex,  
                                    id,  W_mage, W_mhtcm, 
                                    W_mwtkg, W_mbmi,W_meducyrs,  W_nrooms))


#Subset to included studies
#included_studies <- c("MAL-ED", "TanzaniaChild2", "JiVitA-3", "PROVIDE", "LCNI-5", "iLiNS-DOSE", "iLiNS-DYAD-M")
included_studies <- c("MAL-ED", "NIH-Birth", "JiVitA-3", "PROVIDE", "NIH-Crypto", "iLiNS-DOSE", "iLiNS-DYAD-M")
cov <- cov %>% filter(studyid %in% included_studies)
cov$tr <- as.character(cov$tr)
cov$tr[is.na(cov$tr)] <- ""
cov$tr <- factor(cov$tr)

df <- cov %>% filter(studyid=="iLiNS-DOSE")
summary(df$W_mage)
summary(df$W_mbmi)

# #Subset to relevant predictors
# df <- d %>% subset(., select = c(studyid,region, subjid, stunt, wast,haz,whz))
# cov <- d %>% subset(., select = -c(region, subjid, stunt, wast,haz,whz))
# 
# #Impute missingness by study
# table(is.na(cov$W_gagebrth))

cov <- cov %>% 
  group_by(studyid) %>% 
  do(impute_missing_values(., type = "standard", add_indicators = T, prefix = "miss_")$data) %>%
  as.data.frame()
table(cov$miss_W_mbmi)
table(cov$miss_W_mbmi[cov$studyid=="iLiNS-DOSE"])
table(cov$studyid,cov$miss_W_mbmi)

i="iLiNS-DYAD-M" 
j="miss_W_mbmi"

#Mark missingness for variables fully missing from studies and
#median/mode impute values from full data
for(i in unique(cov$studyid)){
  for(j in colnames(cov)[grepl("miss_",colnames(cov))]){
    var <- cov[cov$studyid==i, colnames(cov)==j]
    cov[cov$studyid==i, colnames(cov)==j] <- ifelse(is.na(var) & !is.na(cov[cov$studyid==i, colnames(cov)==gsub("miss_","",j)]), 0, cov[cov$studyid==i, colnames(cov)==j])
    cov[cov$studyid==i, colnames(cov)==j] <- ifelse(is.na(var) & is.na(cov[cov$studyid==i, colnames(cov)==gsub("miss_","",j)]), 1, cov[cov$studyid==i, colnames(cov)==j])
    if(sum(is.na(var)) == length(var)){
      #cov[cov$studyid==i, colnames(cov)==j] <- 1
      if(j %in% colnames(cov)){
        if(class(cov[, colnames(cov)==gsub("miss_","",j)])=="factor"){
          cov[cov$studyid==i, colnames(cov)==gsub("miss_","",j)] <- Mode(cov[, colnames(cov)==gsub("miss_","",j)])
        }else{
          cov[cov$studyid==i, colnames(cov)==gsub("miss_","",j)] <- median(cov[, colnames(cov)==gsub("miss_","",j)], na.rm=T)
        }
      }
    }
  }
}


#Replace factor imputation with "missing" category
for(i in colnames(cov)){
  if(class(cov[,i])=="factor"){
    if(paste0("miss_",i) %in% colnames(cov)){
      levels(cov[,i]) <- c(levels(cov[,i]),"missing")
      cov[cov[,paste0("miss_",i)]==1,i] <- "missing"
    }
  }
}


##Merge outcome and covariates back together
#d <- cbind(df, cov %>% subset(., select = -c(studyid)))

table(is.na(cov))



#save covariates
#saveRDS(d, here("data/clean_data.RDS"))
saveRDS(cov, file=paste0(data_dir,"covariates.RDS"))

