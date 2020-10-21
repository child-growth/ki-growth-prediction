


source(paste0(here::here(), "/0-config.R"))


#load data
d <- readRDS(paste0(ghapdata_dir, "ki-manuscript-dataset.rds"))
d <- d %>% subset(., select=c(studyid, subjid, country, agedays, haz, whz, waz))



#load covariate dataset (one row per child)
cov <- readRDS(paste0(data_dir,"covariates.RDS"))



# included_studies <- c("MAL-ED", "NIH-Birth", "JiVitA-3", "PROVIDE", "NIH-Crypto", "iLiNS-DOSE", "iLiNS-DYAD-M")
# 
# d <- d %>% filter(agedays < 25 * 30.4167) %>%
#   filter(studyid %in% included_studies)

#Temp use only a few mal-ed cohorts
# included_studies <- c("MAL-ED")
# d <- d %>% filter(agedays < 25 * 30.4167) %>%
#   filter(studyid %in% included_studies, country %in% c("BANGLADESH", "BRAZIL", "INDIA"))




#get id dataset
id <- d %>% subset(., select = c(studyid, country, subjid)) %>%
      group_by(studyid, country, subjid) %>% slice(1)

#--------------------------------------------
# Label age by month
#--------------------------------------------

calc.monthly.agecat <-function(d){
  d$agecat <- cut(d$agedays,
                  breaks = c(0:26) * 30.4167 - 30.4167 / 2, include.lowest = F,
                  labels = paste0(0:25, " months")
  )
  levels(d$agecat)[1] <- "Two weeks"
  levels(d$agecat)[2] <- "One month"
  table(d$agecat)
  return(d)
}

d <- calc.monthly.agecat(d)
table(d$agecat)


#--------------------------------------------
# Split data by outcome and 
# drop unrealistic measures depending on 
# anthropometry measure
#--------------------------------------------

stunt <- d %>% filter(haz >= -6 & haz <=6, !is.na(haz)) %>%
  subset(., select = - c(whz, waz)) %>%
  mutate(stunt = 1*(haz < -2), sstunt =  1*(haz < -3))

wast <- d %>% filter(whz >= -5 & whz <=5, !is.na(whz)) %>%
  subset(., select = - c(haz, waz )) %>%
  mutate(wast = 1*(whz < -2), swast =  1*(whz < -3))

waz <- d %>% filter(waz >= -6 & waz <=5, !is.na(waz)) %>%
  subset(., select = - c(haz, whz)) %>%
  mutate(underwt = 1*(waz < -2), sunderwt =  1*(waz < -3))

#--------------------------------------------
# Get monthly Z-score measurements for 1,3,6,9,12,18, and 24 months
#--------------------------------------------

get_monthly_meas <- function(d , month=3, birth=F){
  if(birth){
    df <- d %>% group_by(studyid,country,subjid) %>% filter(agedays < 8, agedays==min(agedays)) %>%
      mutate(agecat="birth")
  }else{
    df <- d %>% group_by(studyid,country,subjid) %>% 
      filter(agedays >= (month-1)*30.4167 & agedays <= (month+1)*30.4167) %>%
      filter(abs(agedays-month*30.4167)==min(abs(agedays-month*30.4167)))
    df$agecat <- as.character(month)
  }
  return(df)
}

stunt_birth <- get_monthly_meas(stunt, birth=T)
stunt3 <- get_monthly_meas(stunt, 3)
stunt6 <- get_monthly_meas(stunt, 6)
stunt9 <- get_monthly_meas(stunt, 9)
stunt12 <- get_monthly_meas(stunt, 12)
stunt18 <- get_monthly_meas(stunt, 18)
stunt24 <- get_monthly_meas(stunt, 24)

#compile and transform to wide
stunt_long <- rbind(stunt_birth, stunt3, stunt6, stunt9, stunt12, stunt18, stunt24)
stunt_wide <- stunt_long %>% subset(., select = -c(agedays)) %>%
  pivot_wider(names_from = agecat, values_from = c(haz, stunt, sstunt))
head(stunt_wide)
dim(stunt_long)
dim(stunt_wide)

wast_birth <- get_monthly_meas(wast, birth=T)
wast3 <- get_monthly_meas(wast, 3)
wast6 <- get_monthly_meas(wast, 6)
wast9 <- get_monthly_meas(wast, 9)
wast12 <- get_monthly_meas(wast, 12)
wast18 <- get_monthly_meas(wast, 18)
wast24 <- get_monthly_meas(wast, 24)

#compile and transform to wide
wast_long <- rbind(wast_birth, wast3, wast6, wast9, wast12, wast18, wast24)
wast_wide <- wast_long %>% subset(., select = -c(agedays)) %>%
  pivot_wider(names_from = agecat, values_from = c(whz, wast, swast))

waz_birth <- get_monthly_meas(waz, birth=T)
waz3 <- get_monthly_meas(waz, 3)
waz6 <- get_monthly_meas(waz, 6)
waz9 <- get_monthly_meas(waz, 9)
waz12 <- get_monthly_meas(waz, 12)
waz18 <- get_monthly_meas(waz, 18)
waz24 <- get_monthly_meas(waz, 24)

#compile and transform to wide
waz_long <- rbind(waz_birth, waz3, waz6, waz9, waz12, waz18, waz24)
waz_wide <- waz_long %>% subset(., select = -c(agedays)) %>%
  pivot_wider(names_from = agecat, values_from = c(waz, underwt, sunderwt))


#merge together datasets and merge in covariates
dim(id)
dim(stunt_wide)
d <- left_join(id, stunt_wide, by = c("studyid", "subjid", "country"))
dim(d)

d <- left_join(d, wast_wide, by = c("studyid", "subjid", "country"))
dim(d)

d <- left_join(d, waz_wide, by = c("studyid", "subjid", "country"))
dim(d)

d <- left_join(d, cov, by = c("studyid", "subjid", "country"))
dim(d)

colnames(d)

saveRDS(d, file=paste0(data_dir,"analysis_data.RDS"))




#save study metadata


#N's for manuscript
# haz %>% ungroup() %>%
#   summarize(Nstudies=length(unique(paste0(studyid, country))),
#             Nchildren=length(unique(paste0(studyid, subjid))))
# whz %>% ungroup() %>%
#   summarize(Nstudies=length(unique(paste0(studyid, country))),
#             Nchildren=length(unique(paste0(studyid, subjid))))

table(d$studyid, !is.na(d$haz_birth))
table(d$studyid, !is.na(d$haz_3))
table(d$studyid, !is.na(d$haz_6))
table(d$studyid, !is.na(d$haz_9))
table(d$studyid, !is.na(d$haz_12))
table(d$studyid, !is.na(d$haz_18))
table(d$studyid, !is.na(d$haz_24))

