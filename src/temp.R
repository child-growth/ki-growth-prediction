

source(paste0(here::here(), "/0-config.R"))

#Load data
d <- readRDS(paste0(ghapdata_dir, "ki-manuscript-dataset.rds")) %>% filter(haz >= -6 & haz <=6, !is.na(haz), agedays < 731) 


ggplot(d, aes(x=agedays)) + geom_histogram()

