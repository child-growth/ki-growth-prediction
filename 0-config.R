

#-------------------------------------
# ki longitudinal analysis manuscripts

# configure data directories
# source base functions
# load libraries
#-------------------------------------
# kiPath <- c("/data/KI/R/x86_64-pc-linux-gnu-library/3.6/" , .libPaths())
# .libPaths(kiPath)


library(here)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(metafor)
library(data.table)
library(mada)
library(pROC)

library(assertthat)
# library(haven)
library(tidyselect)
# options(repos = c(CRAN = "http://cran.rstudio.com/",
#                   deltarho = "http://packages.deltarho.org"))
library(stringr)
library(mgcv)
library(grid)
library(lazyeval)
library(rlang)
library(scales)


# for parallel computing 
# (will need to configure in each script)
library(foreach)
library(doParallel)


library(tlverse)
library(sl3)
library(tmle3)
library(skimr)
library(data.table)
library(SuperLearner)
library(ck37r)
library(caret)
library(ROCR)    
library(origami)
library(future)
library(cvAUC)

# Source base functions  
source(here("/functions/0_pooling_functions.R"))

source(here("/functions/0_prediction_functions.R"))

#Data directories
ghapdata_dir = "/data/KI/UCB-SuperLearner/Manuscript analysis data/"
data_dir = "/data/KI/SLprediction/"


# Set plot themes
theme_ki <- function() {
  theme_bw() %+replace%
    theme(
      strip.background = element_blank(),
      legend.position="none",
      #panel.grid.major = element_line(),
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size=14),
      axis.title = element_text(size=12),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=.1)
    )
}

#hbgdki pallets
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
tableau11 <- c("Black","#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

theme_set(theme_ki())



#Pick prediction variables 
covars <- c("sex", "W_mage", "W_meducyrs", "W_nrooms",
            "W_mwtkg",  "W_mhtcm","W_mbmi", 
            "miss_W_mage", "miss_W_meducyrs", "miss_W_nrooms",
            "miss_W_mwtkg",  "miss_W_mhtcm","miss_W_mbmi")
covars_birth <- c("haz_birth", "whz_birth", "waz_birth", "stunt_birth","wast_birth","underwt_birth", "sstunt_birth","swast_birth","sunderwt_birth", covars)
covars3 <- c("haz_3", "whz_3", "waz_3", "stunt_3","wast_3","underwt_3", "sstunt_3","swast_3","sunderwt_3", covars)
covars6 <- c("haz_6", "whz_6", "waz_6", "stunt_6","wast_6","underwt_6", "sstunt_6","swast_6","sunderwt_6", covars)
covars9 <- c("haz_9", "whz_9", "waz_9", "stunt_9","wast_9","underwt_9", "sstunt_9","swast_9","sunderwt_9", covars)
covars12 <- c("haz_12", "whz_12", "waz_12", "stunt_12","wast_12","underwt_12", "sstunt_12","swast_12","sunderwt_12", covars)

#Set parameters
CV_setting = FALSE
