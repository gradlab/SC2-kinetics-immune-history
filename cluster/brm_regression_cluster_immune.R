## Perform logistic regression on various subsets of the data
## Preamble
library(tidyverse)
library(zoo)
library(lubridate)
library(patchwork)
library(splines)
library(brms)
library(mgcv)
library(data.table)
library(ROCR)
library(tidybayes)
library(Rcpp)
library(future)

options(future.fork.multithreading.enable = FALSE)

#setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")
setwd("~/SC2-kinetics-immune-history")

n_iter <- 2000
rerun_stan <- FALSE
load("data/data_for_regressions.RData")

## For these analyses, we only want to use Ct values after detection
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0)
dat_subset_use <- dat_subset_use %>% filter(!is.na(AgeGroup))

filename_base <- paste0("outputs/immune_models")
if(!file.exists(filename_base)) dir.create(filename_base)

## 24 entries, 12 for each data type
#i <- 10
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(i)

formulas <- list(
    ## Baseline
    bf(low_ct1 ~ s(DaysSinceDetection)),
    
    ## Simple predictors
    bf(low_ct1 ~ VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus)),
    bf(low_ct1 ~ CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=CumulativeExposureNumber)),
    bf(low_ct1 ~ DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=DaysSinceExposureGroup)),
    bf(low_ct1 ~ TiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroup)),
    #bf(low_ct1 ~ TiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroupAlt)),
   # bf(low_ct1 ~ LineageBroad + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)),
    
    ## Add interactions with lineage
    bf(low_ct1 ~ LineageBroad_VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus)),
    bf(low_ct1 ~ LineageBroad_CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_CumulativeExposureNumber)),
    bf(low_ct1 ~ LineageBroad_DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_DaysSinceExposureGroup)),
    #bf(low_ct1 ~ LineageBroad_TiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroup)),
    #bf(low_ct1 ~ LineageBroad_TiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroupAlt)),
    
    ## All of the above but with age
    ## Simple predictors
   
    bf(low_ct1 ~ AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ VaccStatus + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus) + s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ CumulativeExposureNumber + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=CumulativeExposureNumber)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ DaysSinceExposureGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=DaysSinceExposureGroup)+ s(DaysSinceDetection,by=AgeGroup)),
    #bf(low_ct1 ~ TiterGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroup)),
    #bf(low_ct1 ~ TiterGroupAlt + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroupAlt)),
    bf(low_ct1 ~ LineageBroad + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)+ s(DaysSinceDetection,by=AgeGroup)),
    
    ## Add interactions with lineage
    bf(low_ct1 ~ LineageBroad_VaccStatus + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ LineageBroad_CumulativeExposureNumber + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_CumulativeExposureNumber)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ LineageBroad_DaysSinceExposureGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_DaysSinceExposureGroup)+ s(DaysSinceDetection,by=AgeGroup))
    #bf(low_ct1 ~ LineageBroad_TiterGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroup)),
    #bf(low_ct1 ~ LineageBroad_TiterGroupAlt + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroupAlt))
    
    
)

names <- c("baseline","vaccine","cumulative_exposures","days_since_exposure",
           #"titer_group","titer_group_alt",
           "lineage",
           "vaccine_and_lineage","cumulative_exposures_and_lineage","days_since_exposure_and_lineage",
           #"titer_group_and_lineage","titer_group_alt_and_lineage",
           "baseline_age","vaccine_age","cumulative_exposures_age","days_since_exposure_age",
           "lineage_age",
           "vaccine_and_lineage_age","cumulative_exposures_and_lineage_age","days_since_exposure_and_lineage_age"
           )
           
key <- tibble(name=rep(names,2), formula=rep(seq(1, length(formulas),by=1),2), data=rep(1:2, each=length(formulas)))

use_formula <- unlist(formulas[key$formula[i]])
use_data <- key$data[i]
name <- key$name[i]

tmp_dat_frequent <- dat_subset_use %>% filter(DetectionSpeed=="Frequent testing")
tmp_dat_infrequent <- dat_subset_use %>% filter(DetectionSpeed != "Frequent testing")

if(use_data == 1) tmp_dat <- tmp_dat_frequent
if(use_data == 2) tmp_dat <- tmp_dat_infrequent

print(use_formula$formula)

## BASELINE. Just time since detection

if(rerun_stan){
    fit <- brm(data=tmp_dat, use_formula, family=bernoulli(link="logit"),cores=4,
               prior=c(prior_string("normal(0,10)",class="b"),prior_string("normal(0,10)",class="Intercept")),
               iter=n_iter,save_pars=save_pars(all=TRUE))
    save(fit, file=paste0("outputs/immune_models/",name,"_",use_data,".RData"))
} else {
    load(paste0("outputs/immune_models/",name,"_",use_data,".RData"))
}
## Assess performance
pred <- as.data.frame(predict(fit, type = "response"))
pred$pos <- as.numeric(pred$Estimate > 0.5)
dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% summarize(`Proportion correct`=sum(correct)/n())
dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% group_by(low_ct1) %>% 
    summarize(`Proportion correct`=sum(correct)/n()) %>% rename(`Ct<30`=low_ct1)
performance(prediction(pred$Estimate, pull(fit$data, low_ct1)),measure="auc")@y.values[[1]]

print(availableCores())
plan(multicore)
kfold_est <- kfold(fit, chains=1, K=25)
save(kfold_est, file=paste0("outputs/immune_models/",name,"_",use_data,"_kfolds",".RData"))

#loo_est <- loo(fit,reloo=TRUE,chains=1)
#save(loo_est, file=paste0(name,"_",use_data,"_loo_",i,".RData"))

