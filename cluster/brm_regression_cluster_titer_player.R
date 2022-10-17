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
setwd("~/SC2-kinetics-immune-history/")

n_iter <- 2000
rerun_stan <- TRUE
load("~/ct_data/data/data_for_regressions_player_new.RData")


## For these analyses, we only want to use Ct values after detection
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0) %>% ungroup()

## Not players
dat_subset_use <- dat_subset_use %>% filter(Role != "Player")

filename_base <- paste0("outputs/titer_models_nonplayers")
if(!file.exists(filename_base)) dir.create(filename_base)

## 48 options
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(i)

formulas <- list(
    bf(low_ct1 ~ BoostTiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=BoostTiterGroup)),
    bf(low_ct1 ~ BoostTiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=BoostTiterGroupAlt)),
    bf(low_ct1 ~ LineageBroad_BoostTiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_BoostTiterGroup)),
    bf(low_ct1 ~ LineageBroad_BoostTiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_BoostTiterGroupAlt)),
    
    bf(low_ct1 ~ LineageBroad*BoostTiterGroup + 
           s(DaysSinceDetection) + 
           s(DaysSinceDetection,by=BoostTiterGroup) + 
           s(DaysSinceDetection,by=LineageBroad) +
           s(DaysSinceDetection,by=interaction(LineageBroad,BoostTiterGroup))),
    
    bf(low_ct1 ~ LineageBroad*BoostTiterGroupAlt + 
           s(DaysSinceDetection) + 
           s(DaysSinceDetection,by=BoostTiterGroupAlt) + 
           s(DaysSinceDetection,by=LineageBroad) +
           s(DaysSinceDetection,by=interaction(LineageBroad,BoostTiterGroupAlt)))
)

names <- expand_grid(time=c("all","under60","under90","60to90"),freq=c("freq","infreq"),model=seq_along(formulas))
names <- names %>% mutate(name=paste(time,freq,model,sep="_"))

filename <- names$name[i]
use_formula <- unlist(formulas[names$model[i]])
use_data <- names$freq[i]
use_timerange <- names$time[i]

if(use_timerange == "under60"){
    dat_subset_use <- dat_subset_use %>% filter(UseLessThan60 == TRUE)
} else if(use_timerange == "under90"){
    dat_subset_use <- dat_subset_use %>% filter(UseLessThan90 == TRUE)
} else if(use_timerange == "60to90"){
    dat_subset_use <- dat_subset_use %>% filter(Use60to90 == TRUE)
} else {
    dat_subset_use <- dat_subset_use
}

tmp_dat_frequent <- dat_subset_use %>% filter(DetectionSpeed=="Frequent testing")
tmp_dat_infrequent <- dat_subset_use %>% filter(DetectionSpeed != "Frequent testing")

if(use_data == "freq") tmp_dat <- tmp_dat_frequent
if(use_data == "infreq") tmp_dat <- tmp_dat_infrequent

print(use_formula$formula)

if(rerun_stan){
## BASELINE. Just time since detection
    fit <- brm(data=tmp_dat, use_formula, family=bernoulli(link="logit"),cores=4,
               prior=c(prior_string("normal(0,10)",class="b"),prior_string("normal(0,10)",class="Intercept")),
               iter=n_iter,save_pars=save_pars(all=TRUE))

    save(fit, file=paste0("outputs/titer_models_nonplayers/",filename,".RData"))
} else {
    load(paste0("outputs/titer_models_nonplayers/",filename,".RData"))
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
save(kfold_est, file=paste0("outputs/titer_models_nonplayers/",filename,"_kfolds",".RData"))

#loo_est <- loo(fit,reloo=TRUE,chains=1)
#save(loo_est, file=paste0(name,"_",use_data,"_loo_",i,".RData"))

