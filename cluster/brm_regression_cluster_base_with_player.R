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

setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")
#setwd("~/ct_nba")

n_iter <- 2000

load("data/data_for_regressions.RData")

dat_subset_use <- dat_subset_use %>% filter(Role != "Player")

filename_base <- paste0("outputs/")
if(!file.exists(filename_base)) dir.create(filename_base)

## Shouild be 8 entries -- 4 models, 2 datasets for each
i <- 5
#i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(i)

formulas <- list(
    ## Baseline
    bf(low_ct1 ~ s(DaysSinceDetection)),
    
    ## Single predictor
    bf(low_ct1 ~ VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus)),
    bf(low_ct1 ~ LineageBroad + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)),
    bf(low_ct1 ~ LineageBroad_VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus))
)

names <- c("baseline","vaccine","lineage","vaccine_and_lineage")

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
fit <- brm(data=tmp_dat, use_formula, family=bernoulli(link="logit"),cores=4,prior=c(set_prior("normal(0,10)",class="b")),
           iter=n_iter,save_pars=save_pars(all=TRUE))

save(fit, file=paste0("outputs/player_sensitivity/",name,"_",use_data,".RData"))

## Assess performance
pred <- as.data.frame(predict(fit, type = "response"))
pred$pos <- as.numeric(pred$Estimate > 0.5)
dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% summarize(`Proportion correct`=sum(correct)/n())
dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% group_by(low_ct1) %>% 
    summarize(`Proportion correct`=sum(correct)/n()) %>% rename(`Ct<30`=low_ct1)
performance(prediction(pred$Estimate, pull(fit$data, low_ct1)),measure="auc")@y.values[[1]]

