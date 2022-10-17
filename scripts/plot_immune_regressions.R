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
setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")
savewd <- getwd()

## Assess performance
print_classification_accuracy_old <- function(fit){
    pred <- as.data.frame(predict(fit, type = "response"))
    pred$pos <- as.numeric(pred$Estimate > 0.5)
    overall_correct <- dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% summarize(`Proportion correct`=sum(correct)/n())
    correct_by_group <- dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% group_by(low_ct1) %>% 
        summarize(`Proportion correct`=sum(correct)/n()) %>% rename(`Ct<30`=low_ct1)
    auc <- performance(prediction(pred$Estimate, pull(fit$data, low_ct1)),measure="auc")@y.values[[1]]
    return(list(overall_correct, correct_by_group, auc))
}

print_classification_accuracy <- function(fit){
    pred <- kfold_predict(fit, method = "predict")
    pos <- as.numeric(colMeans(pred$yrep) > 0.5)
    estimated <- dplyr::bind_cols(fit$data, pos=pos) %>% mutate(correct= pos == low_ct1)
    
    overall_correct <- estimated %>% summarize(`Proportion correct`=sum(correct)/n())
    correct_by_group <- estimated %>% group_by(low_ct1) %>% 
        summarize(`Proportion correct`=sum(correct)/n()) %>% rename(`Ct<30`=low_ct1)
    auc <- performance(prediction(colMeans(pred$yrep), fit$data$low_ct1),measure="auc")@y.values[[1]]
    return(list(overall_correct, correct_by_group, auc))
}

load("data/data_for_regressions.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0)
filename_base <- paste0("outputs/immune_models")
if(!file.exists(filename_base)) dir.create(filename_base)
i <- 3

names <- c("baseline","vaccine","cumulative_exposures","days_since_exposure",
           "lineage",
           "vaccine_and_lineage","cumulative_exposures_and_lineage","days_since_exposure_and_lineage",
           "baseline_age","vaccine_age","cumulative_exposures_age","days_since_exposure_age",
           "lineage_age",
           "vaccine_and_lineage_age","cumulative_exposures_and_lineage_age","days_since_exposure_and_lineage_age"
)

key <- tibble(name=rep(names,2), formula=rep(seq(1, length(names),by=1),2), data=rep(1:2, each=length(names)))

formulas <- list(
    ## Baseline
    bf(low_ct1 ~ s(DaysSinceDetection)),
    
    ## Simple predictors
    bf(low_ct1 ~ VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus)),
    bf(low_ct1 ~ CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=CumulativeExposureNumber)),
    bf(low_ct1 ~ DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=DaysSinceExposureGroup)),
    bf(low_ct1 ~ LineageBroad + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)),
    
    ## Add interactions with lineage
    bf(low_ct1 ~ LineageBroad_VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus)),
    bf(low_ct1 ~ LineageBroad_CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_CumulativeExposureNumber)),
    bf(low_ct1 ~ LineageBroad_DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_DaysSinceExposureGroup)),
    
    ## All of the above but with age
    ## Simple predictors
    bf(low_ct1 ~ AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ VaccStatus + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus) + s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ CumulativeExposureNumber + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=CumulativeExposureNumber)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ DaysSinceExposureGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=DaysSinceExposureGroup)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ LineageBroad + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)+ s(DaysSinceDetection,by=AgeGroup)),
    
    ## Add interactions with lineage
    bf(low_ct1 ~ LineageBroad_VaccStatus + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ LineageBroad_CumulativeExposureNumber + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_CumulativeExposureNumber)+ s(DaysSinceDetection,by=AgeGroup)),
    bf(low_ct1 ~ LineageBroad_DaysSinceExposureGroup + AgeGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_DaysSinceExposureGroup)+ s(DaysSinceDetection,by=AgeGroup))
)

use_formula <- unlist(formulas[key$formula[i]])
use_data <- key$data[i]
name <- key$name[i]

all_fits <- NULL
all_kfolds <- NULL
all_accuracies <- NULL
all_conditional_effects <- NULL

setwd("outputs/immune_models/")
for(i in 1:nrow(key)){
    print(i)
    
    try({
        filename <- key$name[i]
        use_data <- key$data[i]
        load(paste0(filename,"_",use_data,".RData"))
        
        gdata::mv("fit",filename)
        
        all_fits[[i]] <- get(filename)
        
       #all_accuracies[[i]] <- print_classification_accuracy(get(filename))
        all_conditional_effects[[i]] <- conditional_effects(get(filename))
        load(paste0(filename,"_",use_data,"_kfolds.RData"))
        all_kfolds[[i]] <- kfold_est
        all_accuracies[[i]] <- print_classification_accuracy(all_kfolds[[i]])
    })    
}



key$AUC <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[3]])){return(NA)
        } else {return(x[[3]])}})
)
key$classification_accuracy <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[1]]) }})
)
key$classification_accuracy_over30 <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[2]][1,2]) }})
)
key$classification_accuracy_under30 <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[2]][2,2]) }})
)

key$elpd <- unlist(
    lapply(all_kfolds, function(x){
        if(is.null(x)){ return(NA) } else {return(x$estimates[1,1]) }})
)
key$kfoldsic <- unlist(
    lapply(all_kfolds, function(x){
        if(is.null(x)){ return(NA) } else {return(x$estimates[3,1]) }})
)
key %>% arrange(data, -elpd) %>% filter(!(name %like% "titer"))
key %>% arrange(data, -AUC) %>% filter(!(name %like% "titer"))
key %>% arrange(data, -classification_accuracy_under30) %>% filter(!(name %like% "titer"))

ests_freq <- loo_compare(all_kfolds[which(key$data == 1 & !(key$name %like% "titer"))])
ests_freq <- bind_cols(ests_freq, name= key %>% arrange(data, -elpd) %>% filter(!(name %like% "titer"),data==1) %>% pull(name)) 
ests_freq_res <- ests_freq %>% 
    select(name, elpd_diff, se_diff) %>% 
    mutate(data = 1) %>%
    left_join(key) %>% 
    select(-c("data","formula","elpd"))

ests_infreq <- loo_compare(all_kfolds[which(key$data == 2 & !(key$name %like% "titer"))])
ests_infreq <- bind_cols(ests_infreq, name= key %>% arrange(data, -elpd) %>% filter(!(name %like% "titer"),data==2) %>% pull(name))
ests_infreq_res <- ests_infreq %>% 
    select(name, elpd_diff, se_diff) %>% 
    mutate(data = 2) %>%
    left_join(key) %>% 
    select(-c("data","formula","elpd"))

formula_key <- c("cumulative_exposures_and_lineage"="Cumulative number of exposures and lineage",
                 "vaccine_and_lineage"="Vaccination status and lineage", 
                 "cumulative_exposures"="Cumulative number of exposures", 
                 "days_since_exposure_and_lineage"="Days since previous exposure and lineage", 
                 "vaccine"="Vaccination status", 
                 "days_since_exposure"="Days since previous exposure", 
                 "lineage"="Lineage", 
                 "baseline"="Baseline",
                 
                 "cumulative_exposures_and_lineage_age"="Cumulative number of exposures, lineage and age",
                 "vaccine_and_lineage_age"="Vaccination status, lineage and age", 
                 "cumulative_exposures_age"="Cumulative number of exposures and age", 
                 "days_since_exposure_and_lineage_age"="Days since previous exposure, lineage and age", 
                 "vaccine_age"="Vaccination status and age", 
                 "days_since_exposure_age"="Days since previous exposure and age", 
                 "lineage_age"="Lineage and age", 
                 "baseline_age"="Age"
                 ) 
ests_freq_res$name <- formula_key[ests_freq_res$name]
ests_infreq_res$name <- formula_key[ests_infreq_res$name]

colnames_use <- c("Model","ELPD difference","SE difference","AUC","Classification accuracy","Accuracy (>=30)","Accuracy (<30)")
colnames(ests_freq_res) <- colnames_use
colnames(ests_infreq_res) <- colnames_use


lpd_point_freq  <- cbind(sapply(1:16, function(x) all_kfolds[[x]]$pointwise[,"elpd_kfold"]))

model_weights_freq <- loo::stacking_weights(lpd_point_freq)
#loo::pseudobma_weights(lpd_point)
lpd_point_infreq  <- cbind(sapply(17:32, function(x) all_kfolds[[x]]$pointwise[,"elpd_kfold"]))

model_weights_freq <- loo::stacking_weights(lpd_point_freq)
model_weights_freq
model_weights_freq_res <- bind_cols(key$name[1:16], model_weights_freq)
colnames(model_weights_freq_res) <- c("Model","Weight")
model_weights_freq_res$Weight <- signif(model_weights_freq_res$Weight, 3)


model_weights_infreq <- loo::stacking_weights(lpd_point_infreq)
model_weights_infreq
model_weights_infreq_res <- bind_cols(key$name[17:32], model_weights_infreq)
colnames(model_weights_infreq_res) <- c("Model","Weight")
model_weights_infreq_res$Weight <- signif(model_weights_infreq_res$Weight, 3)


model_weights_freq_res$Model <- formula_key[model_weights_freq_res$Model]
model_weights_infreq_res$Model <- formula_key[model_weights_infreq_res$Model]

ests_freq_res1 <- ests_freq_res %>% left_join(model_weights_freq_res)
ests_freq_res1 <- ests_freq_res1[,c(1,2,3,9,4,5,6,7)]

ests_infreq_res1 <- ests_infreq_res %>% left_join(model_weights_infreq_res)
ests_infreq_res1 <- ests_infreq_res1[,c(1,2,3,9,4,5,6,7)]

write_csv(ests_freq_res1, file="../../figures/immune_models_freq_compare.csv")
write_csv(ests_infreq_res1, file="../../figures/immune_models_infreq_compare.csv")


## Best model -- cumulative exposures and age
tmp <- all_conditional_effects[[which(key$name == "cumulative_exposures_and_lineage_age")[1]]][[4]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_CumulativeExposureNumber %like% "Delta","Delta",
                             ifelse(LineageBroad_CumulativeExposureNumber %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(ExposureNumber = ifelse(LineageBroad_CumulativeExposureNumber %like% "1",1,
                                               ifelse(LineageBroad_CumulativeExposureNumber %like% "2",2,
                                                      ifelse(LineageBroad_CumulativeExposureNumber %like% "3",3,
                                                             ifelse(LineageBroad_CumulativeExposureNumber %like% "4",4,
                                                                    ifelse(LineageBroad_CumulativeExposureNumber %like% "5",5,6  )))))) %>%
    mutate(ExposureNumber=as.factor(ExposureNumber))
p_best_model_freq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=ExposureNumber),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=ExposureNumber)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    ggtitle("Frequent testing") +
    theme_minimal() +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)

tmp <- all_conditional_effects[[which(key$name == "cumulative_exposures_and_lineage_age")[2]]][[4]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_CumulativeExposureNumber %like% "Delta","Delta",
                                            ifelse(LineageBroad_CumulativeExposureNumber %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(ExposureNumber = ifelse(LineageBroad_CumulativeExposureNumber %like% "1",1,
                                               ifelse(LineageBroad_CumulativeExposureNumber %like% "2",2,
                                                      ifelse(LineageBroad_CumulativeExposureNumber %like% "3",3,
                                                             ifelse(LineageBroad_CumulativeExposureNumber %like% "4",4,
                                                                    ifelse(LineageBroad_CumulativeExposureNumber %like% "5",5,6  )))))) %>%
    mutate(ExposureNumber=as.factor(ExposureNumber))
p_best_model_infreq <- ggplot(tmp %>% filter(ExposureNumber != 6)) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=ExposureNumber),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=ExposureNumber)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    theme_minimal() +
    ggtitle("Delayed detection") +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)


## Second best model -- vaccination status
tmp <- all_conditional_effects[[which(key$name == "vaccine_and_lineage_age")[1]]][[4]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_VaccStatus %like% "Delta","Delta",
                                            ifelse(LineageBroad_VaccStatus %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(VaccStatus = ifelse(LineageBroad_VaccStatus %like% "Boosted","Boosted",
                                            ifelse(LineageBroad_VaccStatus %like% "Second dose","Second dose",
                                                   ifelse(LineageBroad_VaccStatus %like% "First dose","First dose",
                                                          ifelse(LineageBroad_VaccStatus %like% "Unvaccinated","Unvaccinated",NA)))))

p_vacc_model_freq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=VaccStatus),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=VaccStatus)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    ggtitle("Frequent testing") +
    theme_minimal() +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)

tmp <- all_conditional_effects[[which(key$name == "vaccine_and_lineage_age")[2]]][[4]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_VaccStatus %like% "Delta","Delta",
                                            ifelse(LineageBroad_VaccStatus %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(VaccStatus = ifelse(LineageBroad_VaccStatus %like% "Boosted","Boosted",
                                          ifelse(LineageBroad_VaccStatus %like% "Second dose","Second dose",
                                                 ifelse(LineageBroad_VaccStatus %like% "First dose","First dose",
                                                        ifelse(LineageBroad_VaccStatus %like% "Unvaccinated","Unvaccinated",NA)))))

p_vacc_model_infreq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=VaccStatus),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=VaccStatus)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    theme_minimal() +
    ggtitle("Delayed detection") +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)

ggsave(p_vacc_model_freq + p_vacc_model_infreq,filename = paste0(savewd,"/figures/immune_models/vacc_lineage_fits.png"),width=10,height=10,dpi=300)
ggsave(p_best_model_freq + p_best_model_infreq,filename = paste0(savewd,"/figures/immune_models/cumu_exposures_fits.png"),width=10,height=10,dpi=300)


