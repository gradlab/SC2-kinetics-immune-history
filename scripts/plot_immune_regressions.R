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
setwd("~/Documents/GitHub/ct_nba")


## Assess performance
print_classification_accuracy <- function(fit){
    pred <- as.data.frame(predict(fit, type = "response"))
    pred$pos <- as.numeric(pred$Estimate > 0.5)
    overall_correct <- dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% summarize(`Proportion correct`=sum(correct)/n())
    correct_by_group <- dplyr::bind_cols(fit$data, pred) %>% mutate(correct= pos == low_ct1) %>% group_by(low_ct1) %>% 
        summarize(`Proportion correct`=sum(correct)/n()) %>% rename(`Ct<30`=low_ct1)
    auc <- performance(prediction(pred$Estimate, pull(fit$data, low_ct1)),measure="auc")@y.values[[1]]
    return(list(overall_correct, correct_by_group, auc))
}

load("data/data_for_regressions.RData")

filename_base <- paste0("outputs/immune_models")
if(!file.exists(filename_base)) dir.create(filename_base)

print(i)

names <- c("baseline","vaccine","cumulative_exposures","days_since_exposure","titer_group","titer_group_alt","lineage",
           "vaccine_and_lineage","cumulative_exposures_and_lineage","days_since_exposure_and_lineage",
           "titer_group_and_lineage","titer_group_alt_and_lineage")

key <- tibble(name=rep(names,2), formula=rep(seq(1, 12,by=1),2), data=rep(1:2, each=12))

formulas <- list(
    ## Baseline
    bf(low_ct1 ~ s(DaysSinceDetection)),
    
    ## Simple predictors
    bf(low_ct1 ~ VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=VaccStatus)),
    bf(low_ct1 ~ CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=CumulativeExposureNumber)),
    bf(low_ct1 ~ DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=DaysSinceExposureGroup)),
    bf(low_ct1 ~ TiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroup)),
    bf(low_ct1 ~ TiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=TiterGroupAlt)),
    bf(low_ct1 ~ LineageBroad + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad)),
    
    ## Add interactions with lineage
    bf(low_ct1 ~ LineageBroad_VaccStatus + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_VaccStatus)),
    bf(low_ct1 ~ LineageBroad_CumulativeExposureNumber + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_CumulativeExposureNumber)),
    bf(low_ct1 ~ LineageBroad_DaysSinceExposureGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_DaysSinceExposureGroup)),
    bf(low_ct1 ~ LineageBroad_TiterGroup + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroup)),
    bf(low_ct1 ~ LineageBroad_TiterGroupAlt + s(DaysSinceDetection) + s(DaysSinceDetection,by=LineageBroad_TiterGroupAlt))
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
        
        all_accuracies[[i]] <- print_classification_accuracy(get(filename))
        all_conditional_effects[[i]] <- conditional_effects(get(filename))
        load(paste0(filename,"_",use_data,"_kfolds.RData"))
        all_kfolds[[i]] <- kfold_est
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
                 "days_since_exposure"="Days since previous exposure", "lineage"="Lineage", "baseline"="Baseline") 
ests_freq_res$name <- formula_key[ests_freq_res$name]
ests_infreq_res$name <- formula_key[ests_infreq_res$name]

colnames_use <- c("Model","ELPD difference","SE difference","AUC","Classification accuracy","Accuracy (>=30)","Accuracy (<30)")
colnames(ests_freq_res) <- colnames_use
colnames(ests_infreq_res) <- colnames_use


freq_models <- c(1,2,3,4,7,8,9,10)
lpd_point_freq  <- cbind(
    all_kfolds[[1]]$pointwise[,"elpd_kfold"],
    all_kfolds[[2]]$pointwise[,"elpd_kfold"],
    all_kfolds[[3]]$pointwise[,"elpd_kfold"],
    all_kfolds[[4]]$pointwise[,"elpd_kfold"],
    all_kfolds[[7]]$pointwise[,"elpd_kfold"],
    all_kfolds[[8]]$pointwise[,"elpd_kfold"],
    all_kfolds[[9]]$pointwise[,"elpd_kfold"],
    all_kfolds[[10]]$pointwise[,"elpd_kfold"])

model_weights_freq <- loo::stacking_weights(lpd_point_freq)
#loo::pseudobma_weights(lpd_point)

infreq_models <- c(13,14,15,16,19,20,21,22)
lpd_point_infreq  <- cbind(
    all_kfolds[[13]]$pointwise[,"elpd_kfold"],
    all_kfolds[[14]]$pointwise[,"elpd_kfold"],
    all_kfolds[[15]]$pointwise[,"elpd_kfold"],
    all_kfolds[[16]]$pointwise[,"elpd_kfold"],
    all_kfolds[[19]]$pointwise[,"elpd_kfold"],
    all_kfolds[[20]]$pointwise[,"elpd_kfold"],
    all_kfolds[[21]]$pointwise[,"elpd_kfold"],
    all_kfolds[[22]]$pointwise[,"elpd_kfold"])


model_weights_freq <- loo::stacking_weights(lpd_point_freq)
model_weights_freq
model_weights_freq_res <- bind_cols(key$name[freq_models], model_weights_freq)
colnames(model_weights_freq_res) <- c("Model","Weight")
model_weights_freq_res$Weight <- signif(model_weights_freq_res$Weight, 3)


model_weights_infreq <- loo::stacking_weights(lpd_point_infreq)
model_weights_infreq
model_weights_infreq_res <- bind_cols(key$name[infreq_models], model_weights_infreq)
colnames(model_weights_infreq_res) <- c("Model","Weight")
model_weights_infreq_res$Weight <- signif(model_weights_infreq_res$Weight, 3)


model_weights_freq_res$Model <- formula_key[model_weights_freq_res$Model]
model_weights_infreq_res$Model <- formula_key[model_weights_infreq_res$Model]

write_csv(ests_freq_res, file="../../figures/immune_models_freq_compare.csv")
write_csv(ests_infreq_res, file="../../figures/immune_models_infreq_compare.csv")
write_csv(model_weights_freq_res %>% arrange(-Weight), file="../../figures/immune_models_freq_weights.csv")
write_csv(model_weights_infreq_res %>% arrange(-Weight), file="../../figures/immune_models_infreq_weights.csv")


## Best model -- cumulative exposures
tmp <- all_conditional_effects[[which(key$name == "cumulative_exposures_and_lineage")[1]]][[3]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_CumulativeExposureNumber %like% "Delta","Delta",
                             ifelse(LineageBroad_CumulativeExposureNumber %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(InfectionNumber = ifelse(LineageBroad_CumulativeExposureNumber %like% "1",1,
                                               ifelse(LineageBroad_CumulativeExposureNumber %like% "2",2,
                                                      ifelse(LineageBroad_CumulativeExposureNumber %like% "3",3,
                                                             ifelse(LineageBroad_CumulativeExposureNumber %like% "4",4,
                                                                    ifelse(LineageBroad_CumulativeExposureNumber %like% "5",5,6  )))))) %>%
    mutate(InfectionNumber=as.factor(InfectionNumber))
p_best_model_freq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=InfectionNumber),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=InfectionNumber)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    ggtitle("Frequent testing") +
    theme_minimal() +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)

tmp <- all_conditional_effects[[which(key$name == "cumulative_exposures_and_lineage")[2]]][[3]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_CumulativeExposureNumber %like% "Delta","Delta",
                                            ifelse(LineageBroad_CumulativeExposureNumber %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(InfectionNumber = ifelse(LineageBroad_CumulativeExposureNumber %like% "1",1,
                                               ifelse(LineageBroad_CumulativeExposureNumber %like% "2",2,
                                                      ifelse(LineageBroad_CumulativeExposureNumber %like% "3",3,
                                                             ifelse(LineageBroad_CumulativeExposureNumber %like% "4",4,
                                                                    ifelse(LineageBroad_CumulativeExposureNumber %like% "5",5,6  )))))) %>%
    mutate(InfectionNumber=as.factor(InfectionNumber))
p_best_model_infreq <- ggplot(tmp %>% filter(InfectionNumber != 6)) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=InfectionNumber),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=InfectionNumber)) +
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
tmp <- all_conditional_effects[[which(key$name == "vaccine_and_lineage")[1]]][[3]]
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

tmp <- all_conditional_effects[[which(key$name == "vaccine_and_lineage")[2]]][[3]]
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

ggsave(p_vacc_model_freq + p_vacc_model_infreq,filename = "figures/immune_models/vacc_lineage_fits.png",width=10,height=10,dpi=300)
ggsave(p_best_model_freq + p_best_model_infreq,filename = "figures/immune_models/cumu_exposures_fits.png",width=10,height=10,dpi=300)



## Antibody titer model
tmp <- all_conditional_effects[[which(key$name == "titer_group_and_lineage")[1]]][[3]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_TiterGroup %like% "Delta","Delta",
                                            ifelse(LineageBroad_TiterGroup %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(TiterGroup = ifelse(LineageBroad_TiterGroup %like% "-1,125","(-1,125]",
                                          ifelse(LineageBroad_TiterGroup %like% "125,250","(125,250]",
                                                 ifelse(LineageBroad_TiterGroup %like% "250,","(250,1e+03]",NA))))

p_titer_model_freq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=TiterGroup),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=TiterGroup)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    ggtitle("Frequent testing") +
    theme_minimal() +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)

tmp <- all_conditional_effects[[which(key$name == "titer_group_and_lineage")[2]]][[3]]
tmp <- tmp %>% mutate(LineageBroad = ifelse(LineageBroad_TiterGroup %like% "Delta","Delta",
                                            ifelse(LineageBroad_TiterGroup %like% "Omicron","Omicron","Other")))
tmp <- tmp %>% mutate(TiterGroup = ifelse(LineageBroad_TiterGroup %like% "-1,125","(-1,125]",
                                          ifelse(LineageBroad_TiterGroup %like% "125,250","(125,250]",
                                                 ifelse(LineageBroad_TiterGroup %like% "250,","(250,1e+03]",NA))))


p_titer_model_infreq <- ggplot(tmp) + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=TiterGroup),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=TiterGroup)) +
    scale_fill_viridis_d() + scale_color_viridis_d() +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    ggtitle("Delayed detection") +
    theme_minimal() +
    coord_cartesian(xlim=c(0,20)) +
    theme(legend.position="bottom",
          plot.background=element_rect(fill="white",color="white"))+
    facet_wrap(~LineageBroad,ncol=1)


ggsave(p_titer_model_freq + p_titer_model_infreq,filename = "figures/immune_models/titer_lineage_fits.png",width=10,height=10,dpi=300)
ggsave(p_vacc_model_freq + p_vacc_model_infreq,filename = "figures/immune_models/vacc_lineage_fits.png",width=10,height=10,dpi=300)
ggsave(p_best_model_freq + p_best_model_infreq,filename = "figures/immune_models/cumu_exposures_fits.png",width=10,height=10,dpi=300)

