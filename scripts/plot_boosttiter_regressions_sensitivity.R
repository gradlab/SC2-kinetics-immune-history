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

filename_base <- paste0("outputs/titer_models_sensitivity")
if(!file.exists(filename_base)) dir.create(filename_base)

## 48 options
#i <- 3
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
           s(DaysSinceDetection,by=LineageBroad_BoostTiterGroup)),
    
    bf(low_ct1 ~ LineageBroad*BoostTiterGroupAlt + 
           s(DaysSinceDetection) + 
           s(DaysSinceDetection,by=BoostTiterGroupAlt) + 
           s(DaysSinceDetection,by=LineageBroad) +
           s(DaysSinceDetection,by=LineageBroad_BoostTiterGroupAlt))
)

names <- expand_grid(time=c("all","under60","under90","60to90"),freq=c("freq","infreq"),model=seq_along(formulas))
names <- names %>% mutate(name=paste(time,freq,model,sep="_"))

all_fits <- NULL
all_accuracies <- NULL
all_conditional_effects <- NULL

newdata_all <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1),
                           BoostTiterGroup=unique(dat_subset_use$BoostTiterGroup),LineageBroad=unique(dat_subset_use$LineageBroad))
newdata_all$LineageBroad_BoostTiterGroup <- paste0(newdata_all$LineageBroad,"_",newdata_all$BoostTiterGroup)


names_tmp <- names

names <- names[c(3, 9),]
#names <- names[c(39, 45),]



setwd("outputs/titer_models_sensitivity//")
for(i in 1:nrow(names)){
    print(i)
    
    try({
        filename <- names$name[i]
        use_formula <- unlist(formulas[names$model[i]])
        use_data <- names$freq[i]
        use_timerange <- names$time[i]
        load(paste0(filename,".RData"))
        
        all_fits[[i]] <- fit
    
        all_accuracies[[i]] <- print_classification_accuracy(fit)
        all_conditional_effects[[i]] <- conditional_effects(fit)
    })
}

names$AUC <- unlist(
    lapply(all_accuracies, function(x){
    if(is.null(x[[3]])){return(NA)
    } else {return(x[[3]])}})
    )
names$classification_accuracy <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[1]]) }})
)
names$classification_accuracy_over30 <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[2]][1,2]) }})
)
names$classification_accuracy_under30 <- unlist(
    lapply(all_accuracies, function(x){
        if(is.null(x[[1]])){ return(NA) } else {return(x[[2]][2,2]) }})
)


convert_to_lineage_and_titer <- function(dat){
    library(data.table)
    
    if("LineageBroad_BoostTiterGroup" %in% colnames(dat)){
        dat <- dat %>% mutate(LineageBroad=ifelse(LineageBroad_BoostTiterGroup %like% "Omicron","Omicron",
                                                  ifelse(LineageBroad_BoostTiterGroup %like% "Delta","Delta",
                                                         ifelse(LineageBroad_BoostTiterGroup %like% "Other","Other",NA))))
        dat <- dat %>% mutate(BoostTiterGroup=ifelse(LineageBroad_BoostTiterGroup %like% "LowNotBoosted","LowNotBoosted",
                                     ifelse(LineageBroad_BoostTiterGroup %like% "HighNotBoosted","HighNotBoosted",
                                            ifelse(LineageBroad_BoostTiterGroup %like% "LowBoosted","LowBoosted",
                                                   ifelse(LineageBroad_BoostTiterGroup %like% "HighBoosted","HighBoosted",NA
                                                   )))))
    } else {
        dat <- dat %>% mutate(LineageBroad=ifelse(LineageBroad_BoostTiterGroupAlt %like% "Omicron","Omicron",
                                                  ifelse(LineageBroad_BoostTiterGroupAlt %like% "Delta","Delta",
                                                         ifelse(LineageBroad_BoostTiterGroupAlt %like% "Other","Other",NA))))
        dat <- dat %>% mutate(BoostTiterGroupAlt=ifelse(LineageBroad_BoostTiterGroupAlt %like% "LowNotBoosted","LowNotBoosted",
                                                     ifelse(LineageBroad_BoostTiterGroupAlt %like% "HighNotBoosted","HighNotBoosted",
                                                            ifelse(LineageBroad_BoostTiterGroupAlt %like% "LowBoosted","LowBoosted",
                                                                   ifelse(LineageBroad_BoostTiterGroupAlt %like% "HighBoosted","HighBoosted",NA
                                                                   )))))
    }
    dat
}

plot_models <- which(names$model %in% c(3,4))
all_plots <- NULL
setwd("~/Documents/GitHub/ct_nba/")

all_dat <- NULL

for(i in seq_along(plot_models)){
    use_i <- plot_models[i]
    names[use_i,]
    print(use_i)
    try({
    use_data <- names$freq[use_i]
    use_timerange <- names$time[use_i]
    use_model <- names$model[use_i]
    filename <- names$name[use_i]
    
    use_data <- ifelse(use_data == "freq","Frequent testing","Delayed detection")
    use_timerange <- ifelse(use_timerange == "all","All infections",
                            ifelse(use_timerange=="under60","Within 60 days of draw",
                             ifelse(use_timerange=="under90","Within 90 days of draw","60 to 90 days after draw")))
    
    use_model <- ifelse(use_model == 3,"250 threshold","400 threshold")
    
    title <- paste(use_data, use_timerange,use_model,sep="; ")
    if(use_model == "250 threshold"){
        tmp_dat <- convert_to_lineage_and_titer(all_conditional_effects[[use_i]]$`DaysSinceDetection:LineageBroad_BoostTiterGroup`)
        all_plots[[i]] <- tmp_dat %>% drop_na() %>% ggplot() + 
            geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=BoostTiterGroup),alpha=0.25) +
            geom_line(aes(x=DaysSinceDetection,y=estimate__,col=BoostTiterGroup)) +
            scale_fill_viridis_d() + scale_color_viridis_d() +
            ggtitle(title) +
            ylab("Probability of Ct value <30") +
            xlab("Days since detection") +
            theme_minimal() +
            coord_cartesian(xlim=c(0,20)) +
            theme(legend.position="bottom",
                  plot.background=element_rect(fill="white",color="white"))+
            facet_wrap(~LineageBroad,ncol=1)
    } else {
        tmp_dat <- convert_to_lineage_and_titer(all_conditional_effects[[use_i]]$`DaysSinceDetection:LineageBroad_BoostTiterGroupAlt`)
        all_plots[[i]] <- tmp_dat %>% drop_na() %>% ggplot() + 
            geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=BoostTiterGroupAlt),alpha=0.25) +
            geom_line(aes(x=DaysSinceDetection,y=estimate__,col=BoostTiterGroupAlt)) +
            scale_fill_viridis_d() + scale_color_viridis_d() +
            ggtitle(title) +
            ylab("Probability of Ct value <30") +
            xlab("Days since detection") +
            theme_minimal() +
            coord_cartesian(xlim=c(0,20)) +
            theme(legend.position="bottom",
                  plot.background=element_rect(fill="white",color="white"))+
            facet_wrap(~LineageBroad,ncol=1)
        
    }
    all_dat[[i]] <- tmp_dat %>% mutate(`Detection group`=use_data)
    ggsave(filename=paste0(getwd(),"/figures/boosttiters_sensitivity/",filename,".png"),all_plots[[i]],width=8,height=10,dpi=300)
    })
    
}
all_dat <- do.call("bind_rows",all_dat)
all_dat <- all_dat %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))
all_dat$`Detection group` <- factor(all_dat$`Detection group`, levels=c("Frequent testing","Delayed detection"))

dat_raw <- dat_subset_use %>% group_by(BoostTiterGroup,LineageBroad, DetectionSpeed,DaysSinceDetection) %>% 
    summarize(prop_low = sum(low_ct1)/n()) %>% rename(`Detection group`=DetectionSpeed)


dat_raw <- dat_raw %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))
dat_raw$`Detection group` <- factor(dat_raw$`Detection group`, levels=c("Frequent testing","Delayed detection"))



figS10 <- all_dat %>% ggplot() + 
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=BoostTiterGroup),alpha=0.25) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=BoostTiterGroup)) +
    #geom_line(data=dat_raw ,
    #          aes(x=DaysSinceDetection, y=prop_low,col=BoostTiterGroup)) +
    scale_fill_viridis_d(name="Immune status") + scale_color_viridis_d(name="Immune status") +
    ylab("Probability of Ct value <30") +
    xlab("Days since detection") +
    theme_classic() +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
    theme( legend.position="bottom",
           strip.background=element_blank(),
           panel.spacing=unit(2,"lines"),
                 strip.text=element_text(face="bold"),
                 axis.text=element_text(size=8),axis.title=element_text(size=8),
                 legend.title=element_text(size=8),legend.text=element_text(size=8),
                 plot.background = element_rect(fill="white",color="white")) +
    facet_grid(LineageBroad~`Detection group`)

write_csv(all_dat, "figures/supplement_titer/figS13_data.csv")

    #ggsave(filename="figures/supplement_titer/figS10.png",plot=figS10,width=8,height=7,dpi=300)
    #ggsave(filename="figures/supplement_titer/figS10.pdf",plot=figS10,width=8,height=7)
