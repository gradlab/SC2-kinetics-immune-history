######################################################
## SCRIPT 2 -- analyze baseline regression results
######################################################

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
library(data.table)

lineage_colors <- c("Delta"="blue","Omicron"="red","Other"="orange","None"="grey70")
lineage_colors1 <- c("Delta"="blue","Omicron"="red","Other"="orange")
vacc_status_colors <- c("Boosted"="darkgreen","Second dose"="purple")
immunecolors <- c("Other: Unvaccinated"="black","Delta: 1-2 doses"="orange3","Omicron: Boosted"="purple3","Omicron: Not Boosted" = "mediumpurple1")

colors <- c("black","tomato","red3","dodgerblue","blue")
names(colors) <- c("Other","Delta1","Delta2","Omicron1","Omicron2")

setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")

load("~/Documents/GitHub/ct_nba/data/data_for_regressions_nonplayers.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0) %>% ungroup()

## Read in the two baseline regressions
load("outputs/immune_models_nonplayers/baseline_1.rdata")
basemodel_freq <- fit
load("outputs/immune_models_nonplayers/baseline_2.rdata")
basemodel_infreq <- fit

## Read in the two Lineage status 
## Read in the vaccine status and lineage regressions
load("outputs/immune_models_nonplayers/vaccine_and_lineage_age_1.rdata")
vacclineagemodel_freq <- fit
load("outputs/immune_models_nonplayers/vaccine_and_lineage_age_2.rdata")
vacclineagemodel_infreq <- fit


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

## Prediction accuracies
base_freq_res <- print_classification_accuracy(basemodel_freq)
base_infreq_res <- print_classification_accuracy(basemodel_infreq)
vacclineage_freq_res <- print_classification_accuracy(vacclineagemodel_freq)
vacclineage_infreq_res <- print_classification_accuracy(vacclineagemodel_infreq)

## Baseline model
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1))
base_freq_draws <- basemodel_freq %>% epred_draws(newdata=newdata)
base_infreq_draws <- basemodel_infreq %>% epred_draws(newdata=newdata)

base_p_dat <- bind_rows(base_freq_draws %>% mutate(Protocol="Frequent testing"), base_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection) %>% summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred))

p_base <- ggplot(base_p_dat, aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=Protocol,y=med),col="None") +
    geom_ribbon(alpha=0.5) +
    geom_line(aes(col=Protocol))+
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(breaks=seq(0,20,by=1)) +
    labs(y="Probability of Ct value <30",x="Days since detection") +
    theme_minimal() +
    theme(legend.position=c(0.7,0.7),
          plot.background = element_rect(fill="white",color="white"))


## Lineage-specific models

## Vaccine and lineage-specific models
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_VaccStatus=unique(vacclineagemodel_freq$data$LineageBroad_VaccStatus),
                       AgeGroup=unique(vacclineagemodel_freq$data$AgeGroup))
vacclineage_freq_draws <- vacclineagemodel_freq %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_VaccStatus=unique(vacclineagemodel_infreq$data$LineageBroad_VaccStatus),
                       AgeGroup=unique(vacclineagemodel_infreq$data$AgeGroup))
vacclineage_infreq_draws <- vacclineagemodel_infreq %>% epred_draws(newdata=newdata)

vacclineage_p_dat <- bind_rows(vacclineage_freq_draws %>% mutate(Protocol="Frequent testing"), 
                               vacclineage_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection,LineageBroad_VaccStatus,AgeGroup) %>% 
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred)) %>%
    mutate(LineageBroad = ifelse(LineageBroad_VaccStatus %like% "Other","Other",
                                 ifelse(LineageBroad_VaccStatus %like% "Delta","Delta","Omicron"))) %>%
    mutate(VaccStatus = ifelse(LineageBroad_VaccStatus %like% "Boosted","Boosted",
                               ifelse(LineageBroad_VaccStatus %like% "Second dose","Second dose",
                                      ifelse(LineageBroad_VaccStatus %like% "First dose","First dose","Unvaccinated"))))
vacclineage_p_dat$Protocol <- factor(vacclineage_p_dat$Protocol,levels=c("Frequent testing","Delayed detection"))

tmp_dat <- vacclineage_p_dat %>% filter(LineageBroad == "Omicron") %>% 
    filter(VaccStatus %in% c("Boosted","Second dose")) %>%
    filter(AgeGroup == "(30,50]") %>%
    mutate(VaccStatus=ifelse(VaccStatus=="Boosted","Omicron: Boosted","Omicron: Not Boosted"))
tmp_dat$VaccStatus <- factor(tmp_dat$VaccStatus,levels=c("Omicron: Boosted","Omicron: Not Boosted"))



p_vacclineage <-  ggplot(tmp_dat, 
                         aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=VaccStatus,y=med),col="None") +
    facet_wrap(~Protocol) +
    geom_ribbon(alpha=0.5) +
    geom_line(aes(col=VaccStatus))+
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Vaccination status",color="Vaccination status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_color_manual(values=immunecolors[3:4]) +
    scale_fill_manual(values=immunecolors[3:4])+
    theme(legend.position=c(0.85,0.75),
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white"))

all_colors <- c("Delta: Boosted"="orange", "Delta: Second dose"="yellow", "Delta: Unvaccinated"="tomato", 
                "Omicron: Boosted"="purple3", "Omicron: First dose"="lightblue", "Omicron: Second dose"="mediumpurple1", 
                "Omicron: Unvaccinated"="darkblue", "Other: First dose"="grey80", "Other: Second dose"="grey50", 
                "Other: Unvaccinated"="black")

samp_sizes <- dat_subset_use %>% ungroup() %>% select(PersonID, VaccStatus,LineageBroad,DetectionSpeed,AgeGroup) %>% distinct() %>% group_by(VaccStatus,LineageBroad,DetectionSpeed) %>% tally() %>% mutate(label=paste0("N=",n)) %>% 
    mutate(VaccStatusLineage = paste0(LineageBroad, ": ", VaccStatus)) %>%
    mutate(y = ifelse(VaccStatus == "Unvaccinated",0.5, 
                      ifelse(VaccStatus == "First dose", 0.6,
                             ifelse(VaccStatus == "Second dose", 0.7,
                                    0.8)))) %>%
    rename(Protocol=DetectionSpeed)
vacclineage_p_dat$Protocol <- factor(vacclineage_p_dat$Protocol, levels=c("Frequent testing","Delayed detection"))
samp_sizes$Protocol <- factor(samp_sizes$Protocol, levels=c("Frequent testing","Delayed detection"))
p_vacclineage_all <-  ggplot(vacclineage_p_dat %>% 
                                 filter(AgeGroup == "(0,30]") %>%
                                 mutate(
                                     VaccStatusLineage = paste0(LineageBroad, ": ", VaccStatus)),col="None") +
    facet_grid(LineageBroad~Protocol) +
    geom_ribbon(alpha=0.25,aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=VaccStatusLineage,y=med)) +
    geom_line(aes(col=VaccStatusLineage,x=DaysSinceDetection,ymin=lower,ymax=upper,fill=VaccStatusLineage,y=med))+
    geom_text(data=samp_sizes,aes(x=22, y=y,label=label,col=VaccStatusLineage),show.legend=FALSE) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Vaccination status",color="Vaccination status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    theme(legend.position="bottom",
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white")) +
    scale_color_manual(values=all_colors) +
    scale_fill_manual(values=all_colors) +
    guides(fill=guide_legend(nrow=3),color=guide_legend(nrow=3))

## Numbers for manuscript
## Baseline model
base_freq_res
base_infreq_res
base_p_dat %>% filter(DaysSinceDetection %in% c(5,10))


## Vaccinelineage interaction model
## Sample sizes
vacclineagemodel_freq$data %>% group_by(LineageBroad_VaccStatus) %>% tally()
vacclineage_freq_res
vacclineage_infreq_res

ggsave(filename="figures/supplement_player/baseline_regression.png",plot=p_base,width=6,height=3,dpi=300)
ggsave(filename="figures/supplement_player/vacclineage_regression.png",plot=p_vacclineage,width=8,height=3,dpi=300)
ggsave(filename="figures/supplement_player/vacclineage_regression_all.png",plot=p_vacclineage_all,width=8,height=8,dpi=300)

ggsave(filename="figures/supplement/figS7.png",plot=p_vacclineage_all,width=8,height=8,dpi=300)
ggsave(filename="figures/supplement/figS7.pdf",plot=p_vacclineage_all,width=8,height=8)
