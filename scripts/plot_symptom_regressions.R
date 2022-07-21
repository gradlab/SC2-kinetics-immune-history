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

load("data/data_for_regressions.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0)

## Read in the two baseline regressions
load("outputs/base_models/symptoms_1.rdata")
sympt_freq <- fit
load("outputs/base_models/symptoms_2.rdata")
sympt_infreq <- fit

## Read in the two Lineage status regressions
load("outputs/base_models/symptoms_lineage_vacc_1.rdata")
symptlineage_freq <- fit
load("outputs/base_models/symptoms_lineage_vacc_2.rdata")
symptlineage_infreq <- fit
 

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
sympt_freq_res <- print_classification_accuracy(sympt_freq)
sympt_infreq_res <- print_classification_accuracy(sympt_infreq)
symptlineage_freq_res <- print_classification_accuracy(symptlineage_freq)
symptlineage_infreq_res <- print_classification_accuracy(symptlineage_infreq)

## Symptoms model
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1),Symptomatic=c("Yes","No"))
sympt_freq_draws <- sympt_freq %>% epred_draws(newdata=newdata)
sympt_infreq_draws <- sympt_infreq %>% epred_draws(newdata=newdata)

base_p_dat <- bind_rows(sympt_freq_draws %>% mutate(Protocol="Frequent testing"), 
                        sympt_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection,Symptomatic) %>% summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred))

p_base <- ggplot(base_p_dat, aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=Symptomatic,y=med),col="None") +
    geom_ribbon(alpha=0.5) +
    geom_line(aes(col=Symptomatic))+
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(breaks=seq(0,20,by=1)) +
    labs(y="Probability of Ct value <30",x="Days since detection") +
    facet_wrap(~Protocol) +
    theme_minimal() +
    theme(legend.position=c(0.7,0.7),
          plot.background = element_rect(fill="white",color="white"))


## Lineage-specific models
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), 
                       LineageBroad_Symptomatic_VaccStatus=unique(symptlineage_freq$data$LineageBroad_Symptomatic_VaccStatus))
lineage_freq_draws <- symptlineage_freq %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), 
                       LineageBroad_Symptomatic_VaccStatus=unique(symptlineage_infreq$data$LineageBroad_Symptomatic_VaccStatus))
lineage_infreq_draws <- symptlineage_infreq %>% epred_draws(newdata=newdata)

vacclineage_p_dat <- bind_rows(lineage_freq_draws %>% mutate(Protocol="Frequent testing"), 
                               lineage_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection,LineageBroad_Symptomatic_VaccStatus) %>% 
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred)) %>%
    mutate(LineageBroad = ifelse(LineageBroad_Symptomatic_VaccStatus %like% "Other","Other",
                                 ifelse(LineageBroad_Symptomatic_VaccStatus %like% "Delta","Delta","Omicron"))) %>%
    mutate(VaccStatus = ifelse(LineageBroad_Symptomatic_VaccStatus %like% "Boosted","Boosted",
                               ifelse(LineageBroad_Symptomatic_VaccStatus %like% "Second dose","Second dose",
                                      ifelse(LineageBroad_Symptomatic_VaccStatus %like% "First dose","First dose","Unvaccinated")))) %>%
    mutate(Symptomatic = ifelse(LineageBroad_Symptomatic_VaccStatus %like% "Yes","Yes",
                                ifelse(LineageBroad_Symptomatic_VaccStatus %like% "No","No","Unknown")))
vacclineage_p_dat$Protocol <- factor(vacclineage_p_dat$Protocol,levels=c("Frequent testing","Delayed detection"))

tmp_dat <- vacclineage_p_dat %>% filter(LineageBroad == "Omicron") %>% filter(VaccStatus %in% c("Boosted","Second dose")) %>%
    mutate(VaccStatus=ifelse(VaccStatus=="Boosted","Omicron: Boosted","Omicron: Not Boosted"))
tmp_dat$VaccStatus <- factor(tmp_dat$VaccStatus,levels=c("Omicron: Boosted","Omicron: Not Boosted"))

samp_sizes <- dat_subset_use %>% filter(LineageBroad == "Omicron") %>% 
    select(PersonID, VaccStatus,DetectionSpeed,Symptomatic) %>% distinct() %>% mutate(VaccStatus=ifelse(VaccStatus=="Boosted","Boosted","Not Boosted")) %>%
    group_by(VaccStatus,DetectionSpeed,Symptomatic) %>% tally() %>% mutate(label=paste0(VaccStatus,"; N=",n)) %>%
    mutate(Symptomatic = as.character(Symptomatic)) %>%
    mutate(Symptomatic = ifelse(is.na(Symptomatic),"Unknown",Symptomatic)) %>% 
    mutate(Symptomatic = paste0("Symptom status: ", Symptomatic)) %>%
    mutate(y=ifelse(VaccStatus == "Boosted",0.92,0.8)) %>%
    rename(Protocol=DetectionSpeed)

samp_sizes$Protocol <- factor(samp_sizes$Protocol, levels=c("Frequent testing","Delayed detection"))

p_vacclineage <-  ggplot(tmp_dat %>% mutate(Symptomatic = paste0("Symptom status: ", Symptomatic)),col="None") +
    facet_grid(Symptomatic~Protocol) +
    geom_ribbon(alpha=0.5, 
                aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=VaccStatus,y=med)) +
    geom_line(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=VaccStatus,y=med,col=VaccStatus))+
    geom_text(data=samp_sizes,aes(x=20,y=y,label=label),size=3) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Vaccination status",color="Vaccination status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_color_manual(values=immunecolors[3:4]) +
    scale_fill_manual(values=immunecolors[3:4])+
    theme(legend.position=c(0.85,0.85),
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white"))



ggsave(filename="figures/supplement/FigS8.png",plot=p_vacclineage,width=8,height=8,dpi=300)
ggsave(filename="figures/supplement/FigS8.pdf",plot=p_vacclineage,width=8,height=8)


