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

load("data/data_for_regressions.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0) %>% ungroup()
dat_subset_use <- dat_subset_use %>% filter(UseLessThan90 == TRUE)

## Load 60-90 day filter sensitivity. Note that this fit does not include age group
load("outputs/titer_models/60to90_freq_2.RData")
freq_60to90 <- fit
load("outputs/titer_models/60to90_infreq_2.RData")
infreq_60to90 <- fit

## Get posterior draws
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(freq_60to90$data$LineageBroad_BoostTiterGroup))
freq_60to90_draws <- freq_60to90 %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(infreq_60to90$data$LineageBroad_BoostTiterGroup))
infreq_60to90_draws <- infreq_60to90 %>% epred_draws(newdata=newdata)

## Get posterior summaries
freq_60to90_res <- freq_60to90_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred)) %>% 
    mutate(Protocol="Frequent testing")
infreq_60to90_res <- infreq_60to90_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred))%>% 
    mutate(Protocol="Delayed detection")

res_60to90 <- bind_rows(freq_60to90_res, infreq_60to90_res)

convert_to_lineage_and_titer <- function(dat){
    library(data.table)
    
    dat <- dat %>% mutate(LineageBroad=ifelse(LineageBroad_BoostTiterGroup %like% "Omicron","Omicron",
                                              ifelse(LineageBroad_BoostTiterGroup %like% "Delta","Delta",
                                                     ifelse(LineageBroad_BoostTiterGroup %like% "Other","Other",NA))))
    dat <- dat %>% mutate(BoostTiterGroup=ifelse(LineageBroad_BoostTiterGroup %like% "LowNotBoosted","LowNotBoosted",
                                                 ifelse(LineageBroad_BoostTiterGroup %like% "HighNotBoosted","HighNotBoosted",
                                                        ifelse(LineageBroad_BoostTiterGroup %like% "LowBoosted","LowBoosted",
                                                               ifelse(LineageBroad_BoostTiterGroup %like% "HighBoosted","HighBoosted",NA
                                                               )))))
    dat
}
res_60to90 <- convert_to_lineage_and_titer(res_60to90)
res_60to90$Protocol <- factor(res_60to90$Protocol, levels=c("Frequent testing","Delayed detection"))

## Get better labels
res_60to90 <- res_60to90 %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Omicron; Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Omicron; Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Omicron; Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Omicron; Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))

res_60to90$BoostTiterGroup <- factor(res_60to90$BoostTiterGroup,
                                  levels=c( "Omicron; Boosted: ≤250 AU/ml", 
                                            "Omicron; Boosted: >250 AU/ml", 
                                            "Omicron; Not Boosted: ≤250 AU/ml", 
                                            "Omicron; Not Boosted: >250 AU/ml"))
res_60to90 <- res_60to90 %>% rename(`Immune status`=BoostTiterGroup)

## Get sample sizes for plot
samp_sizes <- expand_grid(LineageBroad="Omicron",BoostTiterGroup=unique(dat_subset_use$BoostTiterGroup),DetectionSpeed=unique(dat_subset_use$DetectionSpeed)) %>% drop_na() %>%
    left_join(
            dat_subset_use %>% filter(LineageBroad=="Omicron") %>% select(PersonID, BoostTiterGroup,DetectionSpeed) %>% distinct() %>% drop_na() %>%
    group_by(BoostTiterGroup,DetectionSpeed) %>% tally()) %>%
    mutate(n=ifelse(is.na(n),0,n))
## Create label and better factor names
samp_sizes <- samp_sizes %>% rename(Protocol=DetectionSpeed)
samp_sizes$Protocol <- factor(samp_sizes$Protocol, levels=c("Frequent testing","Delayed detection"))
samp_sizes <- samp_sizes %>% mutate(label=paste0("N=",n))
samp_sizes <- samp_sizes %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Omicron; Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Omicron; Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Omicron; Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Omicron; Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))
samp_sizes <- samp_sizes %>% 
    mutate(y = ifelse(
    BoostTiterGroup == "Omicron; Boosted: ≤250 AU/ml", 0.8,
        ifelse(
            BoostTiterGroup == "Omicron; Boosted: >250 AU/ml", 0.7,
            ifelse(
                BoostTiterGroup == "Omicron; Not Boosted: ≤250 AU/ml", 0.6,0.5))))

samp_sizes$BoostTiterGroup <- factor(samp_sizes$BoostTiterGroup,
                                    levels=c( "Omicron; Boosted: ≤250 AU/ml", 
                                              "Omicron; Boosted: >250 AU/ml", 
                                              "Omicron; Not Boosted: ≤250 AU/ml", 
                                              "Omicron; Not Boosted: >250 AU/ml"))
samp_sizes <- samp_sizes %>% rename(`Immune status`=BoostTiterGroup)


pA <- ggplot(data=res_60to90 %>% filter(LineageBroad=="Omicron"),col="None") +
    facet_wrap(~Protocol) +
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                    fill=`Immune status`,y=mean),alpha=0.5) +
    geom_line(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                  fill=`Immune status`,y=mean,col=`Immune status`))+
    geom_text(data=samp_sizes,aes(x=18, y=y,label=label,col=`Immune status`),show.legend=FALSE,size=2.5) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,20),breaks=seq(0,20,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Immune status",color="Immune status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_color_viridis_d(option="D") +
    scale_fill_viridis_d(option="D") +
    theme(legend.position="bottom",
          panel.spacing=unit(1,"lines"),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          plot.tag=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white")) +
    labs(tag="A")


# 100-200 sensitivity -----------------------------------------------------

load("data/data_for_regressions.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0) %>% ungroup()
dat_subset_use <- dat_subset_use %>% filter(TiterSensitivity==TRUE)

## Load 60-90 day filter sensitivity. Note that this fit does not include age group
load("outputs/titer_models_sensitivity/all_freq_2.RData")
freq_sens <- fit
load("outputs/titer_models_sensitivity/all_infreq_2.RData")
infreq_sens <- fit

## Get posterior draws
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(freq_sens$data$LineageBroad_BoostTiterGroup))
freq_sens_draws <- freq_sens %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(infreq_sens$data$LineageBroad_BoostTiterGroup))
infreq_sens_draws <- infreq_sens %>% epred_draws(newdata=newdata)

## Get posterior summaries
freq_sens_res <- freq_sens_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred)) %>% 
    mutate(Protocol="Frequent testing")
infreq_sens_res <- infreq_sens_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred))%>% 
    mutate(Protocol="Delayed detection")

res_sens <- bind_rows(freq_sens_res, infreq_sens_res)
res_sens <- convert_to_lineage_and_titer(res_sens)
res_sens$Protocol <- factor(res_sens$Protocol, levels=c("Frequent testing","Delayed detection"))

## Get better labels
res_sens <- res_sens %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Omicron; Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Omicron; Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Omicron; Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Omicron; Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))

res_sens$BoostTiterGroup <- factor(res_sens$BoostTiterGroup,
                                     levels=c( "Omicron; Boosted: ≤250 AU/ml", 
                                               "Omicron; Boosted: >250 AU/ml", 
                                               "Omicron; Not Boosted: ≤250 AU/ml", 
                                               "Omicron; Not Boosted: >250 AU/ml"))
res_sens <- res_sens %>% rename(`Immune status`=BoostTiterGroup)

## Get sample sizes for plot
samp_sizes <- expand_grid(LineageBroad="Omicron",BoostTiterGroup=unique(dat_subset_use$BoostTiterGroup),DetectionSpeed=unique(dat_subset_use$DetectionSpeed)) %>% drop_na() %>%
    left_join(
        dat_subset_use %>% filter(LineageBroad=="Omicron") %>% select(PersonID, BoostTiterGroup,DetectionSpeed) %>% distinct() %>% drop_na() %>%
            group_by(BoostTiterGroup,DetectionSpeed) %>% tally()) %>%
    mutate(n=ifelse(is.na(n),0,n))
## Create label and better factor names
samp_sizes <- samp_sizes %>% rename(Protocol=DetectionSpeed)
samp_sizes$Protocol <- factor(samp_sizes$Protocol, levels=c("Frequent testing","Delayed detection"))
samp_sizes <- samp_sizes %>% mutate(label=paste0("N=",n))
samp_sizes <- samp_sizes %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Omicron; Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Omicron; Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Omicron; Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Omicron; Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))
samp_sizes <- samp_sizes %>% 
    mutate(y = ifelse(
        BoostTiterGroup == "Omicron; Boosted: ≤250 AU/ml", 0.8,
        ifelse(
            BoostTiterGroup == "Omicron; Boosted: >250 AU/ml", 0.7,
            ifelse(
                BoostTiterGroup == "Omicron; Not Boosted: ≤250 AU/ml", 0.6,0.5))))

samp_sizes$BoostTiterGroup <- factor(samp_sizes$BoostTiterGroup,
                                     levels=c( "Omicron; Boosted: ≤250 AU/ml", 
                                               "Omicron; Boosted: >250 AU/ml", 
                                               "Omicron; Not Boosted: ≤250 AU/ml", 
                                               "Omicron; Not Boosted: >250 AU/ml"))
samp_sizes <- samp_sizes %>% rename(`Immune status`=BoostTiterGroup)


pB <- ggplot(data=res_sens %>% filter(LineageBroad=="Omicron"),col="None") +
    facet_wrap(~Protocol) +
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                    fill=`Immune status`,y=mean),alpha=0.5) +
    geom_line(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                  fill=`Immune status`,y=mean,col=`Immune status`))+
    geom_text(data=samp_sizes,aes(x=18, y=y,label=label,col=`Immune status`),show.legend=FALSE,size=2.5) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,20),breaks=seq(0,20,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Immune status",color="Immune status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_color_viridis_d(option="D") +
    scale_fill_viridis_d(option="D") +
    theme(legend.position="bottom",
          panel.spacing=unit(1,"lines"),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          plot.tag=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white")) +
    labs(tag="B")

pA/pB

ggsave(pA/pB,filename = "figures/supplement/figS13.png",width=8,height=7,dpi=300)
ggsave(pA/pB,filename = "figures/supplement/figS13.pdf",width=8,height=7)


