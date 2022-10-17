######################################################
## Analyze baseline regression results
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

age_key <- c("(0,30]"="<30","(30,50]"="30-50","(50,100]"=">50")
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
load("outputs/immune_models/baseline_age_1.rdata")
basemodel_freq <- fit
load("outputs/immune_models/baseline_age_2.rdata")
basemodel_infreq <- fit

## Read in the vaccine status and lineage regressions
load("outputs/immune_models/vaccine_and_lineage_age_1.rdata")
vacclineagemodel_freq <- fit
load("outputs/immune_models/vaccine_and_lineage_age_2.rdata")
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
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1),AgeGroup=unique(basemodel_freq$data$AgeGroup))
base_freq_draws <- basemodel_freq %>% epred_draws(newdata=newdata)
base_infreq_draws <- basemodel_infreq %>% epred_draws(newdata=newdata)

base_p_dat <- bind_rows(base_freq_draws %>% mutate(Protocol="Frequent testing"), base_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection,AgeGroup) %>% summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred))

p_base <- ggplot(base_p_dat, aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=AgeGroup,y=med),col="None") +
    geom_ribbon(alpha=0.5) +
    geom_line(aes(col=AgeGroup))+
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(breaks=seq(0,20,by=1)) +
    labs(y="Probability of Ct value <30",x="Days since detection") +
    theme_minimal() +
    theme(legend.position="bottom",
          plot.background = element_rect(fill="white",color="white")) +
  facet_wrap(~Protocol,ncol=1)
p_base
ggsave("plots/age/base_regression.png",p_base,width=5,height=7)


## Vaccine and lineage-specific models
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_VaccStatus=unique(vacclineagemodel_freq$data$LineageBroad_VaccStatus),
                       AgeGroup=unique(vacclineagemodel_freq$data$AgeGroup))
vacclineage_freq_effects <- conditional_effects(vacclineagemodel_freq)
vacclineage_infreq_effects <- conditional_effects(vacclineagemodel_infreq)


## Effect of booster status
tmp_effect_vacc <- conditional_effects(vacclineagemodel_freq,effects="DaysSinceDetection:LineageBroad_VaccStatus",conditions=data.frame("AgeGroup"="(0,30]"),re_formula=NA,int_conditions=list(LineageBroad_VaccStatus=c("OmicronBoosted","OmicronSecond dose")))
tmp_effect_vacc <- tmp_effect_vacc$`DaysSinceDetection:LineageBroad_VaccStatus`
tmp_effect_vacc <- tmp_effect_vacc %>% mutate(Protocol="Frequent testing")

tmp_effect_vacc2 <- conditional_effects(vacclineagemodel_infreq,effects="DaysSinceDetection:LineageBroad_VaccStatus",conditions=data.frame("AgeGroup"="(0,30]"),re_formula=NA,int_conditions=list(LineageBroad_VaccStatus=c("OmicronBoosted","OmicronSecond dose"))) 
tmp_effect_vacc2 <- tmp_effect_vacc2$`DaysSinceDetection:LineageBroad_VaccStatus`
tmp_effect_vacc2 <- tmp_effect_vacc2 %>% mutate(Protocol="Delayed detection")

tmp_dat <- bind_rows(tmp_effect_vacc, tmp_effect_vacc2)

tmp_dat <- tmp_dat %>% filter(LineageBroad_VaccStatus %like% "Omicron") %>% filter(LineageBroad_VaccStatus %in% c("OmicronBoosted","OmicronSecond dose")) %>%
  mutate(VaccStatus=ifelse(LineageBroad_VaccStatus=="OmicronBoosted","Omicron: Boosted","Omicron: Not Boosted"))
tmp_dat$VaccStatus <- factor(tmp_dat$VaccStatus,levels=c("Omicron: Boosted","Omicron: Not Boosted"))

tmp_dat$Protocol <- factor(tmp_dat$Protocol,levels=c("Frequent testing","Delayed detection"))


samp_sizes_age <- dat_subset_use %>% 
  ungroup() %>%
  select(PersonID,LineageBroad,AgeGroup,DetectionSpeed) %>% 
  distinct() %>% 
  group_by(AgeGroup,LineageBroad,DetectionSpeed) %>% tally() %>%
  mutate(label=paste0("N=",n)) %>% 
  drop_na() %>%
  mutate(y = ifelse(AgeGroup == "(0,30]",0.3, 
                    ifelse(AgeGroup == "(30,50]", 0.4,0.5))) %>%
  rename(Protocol=DetectionSpeed) %>%
  filter(LineageBroad == "Omicron")

samp_sizes_vacc<- dat_subset_use %>% select(PersonID, VaccStatus,LineageBroad,DetectionSpeed) %>% distinct() %>% group_by(VaccStatus,LineageBroad,DetectionSpeed) %>% tally() %>% mutate(label=paste0("N=",n)) %>% 
  mutate(VaccStatusLineage = paste0(LineageBroad, ": ", VaccStatus)) %>%
  mutate(y = ifelse(VaccStatus == "Unvaccinated",0.5, 
                    ifelse(VaccStatus == "First dose", 0.6,
                           ifelse(VaccStatus == "Second dose", 0.4,
                                  0.5)))) %>%
  rename(Protocol=DetectionSpeed) %>%
  filter(LineageBroad == "Omicron") %>%
  filter(VaccStatus %in% c("Boosted","Second dose")) %>%
  mutate(VaccStatus=ifelse(VaccStatus=="Boosted","Omicron: Boosted","Omicron: Not Boosted"))

samp_sizes_vacc$VaccStatus <- factor(samp_sizes_vacc$VaccStatus,levels=c("Omicron: Boosted","Omicron: Not Boosted"))
  

tmp_dat$Protocol <- factor(tmp_dat$Protocol, levels=c("Frequent testing","Delayed detection"))
samp_sizes_age$Protocol <- factor(samp_sizes_age$Protocol, levels=c("Frequent testing","Delayed detection"))
samp_sizes_vacc$Protocol <- factor(samp_sizes_vacc$Protocol, levels=c("Frequent testing","Delayed detection"))



p_vacclineage <-  ggplot(tmp_dat, col="None") +
    facet_wrap(~Protocol) +
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=VaccStatus),alpha=0.5) +
    geom_line(aes(x=DaysSinceDetection,y=estimate__,col=VaccStatus))+
  geom_text(data=samp_sizes_vacc,aes(x=22, y=y,label=label,col=VaccStatus),show.legend=FALSE) +
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
          plot.background = element_rect(fill="white",color="white"),
          plot.tag=element_text(face="bold"))+
  labs(tag="B")


## Effect of age
## Effect of booster status
tmp_effect_age <- conditional_effects(vacclineagemodel_freq,effects="DaysSinceDetection:AgeGroup",conditions=data.frame("LineageBroad_VaccStatus"="OmicronBoosted"),re_formula=NA)
tmp_effect_age <- tmp_effect_age$`DaysSinceDetection:AgeGroup`
tmp_effect_age <- tmp_effect_age %>% mutate(Protocol="Frequent testing")

tmp_effect_age2 <- conditional_effects(vacclineagemodel_infreq,effects="DaysSinceDetection:AgeGroup",conditions=data.frame("LineageBroad_VaccStatus"="OmicronBoosted"),re_formula=NA)
tmp_effect_age2 <- tmp_effect_age2$`DaysSinceDetection:AgeGroup`
tmp_effect_age2 <- tmp_effect_age2 %>% mutate(Protocol="Delayed detection")

tmp_dat <- bind_rows(tmp_effect_age,tmp_effect_age2)

tmp_dat$Protocol <- factor(tmp_dat$Protocol,levels=c("Frequent testing","Delayed detection"))


p_vacclineage_alt <-  ggplot(tmp_dat,col="None") +
  facet_wrap(~Protocol) +
  geom_ribbon(aes(x=DaysSinceDetection,ymin=lower__,ymax=upper__,fill=AgeGroup),alpha=0.5) +
  geom_line(aes(x=DaysSinceDetection,y=estimate__,col=AgeGroup))+
  geom_text(data=samp_sizes_age,aes(x=22, y=y,label=label,col=AgeGroup),show.legend=FALSE) +
  scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
  scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
  labs(y="Probability of Ct value <30",x="Days since detection",fill="Age group",color="Age group") +
  theme_classic() +
  geom_vline(xintercept=5,linetype="dashed") +
  geom_hline(yintercept=0.05,linetype="dashed") +
  scale_color_viridis_d(option="C") +
  scale_fill_viridis_d(option="C") +
  #scale_color_manual(values=immunecolors[3:4]) +
  #scale_fill_manual(values=immunecolors[3:4])+
  theme(legend.position=c(0.85,0.8),
        axis.text=element_text(size=8),axis.title=element_text(size=8),
        legend.title=element_text(size=6),legend.text=element_text(size=6),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        plot.background = element_rect(fill="white",color="white"),
        plot.tag=element_text(face="bold")) +
  labs(tag="A")

ggsave("plots/age/p_vacc_lineage_age.png",p_vacclineage_alt/p_vacclineage,width=7,height=7)


all_colors <- c("Delta: Boosted"="orange", "Delta: Second dose"="yellow", "Delta: Unvaccinated"="tomato", 
  "Omicron: Boosted"="purple3", "Omicron: First dose"="lightblue", "Omicron: Second dose"="mediumpurple1", 
  "Omicron: Unvaccinated"="darkblue", "Other: First dose"="grey80", "Other: Second dose"="grey50", 
  "Other: Unvaccinated"="black")



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

load("outputs/titer_models_age/all_freq_1.RData")
fit1 <- fit
load("outputs/titer_models_age/all_infreq_1.RData")
fit2 <- fit

newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(fit1$data$LineageBroad_BoostTiterGroup),
                       AgeGroup=unique(fit1$data$AgeGroup))
titer_freq_draws <- fit1 %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(fit2$data$LineageBroad_BoostTiterGroup),
                       AgeGroup=unique(fit2$data$AgeGroup))
titer_infreq_draws <- fit2 %>% epred_draws(newdata=newdata)



titer_p_dat <- bind_rows(titer_freq_draws %>% mutate(Protocol="Frequent testing"), 
                         titer_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
  group_by(Protocol, DaysSinceDetection,LineageBroad_BoostTiterGroup,AgeGroup) %>% 
  summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred)) %>%
  mutate(LineageBroad = ifelse(LineageBroad_BoostTiterGroup %like% "Other","Other",
                               ifelse(LineageBroad_BoostTiterGroup %like% "Delta","Delta","Omicron")))

titer_p_dat$Protocol <- factor(titer_p_dat$Protocol,levels=c("Frequent testing","Delayed detection"))

tmp_dat1 <- titer_p_dat %>% filter(LineageBroad == "Omicron") 


boosttitergroup_key <- c("OmicronHighBoosted" = "Omicron: Boosted >250 AU/ml",
                         "OmicronLowBoosted" = "Omicron: Boosted ≤250 AU/ml",
                         "OmicronHighNotBoosted"="Omicron: Not Boosted >250 AU/ml",
                         "OmicronLowNotBoosted"="Omicron: Not Boosted ≤250 AU/ml")

tmp_dat1$`Immune status` <- boosttitergroup_key[as.character(tmp_dat1$LineageBroad_BoostTiterGroup)]
tmp_dat1$`Detection group` <- factor(tmp_dat1$`Protocol`, levels=c("Frequent testing","Delayed detection"))
tmp_dat1$`Immune status` <- factor(tmp_dat1$`Immune status`, levels=c("Omicron: Boosted ≤250 AU/ml","Omicron: Boosted >250 AU/ml",
                                                                      "Omicron: Not Boosted ≤250 AU/ml","Omicron: Not Boosted >250 AU/ml"))

tmp_dat1$AgeGroup <- age_key[as.character(tmp_dat1$AgeGroup)]
tmp_dat1$AgeGroup <- factor(tmp_dat1$AgeGroup, levels=c("<30","30-50",">50"))


all_samp <- expand_grid(AgeGroup=unique(dat_subset_use$AgeGroup),
            LineageBroad_BoostTiterGroup=unique(dat_subset_use$LineageBroad_BoostTiterGroup),
            DetectionSpeed=unique(dat_subset_use$DetectionSpeed)) %>%
  drop_na()

samp_sizes_titer <- dat_subset_use %>% 
  ungroup() %>%
  select(PersonID, LineageBroad_BoostTiterGroup,AgeGroup,DetectionSpeed) %>% 
  distinct() %>% 
  group_by(AgeGroup,LineageBroad_BoostTiterGroup,DetectionSpeed) %>% 
  tally() 

samp_sizes_titer <- samp_sizes_titer %>% bind_rows(setdiff(all_samp, samp_sizes_titer %>% select(-n)) %>% mutate(n = 0))

samp_sizes_titer <- samp_sizes_titer %>% 
  mutate(n=ifelse(is.na(n),0,n)) %>%
  mutate(label=paste0("N=",n)) %>%
  rename(Protocol=DetectionSpeed) %>% filter(LineageBroad_BoostTiterGroup %like% "Omicron") %>%
  mutate(y = ifelse(LineageBroad_BoostTiterGroup == "OmicronHighBoosted",0.8, 
                    ifelse(LineageBroad_BoostTiterGroup == "OmicronLowBoosted", 0.7,
                           ifelse(LineageBroad_BoostTiterGroup == "OmicronHighNotBoosted", 0.6,
                                  0.5))))  %>%
  drop_na()

samp_sizes_titer$`Immune status` <- boosttitergroup_key[as.character(samp_sizes_titer$LineageBroad_BoostTiterGroup)]
samp_sizes_titer$AgeGroup <- age_key[as.character(samp_sizes_titer$AgeGroup)]
samp_sizes_titer$AgeGroup <- factor(samp_sizes_titer$AgeGroup, levels=c("<30","30-50",">50"))
samp_sizes_titer$Protocol <- factor(samp_sizes_titer$Protocol, levels=c("Frequent testing","Delayed detection"))

p_titerlineage <-  ggplot(tmp_dat1, 
                          col="None") +
  facet_grid(AgeGroup~Protocol) +
  geom_ribbon(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=`Immune status`,group=interaction(`Immune status`,AgeGroup),y=med), alpha=0.5) +
  geom_line(aes(col=`Immune status`,y=med,x=DaysSinceDetection))+
  geom_text(data=samp_sizes_titer,aes(x=18, y=y,label=label,col=`Immune status`),show.legend=FALSE) +
  scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
  scale_x_continuous(limits=c(0,20),breaks=seq(0,20,by=5)) +
  labs(y="Probability of Ct value <30",x="Days since detection",fill="Vaccination status",color="Vaccination status") +
  theme_classic() +
  geom_vline(xintercept=5,linetype="dashed") +
  geom_hline(yintercept=0.05,linetype="dashed") +
  scale_color_viridis_d(option="D") +
  scale_fill_viridis_d(option="D") +
  theme(legend.position="bottom",
        axis.text=element_text(size=8),axis.title=element_text(size=8),
        legend.title=element_text(size=6),legend.text=element_text(size=6),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        plot.background = element_rect(fill="white",color="white"))
p_titerlineage
ggsave("plots/age/p_vacc_titer_age.png",p_titerlineage,width=7,height=7)

