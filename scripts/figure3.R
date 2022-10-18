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

supplementS14 <- TRUE

load("data/data_for_regressions.RData")
dat_subset_use <- dat_subset_use %>% filter(DaysSinceDetection >= 0) %>% ungroup()

## Load main models
load("outputs/titer_models_age/all_freq_1.RData")
freq_age_titer <- fit
load("outputs/titer_models_age/all_infreq_1.RData")
infreq_age_titer <- fit

## Get posterior draws
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(freq_age_titer$data$LineageBroad_BoostTiterGroup),AgeGroup=unique(freq_age_titer$data$AgeGroup))
titer_freq_draws <- freq_age_titer %>% epred_draws(newdata=newdata)
newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(infreq_age_titer$data$LineageBroad_BoostTiterGroup),AgeGroup=unique(infreq_age_titer$data$AgeGroup))
titer_infreq_draws <- infreq_age_titer %>% epred_draws(newdata=newdata)

## Get posterior summaries
titer_freq_res <- titer_freq_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup, AgeGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred)) %>% 
    mutate(Protocol="Frequent testing")
titer_infreq_res <- titer_infreq_draws %>% group_by(DaysSinceDetection,LineageBroad_BoostTiterGroup, AgeGroup) %>%
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),mean=mean(.epred))%>% 
    mutate(Protocol="Delayed detection")

titer_res <- bind_rows(titer_freq_res, titer_infreq_res)

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
titer_res <- convert_to_lineage_and_titer(titer_res)
titer_res$Protocol <- factor(titer_res$Protocol, levels=c("Frequent testing","Delayed detection"))

## Get better labels
titer_res <- titer_res %>% mutate(BoostTiterGroup = case_when(
    BoostTiterGroup=="HighBoosted"~"Omicron; Boosted: >250 AU/ml",
    BoostTiterGroup=="LowBoosted"~"Omicron; Boosted: ≤250 AU/ml",
    BoostTiterGroup=="HighNotBoosted"~"Omicron; Not Boosted: >250 AU/ml",
    BoostTiterGroup=="LowNotBoosted"~"Omicron; Not Boosted: ≤250 AU/ml",
    TRUE~NA_character_))

titer_res <- titer_res %>%
    mutate(AgeGroup=case_when(
    AgeGroup=="(0,30]"~"<30 years",AgeGroup=="(30,50]"~"30-50 years",AgeGroup=="(50,100]"~"50+ years",TRUE~NA_character_
))

titer_res$BoostTiterGroup <- factor(titer_res$BoostTiterGroup,
                                  levels=c( "Omicron; Boosted: ≤250 AU/ml", 
                                            "Omicron; Boosted: >250 AU/ml", 
                                            "Omicron; Not Boosted: ≤250 AU/ml", 
                                            "Omicron; Not Boosted: >250 AU/ml"))
titer_res <- titer_res %>% rename(`Immune status`=BoostTiterGroup)

## Get sample sizes for plot
samp_sizes <- expand_grid(AgeGroup=unique(dat_subset_use$AgeGroup),LineageBroad="Omicron",BoostTiterGroup=unique(dat_subset_use$BoostTiterGroup),DetectionSpeed=unique(dat_subset_use$DetectionSpeed)) %>% drop_na() %>%
    
    left_join(
            dat_subset_use %>% filter(LineageBroad=="Omicron") %>% select(PersonID, AgeGroup,BoostTiterGroup,DetectionSpeed) %>% distinct() %>% drop_na() %>%
    group_by(AgeGroup,BoostTiterGroup,DetectionSpeed) %>% tally()) %>%
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
    mutate(AgeGroup=case_when(
        AgeGroup=="(0,30]"~"<30 years",AgeGroup=="(30,50]"~"30-50 years",AgeGroup=="(50,100]"~"50+ years",TRUE~NA_character_
    ))

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


p_bot <- ggplot(data=titer_res %>% filter(LineageBroad=="Omicron") %>%
           ## Remove 50+ years not boosted high and low
           filter(!(AgeGroup=="50+ years" & Protocol=="Frequent testing" &
                                  LineageBroad_BoostTiterGroup %in% c("OmicronLowNotBoosted","OmicronHighNotBoosted"))),
       col="None") +
    facet_grid(Protocol~AgeGroup) +
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


figSXX <- ggplot(data=titer_res %>% filter(AgeGroup=="30-50 years"),
                col="None") +
    facet_grid(LineageBroad~Protocol) +
    geom_ribbon(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                    fill=`Immune status`,y=mean),alpha=0.5) +
    geom_line(aes(x=DaysSinceDetection,ymin=lower,ymax=upper,
                  fill=`Immune status`,y=mean,col=`Immune status`))+
    #geom_text(data=samp_sizes,aes(x=18, y=y,label=label,col=`Immune status`),show.legend=FALSE,size=2.5) +
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
          plot.background = element_rect(fill="white",color="white"))


load(file="plots/titer_plot.RData")
load(file="plots/titer_plot_hists.RData")
p_titers1 <- p_titers1 + xlab("Antibody titer (AU/ml") + theme(axis.title.y=element_text(size=8,angle=90,vjust=1.5)) +
    labs(tag="A")
p_titers1
fig3 <- p_titers1 + p_bot + plot_layout(heights=c(1,1.75))
fig3

ggsave(filename="figures/figure3_new.png",plot=fig3,width=8,height=7,dpi=300)
ggsave(filename="figures/figure3_new.pdf",plot=fig3,width=8,height=7)
ggsave(filename="figures/figure3_back.pdf",plot=p_titer_hists,width=5,height=2)
