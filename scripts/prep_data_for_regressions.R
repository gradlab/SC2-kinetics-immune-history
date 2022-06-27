############################################################
## CLEAN DATA FOR REGRESSION MODELS
## James Hay <jameshay218@gmail.com> 
## 2022-06-22
## - This script reads in the raw data and generates a bunch of additional variables and filters in preparation for the regression analyses.
############################################################
library(tidyverse)
library(lubridate)
library(data.table)
library(ggbeeswarm)

colors <- c("black","tomato","red3","dodgerblue","blue")
names(colors) <- c("Other","Delta1","Delta2","Omicron1","Omicron2")

setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")

## Ct threshold to count as "low"
low_ct_threshold <- 30
## Antibody titer threshold to count as "high"
high_threshold <- 250

# Clean data --------------------------------------------------------------
## Read in cleaned data
dat <- read_csv("data/ct_data_cleaned.csv")

## Remove BA.2
dat <- dat %>% filter(LineageBroad != "BA.2" | is.na(LineageBroad))

## RAW DATA NUMBERS
dat %>% group_by(PersonID) %>% filter(CumulativeInfectionNumber == max(CumulativeInfectionNumber)) %>% select(PersonID, CumulativeInfectionNumber) %>% distinct()  %>% group_by(CumulativeInfectionNumber) %>% tally()
dat %>% filter(NewInfectionIdentified == 1) %>% group_by(LineageBroad) %>% tally()
dat %>% select(PersonID, CumulativeInfectionNumber, Symp_Ever) %>% distinct() %>% group_by(Symp_Ever) %>% tally()

## Find when individuals were detected relative to their last negative PCR result
dat_detections <- dat %>% 
    filter(DaysSinceDetection == 0) %>%
    group_by(PersonID, CumulativeInfectionNumber) %>%
    mutate(DetectionSpeed = ifelse(DaysSinceNegative >= 2 | is.na(DaysSinceNegative),
                                   "≥2 days",ifelse(DaysSinceNegative < 2, "≤1 days",NA))) %>% 
    mutate(DetectionSpeed=ifelse(is.na(DetectionSpeed),
                                 "≥2 days",DetectionSpeed)) %>%
    select(PersonID, TestDate, DetectionSpeed, DaysSinceNegative)
dat <- dat %>% left_join(dat_detections)%>% 
    group_by(PersonID, CumulativeInfectionNumber) %>%
    fill(DetectionSpeed, .direction="updown")
dat <- dat %>%
    mutate(DetectionSpeed = ifelse(NewInfectionIdentified == 1 & is.na(DetectionSpeed), "≥2 days", DetectionSpeed))
dat %>% filter(NewInfectionIdentified == 1 | DaysSinceDetection == 0) %>% select(PersonID, CumulativeInfectionNumber, DetectionSpeed) %>% distinct() %>% group_by(DetectionSpeed) %>% tally()

## Make the -1 Ct detections NA so we don't use them
dat<- dat %>% mutate(CtT1 = ifelse(!is.na(CtT1) & CtT1 < 0, NA, CtT1),
                     CtT2 = ifelse(!is.na(CtT2) & CtT2 < 0, NA, CtT2))

## Use negatives and assumed negatives as 40
dat <- dat %>% 
    mutate(CtT1=ifelse(TestResult %in% c("Assumed Negative","Negative"),40,CtT1),CtT2=ifelse(TestResult %in% c("Assumed Negative","Negative"),40,CtT2))


# Titer plots -------------------------------------------------------------

## Figure 2C plot
dat_titer_tmp <- dat %>% filter(NewInfectionIdentified == 1,!is.na(Titer)) %>% 
    select(PersonID, LineageBroad,VaccStatus,Titer,Symp_Ever) %>% distinct() %>%
    mutate(category=case_when(
        LineageBroad=="Other" & VaccStatus=="Unvaccinated" ~ "Other/Unvaccinated",
        LineageBroad=="Delta" & VaccStatus=="First dose" ~ "Delta/1-2 doses",
        LineageBroad=="Delta" & VaccStatus=="Second dose" ~ "Delta/1-2 doses",
        LineageBroad=="Omicron" & VaccStatus=="First dose" ~ "Omicron/1-2 doses",
        LineageBroad=="Omicron" & VaccStatus=="Second dose" ~ "Omicron/1-2 doses",
        LineageBroad=="Omicron" & VaccStatus=="Boosted" ~ "Omicron/Boosted",
        TRUE~NA_character_)) %>% filter(!is.na(category))

straptheboot <- function(df, goods, nstraps){
    
    out <- tibble()
    
    for(indexA in 1:nstraps){
        out <- bind_rows(out, (df %>% 
                                   mutate(n=1:n()) %>% 
                                   mutate(cut=ceiling(runif(1)*n())) %>% 
                                   filter(n!=cut) %>% 
                                   select(-n, -cut) %>% 
                                   summarise(across(all_of(goods), ~ mean(.x, na.rm = TRUE)))))
    }
    
    return(out)
    
}

dat_titer_tmp_summ <- dat_titer_tmp %>% group_by(category) %>% 
    straptheboot(goods="Titer",nstraps=1000) %>% 
    group_by(category) %>% 
    summarise(titer_mean=mean(Titer),titer_lwr=quantile(Titer,0.025),titer_upr=quantile(Titer,0.975))
dat_titer_tmp_summ$category <- factor(dat_titer_tmp_summ$category,levels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted"))
dat_titer_tmp$category <- factor(dat_titer_tmp$category,levels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted"))

dat_titer_tmp_summ <- dat_titer_tmp_summ %>% 
    mutate(label=paste0(signif(titer_mean,3)," AU \n(95% CI: ",signif(titer_lwr,3),"-",signif(titer_upr,3),")"))


p_titer_hists <- ggplot(data=dat_titer_tmp_summ) + 
    geom_histogram(data=dat_titer_tmp,aes(y=Titer),alpha=1,fill="grey70",binwidth=25) + 
    facet_wrap(~category,nrow=1) +
    scale_x_continuous(position="top",sec.axis=dup_axis(),expand=c(0,0),breaks=seq(0,100,by=25)) +
    theme_classic() +
    theme(panel.grid=element_blank(),
          axis.line = element_blank(),
          axis.title.x=element_blank(),
          strip.text=element_blank())
p_titer_hists

p_titers1 <- ggplot(data=dat_titer_tmp_summ) + 
    geom_beeswarm(data=dat_titer_tmp%>% mutate(Titer1 = ifelse(Titer == 800, Titer + rnorm(n(),0,5),Titer)),
                  aes(x=Titer1,y=as.numeric(category),fill=category,col=category),groupOnX=FALSE,cex=0.5,size=0.75,alpha=0.5) +
    geom_text(aes(y=as.numeric(category)+0.25,x=titer_mean + 70,label=label),size=2.2)+
    geom_point(aes(y=as.numeric(category), x=titer_mean),size=3) + 
    geom_segment(aes(y=as.numeric(category), yend=as.numeric(category), x=titer_lwr, xend=titer_upr)) + 
    geom_segment(aes(y=as.numeric(category)-0.2, yend=as.numeric(category)+0.2, x=titer_lwr, xend=titer_lwr)) + 
    geom_segment(aes(y=as.numeric(category)-0.2, yend=as.numeric(category)+0.2, x=titer_upr, xend=titer_upr)) + 
    theme_classic() +
    #geom_jitter(data=dat_titer_tmp%>% mutate(Titer1 = ifelse(Titer == 800, Titer + rnorm(n(),0,5),Titer)),
    #           aes(x=Titer1,y=as.numeric(category),fill=category,col=category),width=0,height=0.25)
    
    scale_color_manual(values=c("Other/Unvaccinated"="black","Delta/1-2 doses"="orange3","Omicron/1-2 doses"="mediumpurple1","Omicron/Boosted"="purple3")) +
    scale_fill_manual(values=c("Other/Unvaccinated"="black","Delta/1-2 doses"="orange3","Omicron/1-2 doses"="mediumpurple1","Omicron/Boosted"="purple3")) +
    xlab("Titer (AU/ml)") +
    ylab("Category") +
    scale_y_continuous(breaks=seq(1,4,by=1),labels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted")) +
    theme(legend.position = "none",plot.tag=element_text(face="bold"),axis.title.y=element_blank(),
          axis.text=element_text(size=8),axis.title=element_text(size=8)) +
    coord_flip() +
    labs(tag="C")

save(p_titer_hists,file="plots/titer_plot_hists.RData")
save(p_titers1,file="plots/titer_plot.RData")

## Stratify by symptom status
dat_titer_tmp_summ1 <- dat_titer_tmp %>% group_by(category,Symp_Ever) %>% 
    straptheboot(goods="Titer",nstraps=100) %>% 
    group_by(category,Symp_Ever) %>% 
    summarise(titer_mean=mean(Titer),titer_lwr=quantile(Titer,0.025),titer_upr=quantile(Titer,0.975))


dat_titer_tmp_summ1$category <- factor(dat_titer_tmp_summ1$category,levels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted"))
dat_titer_tmp$category <- factor(dat_titer_tmp$category,levels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted"))

dat_titer_tmp_summ1 <- dat_titer_tmp_summ1 %>% 
    mutate(label=paste0(signif(titer_mean,3)," AU \n(95% CI: ",signif(titer_lwr,3),"-",signif(titer_upr,3),")"))

p_titers2 <- ggplot(data=dat_titer_tmp_summ1) + 
    geom_beeswarm(data=dat_titer_tmp%>% mutate(Titer1 = ifelse(Titer == 800, Titer + rnorm(n(),0,5),Titer)),
                  aes(x=Titer1,y=as.numeric(Symp_Ever),fill=category,col=category),groupOnX=FALSE,cex=0.5,size=0.75,alpha=0.5) +
    geom_text(aes(y=as.numeric(Symp_Ever)+0.25,x=titer_mean + 70,label=label),size=2.2)+
    geom_point(aes(y=as.numeric(Symp_Ever), x=titer_mean),size=3) + 
    geom_segment(aes(y=as.numeric(Symp_Ever), yend=as.numeric(Symp_Ever), x=titer_lwr, xend=titer_upr)) + 
    geom_segment(aes(y=as.numeric(Symp_Ever)-0.2, yend=as.numeric(Symp_Ever)+0.2, x=titer_lwr, xend=titer_lwr)) + 
    geom_segment(aes(y=as.numeric(Symp_Ever)-0.2, yend=as.numeric(Symp_Ever)+0.2, x=titer_upr, xend=titer_upr)) + 
    theme_classic() +
    #geom_jitter(data=dat_titer_tmp%>% mutate(Titer1 = ifelse(Titer == 800, Titer + rnorm(n(),0,5),Titer)),
    #           aes(x=Titer1,y=as.numeric(category),fill=category,col=category),width=0,height=0.25)
    
    scale_color_manual(values=c("Other/Unvaccinated"="black","Delta/1-2 doses"="orange3","Omicron/1-2 doses"="mediumpurple1","Omicron/Boosted"="purple3")) +
    xlab("Titer (AU/ml)") +
    ylab("Category") +
    scale_y_continuous(breaks=seq(1,4,by=1),labels=c("Other/Unvaccinated","Delta/1-2 doses","Omicron/1-2 doses","Omicron/Boosted")) +
    theme(legend.position = "none",plot.tag=element_text(face="bold"),axis.title.y=element_blank(),
          axis.text=element_text(size=8),axis.title=element_text(size=8)) +
    coord_flip() +
    facet_wrap(~category,ncol=1) +
    labs(tag="C")



# Symptom data ------------------------------------------------------------
## How many symptomatic individuals?
dat %>% 
    filter(NewInfectionIdentified == 1) %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Ever) %>% distinct() %>% 
    group_by(Symp_Ever) %>% tally() %>%
    pivot_wider(names_from=Symp_Ever,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="Unknown") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic))

## By variant
dat %>% 
    filter(NewInfectionIdentified == 1, TestResult == "Positive") %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Ever,LineageBroad) %>% 
    distinct() %>% group_by(Symp_Ever, LineageBroad) %>% tally() %>%
    pivot_wider(names_from=Symp_Ever,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="Unknown") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic))


## How many symptomatic at detection?
## 503 detection as symptomatic, 605 detected as asymptomatic, rest unknown
dat %>% 
    filter(NewInfectionIdentified == 1, TestResult == "Positive") %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Init) %>% distinct() %>% 
    group_by(Symp_Init) %>% tally() %>%
    pivot_wider(names_from=Symp_Init,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="NA") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic))

## By variant
dat %>% 
    filter(NewInfectionIdentified == 1, TestResult == "Positive") %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Init,LineageBroad) %>% 
    distinct() %>% group_by(Symp_Init, LineageBroad) %>% tally() %>%
    pivot_wider(names_from=Symp_Init,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="NA") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic))

## A far higher proportion were symptomatic initially with Delta -- this is likely an ascertainment issue.
## Actually, still seems higher when stratifying by detection speed
dat %>% 
    filter(NewInfectionIdentified == 1, TestResult == "Positive") %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Init,LineageBroad,DetectionSpeed) %>% 
    distinct() %>% group_by(Symp_Init, LineageBroad,DetectionSpeed) %>% tally() %>%
    pivot_wider(names_from=Symp_Init,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="NA") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic)) %>%
    filter(!is.na(DetectionSpeed))

## Look at the incubation period distribution -- at least, time from detection to symptom onset
dat_incubation_periods <- dat %>% 
    filter(DaysSinceDetection == 0) %>%
    select(PersonID, TestDate, CumulativeInfectionNumber,Symp_Init,Symp_OnsetDate,LineageBroad,DetectionSpeed) %>% 
    mutate(IncubationPeriod = Symp_OnsetDate - TestDate) %>%
    distinct() 
dat_incubation_periods %>% group_by(LineageBroad,DetectionSpeed) %>%
    filter(!is.na(IncubationPeriod)) %>%
    summarize(MeanIncubation=mean(IncubationPeriod,na.rm=TRUE), N=n())

dat_incubation_periods %>% group_by(DetectionSpeed) %>%
    filter(!is.na(IncubationPeriod)) %>%
    summarize(MeanIncubation=mean(IncubationPeriod,na.rm=TRUE), N=n(),
              MedianDelay = median(IncubationPeriod,na.rm=TRUE))

p1 <- ggplot(dat_incubation_periods %>% filter(!is.na(IncubationPeriod),LineageBroad !="None")) +
    geom_histogram(aes(x=IncubationPeriod,fill=LineageBroad),binwidth=1,col="grey10",alpha=0.5) +
    geom_vline(data=.%>%group_by(LineageBroad,DetectionSpeed) %>% 
                   summarize(MeanIncubation=mean(IncubationPeriod,na.rm=TRUE)),aes(xintercept=MeanIncubation),
               linetype="dashed") +
    geom_vline(xintercept=0,size=1) +
    scale_fill_manual(name="",values=c("Delta"="red3","Omicron"="blue","Other"="black")) +
    coord_cartesian(xlim=c(-10,10)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(-10,10,by=1)) +
    ylab("Count") +
    xlab("Delay from detection to symptom onset") +
    theme_classic() +
    facet_grid(LineageBroad~DetectionSpeed) +
    theme(
        legend.title=element_text(size=8),legend.text=element_text(size=8),
        strip.background = element_blank(),strip.text=element_text(face="bold"))


## Look at gap between peak viral load and onset
dat_peak <- dat %>% 
    filter(TimeRelToPeak == 0) %>%
    select(PersonID, TestDate, CumulativeInfectionNumber,Symp_Init,Symp_OnsetDate,LineageBroad,DetectionSpeed) %>% 
    mutate(Delay = TestDate - Symp_OnsetDate) %>%
    distinct() 
dat_peak %>% group_by(LineageBroad,DetectionSpeed) %>%
    filter(!is.na(Delay)) %>%
    summarize(MeanDelay=mean(Delay,na.rm=TRUE), N=n(),
              MedianDelay = median(Delay,na.rm=TRUE))

dat_peak %>% group_by(DetectionSpeed) %>%
    filter(!is.na(Delay)) %>%
    summarize(MeanDelay=mean(Delay,na.rm=TRUE), N=n(),
              MedianDelay = median(Delay,na.rm=TRUE))

## Most symptom onsets are before peak viral load by a few days
p2 <- ggplot(dat_peak %>% filter(!is.na(Delay),LineageBroad !="None")) +
    geom_histogram(aes(x=Delay,fill=LineageBroad),binwidth=1,col="grey10",alpha=0.5) +
    geom_vline(data=.%>%group_by(LineageBroad,DetectionSpeed) %>% 
                   summarize(MeanDelay=mean(Delay,na.rm=TRUE)),aes(xintercept=MeanDelay),
               linetype="dashed") +
    geom_vline(xintercept=0,size=1) +
    scale_fill_manual(name="",values=c("Delta"="red3","Omicron"="blue","Other"="black")) +
    coord_cartesian(xlim=c(-10,10)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(-10,10,by=1)) +
    ylab("Count") +
    xlab("Time from symptom onset to peak Ct") +
    theme_classic() +
    facet_grid(LineageBroad~DetectionSpeed) +
    theme(
        legend.title=element_text(size=8),legend.text=element_text(size=8),
        strip.background = element_blank(),strip.text=element_text(face="bold"))

ggsave("figures/supplement/incubation_period_detection.png",p1,height=6,width=8,units="in",dpi=300)
ggsave("figures/supplement/incubation_period_peak.png",p2,height=6,width=8,units="in",dpi=300)

## More symptomatics in delayed detection
## Yes, slightly higher proportion
tab_for_chisq <- dat %>% 
    filter(NewInfectionIdentified == 1, TestResult == "Positive") %>%
    select(PersonID, CumulativeInfectionNumber,Symp_Ever,DetectionSpeed) %>% 
    distinct() %>% group_by(Symp_Ever, DetectionSpeed) %>% tally() %>%
    pivot_wider(names_from=Symp_Ever,values_from=n) %>%
    rename("Asymptomatic"="No","Symptomatic"="Yes","Unknown"="Unknown") %>%
    mutate(prop=Symptomatic/(Asymptomatic+Symptomatic)) %>%
    filter(!is.na(DetectionSpeed))

to_test <- as.table(rbind(as.numeric(tab_for_chisq[2,c(2,4)]),as.numeric(tab_for_chisq[1,c(2,4)])))
dimnames(to_test) <- list(DetectionSpeed=c("Delayed","Frequent"),SymptomStatus=c("Asymptomatic","Symptomatic"))
chisq.test(to_test)

# Plot raw trajectories ---------------------------------------------------
## Remove any other unknown lineages, removing 291 "None" infections
dat <- dat %>% filter(LineageBroad != "None")
dat$LineageBroad <- factor(dat$LineageBroad, levels=c("Other","Delta","Omicron"))
dat %>% filter(NewInfectionIdentified==1) %>% group_by(LineageBroad) %>% tally()

## Don't use dummy data for plot
dat_subset <- dat %>% filter(TestResult != "No sample", TestResult != "Positive (no CT)")
dat_subset <- dat_subset %>% filter(DaysSinceDetection <= 25) %>% filter(CumulativeInfectionNumber > 0)
p1_key <- c("≤1 days"="Frequent testing \n(≤1 days since last non-positive PCR)",
            "≥2 days"="Delayed detection \n(≥2 days since last non-positive PCR)")
dat_subset$DetectionSpeed <- p1_key[dat_subset$DetectionSpeed]
dat_subset$DetectionSpeed <- factor(dat_subset$DetectionSpeed, 
                                    levels=c("Frequent testing \n(≤1 days since last non-positive PCR)",
                                             "Delayed detection \n(≥2 days since last non-positive PCR)"))

dat_summary <- dat_subset %>% group_by(DetectionSpeed, DaysSinceDetection, LineageBroad) %>% 
    summarize(med_ct=mean(CtT1,na.rm=TRUE),
              lower_50 = quantile(CtT1, 0.25,na.rm=TRUE),
              upper_50 = quantile(CtT1, 0.75,na.rm=TRUE),
              lower_95 = quantile(CtT1, 0.025,na.rm=TRUE),
              upper_95 = quantile(CtT1, 0.975,na.rm=TRUE)) 

dat_subset$LineageBroad <- factor(dat_subset$LineageBroad, levels=c("Delta","Omicron","Other","None"))
dat_summary$LineageBroad <- factor(dat_summary$LineageBroad, levels=c("Delta","Omicron","Other","None"))

## Plot trajectories with mean Ct trajectory
samp_sizes <- dat_subset %>% select(PersonID, LineageBroad,DetectionSpeed) %>% distinct() %>% group_by(LineageBroad,DetectionSpeed) %>% tally() %>% mutate(label=paste0("N=",n))
p_traj <-  ggplot()  + 
    geom_line(data=dat_subset, aes(x=DaysSinceDetection,y=CtT1,group=PersonID,col=LineageBroad),alpha=0.25,size=0.1) +
    geom_line(data=dat_summary, aes(x=DaysSinceDetection,y=med_ct,col=LineageBroad),size=0.75) + 
    geom_text(data=samp_sizes,aes(x=20,y=20,label=label),size=3) +
    #geom_line(data=dat_summary2, aes(x=DaysSinceDetection,y=mean_ct,linetype="Only positive",col=LineageBroad),size=0.75,linetype="dashed") + 
    scale_y_continuous(trans="reverse",expand=c(0,0),breaks=seq(10,40,by=5)) + 
    scale_x_continuous(limits=c(-5,25),breaks=seq(-5,25,by=5)) +
    #facet_grid(DetectionSpeed~LineageBroad) + 
    facet_grid(LineageBroad~DetectionSpeed) +
    theme_classic() + 
    xlab("Days since detection") + 
    ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="dotted",col="black")+
    #scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    scale_color_manual(name="Lineage",values=c("Omicron"=unname(colors["Omicron2"]),"Delta"=unname(colors["Delta2"]),"Other"=unname(colors["Other"]))) +
    theme(legend.position="none", panel.grid.minor = element_blank(),
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background = element_blank(),strip.text=element_text(face="bold"))
p_traj



## Same thing but stratified by symptom status

dat_summary_symp <- dat_subset %>% group_by(DetectionSpeed, DaysSinceDetection, LineageBroad,Symp_Ever) %>% 
    summarize(med_ct=mean(CtT1,na.rm=TRUE),
              lower_50 = quantile(CtT1, 0.25,na.rm=TRUE),
              upper_50 = quantile(CtT1, 0.75,na.rm=TRUE),
              lower_95 = quantile(CtT1, 0.025,na.rm=TRUE),
              upper_95 = quantile(CtT1, 0.975,na.rm=TRUE)) 
dat_subset$LineageBroad <- factor(dat_subset$LineageBroad, levels=c("Delta","Omicron","Other","None"))
dat_summary_symp$LineageBroad <- factor(dat_summary_symp$LineageBroad, levels=c("Delta","Omicron","Other","None"))
## Plot trajectories with mean Ct trajectory
dat_subset <- dat_subset %>% mutate(Symptomatic = ifelse(Symp_Ever == "Unknown", "Unknown", ifelse(Symp_Ever == "Yes","Symptomatic","Asymptomatic")))
dat_summary_symp <- dat_summary_symp %>% mutate(Symptomatic = ifelse(Symp_Ever == "Unknown", "Unknown", ifelse(Symp_Ever == "Yes","Symptomatic","Asymptomatic")))

dat_subset$group <- paste0(dat_subset$LineageBroad,", ", dat_subset$Symptomatic)
dat_summary_symp$group <- paste0(dat_summary_symp$LineageBroad,", ", dat_summary_symp$Symptomatic)

p_traj_symp <-  ggplot()  + 
    geom_line(data=dat_subset, aes(x=DaysSinceDetection,y=CtT1,group=PersonID,col=group),alpha=0.25,size=0.1) +
    geom_line(data=dat_summary_symp, aes(x=DaysSinceDetection,y=med_ct,col=group),size=0.75) + 
    #geom_line(data=dat_summary2, aes(x=DaysSinceDetection,y=mean_ct,linetype="Only positive",col=LineageBroad),size=0.75,linetype="dashed") + 
    scale_y_continuous(trans="reverse",expand=c(0,0),breaks=seq(10,40,by=5)) + 
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,by=5)) +
    #facet_grid(DetectionSpeed~LineageBroad) + 
    facet_grid(LineageBroad~DetectionSpeed) +
    theme_classic() + 
    xlab("Days since detection") + 
    ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="dotted",col="black")+
    #scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    scale_color_manual(name="Lineage",values=c("Omicron, Unknown"="blue3","Delta, Unknown" = "red3","Other, Unknown"="black",
                                               
                                               "Omicron, Symptomatic"="purple","Omicron, Asymptomatic"="dodgerblue",
                                               "Delta, Symptomatic"="orange","Delta, Asymptomatic"="tomato",
                                               "Other, Symptomatic"="grey40","Other, Asymptomatic"="grey60")) +
    theme(legend.position="bottom", panel.grid.minor = element_blank(),
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background = element_blank(),strip.text=element_text(face="bold"))
p_traj_symp





## Don't use dummy data
dat <- dat %>% filter(TestResult != "No sample", TestResult != "Positive (no CT)", TestResult != "Positive (external)")

## Only use actual infections and between days 0 and 25 post detection
## There is one person with 6 cumulative exposures -- report on this person separately but don't include here
dat_subset <- dat %>% filter(DaysSinceDetection <= 25) %>% filter(CumulativeInfectionNumber > 0) %>% 
    mutate(PriorExposures=CumulativeExposureNumber-1)
dat_subset$PriorExposures <- as.factor(dat_subset$PriorExposures)

## Re-level vaccinations for regression
dat_subset <- dat_subset %>% filter(VaccStatus != "No record")
dat_subset$VaccStatus <- factor(dat_subset$VaccStatus, levels=c("Unvaccinated","First dose","Second dose","Boosted"))


## Find proportion of Cts < 30
dat_subset <- dat_subset %>% mutate(low_ct1 = !is.na(CtT1) & CtT1 < low_ct_threshold, low_ct2 = !is.na(CtT2) & CtT2 < low_ct_threshold)
dat_subset <- dat_subset %>% mutate(low_ct1 = as.numeric(low_ct1), low_ct2=as.numeric(low_ct2))

## Sample sizes
dat_subset <- dat_subset %>% mutate(DetectionSpeed = ifelse(DetectionSpeed=="≤1 days",
                                                            "Frequent testing","Delayed detection"))

## Factor levels
dat_subset$VaccStatus <- factor(dat_subset$VaccStatus,levels=c("Unvaccinated","First dose","Second dose","Boosted"))
dat_subset$CumulativeExposureNumber <- as.factor(dat_subset$CumulativeExposureNumber)
dat_subset$DaysSinceExposureGroup <- as.factor(dat_subset$DaysSinceExposureGroup)
dat_subset$LineageBroad <- factor(dat_subset$LineageBroad,levels=c("Other","Delta","Omicron"))
dat_subset$TiterGroup <- as.factor(dat_subset$TiterGroup)
dat_subset$TiterGroupAlt <- as.factor(dat_subset$TiterGroupAlt)
dat_subset$Symptomatic <- factor(dat_subset$Symp_Ever,levels=c("No","Yes"))

dat_subset$LineageBroad_VaccStatus <- paste0(dat_subset$LineageBroad, dat_subset$VaccStatus)
dat_subset$LineageBroad_CumulativeExposureNumber <- paste0(dat_subset$LineageBroad, dat_subset$CumulativeExposureNumber)
dat_subset$LineageBroad_DaysSinceExposureGroup <- paste0(dat_subset$LineageBroad, dat_subset$DaysSinceExposureGroup)
dat_subset$LineageBroad_TiterGroup <- paste0(dat_subset$LineageBroad, dat_subset$TiterGroup)
dat_subset$LineageBroad_TiterGroupAlt <- paste0(dat_subset$LineageBroad, dat_subset$TiterGroupAlt)
dat_subset$LineageBroad_BoostTiterGroup <- paste0(dat_subset$LineageBroad, dat_subset$BoostTiterGroup)
dat_subset$LineageBroad_BoostTiterGroupAlt <- paste0(dat_subset$LineageBroad, dat_subset$BoostTiterGroupAlt)
dat_subset$LineageBroad_Symptomatic <- paste0(dat_subset$LineageBroad, dat_subset$Symptomatic)
dat_subset$LineageBroad_Symptomatic_VaccStatus <- paste0(dat_subset$LineageBroad, dat_subset$Symptomatic, dat_subset$VaccStatus)

## NA titer groups should all be NA
dat_subset <- dat_subset %>% mutate(LineageBroad_BoostTiterGroup = as.character(ifelse(is.na(BoostTiterGroup),NA, LineageBroad_BoostTiterGroup)),
                        LineageBroad_BoostTiterGroupAlt = as.character(ifelse(is.na(BoostTiterGroupAlt),NA, LineageBroad_BoostTiterGroupAlt)),
                                    LineageBroad_TiterGroup = as.character(ifelse(is.na(TiterGroup), NA, LineageBroad_TiterGroup)),
                                    LineageBroad_TiterGroupAlt = as.character(ifelse(is.na(TiterGroupAlt), NA, LineageBroad_TiterGroupAlt)))


dat_subset$LineageBroad_VaccStatus <- as.factor(dat_subset$LineageBroad_VaccStatus)
dat_subset$LineageBroad_CumulativeExposureNumber <- as.factor(dat_subset$LineageBroad_CumulativeExposureNumber)
dat_subset$LineageBroad_DaysSinceExposureGroup <- as.factor(dat_subset$LineageBroad_DaysSinceExposureGroup)
dat_subset$LineageBroad_TiterGroup <- as.factor(dat_subset$LineageBroad_TiterGroup)
dat_subset$LineageBroad_TiterGroupAlt <- as.factor(dat_subset$LineageBroad_TiterGroupAlt)
dat_subset$LineageBroad_BoostTiterGroup <- as.factor(dat_subset$LineageBroad_BoostTiterGroup)
dat_subset$LineageBroad_BoostTiterGroupAlt <- as.factor(dat_subset$LineageBroad_BoostTiterGroupAlt)


boosttiter_sample_sizes <- dat_subset %>% filter(NewInfectionIdentified == 1, !is.na(BoostTiterGroup))  %>% 
    group_by(LineageBroad, BoostTiterGroup,DetectionSpeed) %>%
    tally()
write_csv(boosttiter_sample_sizes %>% arrange(DetectionSpeed, LineageBroad, BoostTiterGroup, n), "figures/boosttiter_sample_sizes.csv")


dat_subset_use <- dat_subset %>% 
    select(low_ct1, 
           PersonID,
           CtT1,
           VaccStatus, 
           DaysSinceDetection, 
           LineageBroad, 
           Symptomatic,
           DetectionSpeed,
           LastNegative,
           DaysSinceExposureGroup,
           CumulativeExposureNumber,
           CumulativeInfectionNumber,
           TiterGroup,
           TiterGroupAlt,
           LineageBroad_VaccStatus,
           LineageBroad_CumulativeExposureNumber,
           LineageBroad_DaysSinceExposureGroup,
           LineageBroad_TiterGroup,
           LineageBroad_TiterGroupAlt,
           BoostTiterGroup,
           BoostTiterGroupAlt,
           LineageBroad_BoostTiterGroup,
           LineageBroad_BoostTiterGroupAlt,
           UseLessThan60,
           UseLessThan90,
           Use60to90,
           LineageBroad_Symptomatic,
           LineageBroad_Symptomatic_VaccStatus,
           TiterSensitivity
    ) %>% distinct()

dat_subset_use %>% select(PersonID, TiterGroup) %>% ungroup() %>% distinct() %>% group_by(TiterGroup) %>% tally()
dat_subset_use %>% select(PersonID, BoostTiterGroup, LineageBroad) %>% ungroup() %>% distinct() %>% group_by(BoostTiterGroup,LineageBroad) %>% tally()
dat_subset_use %>% select(PersonID, VaccStatus, LineageBroad) %>% ungroup() %>% distinct() %>% group_by(VaccStatus,LineageBroad) %>% tally()


save(p_traj,file="plots/traj_plot.RData")
ggsave(filename="figures/all_trajectories.png",p_traj,height=8,width=8,units="in",dpi=300)
ggsave(filename="figures/all_trajectories.pdf",p_traj,height=8,width=8)

dat_subset_use %>% filter(!is.na(BoostTiterGroup)) %>%group_by(DaysSinceDetection,BoostTiterGroup,DetectionSpeed) %>% 
    summarize(mean_ct=mean(CtT1)) %>% ggplot() + geom_line(aes(x=DaysSinceDetection,y=mean_ct,col=BoostTiterGroup)) + facet_wrap(~DetectionSpeed) + scale_y_continuous(trans="reverse") + 
    scale_x_continuous(limits=c(-5,10))

dat_subset_use %>% group_by(DaysSinceDetection,DetectionSpeed,LineageBroad) %>% 
    summarize(mean_ct=mean(CtT1)) %>% ggplot() + geom_line(aes(x=DaysSinceDetection,y=mean_ct,col=LineageBroad)) + facet_wrap(~DetectionSpeed) + scale_y_continuous(trans="reverse") + 
    scale_x_continuous(limits=c(-3,5))

save(dat_subset_use, file="data/data_for_regressions.RData")

tmp <- dat_subset %>% left_join(dat_subset %>% 
                                    select(PersonID, CumulativeInfectionNumber,LineageBroad, DetectionSpeed) %>%
                                    group_by(PersonID, CumulativeInfectionNumber,LineageBroad, DetectionSpeed) %>%
                                    tally() %>% 
                                    filter(n > 3) %>%
                                    select(-n) %>%
                                    distinct() %>% 
                                    group_by(LineageBroad, DetectionSpeed) %>%sample_n(3) %>% mutate(Bolden=TRUE)) %>% mutate(Bolden=ifelse(is.na(Bolden),FALSE,TRUE))
ggplot()  + 
    geom_line(data= tmp %>% filter(Bolden==FALSE),aes(x=DaysSinceDetection,y=CtT1,group=PersonID,col=LineageBroad),alpha=0.25,size=0.1) +
    geom_line(data= tmp %>% filter(Bolden==TRUE),aes(x=DaysSinceDetection,y=CtT1,group=PersonID,col=LineageBroad),alpha=1,size=0.75) +
    scale_y_continuous(trans="reverse",expand=c(0,0),breaks=seq(10,40,by=5)) + 
    scale_x_continuous(limits=c(-5,25),breaks=seq(-5,25,by=5)) +
    #facet_grid(DetectionSpeed~LineageBroad) + 
    facet_grid(LineageBroad~DetectionSpeed) +
    theme_classic() + 
    xlab("Days since detection") + 
    ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="dotted",col="black")+
    #scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    scale_color_manual(name="Lineage",values=c("Omicron"=unname(colors["Omicron2"]),"Delta"=unname(colors["Delta2"]),"Other"=unname(colors["Other"]))) +
    theme(legend.position="none", panel.grid.minor = element_blank(),
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background = element_blank(),strip.text=element_text(face="bold"))

tmp_summary <- tmp %>% group_by(DetectionSpeed, TimeRelToPeak, LineageBroad) %>% 
    summarize(med_ct=mean(CtT1,na.rm=TRUE),
              lower_50 = quantile(CtT1, 0.25,na.rm=TRUE),
              upper_50 = quantile(CtT1, 0.75,na.rm=TRUE),
              lower_95 = quantile(CtT1, 0.025,na.rm=TRUE),
              upper_95 = quantile(CtT1, 0.975,na.rm=TRUE)) 

ggplot()  + 
    geom_line(data= tmp %>% filter(Bolden==FALSE),aes(x=TimeRelToPeak,y=CtT1,group=PersonID,col=LineageBroad),alpha=0.25,size=0.1) +
    #geom_line(data= tmp %>% filter(Bolden==TRUE),aes(x=TimeRelToPeak,y=CtT1,group=PersonID,col=LineageBroad),alpha=1,size=0.75) +
    geom_line(data=tmp_summary, aes(x=TimeRelToPeak,y=med_ct,col=LineageBroad),size=0.75) + 
    scale_y_continuous(trans="reverse",expand=c(0,0),breaks=seq(10,40,by=5)) + 
    scale_x_continuous(limits=c(-20,20),breaks=seq(-20,20,by=5)) +
    #facet_grid(DetectionSpeed~LineageBroad) + 
    facet_grid(LineageBroad~DetectionSpeed) +
    theme_classic() + 
    xlab("Days since detection") + 
    ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="dotted",col="black")+
    #scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    scale_color_manual(name="Lineage",values=c("Omicron"=unname(colors["Omicron2"]),"Delta"=unname(colors["Delta2"]),"Other"=unname(colors["Other"]))) +
    theme(legend.position="none", panel.grid.minor = element_blank(),
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background = element_blank(),strip.text=element_text(face="bold"))

dat_subset %>% group_by(PersonID, CumulativeInfectionNumber,DetectionSpeed,LineageBroad) %>% filter(TimeRelToPeak == min(TimeRelToPeak)) %>% group_by(DetectionSpeed,LineageBroad) %>% summarize(mean_wait = mean(TimeRelToPeak),sd_wait = sd(TimeRelToPeak))
