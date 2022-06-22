############################################################
## COUNT NUMBER OF REBOUNDS
## James Hay <jameshay218@gmail.com> 
## 2022-06-22
## - Classifies as extracts trajectories which count as rebound. Change the two variables (L20 and L22) to identify rebounds under different criteria.
############################################################
setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")
library(tidyverse)
library(data.table)
library(patchwork)


lineage_colors <- c("Omicron"="blue","Delta"="tomato","Other"="black")
lineage_colors1 <- c("Omicron"="dodgerblue","Delta"="red3","Other"="black")
vacc_status_colors <- c("Boosted"="purple3","Second dose"="mediumpurple1")

dat_to_save <- read_csv("data/ct_data_cleaned.csv")

## To what Ct value does clearance need to occur to count as initially cleared?
ct_rebound_threshold <- 30
## How many consecutive days does the trajectory need to be above this threshold to be counted as cleared?
days_clearance <- 2

## Find when individuals were detected relative to their last negative PCR result
dat_detections <- dat_to_save %>% 
    filter(DaysSinceDetection == 0) %>%
    group_by(PersonID, CumulativeInfectionNumber) %>%
    mutate(DetectionSpeed = ifelse(DaysSinceNegative >= 2 | is.na(DaysSinceNegative),
                                   ">=2 days",ifelse(DaysSinceNegative < 2, "<=1 days",NA))) %>% 
    mutate(DetectionSpeed=ifelse(is.na(DetectionSpeed),
                                 ">=2 days",DetectionSpeed)) %>%
    select(PersonID, CumulativeInfectionNumber,TestDate, DetectionSpeed, DaysSinceNegative)
dat_to_save <- dat_to_save %>% left_join(dat_detections)%>% 
    fill(DetectionSpeed, .direction="down")
dat_to_save <- dat_to_save %>% mutate(DetectionSpeed = ifelse(is.na(DetectionSpeed),">=2 days",DetectionSpeed))

## Make the -1 Ct detections NA so we don't use them
dat_to_save<- dat_to_save %>% mutate(CtT1 = ifelse(!is.na(CtT1) & CtT1 < 0, NA, CtT1),
                     CtT2 = ifelse(!is.na(CtT2) & CtT2 < 0, NA, CtT2))


## Don't use dummy data
dat_subset <- dat_to_save %>% filter(TestResult != "No sample", TestResult != "Positive (no CT)", TestResult != "Positive (external)")
dat_subset <- dat_subset %>% #filter(DaysSinceDetection <= 25) %>% 
    filter(CumulativeInfectionNumber > 0)
p1_key <- c("<=1 days"="Frequent testing",
            ">=2 days"="Delayed detection")
dat_subset$DetectionSpeed <- p1_key[dat_subset$DetectionSpeed]
dat_subset$DetectionSpeed <- factor(dat_subset$DetectionSpeed, levels=c("Frequent testing","Delayed detection"))

## Remove any other unknown lineages, removing 291 "None" infections
dat_subset <- dat_subset %>% filter(LineageBroad != "None")
dat_subset$LineageBroad <- factor(dat_subset$LineageBroad, levels=c("Other","Delta","Omicron","None"))


## Stat 1 -- how many tests do we have on each day post detection?
p_samples <- dat_subset %>% filter(DaysSinceDetection <= 25,!is.na(CtT1)) %>% group_by(DetectionSpeed,DaysSinceDetection) %>% tally() %>% ggplot() + geom_line(aes(x=DaysSinceDetection,y=n)) + 
    ylab("Number of samples available") +
    scale_y_continuous(breaks=seq(0,2200,by=500),limits=c(0,2200)) +
    facet_wrap(~DetectionSpeed,ncol=1) +
    xlab("Days since initial detection") +
    theme_linedraw()

ggsave(paste0("plots/rebounds/number_of_tests",days_clearance,".png"),p_samples,height=6,width=7,units="in",dpi=300)


## Stat 2 -- how long after first being Ct > 30 do we keep testing people?
## Find the first day with CtT1 > 30 after the peak
if(days_clearance == 4){
    end_of_infection <- dat_subset %>% #filter(DaysSinceDetection <= 25) %>%
        group_by(PersonID,CumulativeInfectionNumber) %>%
        ## Only times where we've had two consecutive Ct > 30 and before that was Ct < 30
        filter(CtT1 >= ct_rebound_threshold & 
                   lag(CtT1,1) >= ct_rebound_threshold & 
                   lag(CtT1,2) >= ct_rebound_threshold & 
                   lag(CtT1,3) >= ct_rebound_threshold & 
                   lag(CtT1,4) < ct_rebound_threshold) %>% 
        ## Find first day after two consecutive Ct > 30
        filter(DaysSinceDetection == min(DaysSinceDetection)) %>%
        select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>%
        rename(EndOfInfection = DaysSinceDetection)
} else if (days_clearance == 3){
    end_of_infection <- dat_subset %>% #filter(DaysSinceDetection <= 25) %>%
        group_by(PersonID,CumulativeInfectionNumber) %>%
        ## Only times where we've had two consecutive Ct > 30 and before that was Ct < 30
        filter(CtT1 >= ct_rebound_threshold & 
                   lag(CtT1,1) >= ct_rebound_threshold & 
                   lag(CtT1,2) >= ct_rebound_threshold & 
                   lag(CtT1,3) < ct_rebound_threshold) %>% 
        ## Find first day after two consecutive Ct > 30
        filter(DaysSinceDetection == min(DaysSinceDetection)) %>%
        select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>%
        rename(EndOfInfection = DaysSinceDetection)
} else {
    end_of_infection <- dat_subset %>% #filter(DaysSinceDetection <= 25) %>%
        group_by(PersonID,CumulativeInfectionNumber) %>%
        ## Only times where we've had two consecutive Ct > 30 and before that was Ct < 30
        filter(CtT1 >= ct_rebound_threshold & 
                   lag(CtT1,1) >= ct_rebound_threshold & 
                   lag(CtT1,2) < ct_rebound_threshold) %>% 
        ## Find first day after two consecutive Ct > 30
        filter(DaysSinceDetection == min(DaysSinceDetection)) %>%
        select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>%
        rename(EndOfInfection = DaysSinceDetection)
}

## Denominator
end_of_infection %>% select(PersonID, CumulativeInfectionNumber) %>% distinct() %>% nrow()

end_of_infection_blank <- end_of_infection %>% left_join(dat_subset) %>%  
    group_by(PersonID,CumulativeInfectionNumber) %>%
    ## Find days after the two consecutive Ct > 30
    filter(DaysSinceDetection <=50,TestResult %in% c("Negative","Positive","Inconclusive")) %>% mutate(DaysExtra = DaysSinceDetection - EndOfInfection)

## How many samples X days after assumed clearance?
p_samples_2 <- end_of_infection_blank %>% group_by(DetectionSpeed,DaysExtra) %>% summarize(n=n(),n_pos=sum(CtT1 < ct_rebound_threshold),prop_low=n_pos/n) %>%  
    filter(DaysExtra >= 0) %>% 
    ggplot() + 
    geom_line(aes(x=DaysExtra,y=n)) +
    facet_wrap(~DetectionSpeed,ncol=1) +
    theme_linedraw() +
    xlab("Days since two consecutive Ct > 30") +
    ylab("Number of samples")

ggsave(paste0("plots/rebounds/number_of_tests_after_clearance",days_clearance,".png"),p_samples_2,height=6,width=7,units="in",dpi=300)


## Rebounds are classified as having two consecutive Ct < 30 after being "cleared"
## 2 days of Ct < 30 -- easiest rebound
rebounds1 <- end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & DaysExtra > 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber,DaysSinceDetection) %>% distinct() %>% 
    mutate(Rebound="2 x Ct<30") %>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 -- second easiest rebound
rebounds2 <-  end_of_infection_blank %>% filter(CtT1 < 29 & lag(CtT1, 1) < 29 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="2 x Ct<29")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 and 2 Ct difference between each consecutive measurement -- difficult rebound
rebounds3 <-  end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & DaysExtra >= 0 & (CtT1 < (lag(CtT1,1) + 1))) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="2 x Ct<30 and decrease by 2+ cycles")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 25 -- most stringent rebound
rebounds4 <-  end_of_infection_blank %>% filter(CtT1 < 25 & lag(CtT1, 1) < 25 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="2 x Ct<25")%>%
    rename(ReboundDate = DaysSinceDetection)

rebounds <- bind_rows(rebounds1, rebounds3,rebounds4)
rebounds$Rebound <- factor(rebounds$Rebound, levels=c("2 x Ct<30","2 x Ct<29", "2 x Ct<30 and decrease by 2+ cycles","2 x Ct<25"))
rebounds$r_numeric <- as.numeric(rebounds$Rebound)

end_of_infection <- end_of_infection_blank %>% left_join(rebounds)
end_of_infection <- end_of_infection %>% mutate(ReboundLog = ifelse(!is.na(Rebound),TRUE,FALSE))
end_of_infection1 <- end_of_infection %>% group_by(PersonID, CumulativeInfectionNumber) %>% filter(r_numeric == max(r_numeric))


p1 <- end_of_infection1 %>% ggplot() + 
    geom_line(aes(y=CtT1, x=DaysExtra,group=interaction(PersonID, CumulativeInfectionNumber),col=ReboundLog,alpha=ReboundLog,size=ReboundLog)) + 
    scale_y_continuous(trans="reverse",breaks=seq(10,40,by=5)) + 
    scale_x_continuous(breaks=seq(-25,25,by=5)) +
    coord_cartesian(xlim=c(-25,25)) +
    geom_vline(xintercept=0,col="red",linetype="dashed") + 
    geom_hline(yintercept=30,col="red",linetype="dashed") + 
    xlab("Days since two consecutive Ct > 30") + 
    scale_color_manual(name="Rebound",values=c("FALSE"="black","TRUE"="blue")) +
    scale_alpha_manual(name="Rebound",values=c("FALSE"=0.1,"TRUE"=1)) +
    scale_size_manual(name="Rebound",values=c("FALSE"=0.5,"TRUE"=0.5)) +
    ylab("Ct value") +
    facet_wrap(~DetectionSpeed,ncol=1) +
    theme_linedraw() +
    theme(legend.position="bottom")

samp_sizes <- end_of_infection1 %>% filter(ReboundLog == TRUE) %>% select(PersonID, CumulativeInfectionNumber,DetectionSpeed) %>% distinct() %>% 
  group_by(DetectionSpeed) %>% tally()%>% mutate(label=paste0("N=",n))

p2 <- end_of_infection1 %>% filter(ReboundLog == TRUE) %>%
    ggplot() + 
    geom_line(aes(y=CtT1, x=TimeRelToPeak,group=interaction(PersonID, CumulativeInfectionNumber,LineageBroad),col=LineageBroad),alpha=1,size=0.4) + 
    geom_text(data=samp_sizes,aes(x=15,y=20,label=label),size=3,show.legend=FALSE) +
    scale_y_continuous(trans="reverse",expand=c(0,0),breaks=seq(10,40,by=5)) + 
    scale_x_continuous(breaks=seq(-25,25,by=5)) +
    coord_cartesian(xlim=c(-10,25)) +
    geom_hline(yintercept=30,col="black",linetype="dotted") + 
    geom_vline(xintercept=0,linetype="dashed",col="black") +
    xlab("Day relative to peak Ct") + 
    ylab("Ct value") +
    scale_color_manual(name="Lineage",values=lineage_colors) +
    facet_wrap(~DetectionSpeed,ncol=2) +
    theme_classic() +
    theme(legend.position="bottom",
          legend.title=element_text(size=8),legend.text=element_text(size=8),
          strip.background = element_blank(),strip.text=element_text(face="bold"))
p2


p3 <- end_of_infection1 %>% 
    filter(ReboundLog == TRUE)  %>%
    mutate(key = paste0("ID: ",PersonID,", infection: ",CumulativeInfectionNumber)) %>%
    ggplot() + 
    geom_rect(data=end_of_infection1 %>% filter(ReboundLog == TRUE)  %>%
                  mutate(key = paste0("ID: ",PersonID,", infection: ",CumulativeInfectionNumber)) %>%
                  select(key, Rebound) %>% distinct(), 
              aes(fill=Rebound),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.25) +
    
    geom_hline(yintercept=30,linetype="dashed",col="#440154FF") +
    #geom_hline(yintercept=29,linetype="dashed",col="#3B528BFF") +
    geom_hline(yintercept=25,linetype="dashed",col="#21908CFF") +
    geom_vline(aes(xintercept=EndOfInfection),col="red",linetype="dotted") +
    geom_vline(aes(xintercept=ReboundDate),col="red") +
    scale_fill_viridis_d(name="Rebound definition",guide="legend") +
    geom_line(aes(y=CtT1, x=DaysSinceDetection)) + 
    geom_point(aes(y=CtT1, x=DaysSinceDetection),size=1,shape=4) + 
    scale_y_continuous(trans="reverse",breaks=seq(10,40,by=5)) + 
    scale_x_continuous(breaks=seq(-25,50,by=5)) +
    coord_cartesian(xlim=c(0,25)) +
    xlab("Days since initial detection") + 
    ylab("Ct value") +
    facet_wrap(~key,ncol=6) +
    theme_linedraw() +
    theme(legend.position="bottom")

N <- end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber,LineageBroad) %>% group_by(LineageBroad) %>% distinct() %>% tally() %>% rename(Total=n)

table1 <- end_of_infection %>% ungroup() %>%  select(PersonID,CumulativeInfectionNumber, Rebound,LineageBroad) %>% distinct() %>% group_by(Rebound,LineageBroad) %>% tally() %>%
    mutate(Rebound=ifelse(is.na(Rebound),"No rebound",as.character(Rebound))) %>%
    pivot_wider(names_from=Rebound,values_from=n) %>%
    left_join(N)
if(!("2 x Ct<30" %in% colnames(table1))) {
    table1$`2 x Ct<30` <- NA
}
if(!("2 x Ct<25" %in% colnames(table1))) {
    table1$`2 x Ct<25` <- NA
}
if(!("2 x Ct<29" %in% colnames(table1))) {
    table1$`2 x Ct<29` <- NA
}
if(!("2 x Ct<30 and decrease by 2+ cycles" %in% colnames(table1))) {
    table1$`2 x Ct<30 and decrease by 2+ cycles` <- NA
}   
table1 <- table1 %>%  
  mutate(`2 x Ct<30` = ifelse(is.na(`2 x Ct<30`),0,`2 x Ct<30`)) %>%
    mutate(`2 x Ct<29` = ifelse(is.na(`2 x Ct<29`),0,`2 x Ct<29`)) %>%
    mutate(`2 x Ct<30 and decrease by 2+ cycles` = ifelse(is.na(`2 x Ct<30 and decrease by 2+ cycles`),0,`2 x Ct<30 and decrease by 2+ cycles`)) %>%
    mutate(`2 x Ct<25` = ifelse(is.na(`2 x Ct<25`),0,`2 x Ct<25`)) %>%
    
    mutate(`2 x Ct<30`=paste0(`2 x Ct<30`, " (",signif(`2 x Ct<30`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<29`=paste0(`2 x Ct<29`, " (",signif(`2 x Ct<29`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<30 and decrease by 2+ cycles`=paste0(`2 x Ct<30 and decrease by 2+ cycles`, " (",signif(`2 x Ct<30 and decrease by 2+ cycles`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<25`=paste0(`2 x Ct<25`, " (",signif(`2 x Ct<25`*100/Total,3),"%)")) %>%
    select(-`No rebound`) %>%
    rename(Lineage = LineageBroad,`Total infections`=Total)

N <- end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber,VaccStatus) %>% group_by(VaccStatus) %>% distinct() %>% tally() %>% rename(Total=n)
table2 <- end_of_infection %>% ungroup() %>%  select(PersonID,CumulativeInfectionNumber, Rebound,VaccStatus) %>% distinct() %>% group_by(Rebound,VaccStatus) %>% 
    tally() %>%
    mutate(Rebound=ifelse(is.na(Rebound),"No rebound",as.character(Rebound))) %>%
    pivot_wider(names_from=Rebound,values_from=n) %>%
    left_join(N)

if(!("2 x Ct<30" %in% colnames(table2))) {
    table2$`2 x Ct<30` <- NA
}
if(!("2 x Ct<25" %in% colnames(table2))) {
    table2$`2 x Ct<25` <- NA
}
if(!("2 x Ct<29" %in% colnames(table2))) {
    table2$`2 x Ct<29` <- NA
}
if(!("2 x Ct<30 and decrease by 2+ cycles" %in% colnames(table2))) {
    table2$`2 x Ct<30 and decrease by 2+ cycles` <- NA
}   

table2 <- table2 %>%
    mutate(`2 x Ct<30` = ifelse(is.na(`2 x Ct<30`),0,`2 x Ct<30`)) %>%
    mutate(`2 x Ct<29` = ifelse(is.na(`2 x Ct<29`),0,`2 x Ct<29`)) %>%
    mutate(`2 x Ct<30 and decrease by 2+ cycles` = ifelse(is.na(`2 x Ct<30 and decrease by 2+ cycles`),0,`2 x Ct<30 and decrease by 2+ cycles`)) %>%
    mutate(`2 x Ct<25` = ifelse(is.na(`2 x Ct<25`),0,`2 x Ct<25`)) %>%
    
    mutate(`2 x Ct<30`=paste0(`2 x Ct<30`, " (",signif(`2 x Ct<30`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<29`=paste0(`2 x Ct<29`, " (",signif(`2 x Ct<29`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<30 and decrease by 2+ cycles`=paste0(`2 x Ct<30 and decrease by 2+ cycles`, " (",signif(`2 x Ct<30 and decrease by 2+ cycles`*100/Total,3),"%)")) %>%
    mutate(`2 x Ct<25`=paste0(`2 x Ct<25`, " (",signif(`2 x Ct<25`*100/Total,3),"%)")) %>%
    select(-`No rebound`) %>%
    rename(`Vaccination status` = VaccStatus,`Total infections`=Total)



table1_for_chisq <- end_of_infection %>% select(PersonID, CumulativeInfectionNumber, LineageBroad, ReboundLog) %>% mutate(LineageBroadAlt = ifelse(LineageBroad == "Omicron","Omicron","Other")) %>% distinct()
table1_for_chisq %>% ungroup() %>% infer::chisq_test(ReboundLog ~ LineageBroadAlt)

table2_for_chisq <- end_of_infection %>% select(PersonID, CumulativeInfectionNumber, VaccStatus, ReboundLog) %>% mutate(VaccStatus = ifelse(VaccStatus == "Boosted","Boosted","Not boosted")) %>% distinct() 
table2_for_chisq %>% ungroup() %>% infer::chisq_test(ReboundLog ~ VaccStatus)

## Stat 3 -- how long after first being Ct > 30 do we keep testing people?
## Find the first day with CtT1 > 30 after the peak

## Rebounds are classified as having two consecutive Ct < 30 after being "cleared"
## 2 days of Ct < 30 -- easiest rebound
rebounds1 <- end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & lag(CtT1,2) < 30, DaysExtra > 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber,DaysSinceDetection) %>% distinct() %>% 
    mutate(Rebound="3 x Ct<30") %>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 -- second easiest rebound
rebounds2 <-  end_of_infection_blank %>% filter(CtT1 < 29 & lag(CtT1, 1) < 29 & lag(CtT1, 2) < 29 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="3 x Ct<29")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 and 2 Ct difference between each consecutive measurement -- difficult rebound
rebounds3 <-  end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & lag(CtT1, 2) < 30 & DaysExtra >= 0 & (CtT1 < (lag(CtT1,1) + 1)) & (lag(CtT1,1) < (lag(CtT1,2) + 1))) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="3 x Ct<30 and decrease by 2+ cycles")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 25 -- most stringent rebound
rebounds4 <-  end_of_infection_blank %>% filter(CtT1 < 25 & lag(CtT1, 1) < 25  & lag(CtT1, 2) < 25 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="3 x Ct<25")%>%
    rename(ReboundDate = DaysSinceDetection)
rebounds <- bind_rows(rebounds1, rebounds3, rebounds4)
rebounds$Rebound <- factor(rebounds$Rebound, levels=c("3 x Ct<30","3 x Ct<29", "3 x Ct<30 and decrease by 2+ cycles","3 x Ct<25"))
rebounds$r_numeric <- as.numeric(rebounds$Rebound)
end_of_infection_stringent <- end_of_infection_blank %>% left_join(rebounds)
end_of_infection_stringent <- end_of_infection_stringent %>% mutate(ReboundLog = ifelse(!is.na(Rebound),TRUE,FALSE))
end_of_infection_stringent1 <- end_of_infection_stringent %>% group_by(PersonID, CumulativeInfectionNumber) %>% filter(r_numeric == max(r_numeric))


## Stat 4 -- how long after first being Ct > 30 do we keep testing people?
## Find the first day with CtT1 > 30 after the peak
## Rebounds are classified as having two consecutive Ct < 30 after being "cleared"
## 2 days of Ct < 30 -- easiest rebound
rebounds1 <- end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & lag(CtT1,2) < 30 & lag(CtT1, 3) < 30 & DaysExtra > 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber,DaysSinceDetection) %>% distinct() %>% 
    mutate(Rebound="4 x Ct<30") %>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 -- second easiest rebound
rebounds2 <-  end_of_infection_blank %>% filter(CtT1 < 29 & lag(CtT1, 1) < 29 & lag(CtT1, 2) < 29 & lag(CtT1, 3) < 29 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="4 x Ct<29")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 29 and 2 Ct difference between each consecutive measurement -- difficult rebound
rebounds3 <-  end_of_infection_blank %>% filter(CtT1 < 30 & lag(CtT1, 1) < 30 & lag(CtT1, 2) < 30 & lag(CtT1, 3) < 30 & 
                                              DaysExtra >= 0 & 
                                              (CtT1 < (lag(CtT1,1) + 1)) & (lag(CtT1,1) < (lag(CtT1,2) + 1)) & (lag(CtT1,2) < (lag(CtT1,3) + 1))) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="4 x Ct<30 and decrease by 2+ cycles")%>%
    rename(ReboundDate = DaysSinceDetection)
## 2 days of Ct < 25 -- most stringent rebound
rebounds4 <-  end_of_infection_blank %>% filter(CtT1 < 25 & lag(CtT1, 1) < 25  & lag(CtT1, 2) < 25 & lag(CtT1, 3) < 25 & DaysExtra >= 0) %>% 
    group_by(PersonID, CumulativeInfectionNumber ) %>% filter(TestDate == min(TestDate)) %>%
    select(PersonID, CumulativeInfectionNumber, DaysSinceDetection) %>% distinct() %>% mutate(Rebound="4 x Ct<25")%>%
    rename(ReboundDate = DaysSinceDetection)
rebounds <- bind_rows(rebounds1, rebounds3, rebounds4)
rebounds$Rebound <- factor(rebounds$Rebound, levels=c("4 x Ct<30","4 x Ct<29", "4 x Ct<30 and decrease by 2+ cycles","4 x Ct<25"))
rebounds$r_numeric <- as.numeric(rebounds$Rebound)
end_of_infection_more_stringent <- end_of_infection_blank %>% left_join(rebounds)
end_of_infection_more_stringent <- end_of_infection_more_stringent %>% mutate(ReboundLog = ifelse(!is.na(Rebound),TRUE,FALSE))
end_of_infection_more_stringent1 <- end_of_infection_more_stringent %>% group_by(PersonID, CumulativeInfectionNumber) %>% filter(r_numeric == max(r_numeric))


## How many rebounds out of possible detections?
part1 <- end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber, Rebound) %>% distinct() %>% 
    mutate(Rebound=ifelse(is.na(Rebound),"No rebound",as.character(Rebound))) %>%
    group_by(Rebound) %>% tally() %>%
    mutate(Total=end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber) %>% distinct() %>% tally() %>% rename(Total=n) %>% pull(Total)) %>%
    mutate(Prop=n/Total) %>% 
    mutate(`Initial clearance: consecutive days with Ct≥30 following a Ct<30`=paste0("≥",days_clearance),
           `Rebound: Subsequent consecutive days with Ct < 30`="≥2") %>%
    rename(Rebounds=n) %>%
    mutate(Percentage=paste0(signif(Prop*100,3),"%")) %>%
    select(-Prop) %>%
    filter(Rebound != "No rebound")

## How many rebounds out of possible detections?
part2 <- end_of_infection_stringent %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber, Rebound) %>% distinct() %>% 
    mutate(Rebound=ifelse(is.na(Rebound),"No rebound",as.character(Rebound))) %>%
    group_by(Rebound) %>% tally() %>%
    mutate(Total=end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber) %>% distinct() %>% tally() %>% rename(Total=n) %>% pull(Total)) %>%
    mutate(Prop=n/Total) %>% 
    mutate(`Initial clearance: consecutive days with Ct≥30 following a Ct<30`=paste0("≥",days_clearance),
           `Rebound: Subsequent consecutive days with Ct < 30`="≥3") %>%
    rename(Rebounds=n) %>%
    mutate(Percentage=paste0(signif(Prop*100,3),"%")) %>%
    select(-Prop) %>%
    filter(Rebound != "No rebound")

## How many rebounds out of possible detections?
part3 <- end_of_infection_more_stringent %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber, Rebound) %>% distinct() %>% 
    mutate(Rebound=ifelse(is.na(Rebound),"No rebound",as.character(Rebound))) %>%
    group_by(Rebound) %>% tally() %>%
    mutate(Total=end_of_infection %>% ungroup() %>% select(PersonID,CumulativeInfectionNumber) %>% distinct() %>% tally() %>% rename(Total=n) %>% pull(Total)) %>%
    mutate(Prop=n/Total) %>% 
    mutate(`Initial clearance: consecutive days with Ct≥30 following a Ct<30`=paste0("≥",days_clearance),
           `Rebound: Subsequent consecutive days with Ct < 30`="≥4") %>%
    rename(Rebounds=n) %>%
    mutate(Percentage=paste0(signif(Prop*100,3),"%")) %>%
    select(-Prop) %>%
    filter(Rebound != "No rebound")

res <- bind_rows(part1,part2,part3)


load("plots/traj_plot.RData")
fig1 <- (p_traj + labs(tag="A") +
             theme(plot.tag=element_text(face="bold"))) + (p2 + labs(tag="B")+ theme(plot.tag=element_text(face="bold"),legend.position=c(0.95,0.75),legend.background=element_blank())) + plot_layout(ncol=1,heights=c(3,1))

ggsave(filename="figures/figure1.png",plot=fig1,width=8,height=8,dpi=300)
ggsave(filename="figures/figure1.pdf",plot=fig1,width=8,height=8)


ggsave(paste0("plots/rebounds/rebounds_by_clearance_",days_clearance,"days.png"),p1,height=7,width=7,units="in",dpi=300)
ggsave(paste0("plots/rebounds/rebounds_by_detection_",days_clearance,"days.png"),p2,height=7,width=7,units="in",dpi=300)
ggsave(paste0("plots/rebounds/rebounds_individuals_",days_clearance,".png"),p3,height=10,width=8,units="in",dpi=300)

write_csv(table1,paste0("plots/rebounds/table1_lineage_",days_clearance,".csv"))
write_csv(table2,paste0("plots/rebounds/table2_vaccine_",days_clearance,".csv"))
write_csv(res,paste0("plots/rebounds/all_res_",days_clearance,".csv"))
