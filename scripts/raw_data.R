######################################################
## SCRIPT 1 -- summary of data
######################################################
library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)

setwd("~/Documents/GitHub/SC2-kinetics-immune-history/")

lineage_colors <- c("Delta"="red3","Omicron"="blue","Other"="black","None"="grey70")

## Set to FALSE if you have the original data files
anonymous <- TRUE

## Read in cleaned data
dat <- read_csv("data/ct_data_cleaned.csv")
repeat_dat <- read_csv("data/ct_data_cleaned_repeats.csv")

## Number of individuals
bind_rows(dat, repeat_dat) %>% select(PersonID) %>% distinct() %>% nrow()

## Number of tests
bind_rows(dat, repeat_dat) %>% filter(TestResult %in% c("Negative","Inconclusive","Positive")) %>% group_by(TestResult) %>% tally()
repeat_dat %>% filter(TestResult %in% c("Negative","Inconclusive","Positive")) %>% group_by(TestResult) %>% tally()

## Age distribution
if(!anonymous){
    dat %>% select(PersonID, Age) %>% distinct() %>% group_by(PersonID) %>% filter(Age == min(Age,na.rm=TRUE)) %>% ungroup() %>% drop_na() %>% filter(!is.na(Age)) %>% summarize(mean_age=mean(Age),lower_95=quantile(Age,0.025),upper_95=quantile(Age,0.975))
}

## Vaccination status
dat %>% group_by(PersonID) %>% filter(CumulativeInfectionNumber == max(CumulativeInfectionNumber), TestDate == max(TestDate)) %>%
    group_by(VaccStatus) %>% tally()


dat %>% group_by(PersonID) %>% filter(CumulativeInfectionNumber == max(CumulativeInfectionNumber), TestDate == max(TestDate), VaccStatus != "No record") %>% 
    ungroup() %>% tally()



## Date range of actual tests
if(!anonymous){
    ## Role distribution
    role_key <- c("Player"="Player", "Staff"="Staff", "Misc"="Misc", "Team Staff"="Staff", "Household"="Household", 
                  "Game Day Assistant"="Misc", 
                  "Contractor/Consultant"="Contractor/Consultant", "Vendor"="Vendor", "League Protocol Compliance Officer"="Misc", 
                  "vendor"="Vendor", "Referee"="Referee", "NBA Staff"="Staff", "Referee Trainer"="Referee", "All-Star"="Player", 
                  "Referee Household"="Household", "Guest"="Guest", "Drug Free Sport"="Misc", "Combine Player"="Player", 
                  "guest"="Guest", "Summer_league"="Misc", "Broadcast"="Misc")
    dat$Role <- role_key[dat$Role]
    dat %>% filter(TestResult %in% c("Negative","Positive","Inconclusive")) %>% 
        summarize(min_date=min(TestDate,na.rm=TRUE),max_date=max(TestDate,na.rm=TRUE))

    min_date <- dat %>% filter(TestResult %in% c("Negative","Positive","Inconclusive")) %>% summarize(x=min(TestDate,na.rm=TRUE)) %>% pull(x)

    dat %>% filter(PersonID %in% multiple_roles) %>% select(PersonID, Role) %>% distinct()
    
    ## Correct some roles
    dat[dat$PersonID == 969,"Role"] <- "Referee"
    dat[dat$PersonID == 1031,"Role"] <- "Staff"
    dat[dat$PersonID == 1162,"Role"] <- "Player"
    dat[dat$PersonID == 1794,"Role"] <- "Staff"
    dat[dat$PersonID == 1815,"Role"] <- "Staff"
    dat[dat$PersonID == 1822,"Role"] <- "Player"
    dat[dat$PersonID == 1853,"Role"] <- "Household"
    dat[dat$PersonID == 1860,"Role"] <- "Household"
}

bind_rows(repeat_dat, dat) %>% filter(TestResult %in% c("Negative","Inconclusive","Positive")) %>% nrow()
bind_rows(repeat_dat, dat)  %>% filter(TestResult %in% c("Negative","Inconclusive","Positive")) %>% group_by(TestResult) %>% tally()

## Number of positive tests
dat %>% filter(TestResult %in% c("Positive","Inconclusive")) %>% distinct() %>% nrow()
dat %>% group_by(TestResult) %>% tally()

repeat_dat %>% filter(TestResult %in% c("Positive","Inconclusive")) %>% distinct() %>% nrow()
## Number of negative tests
dat %>% filter(TestResult %in% c("Negative")) %>% distinct() %>% nrow()
repeat_dat %>% filter(TestResult %in% c("Negative")) %>% distinct() %>% nrow()

dat %>% filter(NewInfectionIdentified == 1) %>% nrow()
dat %>% filter(NewInfectionIdentified == 1) %>% select(PersonID) %>% distinct() %>% nrow()


## Label Broad Omicron and Delta lineages

lineages <- unique(dat$Lineage)
delta_lineages <-lineages[lineages %like% "AY" | lineages == "B.1.617.2"]
omicron_lineages <- c("Suspected Omicron",lineages[lineages %like% "BA.1"])
ba2 <- lineages[lineages %like% "BA.2"]

dat <- dat %>% mutate(LineageBroadConfirmed = ifelse(Lineage %in% delta_lineages, "Delta",
                                            ifelse(Lineage %in% omicron_lineages, "Omicron",
                                                   ifelse(Lineage %in% ba2,"BA.2",
                                                          ifelse(Lineage == "None", "None","Other")))))
dat %>% filter(NewInfectionIdentified == 1) %>% group_by(LineageBroadConfirmed) %>% tally()

## Number of infections
dat %>% group_by(PersonID) %>% filter(CumulativeInfectionNumber == max(CumulativeInfectionNumber)) %>% select(PersonID, CumulativeInfectionNumber) %>% distinct()  %>% group_by(CumulativeInfectionNumber) %>% tally()

## Number of exposures
dat %>% group_by(PersonID) %>% filter(CumulativeExposureNumber == max(CumulativeExposureNumber)) %>% select(PersonID, CumulativeExposureNumber) %>% distinct()  %>% group_by(CumulativeExposureNumber) %>% tally()

## Vaccination statuses by the end
dat %>% group_by(PersonID) %>% filter(CumulativeInfectionNumber == max(CumulativeInfectionNumber), TestDate == max(TestDate)) %>% group_by(VaccStatus) %>% tally()
 

## Number of infections detected outside of the main protocol
dat %>% filter(NewInfectionIdentified == 1, TestResult %in% c("Positive (no CT)","Positive (external)")) %>% nrow

## Detection speed
dat_detections <- dat %>% 
    ungroup() %>%
    filter(NewInfectionIdentified == 1) %>%
    group_by(PersonID, CumulativeInfectionNumber) %>%
    mutate(DetectionSpeed = ifelse(DaysSinceNegative >= 2 | is.na(DaysSinceNegative),
                                   ">=2 days",ifelse(DaysSinceNegative < 2, "<=1 days",NA))) %>% 
    mutate(DetectionSpeed=ifelse(is.na(DetectionSpeed),
                                 ">=2 days",DetectionSpeed)) %>%
    select(PersonID, TestDate, DetectionSpeed, DaysSinceNegative)

## Tally of detection speeds
dat_detections %>% group_by(DetectionSpeed) %>% tally()


## Imputations
dat_indivs <- dat %>% filter(NewInfectionIdentified > 0) %>% 
    select(PersonID, CumulativeInfectionNumber, Lineage,LineageBroad,VaccStatus,CumulativeExposureNumber,Titer,TestDate)


# Lineages ----------------------------------------------------------------
## Note this cannot run with the anonymized dataset
if(!anonymous){
    proportions <- dat_indivs %>% filter(LineageBroad != "None") %>% mutate(TestDate=lubridate::round_date(TestDate, "week")) %>%  
        group_by(TestDate,LineageBroad) %>% tally() %>% 
        left_join(dat_indivs %>% mutate(TestDate=lubridate::round_date(TestDate, "week")) %>%
                      filter(LineageBroad != "None") %>%  group_by(TestDate) %>% summarize(N=n())) %>%
        mutate(prop=n/N)
    
    
    proportions_day <- dat_indivs %>% filter(LineageBroad != "None") %>%
        group_by(TestDate,LineageBroad) %>% tally() %>% 
        left_join(dat_indivs %>% filter(LineageBroad != "None") %>%  group_by(TestDate) %>% summarize(N=n())) %>%
        mutate(prop=n/N)
    
    ## Assess the end of the Alpha era
    ## End of the pre-delta era -- 29th May
    proportions_day %>% filter(LineageBroad == "Other") %>% filter(n > 0) %>% tail(10)
    ## See 1 "other" July 13th and 1 July 18th, then the rest is all Delta
    proportions_day %>% filter(LineageBroad == "Delta") %>% filter(n > 0) %>% filter(TestDate <="2021-08-01") %>% tail(31)
    ## How many "None" in the time period when Delta was circulating but before Delta was clearly dominant?
    dat %>% filter(NewInfectionIdentified == 1, LineageBroad == "None", TestDate >= "2021-05-29",TestDate <= "2021-07-18", CtT1 < 40, CtT1 > 0) %>% nrow()
    
    ## Assess the end of the Delta era
    ## End of the pre-Omicron era -- 8th December
    proportions_day %>% filter(LineageBroad == "Delta") %>% filter(n > 0) %>% tail(31) %>% View
    ## By 9th December, all Omicron. Omicron first detected 3rd December
    proportions_day %>% filter(LineageBroad == "Omicron") %>% filter(n > 0) %>% filter(TestDate <="2021-12-12") %>% tail(10)
    ## How many "None" in the time period when Omicron was circulating but before Omicron was clearly dominant?
    dat %>% filter(NewInfectionIdentified == 1, LineageBroad == "None", TestDate >= "2021-12-03",TestDate <= "2021-12-08", CtT1 < 40, CtT1 > 0) %>% nrow()
    
    dat %>% filter(NewInfectionIdentified == 1, LineageBroad == "None", TestDate >= "2021-12-03") %>% nrow()
    
    p_lineage_counts <- ggplot(dat_indivs %>% mutate(TestDate=lubridate::round_date(TestDate, "week")) %>%   
                                   group_by(TestDate,LineageBroad) %>% tally()) +
        annotate("rect",xmin=as.Date("2020-01-01"),xmax=as.Date("2021-05-20"),ymin=-10,ymax=700,fill="black",alpha=0.1) +
        annotate("rect",xmin=as.Date("2021-07-22"),xmax=as.Date("2021-11-25"),ymin=-10,ymax=700,fill="red3",alpha=0.25) +
        annotate("rect",xmin=as.Date("2021-12-09"),xmax=as.Date("2022-03-01"),ymin=-10,ymax=700,fill="blue",alpha=0.25) +
        geom_bar(aes(x=TestDate,y=n,fill=LineageBroad),col="black",stat="identity",size=0.25,alpha=0.8)  +
        ylab("Count") + xlab("Detection date") +
        scale_x_date(breaks="1 month") +
        scale_y_continuous(expand=c(0,0),breaks=seq(0,600,by=50)) +
        coord_cartesian(ylim=c(0,625),xlim=range(dat$TestDate)) +
        scale_fill_manual(name="Lineage",values=lineage_colors) +
        theme_classic() +
        theme(legend.position=c(0.2,0.6),
              axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank(),
              plot.tag=element_text(face="bold"),
              legend.background = element_rect(fill="white",color="black")) +
        geom_vline(xintercept=as.Date(c("2021-05-20","2021-07-22")),size=1,linetype="dashed")+
        geom_vline(xintercept=as.Date(c("2021-11-25","2021-12-09")),size=1,linetype="dashed") +
        labs(tag="A")
    
    p_lineage_sequenced_props <- ggplot(proportions) +
        geom_bar(aes(x=TestDate,y=prop,fill=LineageBroad),stat="identity",col="black",size=0.25,alpha=0.8)  +
        ylab("Proportion") + xlab("Detection date") +
        scale_x_date(breaks="1 month",limits=range(dat$TestDate)) +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(name="Lineage",values=lineage_colors) +
        theme_classic() +
        theme(legend.position="none",
              axis.text.x=element_text(angle=45,hjust=1),
              plot.tag=element_text(face="bold"),) +
        geom_vline(xintercept=as.Date(c("2021-05-20","2021-07-22")),size=1,linetype="dashed")+
        geom_vline(xintercept=as.Date(c("2021-11-25","2021-12-09")),size=1,linetype="dashed")+
        labs(tag="B")
    
    
    p_main <- p_lineage_counts/p_lineage_sequenced_props
    ggsave("figures/supplement/lineage_proportions.png",p_main,height=6,width=7)
}
    
## Number of antibody titers
dat %>% filter(!is.na(Titer)) %>% select(PersonID) %>% distinct() %>% nrow()
if(!anonymous){
    dat %>% filter(!is.na(Titer)) %>% summarize(range(TestDate))
    dat %>% filter(!is.na(Titer)) %>% ggplot() + geom_histogram(aes(x=TestDate)) + scale_x_date(breaks="1 week") + theme(axis.text.x=element_text(angle=45,hjust=1))
}
dat %>% filter(!is.na(Titer)) %>% group_by(PersonID) %>% tally() %>% filter(n > 1) %>% nrow()


## Mark when titer was actually measured
dat <- dat %>% mutate(TiterMeasured = !is.na(Titer))

if(!anonymous){
    ## Add titer date to date titer was measured. Then, for each individual, carry this date downward. Note we'll have some individuals with multiple titers
    dat <- dat %>% mutate(TiterDate = ifelse(TiterMeasured == 1, as.Date(TestDate), NA)) %>%
        group_by(PersonID) %>% fill(TiterDate,.direction="down")
    dat$TiterDate <- as.Date(dat$TiterDate,origin="1970-01-01")

## Time between second dose and titer draw
p_lag_1 <- dat %>% select(PersonID,VaccineDose2Date,BoosterDate,TiterDate) %>% distinct() %>% mutate(lag = TiterDate - VaccineDose2Date) %>% 
    ggplot() + geom_histogram(aes(x=lag),fill='white',col="grey40",binwidth=10) + theme_classic() +
    geom_vline(data=.%>% ungroup() %>% summarize(y=median(lag,na.rm=TRUE)),aes(xintercept=y),linetype="dashed") +
    xlab("Time between second vaccine dose and titer draw (days)") +
    scale_x_continuous(limits=c(-125,310)) +
    labs(tag="A")+
    ylab("Count") +
    theme(plot.tag=element_text(face="bold"))
## Time between booster dose and titer draw
p_lag_2 <- dat %>% select(PersonID,VaccineDose2Date,BoosterDate,TiterDate) %>% distinct() %>% mutate(lag = TiterDate - BoosterDate) %>% 
    ggplot() + geom_histogram(aes(x=lag),fill='white',col="grey40",binwidth=10) +
    geom_vline(data=.%>% ungroup() %>% summarize(y=median(lag,na.rm=TRUE)),aes(xintercept=y),linetype="dashed") +
    ylab("Count") +
    scale_x_continuous(limits=c(-125,310)) +
    theme_classic() +
    xlab("Time between booster dose and titer draw (days)") +
    ylab("Count") +
    labs(tag="B") +
    theme(plot.tag=element_text(face="bold"))

supp_p_lags <- p_lag_1/p_lag_2
ggsave("figures/supplement/vaccine_titer_lags.png",supp_p_lags,height=8,width=7)

## How many individuals had an infection inbetween 2nd dose and titer?
dat %>% filter(NewInfectionIdentified== TRUE) %>% filter(TestDate >= VaccineDose2Date, TestDate <= TiterDate) %>% group_by(LineageBroad) %>% tally()
dat %>% filter(NewInfectionIdentified== TRUE) %>% filter(TestDate <= BoosterDate, TestDate >= TiterDate) %>% group_by(LineageBroad) %>% tally()

}

