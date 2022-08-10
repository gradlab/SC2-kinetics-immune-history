load("outputs/titer_models_age/all_freq_1.RData")
fit1 <- fit
load("outputs/titer_models_age/all_infreq_1.RData")
fit2 <- fit

newdata <- expand_grid(DaysSinceDetection=seq(0,20,by=0.1), LineageBroad_BoostTiterGroup=unique(fit1$data$LineageBroad_BoostTiterGroup),
                       AgeGroup=unique(fit1$data$AgeGroup))
titer_freq_draws <- fit1 %>% epred_draws(newdata=newdata)
titer_infreq_draws <- fit2 %>% epred_draws(newdata=newdata)



titer_p_dat <- bind_rows(titer_freq_draws %>% mutate(Protocol="Frequent testing"), 
                         titer_infreq_draws %>% mutate(Protocol = "Delayed detection")) %>%
    group_by(Protocol, DaysSinceDetection,LineageBroad_BoostTiterGroup,AgeGroup) %>% 
    summarize(lower=quantile(.epred,0.025),upper=quantile(.epred,0.975),med=median(.epred)) %>%
    mutate(LineageBroad = ifelse(LineageBroad_BoostTiterGroup %like% "Other","Other",
                                 ifelse(LineageBroad_BoostTiterGroup %like% "Delta","Delta","Omicron")))

titer_p_dat$Protocol <- factor(vacclineage_p_dat$Protocol,levels=c("Frequent testing","Delayed detection"))

tmp_dat1 <- titer_p_dat %>% filter(LineageBroad == "Omicron") 


boosttitergroup_key <- c("OmicronHighBoosted" = "Omicron: Boosted >250 AU/ml",
                         "OmicronLowBoosted" = "Omicron: Boosted ≤250 AU/ml",
                         "OmicronHighNotBoosted"="Omicron: Not Boosted >250 AU/ml",
                         "OmicronLowNotBoosted"="Omicron: Not Boosted ≤250 AU/ml")

tmp_dat1$`Immune status` <- boosttitergroup_key[as.character(tmp_dat1$LineageBroad_BoostTiterGroup)]
tmp_dat1$`Detection group` <- factor(tmp_dat1$`Protocol`, levels=c("Frequent testing","Delayed detection"))
tmp_dat1$`Immune status` <- factor(tmp_dat1$`Immune status`, levels=c("Omicron: Boosted ≤250 AU/ml","Omicron: Boosted >250 AU/ml",
                                                                      "Omicron: Not Boosted ≤250 AU/ml","Omicron: Not Boosted >250 AU/ml"))

age_key <- c("(0,30]"="<30","(30,50]"="30-50","(50,100]"=">50")
tmp_dat1$AgeGroup <- age_key[as.character(tmp_dat1$AgeGroup)]
tmp_dat1$AgeGroup <- factor(tmp_dat1$AgeGroup, levels=c("<30","30-50",">50"))
p_titerlineage <-  ggplot(tmp_dat1, 
                         aes(x=DaysSinceDetection,ymin=lower,ymax=upper,fill=`Immune status`,group=interaction(`Immune status`,AgeGroup),y=med),col="None") +
    facet_grid(AgeGroup~Protocol) +
    geom_ribbon(alpha=0.5) +
    geom_line(aes(col=`Immune status`))+
    scale_y_continuous(limits=c(0,1),expand=c(0,0), breaks=seq(0,1,by=0.2)) +
    scale_x_continuous(limits=c(0,20),breaks=seq(0,20,by=5)) +
    labs(y="Probability of Ct value <30",x="Days since detection",fill="Vaccination status",color="Vaccination status") +
    theme_classic() +
    geom_vline(xintercept=5,linetype="dashed") +
    geom_hline(yintercept=0.05,linetype="dashed") +
    scale_color_viridis_d(option="D") +
    scale_fill_viridis_d(option="D") +
    theme(legend.position=c(0.85,0.85),
          axis.text=element_text(size=8),axis.title=element_text(size=8),
          legend.title=element_text(size=6),legend.text=element_text(size=6),
          strip.background=element_blank(),
          strip.text=element_text(face="bold"),
          plot.background = element_rect(fill="white",color="white"))
p_titerlineage
