library(tidyverse)
library(lazyeval)
library(rstan) 
library(shinystan) 
library(purrr)
library(data.table)
options(mc.cores=parallel::detectCores())
source('code/utils.R')
source("code/set_global_pars.R")

if(run_pars$titeralt){
	titercut <- 400
	ct_dat_refined <- mutate(ct_dat_refined, ThisTiterGroup=TiterGroupAlt)
	titerlow <- c("(-1,13]","(13,400]")
	titerhigh <- c("(400,799]","(799,1e+03]")
} else {
	titercut <- 250
	ct_dat_refined <- mutate(ct_dat_refined, ThisTiterGroup=TiterGroup)
	titerlow <- c("(-1,125]","(125,250]")
	titerhigh <- c("(250,1e+03]")
}

# Define a pared-down dataset for passing to Stan: 
indiv_data <- ct_dat_refined %>% 
	filter(!(RowID %in% run_pars$excluded_rows)) %>%
	mutate(is_exposed=case_when(CumulativeExposureNumber>1~1,TRUE~0)) %>% 
	mutate(titercat=case_when(
		LineageBroad=="Other" & is_exposed==0 ~ 1,
		LineageBroad=="Delta" & is_exposed==1 & ThisTiterGroup%in%titerlow ~ 2,
		LineageBroad=="Delta" & is_exposed==1 & ThisTiterGroup%in%titerhigh ~ 3,
		LineageBroad=="Omicron" & is_exposed==1 & ThisTiterGroup%in%titerlow ~ 4,
		LineageBroad=="Omicron" & is_exposed==1 & ThisTiterGroup%in%titerhigh ~ 5,
		TRUE~NA_real_)) %>% 	
	mutate(vaxstatus=case_when(
		LineageBroad=="Other" & VaccStatus=="Unvaccinated" ~ 1,
		LineageBroad=="Delta" & VaccStatus=="First dose" ~ 2,
		LineageBroad=="Delta" & VaccStatus=="Second dose" ~ 2,
		LineageBroad=="Omicron" & VaccStatus=="First dose" ~ 3,
		LineageBroad=="Omicron" & VaccStatus=="Second dose" ~ 3,
		LineageBroad=="Omicron" & VaccStatus=="Boosted" ~ 4,
		TRUE~NA_real_)) %>% 
	mutate(sympstatus=case_when(
		LineageBroad=="Omicron" & Symptomatic=="No"~1,
		LineageBroad=="Omicron" & Symptomatic=="Yes"~2,
		TRUE~NA_real_)) %>% 
	mutate(titersympstatus=case_when(
		LineageBroad=="Omicron" & ThisTiterGroup%in%titerlow & Symptomatic=="No"~1,
		LineageBroad=="Omicron" & ThisTiterGroup%in%titerlow & Symptomatic=="Yes"~2,
		LineageBroad=="Omicron" & ThisTiterGroup%in%titerhigh & Symptomatic=="No"~3,
		LineageBroad=="Omicron" & ThisTiterGroup%in%titerhigh & Symptomatic=="Yes"~4,
		TRUE~NA_real_)) %>% 
	mutate(vaxsympstatus=case_when(
		LineageBroad=="Omicron" & VaccStatus=="First dose" & Symptomatic=="No"~1,
		LineageBroad=="Omicron" & VaccStatus=="Second dose" & Symptomatic=="No"~1,
		LineageBroad=="Omicron" & VaccStatus=="First dose" & Symptomatic=="Yes"~2,
		LineageBroad=="Omicron" & VaccStatus=="Second dose" & Symptomatic=="Yes"~2,
		LineageBroad=="Omicron" & VaccStatus=="Boosted" & Symptomatic=="No"~3,
		LineageBroad=="Omicron" & VaccStatus=="Boosted" & Symptomatic=="Yes"~4,
		TRUE~NA_real_))

if(run_pars$immunevar=="titer"){
	indiv_data <- filter(indiv_data, !is.na(titercat))
} else if(run_pars$immunevar=="vax"){
	indiv_data <- filter(indiv_data, !is.na(vaxstatus))
} else if(run_pars$immunevar=="symp"){
	indiv_data <- filter(indiv_data, !is.na(sympstatus))
} else if(run_pars$immunevar=="titersymp"){
	indiv_data <- filter(indiv_data, !is.na(titersympstatus))
} else if(run_pars$immunevar=="vaxsymp"){
	indiv_data <- filter(indiv_data, !is.na(vaxsympstatus))
}

indiv_data <- indiv_data	%>% 
	mutate(id=InfectionEvent) %>%
	mutate(t=TestDateIndex) %>%
	mutate(y=CtT1) %>%
	trim_negatives(global_pars) %>% 
	clean_infection_events() %>%
	mutate(id_clean=InfectionEventClean) %>% 
	ungroup() %>% 
	select(id, id_clean, t, y, titercat, vaxstatus, sympstatus, titersympstatus, vaxsympstatus, ThisTiterGroup, VaccStatus, Symptomatic, LineageBroad) %>% 
	arrange(id_clean, t)

# Store the number of infection events we've kept: 
n_indiv <- length(unique(indiv_data$id))

titercat <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, titercat) %>%
	arrange(id) %>%
	pull(titercat)

vaxstatus <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, vaxstatus) %>%
	arrange(id) %>%
	pull(vaxstatus)

sympstatus <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, sympstatus) %>%
	arrange(id) %>%
	pull(sympstatus)

titersympstatus <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, titersympstatus) %>%
	arrange(id) %>%
	pull(titersympstatus)

vaxsympstatus <- indiv_data %>%
	group_by(id) %>%
	slice(1) %>%
	select(id, vaxsympstatus) %>%
	arrange(id) %>%
	pull(vaxsympstatus)

# Useful dataframe for mapping official ID to Stan ID:
id_map <- indiv_data %>% 
	group_by(id) %>%
	summarise(id_clean=first(id_clean)) %>% 
	select(id, id_clean) %>%
	mutate(id_clean=as.character(id_clean))

prior_pars <- list(
	tp_prior=run_pars$tp_prior,
	dp_midpoint=run_pars$dp_midpoint,
	wp_midpoint=run_pars$wp_midpoint,
	wr_midpoint=run_pars$wr_midpoint,
	sigma_prior=run_pars$sigma_prior,
	lambda=run_pars$lambda,
	fpmean=run_pars$fpmean		# so that 90% of mass is <1 and 99% is <2
)	

if(run_pars$immunevar=="titer"){
	prior_pars$immunecat <- titercat
	prior_pars$max_immunecat <- max(titercat)
} else if(run_pars$immunevar=="vax"){
	prior_pars$immunecat <- vaxstatus
	prior_pars$max_immunecat <- max(vaxstatus)
} else if(run_pars$immunevar=="symp"){
	prior_pars$immunecat <- sympstatus
	prior_pars$max_immunecat <- max(sympstatus)
}  else if(run_pars$immunevar=="titersymp"){
	prior_pars$immunecat <- titersympstatus
	prior_pars$max_immunecat <- max(titersympstatus)
}  else if(run_pars$immunevar=="vaxsymp"){
	prior_pars$immunecat <- vaxsympstatus
	prior_pars$max_immunecat <- max(vaxsympstatus)
} else {
	stop("Invalid immunevar in run_pars")
}

