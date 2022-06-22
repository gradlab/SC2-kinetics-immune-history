# ==============================================================================
# Import
# ==============================================================================

library(tidyverse) 
library(purrr) 
library(lubridate) 

source("code/utils.R")
source("code/utils_private.R")

load('data/data_for_regressions.RData')

# Append an index for the infection event (UPDATE with CumulativeInfectionNumber)
dat_subset_use <- dat_subset_use %>% 
	left_join((dat_subset_use %>% 
		group_by(PersonID, CumulativeInfectionNumber) %>% 
		slice(1) %>% 
		select(PersonID, CumulativeInfectionNumber) %>% 
		ungroup() %>% 
		mutate(InfectionEvent=1:n())),
	by=c("PersonID", "CumulativeInfectionNumber"))

dat_subset_use$LineageBroad <- as.character(dat_subset_use$LineageBroad) 
dat_subset_use$CumulativeExposureNumber <- as.character(dat_subset_use$CumulativeExposureNumber) 
dat_subset_use$CumulativeExposureNumber <- as.numeric(dat_subset_use$CumulativeExposureNumber) 

# ==============================================================================
# Refine
# ==============================================================================

ct_dat_refined <- dat_subset_use %>% 
	# -------------------------------------------------------------------------
	# Convert to Yale Ct values, to facilitiate conversion to GE/ml
	mutate(CtT1= -6.25 + 1.34*CtT1) %>%
	# Cut off anything below the limit of detection 
	mutate(CtT1=case_when(CtT1>global_pars[["lod"]] ~ global_pars[["lod"]], TRUE~CtT1)) %>%
	replace_na(list(CtT1=global_pars[["lod"]])) %>%
	# ---------------------------------------------------------------------
	# Keep only infection events whose Ct reaches below 32:
	group_by(PersonID,InfectionEvent) %>%
	mutate(MinCt=min(CtT1, na.rm=TRUE)) %>%
	filter(MinCt <= 32) %>%
	select(-MinCt) %>%
	ungroup() %>% 
	# ---------------------------------------------------------------------
	# Keep only infection events with at least 3 non-lod Ct values:
	mutate(IsUnderLOD=case_when(CtT1<global_pars[["lod"]]~1,TRUE~0)) %>%
	group_by(PersonID,InfectionEvent) %>%
	mutate(TotalUnderLOD=sum(IsUnderLOD)) %>%
	filter(TotalUnderLOD >= 3) %>%
	select(-IsUnderLOD) %>%
	select(-TotalUnderLOD) %>%
	ungroup() %>% 
	# ---------------------------------------------------------------------
	# Clean infection events:
	clean_infection_events() %>% 
	select(-InfectionEvent) %>% 
	rename(InfectionEvent=InfectionEventClean) %>% 
	ungroup() %>% 
	# ---------------------------------------------------------------------
	# Define test date index by infection event: 
	make_test_date_index_infevent() %>% 
	# ---------------------------------------------------------------------
	# Order: 
	arrange(PersonID, InfectionEvent, TestDateIndex) %>% 
	mutate(RowID=1:n())

# # Check vax status
# ct_dat_refined %>% 
# 	group_by(InfectionEvent) %>% 
# 	summarise(NVaccStatus=n_distinct(VaccStatus)) %>% 
# 	arrange(desc(NVaccStatus))

# ct_dat_refined %>% 
	# filter(InfectionEvent%in%c(247,723,836)) %>% 
	# select(PersonID, CumInfNum=CumulativeInfectionNumber, IE=InfectionEvent, CtT1, VaccStatus, TD=TestDateIndex, DDet=DaysSinceDetection) %>% 
	# print(n=Inf) 

# # check symptoms
# ct_dat_refined %>% 
# 	group_by(InfectionEvent) %>% 
# 	summarise(NSymp=n_distinct(Symptomatic)) %>% 
# 	arrange(desc(NSymp))

# ct_dat_refined %>% 
# 	filter(InfectionEvent%in%c(46,59,83,93,95)) %>% 
# 	select(PersonID, IE=InfectionEvent, CtT1, Symptomatic, TD=TestDateIndex, DDet=DaysSinceDetection) %>% 
# 	print(n=Inf) 

# # check titers
# ct_dat_refined %>% 
# 	group_by(PersonID) %>% 
# 	summarise(NTiter=n_distinct(TiterGroupAlt)) %>% 
# 	arrange(desc(NTiter))

# ct_dat_refined %>% 
# 	filter(InfectionEvent%in%c(379)) %>% 
# 	select(PersonID, CumInfNum=CumulativeInfectionNumber, IE=InfectionEvent, CtT1, TiterGroup, TD=TestDateIndex, DDet=DaysSinceDetection) %>% 
# 	print(n=Inf) 

# # check role
# ct_dat_refined %>% 
# 	group_by(InfectionEvent) %>% 
# 	summarise(NRole=n_distinct(Role)) %>% 
# 	arrange(desc(NRole))

# ct_dat_refined %>% 
# 	filter(InfectionEvent%in%c(836)) %>% 
# 	select(PersonID, CumInfNum=CumulativeInfectionNumber, IE=InfectionEvent, CtT1, Role, TD=TestDateIndex, DDet=DaysSinceDetection) %>% 
# 	print(n=Inf) 

# # check lineage
# ct_dat_refined %>% 
# 	group_by(InfectionEvent) %>% 
# 	summarise(NLin=n_distinct(LineageBroad)) %>% 
# 	arrange(desc(NLin))

# Make vaccination status unique by infection event (vax status nearest to day of detection): 
VaccStatusUnique <- ct_dat_refined %>% 
	mutate(DaysSinceDetection=abs(DaysSinceDetection)) %>% 
	group_by(InfectionEvent) %>% 
	filter(DaysSinceDetection==min(DaysSinceDetection)) %>% 
	select(InfectionEvent, VaccStatus) 
ct_dat_refined <- ct_dat_refined %>% 
	select(-VaccStatus) %>% 
	left_join(VaccStatusUnique, by="InfectionEvent")

# Make titer unique by person (max titer across all infections; changes only one person):
ct_dat_refined <- ct_dat_refined %>% 
	mutate(TiterPrecedence=case_when(
		TiterGroup=="(250,1e+03]"~1,
		TiterGroup=="(125,250]"~2,
		TiterGroup=="(-1,125]"~3,
		TRUE~4
		)) %>% 
	group_by(PersonID) %>% 
	mutate(MaxTiterPrecedence=min(TiterPrecedence)) %>% 
	mutate(TiterGroup=case_when(
		MaxTiterPrecedence==1~"(250,1e+03]",
		MaxTiterPrecedence==2~"(125,250]",
		MaxTiterPrecedence==3~"(-1,125]",
		MaxTiterPrecedence==4~NA_character_
		)) %>% 
	select(-TiterPrecedence, -MaxTiterPrecedence)

ct_dat_refined <- ct_dat_refined %>% 
	mutate(TiterPrecedence=case_when(
		TiterGroupAlt=="(799,1e+03]"~1,
		TiterGroupAlt=="(400,799]"~2,
		TiterGroupAlt=="(13,400]"~3,
		TiterGroupAlt=="(-1,13]"~4,
		TRUE~5
		)) %>% 
	group_by(PersonID) %>% 
	mutate(MaxTiterPrecedence=min(TiterPrecedence)) %>% 
	mutate(TiterGroupAlt=case_when(
		MaxTiterPrecedence==1~"(799,1e+03]",
		MaxTiterPrecedence==2~"(400,799]",
		MaxTiterPrecedence==3~"(13,400]",
		MaxTiterPrecedence==4~"(-1,13]",
		MaxTiterPrecedence==5~NA_character_
		)) %>% 
	select(-TiterPrecedence, -MaxTiterPrecedence)

# Make symptom status unique by infection event, preferring symptoms:
ct_dat_refined <- ct_dat_refined %>% 
	mutate(SymptomaticPrecedence=case_when(
		Symptomatic=="Yes"~1,
		Symptomatic=="No"~2,
		TRUE~3
		)) %>% 
	group_by(InfectionEvent) %>% 
	mutate(MaxSymptomaticPrecedence=min(SymptomaticPrecedence)) %>% 
	mutate(Symptomatic=case_when(
		MaxSymptomaticPrecedence==1~"Yes",
		MaxSymptomaticPrecedence==2~"No",
		MaxSymptomaticPrecedence==3~NA_character_
		)) %>% 
	select(-SymptomaticPrecedence, -MaxSymptomaticPrecedence)

# ------------------------------------------------------------------------------

exposurelist_rows <- ct_dat_refined %>% 
	filter(TiterSensitivity) %>% 
	pull(RowID) 
