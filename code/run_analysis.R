library(tidyverse) 
library(scales)

source('code/utils.R')
source('code/utils_private.R')
source("code/set_global_pars.R")
source('code/read_ct_dat.R')
source("code/set_run_pars.R")

# # # For fitting to a subset of the data: 
# ct_dat_refined_raw <- ct_dat_refined
# InfectionEventvec <- sample(unique(ct_dat_refined_raw$InfectionEvent),750)
# ct_dat_refined <- ct_dat_refined_raw %>% 
# 	filter(InfectionEvent %in% InfectionEventvec)

for(run_pars_index in 9:9){

	run_pars <- run_pars_list[[run_pars_index]]

	# source("code/summarise_pop.R")
	source("code/fit_posteriors_preamble.R")
	source("code/fit_posteriors.R")

	save(indiv_data, file=paste0("output/run_pars_",run_pars_index,"/indiv_data.RData"))
	save(ct_fit, file=paste0("output/run_pars_",run_pars_index,"/ct_fit.RData"))
	
	source("code/make_figures.R")
	source("code/save_figures.R")

	source("code/clearbigvars.R")

	print(paste0("Done with index ",run_pars_index," ------------"))

}

# launch_shinystan_nonblocking(ct_fit)

# ==============================================================================
# Figure making: 
# ==============================================================================

# for(run_pars_index in 2:2){

# 	run_pars <- run_pars_list[[run_pars_index]]
# 	load(paste0("output/run_pars_",run_pars_index,"/indiv_data.RData"))
# 	load(paste0("output/run_pars_",run_pars_index,"/ct_fit.RData"))
# 	fitlist <- rstan::extract(ct_fit) 
# 	fitlist_chain <- rstan::extract(ct_fit, pars=c("log_dpadj[2]","log_wpadj[2]"), permuted=FALSE) 

# 	source("code/make_figures.R")
# 	source("code/save_figures.R")

# 	source("code/clearbigvars.R")

# 	print(paste0("Done with index ",run_pars_index," ------------"))

# }


