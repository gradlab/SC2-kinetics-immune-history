postdf_dp <- get_long_output(fitlist,c("log_dp_mean")) %>% 
	rename(log_mean=value) %>% 
	select(-name, -id) %>% 
	left_join(
		(get_long_output(fitlist,c("log_dpadj")) %>% 
		rename(log_adj=value) %>% 
		select(-name)),
		by="iteration") %>% 
	select(iteration, id, log_mean, log_adj) %>% 
	mutate(post_dp=exp(log_mean+log_adj)*(run_pars$dp_midpoint)) %>% 
	select(iteration, id, post_dp)

postdf_wp <- get_long_output(fitlist,c("log_wp_mean")) %>% 
	rename(log_mean=value) %>% 
	select(-name, -id) %>% 
	left_join(
		(get_long_output(fitlist,c("log_wpadj")) %>% 
		rename(log_adj=value) %>% 
		select(-name)),
		by="iteration") %>% 
	select(iteration, id, log_mean, log_adj) %>% 
	mutate(post_wp=exp(log_mean+log_adj)*(run_pars$wp_midpoint)) %>% 
	select(iteration, id, post_wp)

postdf_wr <- get_long_output(fitlist,c("log_wr_mean")) %>% 
	rename(log_mean=value) %>% 
	select(-name, -id) %>% 
	left_join(
		(get_long_output(fitlist,c("log_wradj")) %>% 
		rename(log_adj=value) %>% 
		select(-name)),
		by="iteration") %>% 
	select(iteration, id, log_mean, log_adj) %>% 
	mutate(post_wr=exp(log_mean+log_adj)*(run_pars$wr_midpoint)) %>% 
	select(iteration, id, post_wr)


