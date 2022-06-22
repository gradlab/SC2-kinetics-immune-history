# ==============================================================================
# Import and key definitions
# ==============================================================================

library(tidyverse) 
library(beeswarm) 
library(ggbeeswarm) 

if(run_pars$immunevar=="titer"){
	immunelabs <- c("Other: Unexposed", paste0("Delta: Exposed, <",titercut), paste0("Delta: Exposed, >",titercut), paste0("Omicron: Exposed, ",titercut), paste0("Omicron: Exposed, >",titercut))
	immunelabsdf <- tibble(label=immunelabs, id=1:length(immunelabs)) %>% 
		mutate(label=factor(label,levels=immunelabs))
	immunetitle <- "titer"
	immunecolors <- c("black","tomato","red3","dodgerblue","blue")
	immunetags <- c("A)","B)","C)","D)","E)")
	indiv_data$immunecat <- indiv_data$titercat
} else if(run_pars$immunevar=="vax"){
	immunelabs <- c("Other: Unvaccinated","Delta: 1-2 doses","Omicron: 1-2 doses", "Omicron: Boosted")
	immunelabsdf <- tibble(label=immunelabs, id=1:length(immunelabs)) %>% 
		mutate(label=factor(label,levels=immunelabs))
	immunetitle <- "vaccination status"	
	immunecolors <- c("black","orange3","mediumpurple1","purple3")
	immunetags <- c("A)","B)","C)","D)")
	indiv_data$immunecat <- indiv_data$vaxstatus
}  else if(run_pars$immunevar=="symp"){
	immunelabs <- c("Omicron: No symptoms","Omicron: Symptoms")
	immunelabsdf <- tibble(label=immunelabs, id=1:length(immunelabs)) %>% 
		mutate(label=factor(label,levels=immunelabs))
	immunetitle <- "symptom status"	
	immunecolors <- c("mediumpurple1","purple3")
	immunetags <- c("A)","B)")
	indiv_data$immunecat <- indiv_data$sympstatus
}   else if(run_pars$immunevar=="titersymp"){
	immunelabs <- c("Omicron: Low titer, no symptoms","Omicron: Low titer, symptoms", "Omicron: High tier, no symptoms", "Omicron: High titer, symptoms")
	immunelabsdf <- tibble(label=immunelabs, id=1:length(immunelabs)) %>% 
		mutate(label=factor(label,levels=immunelabs))
	immunetitle <- "symptom status"	
	immunecolors <- c("mediumpurple1","purple3","black","darkblue")
	immunetags <- c("A)","B)","C)","D)")
	indiv_data$immunecat <- indiv_data$titersympstatus
}   else if(run_pars$immunevar=="vaxsymp"){
	immunelabs <- c("Omicron: 1-2 doses, no symptoms","Omicron: 1-2 doses, symptoms", "Omicron: boosted, no symptoms", "Omicron: boosted, symptoms")
	immunelabsdf <- tibble(label=immunelabs, id=1:length(immunelabs)) %>% 
		mutate(label=factor(label,levels=immunelabs))
	immunetitle <- "symptom status"	
	immunecolors <- c("mediumpurple1","purple3","black","darkblue")
	immunetags <- c("A)","B)","C)","D)")
	indiv_data$immunecat <- indiv_data$vaxsympstatus
}  

params_indiv <- get_wide_output(fitlist, c("tp","dp","wp","wr")) %>% 
	left_join((indiv_data %>% 
			group_by(id_clean) %>% 
			slice(1) %>% 
			select(id_clean, titercat, vaxstatus, sympstatus, titersympstatus, vaxsympstatus, immunecat)),
		by=c("id"="id_clean"))

meanvals_indiv <- params_indiv %>% 
	group_by(id) %>% 
	summarise(tp=mean(tp),dp=mean(dp),wp=mean(wp),wr=mean(wr),titercat=first(titercat),vaxstatus=first(vaxstatus),sympstatus=first(sympstatus),titersympstatus=first(titersympstatus),vaxsympstatus=first(vaxsympstatus),immunecat=first(immunecat))

summaryvals_indiv <- params_indiv %>% 
	group_by(id) %>% 
	summarise(
		tp_mean=mean(tp),
		tp_lwr=quantile(tp,0.025),
		tp_upr=quantile(tp,0.975),
		dp_mean=mean(dp),
		dp_lwr=quantile(dp,0.025),
		dp_upr=quantile(dp,0.975),
		wp_mean=mean(wp),
		wp_lwr=quantile(wp,0.025),
		wp_upr=quantile(wp,0.975),
		wr_mean=mean(wr),
		wr_lwr=quantile(wr,0.025),
		wr_upr=quantile(wr,0.975),
		titercat=first(titercat),
		vaxstatus=first(vaxstatus),
		sympstatus=first(sympstatus),
		titersympstatus=first(titersympstatus),
		vaxsympstatus=first(vaxsympstatus),
		immunecat=first(immunecat)) %>% 
	rename(id_clean=id) %>% 
	left_join((indiv_data %>% 
			select(id,id_clean) %>% 
			group_by(id) %>% 
			slice(1)),by="id_clean") %>% 
	left_join((ct_dat_refined %>% 
			group_by(InfectionEvent) %>% 
			slice(1) %>% 
			select(InfectionEvent, PersonID) %>% 
			rename(id=InfectionEvent)),
		by="id")

# ==============================================================================
# Chain plots
# ==============================================================================

chainplot_dpadj <- plot_chains(fitlist_chain, "log_dpadj[2]")
chainplot_wpadj <- plot_chains(fitlist_chain, "log_wpadj[2]")

# ==============================================================================
# Individual trajectories
# ==============================================================================

id_sample <- sample(unique(params_indiv$id), min(length(unique(indiv_data$id)), 100))
iter_sample <- sort(sample(1:max(params_indiv$iteration),50))
fig_indivfits <- indiv_data %>% 
	select(-id) %>% 
	mutate(id=id_clean) %>% 
	filter(id %in% id_sample) %>% 
	ggplot() + 
		geom_point(aes(x=t,y=y), size=0.5)  +
		geom_segment(data=(filter(params_indiv,id%in%id_sample & iteration%in%iter_sample)),
			aes(x=tp-wp, xend=tp, y=40, yend=40-dp),
			size=0.1, alpha=0.2) + 
		geom_segment(data=(filter(params_indiv,id%in%id_sample & iteration%in%iter_sample)),
			aes(x=tp, xend=tp+wr, y=40-dp, yend=40),
			size=0.1, alpha=0.2) + 
		theme_minimal() + 
		scale_y_reverse() + 
		facet_wrap(~id)

# ==============================================================================
# Posterior distributions for the adjustments
# ==============================================================================

fig_dpadjhist <- get_long_output(fitlist, c("log_dpadj")) %>% 
	filter(id>1) %>% 
	left_join(immunelabsdf, by="id") %>% 
	mutate(value=exp(value)) %>% 
	ggplot() + 
		geom_histogram(aes(x=value, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=value, col=label), adjust=2) + 
		geom_vline(aes(xintercept=1)) + 
		scale_color_manual(values=immunecolors[-1]) + 
		scale_fill_manual(values=immunecolors[-1]) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title="dp adjustment factor", subtitle=paste0("by ",immunetitle))

fig_wpadjhist <- get_long_output(fitlist, c("log_wpadj")) %>% 
	filter(id>1) %>% 
	left_join(immunelabsdf, by="id") %>% 
	mutate(value=exp(value)) %>% 
	ggplot() + 
		geom_histogram(aes(x=value, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=value, col=label), adjust=2) + 
		geom_vline(aes(xintercept=1)) + 
		scale_color_manual(values=immunecolors[-1]) + 
		scale_fill_manual(values=immunecolors[-1]) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title="wp adjustment factor", subtitle=paste0("by ",immunetitle))

fig_wradjhist <- get_long_output(fitlist, c("log_wradj")) %>% 
	filter(id>1) %>% 
	left_join(immunelabsdf, by="id") %>% 
	mutate(value=exp(value)) %>% 
	ggplot() + 
		geom_histogram(aes(x=value, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=value, col=label), adjust=2) + 
		geom_vline(aes(xintercept=1)) + 
		scale_color_manual(values=immunecolors[-1]) + 
		scale_fill_manual(values=immunecolors[-1]) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title="wr adjustment factor", subtitle=paste0("by ",immunetitle))

# ==============================================================================
# Posterior distributions for the parameters themselves
# ==============================================================================

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

postdf_rp <- inner_join(postdf_dp, postdf_wp, by=c("iteration","id")) %>% 
	mutate(post_rp=post_dp/post_wp) %>% 
	select(iteration, id, post_rp)

postdf_rr <- inner_join(postdf_dp, postdf_wr, by=c("iteration","id")) %>% 
	mutate(post_rr=post_dp/post_wr) %>% 
	select(iteration, id, post_rr)


fig_dphist <- postdf_dp %>% 
	left_join(immunelabsdf, by="id") %>% 
	ggplot() + 
		geom_histogram(aes(x=global_pars[["lod"]]-post_dp, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=global_pars[["lod"]]-post_dp, col=label), adjust=2) +
		scale_x_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) + 
		scale_color_manual(values=immunecolors) + 
		scale_fill_manual(values=immunecolors) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title=paste0("Peak Ct by ",immunetitle), x="Ct", y="Density")
fig_dphist_nolegend <- fig_dphist + theme(legend.position="none", plot.title=element_blank())

fig_wphist <- postdf_wp %>% 
	left_join(immunelabsdf, by="id") %>% 
	ggplot() + 
		geom_histogram(aes(x=post_wp, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=post_wp, col=label), adjust=2) +
		scale_color_manual(values=immunecolors) + 
		scale_fill_manual(values=immunecolors) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title=paste0("Proliferation time by ",immunetitle), x="Proliferation time (days)", y="Density")

fig_wrhist <- postdf_wr %>% 
	left_join(immunelabsdf, by="id") %>% 
	ggplot() + 
		geom_histogram(aes(x=post_wr, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=post_wr, col=label), adjust=2) +
		scale_color_manual(values=immunecolors) + 
		scale_fill_manual(values=immunecolors) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title=paste0("Clearance time by ",immunetitle), x="Clearance time (days)", y="Density")

fig_rphist <- postdf_rp %>% 
	left_join(immunelabsdf, by="id") %>% 
	ggplot() + 
		geom_histogram(aes(x=post_rp, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=post_rp, col=label), adjust=2) +
		scale_color_manual(values=immunecolors) + 
		scale_fill_manual(values=immunecolors) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title=paste0("Proliferation rate by ",immunetitle), x="Proliferation rate (Ct/day)", y="Density")

fig_rrhist <- postdf_rr %>% 
	left_join(immunelabsdf, by="id") %>% 
	ggplot() + 
		geom_histogram(aes(x=post_rr, fill=label, y=..density..), alpha=0.4, position="identity", bins=50) + 
		geom_density(aes(x=post_rr, col=label), adjust=2) +
		scale_color_manual(values=immunecolors) + 
		scale_fill_manual(values=immunecolors) + 
		theme_classic() + 
		theme(text=element_text(size=9), legend.title=element_blank()) + 
		labs(title=paste0("Clearance rate by ",immunetitle), x="Clearance rate (Ct/day)", y="Density")

# ==============================================================================
# Trajectory summaries 
# ==============================================================================

tpdf <- get_long_output(fitlist, c("tp")) %>% 
	group_by(id) %>% 
	summarise(tp=mean(value)) %>%
	rename(id_clean=id)

postdf_overall <- postdf_dp %>% 
	left_join(postdf_wp, by=c("iteration","id")) %>% 
	left_join(postdf_wr, by=c("iteration","id")) 
postdf_overall_summary <- postdf_overall %>% 
	group_by(id) %>% 
	summarise(post_dp_mean=mean(post_dp), 
		post_dp_lwr=quantile(post_dp,0.025), 
		post_dp_upr=quantile(post_dp,0.975), 
		post_wp_mean=mean(post_wp), 
		post_wp_lwr=quantile(post_wp,0.025), 
		post_wp_upr=quantile(post_wp,0.975),
		post_wr_mean=mean(post_wr),
		post_wr_lwr=quantile(post_wr,0.025), 
		post_wr_upr=quantile(post_wr,0.975))

boundinterp <- function(val,xstart,xend,ystart,yend){
	out <- (yend-ystart)/(xend-xstart)*(val-xstart)+ystart
	return(out)
}

boundvalsup <- seq(from=-10,to=0,by=0.02)
boundvalsdown <- seq(from=0,to=20,by=0.02)

boundsup <- lapply(min(postdf_overall$id):max(postdf_overall$id),
	function(y){
		xstart <- -filter(postdf_overall,id==y)$post_wp
		xend <- rep(0,length(xstart))
		ystart <- rep(0, length(xstart))
		yend <- filter(postdf_overall,id==y)$post_dp
		bounds <- (lapply(boundvalsup, boundinterp, xstart, xend, ystart, yend) %>%
			map(~ tibble(lwr=quantile(.,0.025), mean=mean(.), upr=quantile(.,0.975))) %>%
			bind_rows() %>% 
			mutate(t=boundvalsup))
		return(bounds)
	}
	)

boundsdown <- lapply(min(postdf_overall$id):max(postdf_overall$id),
	function(y){
		xend <- filter(postdf_overall,id==y)$post_wr
		xstart <- rep(0,length(xend))
		ystart <- filter(postdf_overall,id==y)$post_dp
		yend <- rep(0, length(xend))
		bounds <- (lapply(boundvalsdown, boundinterp, xstart, xend, ystart, yend) %>%
			map(~ tibble(lwr=quantile(.,0.025), mean=mean(.), upr=quantile(.,0.975))) %>%
			bind_rows() %>% 
			mutate(t=boundvalsdown))
		return(bounds)
	}
	)

fig_viztrajectories <- indiv_data %>% 
	select(id_clean, t, y, immunecat) %>% 
	filter(y<40) %>% 
	left_join(tpdf, by="id_clean") %>%  
	mutate(t=t-tp) %>% 
	split(.$immunecat) %>% 
	imap(~ ggplot(.x, aes(x=t, y=y)) + 
		geom_point(alpha=0.1, size=0.5, col=immunecolors[as.numeric(.y)]) + 
		coord_cartesian(ylim=c(40,10), xlim=c(-12,25), expand=FALSE) + 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
		labs(x="Days since peak viral load (estimated)", y="Ct", tag=immunetags[as.numeric(.y)]) + 
		theme_classic() + 
		theme(text=element_text(size=9))) 

for(indexA in 1:max(postdf_overall$id)){
	fig_viztrajectories[[indexA]] <- fig_viztrajectories[[indexA]] + 
		geom_ribbon(data=boundsup[[indexA]], aes(x=t, y=mean, ymin=global_pars[["lod"]]-lwr, ymax=global_pars[["lod"]]-upr), alpha=0.5, fill=immunecolors[indexA]) + 
		geom_line(data=boundsup[[indexA]], aes(x=t, y=global_pars[["lod"]]-lwr),size=0.1, col=immunecolors[indexA]) +
		geom_line(data=boundsup[[indexA]], aes(x=t, y=global_pars[["lod"]]-upr),size=0.1, col=immunecolors[indexA]) +
		geom_line(data=boundsup[[indexA]], aes(x=t, y=global_pars[["lod"]]-mean), col=immunecolors[indexA]) + 
		geom_ribbon(data=boundsdown[[indexA]], aes(x=t, y=mean, ymin=global_pars[["lod"]]-lwr, ymax=global_pars[["lod"]]-upr), alpha=0.5, fill=immunecolors[indexA]) + 
		geom_line(data=boundsdown[[indexA]], aes(x=t, y=global_pars[["lod"]]-lwr),size=0.1, col=immunecolors[indexA]) +
		geom_line(data=boundsdown[[indexA]], aes(x=t, y=global_pars[["lod"]]-upr),size=0.1, col=immunecolors[indexA]) +
		geom_line(data=boundsdown[[indexA]], aes(x=t, y=global_pars[["lod"]]-mean), col=immunecolors[indexA])
}


# ==============================================================================
# Whiskers: 
# ==============================================================================

pointsize <- 0.8
alphabee <- 0.2
whiskersize <- 0.3
whiskerwidth <- 0.1
fig_dp_whiskers <- meanvals_indiv %>% 
	select(dp, id=immunecat) %>% 
	ggplot(aes(x=factor(id), col=factor(id), y=global_pars[["lod"]]-dp)) + 
		geom_beeswarm(size=pointsize, alpha=alphabee, cex=1.2) + 
		geom_segment(data=postdf_overall_summary, aes(x=id, xend=id, y=global_pars[["lod"]]-post_dp_lwr, yend=global_pars[["lod"]]-post_dp_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=global_pars[["lod"]]-post_dp_lwr, yend=global_pars[["lod"]]-post_dp_lwr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=global_pars[["lod"]]-post_dp_upr, yend=global_pars[["lod"]]-post_dp_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-1.2*whiskerwidth, xend=id+1.2*whiskerwidth, y=global_pars[["lod"]]-post_dp_mean, yend=global_pars[["lod"]]-post_dp_mean), col="black", size=whiskersize) +
		scale_color_manual(values=immunecolors,guide="none") + 
		scale_x_discrete(labels=immunelabs) + 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name="log10 RNA copies per ml")) + 
		# scale_y_reverse() + 
		theme_classic() + 
		theme(
			text=element_text(size=9),
			axis.text.x=element_text(angle=30,hjust=1)) + 
		labs(x="", y="Peak Ct")


fig_wp_whiskers <- meanvals_indiv %>% 
	select(wp, id=immunecat) %>% 
	ggplot(aes(x=factor(id), col=factor(id), y=wp)) + 
		geom_beeswarm(size=pointsize, alpha=alphabee, cex=0.8) + 
		geom_segment(data=postdf_overall_summary, aes(x=id, xend=id, y=post_wp_lwr, yend=post_wp_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=post_wp_lwr, yend=post_wp_lwr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=post_wp_upr, yend=post_wp_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-1.2*whiskerwidth, xend=id+1.2*whiskerwidth, y=post_wp_mean, yend=post_wp_mean), col="black", size=whiskersize) +  
		scale_color_manual(values=immunecolors,guide="none") + 
		scale_x_discrete(labels=immunelabs) + 
		# scale_y_reverse() + 
		theme_classic() + 
		theme(
			text=element_text(size=9),
			axis.text.x=element_text(angle=30,hjust=1)) + 
		labs(x="", y="Proliferation time (days)")


fig_wr_whiskers <- meanvals_indiv %>% 
	select(wr, id=immunecat) %>% 
	ggplot(aes(x=factor(id), col=factor(id), y=wr)) + 
		geom_beeswarm(size=pointsize, alpha=alphabee, cex=1.0) + 
		geom_segment(data=postdf_overall_summary, aes(x=id, xend=id, y=post_wr_lwr, yend=post_wr_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=post_wr_lwr, yend=post_wr_lwr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-whiskerwidth, xend=id+whiskerwidth, y=post_wr_upr, yend=post_wr_upr), col="black", size=whiskersize) +  
		geom_segment(data=postdf_overall_summary, aes(x=id-1.2*whiskerwidth, xend=id+1.2*whiskerwidth, y=post_wr_mean, yend=post_wr_mean), col="black", size=whiskersize) +  
		scale_color_manual(values=immunecolors,guide="none") + 
		scale_x_discrete(labels=immunelabs) + 
		# scale_y_reverse() + 
		theme_classic() + 
		theme(
			text=element_text(size=9),
			axis.text.x=element_text(angle=30,hjust=1)) + 
		labs(x="", y="Clearance time (days)")


# ==============================================================================
# Bootstrapping: 
# ==============================================================================

# fig_titerboxes_thisvar <- indiv_data %>% 
# 	group_by(id_clean) %>% 
# 	slice(1) %>% 
# 	group_by(immunecat) %>% 
# 	straptheboot(goods="Titer",nstraps=1000) %>% 
# 	group_by(immunecat) %>% 
# 	summarise(titer_mean=mean(Titer),titer_lwr=quantile(Titer,0.025),titer_upr=quantile(Titer,0.975)) %>% 
# 	ggplot() + 
# 		geom_point(aes(x=immunecat, y=titer_mean)) + 
# 		geom_segment(aes(x=as.numeric(immunecat), xend=as.numeric(immunecat), y=titer_lwr, yend=titer_upr)) + 
# 		geom_segment(aes(x=as.numeric(immunecat)-0.1, xend=as.numeric(immunecat)+0.1, y=titer_lwr, yend=titer_lwr)) + 
# 		geom_segment(aes(x=as.numeric(immunecat)-0.1, xend=as.numeric(immunecat)+0.1, y=titer_upr, yend=titer_upr)) + 
# 		theme_classic() + 
# 		scale_x_continuous(breaks=1:length(immunelabs), labels=immunelabs) + 
# 		labs(x="", y="Mean titer (95% CI)") + 
# 		theme(text=element_text(size=9),
# 			axis.text.x=element_text(angle=30,hjust=1))

# ==============================================================================
# Summarise fits
# ==============================================================================

summdf_geml <- postdf_dp %>% 
	group_by(id) %>% 
	summarise(mean=round(convert_Ct_logGEML(global_pars[["lod"]]-mean(post_dp)),1), lwr=round(convert_Ct_logGEML(global_pars[["lod"]]-quantile(post_dp,0.025)),1), upr=round(convert_Ct_logGEML(global_pars[["lod"]]-quantile(post_dp,0.975)),1)) %>% 
	mutate(var="Peak GEML") %>% 
	left_join(tibble(name=immunelabs, id=1:length(immunelabs)), by="id") %>% 
	select(var,name,mean,lwr,upr) %>% 
	mutate(string=paste0(mean," (",lwr,", ",upr,")"))

summdf_ct <- postdf_dp %>% 
	group_by(id) %>% 
	summarise(mean=round(global_pars[["lod"]]-mean(post_dp),1), lwr=round(global_pars[["lod"]]-quantile(post_dp,0.975),1), upr=round(global_pars[["lod"]]-quantile(post_dp,0.025),1)) %>% 
	mutate(var="Peak Ct") %>% 
	left_join(tibble(name=immunelabs, id=1:length(immunelabs)), by="id") %>% 
	select(var,name,mean,lwr,upr) %>% 
	mutate(string=paste0(mean," (",lwr,", ",upr,")"))

summdf_wp <- postdf_wp %>% 
	group_by(id) %>% 
	summarise(mean=round(mean(post_wp),1), lwr=round(quantile(post_wp,0.025),1), upr=round(quantile(post_wp,0.975),1)) %>% 
	mutate(var="Proliferation time") %>% 
	left_join(tibble(name=immunelabs, id=1:length(immunelabs)), by="id") %>% 
	select(var,name,mean,lwr,upr) %>% 
	mutate(string=paste0(mean," (",lwr,", ",upr,")"))

summdf_wr <- postdf_wr %>% 
	group_by(id) %>% 
	summarise(mean=round(mean(post_wr),1), lwr=round(quantile(post_wr,0.025),1), upr=round(quantile(post_wr,0.975),1)) %>% 
	mutate(var="Clearance time") %>% 
	left_join(tibble(name=immunelabs, id=1:length(immunelabs)), by="id") %>% 
	select(var,name,mean,lwr,upr) %>% 
	mutate(string=paste0(mean," (",lwr,", ",upr,")"))

summdf_overall <- bind_rows(
	summdf_ct,
	summdf_geml,
	summdf_wp,
	summdf_wr)

