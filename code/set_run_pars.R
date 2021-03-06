run_pars_list <- list(
	# 1) Basic set, titer
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="titer" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 2) Basic set, vax
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="vax" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 3) Basic set, symptoms
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="symp" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 4) Basic set, titer + symptoms
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="titersymp" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 5) Basic set, vax status + symptoms
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="vaxsymp" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 6) By titer, restrict to people with titer taken 100-200 days from exposure
	list(excluded_rows=setdiff(ct_dat_refined$RowID, exposurelist_rows),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="titer" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 7) By vax status, restrict to people with titer taken 100-200 days from exposure
	list(excluded_rows=setdiff(ct_dat_refined$RowID, exposurelist_rows),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="vax" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 8) Basic set, titer with cut @ 400
	list(excluded_rows=c(),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=TRUE,
		immunevar="titer" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		),
	# 9) Sensitivity: restrict to frequent testers
	list(excluded_rows=(ct_dat_refined %>% filter(DetectionSpeed!="Frequent testing") %>% pull(RowID)),
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="titer" # of: "titer" "vax" "symp" "titersymp" "vaxsymp"
		), 
	# 10) high vs low titer among omicron unboosted
	list(excluded_rows=c(),
		# (ct_dat_refined %>% filter(!(LineageBroad=="Omicron" & VaccStatus!="Boosted")) %>% pull(RowID))
		tp_prior=c(0,2),
		dp_midpoint=20,
		wp_midpoint=5,
		wr_midpoint=12,
		sigma_prior=c(0,0.5),
		lambda=0.01,
		fpmean=1/log(10),
		titeralt=FALSE,
		immunevar="unboosttiter" # of: "titer" "vax" "symp" "titersymp" "vaxsymp" "unboosttiter"
		)
	)
