#### ---- Short-term counterfactual (firstsite or previoussite) ---- ####
#
# For every population (ECOREGION x SPECIES) compares the exposure realised at the
# yearly optimum sites (obs.var, a column of sp.optim.<kind>.shift, e.g.
# shifted.cumtemp.above.mean) with the exposure the population would have had at a
# reference site: ref.var = obs.var - diff.var, where diff.var is a column of the
# lead-time panel sp.optim.<kind>.shift.leadtime (sp_optimloc_ecoreg_shift_lead_time.R)
# giving, for each sampling year, the difference between exposure at the realised
# site and at the first peak site (diff.FirstSite.*) or at the previous year's peak
# site (diff.PreviousSite.*). Both trajectories cover the SAMPLING YEARS only
# ("short-term"). obs and ref are stacked (CONTR = obs | ref) and fit with
#   ENV.VAR ~ YEAR * CONTR (+ (1|group) | ou | ar1)             -> Trend_Environmental
# whose YEAR:CONTR term is Slope.Dev, the counterfactual statistic. The same
# stacked panel is also fit with
#   ABUND ~ ENV.VAR * CONTR                                    -> Abund_Environmental
#   ABUND ~ YEAR * CONTR                                       -> Abund_Trend
# A reduced model without the CONTR terms gives the *.Alone rows.
#
# Args (6 columns, one row per model spec; leadtime.parms.temp.firstsite.txt or
# leadtime.parms.temp.previoussite.txt, chosen in optimloc_leadtime_trend.sh):
#   var.type      Env | Abund
#   obs.env.var   e.g. shifted.cumtemp.above.mean   (column in sp.optim.<kind>.shift)
#   contrast.var  the matching diff.FirstSite.* or diff.PreviousSite.* column
#                 (column in sp.optim.<kind>.shift.leadtime) - this is what selects
#                 the first-site vs previous-site reference
#   autocor       no_autocor | ou | ar1
#   trans         NULL | Standardize
#   mod.family    gaussian | lognormal | nbinom2
#
# Inputs (paths from config.R): sp.optim.<kind>.shift.RData,
# sp.optim.<kind>.shift.leadtime.RData, fish.traits.dat.RData (thermal guild,
# feeding type and realm are attached as metadata columns, not used in the models).
# Data kind (abund | density) = which pair of load() lines is active below.
#
# Output: one file per model spec, <var.type>_<obs.env.var>_<autocor>_<trans>_<family>
# .RData in Optim_leadtime_trend/, a 6-row-per-population effect table (Ref.Int,
# Ref.Slope, Int.Dev, Slope.Dev, Int.Alone, Slope.Alone; with AIC and warnings).
# The folder name is fixed: rename it after each run (Optim_<kind>_leadtime_trend
# for firstsite, Optim_<kind>_leadtime_prevsite_trend for previoussite) so that
# cluster_summaries.R can reassemble it into optim.<kind>.leadtime.trend /
# optim.<kind>.leadtime.prevsite.trend.
#
# Difference from optimloc_ecoregion_leadtime_trend_longterm.R:
#   - reference trajectory: here from the lead-time panel differences, sampling
#     years only; there from a long-term reference panel covering 1993 -> last
#     sampling year (sp_optimloc_ecoreg_longterm_ref_build.R), stacked with its
#     native year range;
#   - first-site vs previous-site: here chosen by the contrast.var column of the
#     parameter file (6 args); there by the ref.mode variable at the top of the
#     script (5 args, ref.env.var = obs.env.var);
#   - models: here Trend_Environmental + Abund_Environmental + Abund_Trend; there
#     Trend_Environmental only;
#   - output folder: fixed name here (rename after each run); there
#     Optim_leadtime_trend_<ref.mode>_longterm, so nothing to rename.
# The main-text counterfactual uses the long-term version; this script provides the
# short-term twin for consistency checks.
#######################################################################

# load libraries
require(tidyr)
require(lubridate)
require(dplyr)
require(forcats)
require(glmmTMB)
require(datawizard)
#
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=20)
source("config.R")
if(!file.exists(file.path(dir_data, "Optim_leadtime_trend")))
	dir.create(file.path(dir_data, "Optim_leadtime_trend"))

# load data

#load(input_file("sp.optim.abund.shift.RData"))
#load(input_file("sp.optim.abund.shift.leadtime.RData"))
load(input_file("sp.optim.density.shift.RData"))
load(input_file("sp.optim.density.shift.leadtime.RData"))

load(input_file("fish.traits.dat.RData"))

# steps for trend of optimloc shift data analysis on HPC - 
shift.df <- sp.optim.density.shift |>
		rename(LAT=SHIFTED.LAT, LON=SHIFTED.LON)
shift.df.leadtime <- sp.optim.density.shift.leadtime |>
		rename(LAT=SHIFTED.LAT, LON=SHIFTED.LON)

ftr <- fish.traits.dat |>
		select(SPECIES, FeedingType, sst_mean_TM) |>
		rename(temp.optim=sst_mean_TM)
attr(ftr$temp.optim, "names") <- NULL

shifted.temp.df <- shift.df |>
		left_join(ftr, by=c("SPECIES.ORIG"="SPECIES")) |>
		left_join(shift.df.leadtime,
				by=c("LAT","LON","SAMPLED_AREA","YEAR","SPECIES","SPECIES.ORIG")) |>
		group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
		mutate(
				ndays=ifelse(leap_year(YEAR)=="TRUE", 366, 365),
				TEMP.OPTIM=mean(temp.optim),
				THERMAL.GUILD=case_when(
						temp.optim<23 ~ "Temperate",
						temp.optim>=23 ~ "Tropical"
				),
				FEEDING.TYPE=case_when(FeedingType%in%c("browsing on substrate") ~ "Microphages",
						FeedingType%in%c("selective plankton feeding",
								"filtering plankton","variable") ~ "Planktivores",
						FeedingType%in%c("grazing on aquatic plants") ~ "Grazers",
						FeedingType%in%c("hunting macrofauna (predator)","picking parasites off a host (cleaner)") ~ "Carnivores",
						FeedingType%in%c("other") ~ "Other",
						TRUE ~ as.character(FeedingType)),
				REALM=case_when(
						(LAT >= 40&LAT <= 60)|(LAT >= -60&LAT <= -40) ~ "Temperate",
						(LAT > 23.5&LAT < 40)|(LAT > -40&LAT < -23.5) ~ "Subtropical",
						(LAT <= 23.5&LAT >= -23.5) ~ "Tropical"),
				
				THERMAL.GUILD=factor(THERMAL.GUILD),
				FEEDING.TYPE=factor(FEEDING.TYPE),
				REALM=factor(REALM),
				REALM=fct_relevel(REALM, c("Tropical","Subtropical","Temperate"))				
		) |>
		relocate(REALM, .before=ECOREGION) |>
		ungroup()

# function to capture errors and warnings; it returns NULL if the model fails with error,
# otherwise it returns the model output list with the warning message as the last element
# if the model runs with a warning.
# inspired by:
# https://stackoverflow.com/questions/68084740/r-trycatch-but-retain-the-expression-result-in-the-case-of-a-warning
myCatch <- function(expr, ...) {
	
	res <- NULL
	diag <- NULL
	res <- tryCatch(
			expr = {
				withCallingHandlers(
						expr = expr,
						# If expression throws a warning, record diagnostics without halting,
						# so as to store the result of the expression.
						warning = function(w){
							parent <- parent.env(environment())
							parent$diag <- w
						}
				)
			},
			error = function(e) NULL
	)
	
	diagnostic <- function(result, w){
		last.message <<- w$message
	}
	
	warning = diagnostic(res, diag)
	
	if(!is.null(res)&!is.null(warning)) {
		res$warning <- warning
	}
	
	return(res)
}

args <- commandArgs(trailingOnly = TRUE)
var.type <- as.character(args[1])
obs.env.var <- as.character(args[2])
contrast.var <- as.character(args[3])
autocor <- as.character(args[4])
trans <- as.character(args[5])
mod.family <- gsub("\r", "", as.character(args[6]))

min.yrs <- 4

sp.id <- unique(shifted.temp.df$SPECIES) 

shifted.trend.analysis.leadtime <- foreach(i=1:length(sp.id), .combine=rbind) %dopar% {
	
	cat('Doing Species ', i, ' of ', length(sp.id), '\n', sep = '')
	
	tmp.df <- shifted.temp.df |> filter(SPECIES%in%sp.id[i])
	
	# check number of years and filter for ecoregions with n.yrs>=min.yrs
	n.yrs <- tmp.df |> group_by(ECOREGION) |> summarise(nyrs=n(), .groups="drop") |>
			filter(nyrs>=min.yrs)
	
	pop.ecoreg.df <- tmp.df |> filter(ECOREGION%in%n.yrs$ECOREGION)
	
	# proceed if there are at least 5 years of observations
	if(nrow(pop.ecoreg.df)>=min.yrs) {
		
		model.out <- foreach(j=1:nrow(n.yrs), .combine="rbind") %do% { 		
			
			pop.df <- pop.ecoreg.df |> filter(ECOREGION%in%unlist(n.yrs[j,"ECOREGION"]))
			
			if(var.type=="Env") {
				
				# temporal trend of environmental predictor (temperature)
				# transform YEAR into decimal values and use function numFactor in glmmTMB t#o 
				# generate factor yrs with irregularly time spaced level#s
				# also generate a fake group (1) to fit a model with temporal autorcorrelation 
				if(trans=="Standardize") {
					
					# standardize YEAR and ENV.VAR
					lm1.df <- pop.df |>
							mutate(
									obs.var=!!sym(obs.env.var),
									diff.var=!!sym(contrast.var),
									ref.var=obs.var-diff.var,
									YEAR=(YEAR-YEAR[1]) + 1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1) |>
							drop_na(ref.var) |>
							pivot_longer(any_of(c("ref.var","obs.var")),
									names_to='CONTR', values_to="ENV.VAR") |>
							arrange(CONTR) |>
							mutate(
									# standardize YEAR and ENV.VAR
									YEAR=as.vector(standardize(YEAR)),
									ENV.VAR=as.vector(standardize(ENV.VAR)),
									CONTR=as.factor(CONTR),
									CONTR=fct_relevel(CONTR, c("ref.var","obs.var"))
							)
					
				} else if(trans=="NULL") {
					
					lm1.df <- pop.df |>
							mutate(
									obs.var=!!sym(obs.env.var),
									diff.var=!!sym(contrast.var),
									ref.var=obs.var-diff.var,
									YEAR=(YEAR-YEAR[1])+1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1) |>
							drop_na(ref.var) |>
							pivot_longer(any_of(c("ref.var","obs.var")),
									names_to='CONTR', values_to="ENV.VAR") |>
							arrange(CONTR) |>
							mutate(
									CONTR=as.factor(CONTR),
									CONTR=fct_relevel(CONTR, c("ref.var","obs.var"))
							)
					
				}
				
				if(autocor=="no_autocor") {
					
					my.fm1 <- formula(ENV.VAR ~ YEAR*CONTR + (1|group))
					
				} else if(autocor=="ou") {
					
					my.fm1 <- formula(ENV.VAR ~ YEAR*CONTR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm1 <- formula(ENV.VAR ~ YEAR*CONTR + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				m1 <- myCatch(glmmTMB(my.fm1,
								data=lm1.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family))
				
				if(any(names(m1)=="warning")) {
					
					mw1 <- m1$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw1)|
							grepl("false convergence", mw1)) {
						
						m1 <- myCatch(glmmTMB(my.fm1,
										data=lm1.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optimizer=optim,
												optArgs=list(method="BFGS"))))
						
						if(any(names(m1)=="warning")) {
							mw1 <- m1$warning
						}
						
					}
					
					if(grepl("limit reached without convergence", mw1)) {
						
						m1 <- myCatch(glmmTMB(my.fm1,
										data=lm1.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				# run model only for observed data to assess whether slope == 0 (for annihilation scenario)
				m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
								data=lm1.df |> filter(CONTR%in%"obs.var"),
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family))
				
				if(any(names(m2)=="warning")) {
					
					mw2 <- m2$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw2)|
							grepl("false convergence", mw2)) {
						
						m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
										data=lm1.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optimizer=optim,
												optArgs=list(method="BFGS"))))
						
						if(any(names(m2)=="warning")) {
							mw2 <- m2$warning
						}
						
					}
					
					if(grepl("limit reached without convergence", mw2)) {
						
						m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
										data=lm1.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				if(!is.null(m1)|!is.null(m2)) {
					
					if(!is.null(m1)&!is.null(m2)) {
						
						if(any(names(m1)=="warning")) {
							
							w1 <- m1$warning} else{ w1 <- NA }
						
						if(any(names(m2)=="warning")) {
							
							w2 <- m2$warning} else{ w2 <- NA }
						
						coefs.out1 <- data.frame(rbind(summary(m1)$coefficients$cond,
										summary(m2)$coefficients$cond), 
										AIC=c(rep(AIC(m1),4),rep(AIC(m2),2)), Warning=c(rep(w1,4),rep(w2,2)))
						
					} else if (!is.null(m1)&is.null(m2)) {
						
						if(any(names(m1)=="warning")) {
							
							w1 <- m1$warning} else{ w1 <- NA }
						
						w2 <- NA
						
						coefs.out1 <- data.frame(rbind(summary(m1)$coefficients$cond,
										rep(NA,4), rep(NA,4)),
										AIC=c(rep(AIC(m1),4),rep(AIC(m2),2)), Warning=c(rep(w1,4), rep(w2,2)))
						
					} else if (is.null(m1)&!is.null(m2)) {
						
						w1 <-  NA
						
						if(any(names(m2)=="warning")) {
							
							w2 <- m2$warning} else{ w2 <- NA }
						
						coefs.out1 <- data.frame(rbind(
										rep(NA,4), rep(NA,4), rep(NA,4), rep(NA,4),
										summary(m2)$coefficients$cond),
										AIC=c(rep(AIC(m1),4),rep(AIC(m2),2)), Warning=c(rep(w1,4),rep(w2,2)))
						
					} 
					
					colnames(coefs.out1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out1) <- NULL
					
					out1 <- data.frame(
							REALM=lm1.df[1,"REALM"],
							ECOREGION=lm1.df[1,"ECOREGION"],
							SPECIES=lm1.df[1,"SPECIES"],
              SPECIES.ORIG=lm1.df[1,"SPECIES.ORIG"],
							FEEDING.TYPE=lm1.df[1,"FEEDING.TYPE"],
							TEMP.OPTIM=lm1.df[1,"TEMP.OPTIM"],
							THERMAL.GUILD=lm1.df[1,"THERMAL.GUILD"],
							FlagAnalysis="Trend_Environmental",
							NYEARS=lm1.df[1,"n.yrs"],
							Effect=c(
									"Ref.Int",
									"Ref.Slope",
									"Int.Dev",
									"Slope.Dev",
									"Int.Alone",
									"Slope.Alone"
							),
							coefs.out1
					)
					
				} else if (is.null(m1)&is.null(m2)) { out1 <- NULL }
				
				if(!is.null(out1)) { out <- out1 }
				
			} else if (var.type=="Abund") {
				
				# abundance as response variable: relationship with environmental predictor
				# and temporal trends
				if(trans=="Standardize") {
					
					# standardize YEAR and environmental predictor (!!sym(obs.env.var))
					lm2.df <- pop.df |>
							mutate(
									obs.var=abund,
									diff.var=!!sym(contrast.var), #diff.FirstSite.abund,
									ref.var=obs.var-diff.var,
									YEAR=(YEAR-YEAR[1]) + 1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1) |>
							drop_na(ref.var) |>
							pivot_longer(c(any_of(c("ref.var","obs.var"))),
									names_to='CONTR', values_to=c("ABUND")) |>
							arrange(CONTR) |>
							mutate(
									# standardize YEAR and ENV.VAR
									#YEAR=as.vector(standardize(YEAR)),
									ENV.VAR=!!sym(obs.env.var),
									ENV.VAR=if_else(CONTR=="ref.var", !!sym(obs.env.var)-!!sym(contrast.var), ENV.VAR),
									# standardize by centering YEAR and ENV.VAR on the ref.var value at
									# the last sampling date, such that the intercept tests for a difference
									# in abundance between ref.var and obs.var at the last sampling date
									cent.val=ENV.VAR[n()],
									ENV.VAR=as.vector(standardize(ENV.VAR, center=cent.val[1])),
									YEAR=as.vector(standardize(YEAR, center=max(YEAR))),
									CONTR=as.factor(CONTR),
									CONTR=fct_relevel(CONTR, c("ref.var","obs.var"))
							)
					
				} else if (trans=="NULL") {
					
					lm2.df <- pop.df |>
							mutate(
									obs.var=abund,
									diff.var=!!sym(contrast.var), #diff.FirstSite.abund,
									ref.var=obs.var-diff.var,
									YEAR=(YEAR-YEAR[1]) + 1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1) |>
							drop_na(ref.var) |>
							pivot_longer(any_of(c("ref.var","obs.var")),
									names_to='CONTR', values_to="ABUND") |>
							arrange(CONTR) |>
							mutate(
									ENV.VAR=!!sym(obs.env.var),
									ENV.VAR=if_else(CONTR=="ref.var", !!sym(obs.env.var)-!!sym(contrast.var), ENV.VAR),
									# center YEAR and ENV.VAR on the ref.var value at the last sampling
									# date, such that the intercept tests for a difference in abundance
									# between ref.var and obs.var at the last sampling date
									cent.val=ENV.VAR[n()],
									ENV.VAR=ENV.VAR-cent.val[1],
									YEAR=YEAR-max(YEAR),
									CONTR=as.factor(CONTR),
									CONTR=fct_relevel(CONTR, c("ref.var","obs.var"))
							)
					
				}
				
				# relationship with environmental predictors
				if(autocor=="no_autocor") {
	
					my.fm2 <- formula(ABUND ~ ENV.VAR*CONTR + (1|group))
	
				} else if(autocor=="ou") {
					
					my.fm2 <- formula(ABUND ~ ENV.VAR*CONTR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm2 <- formula(ABUND ~ ENV.VAR*CONTR + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				# NOTE: abundance values are a constant through time for ref.var (the value estimated
				# at the optimal abundance site in the first year). Thus, the reference slope is zero
				# and the difference in slopes between CONTR levels is actually a test of the hypothesis
				# that optimal fish abundance deviates from zero with varying ENV.VAR.
				
				# Since ENV.VAR is centered on the value observed the last sampling date at the
				# first optimal site (i.e. the value of ENV.VAR that fish would have experienced if
				# the optimal site did not change - ref.env), the difference in intercepts between
				# levels of CONTR provide a contrast of fish abundance estimated at the last optimal
				# site (last peak abundance for obs.var) vs the initial abundance at ref.var - i.e.
				# whether shifting latitude through time is associated with an increase, a decrease
				# or no change in fish abundance at the final optimal site compared to the abundance
				# estimated at the first optimal site in the time series.
				m3 <- myCatch(glmmTMB(my.fm2,
								data=lm2.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family,
								na.action="na.omit"))
				
				if(any(names(m3)=="warning")) {
					
					mw3 <- m3$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw3)|
							grepl("false convergence", mw3)) {
						
						m3 <- myCatch(glmmTMB(my.fm2,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optimizer=optim,
												optArgs=list(method="BFGS"))))
						
						if(any(names(m3)=="warning")) {
							mw3 <- m3$warning
						}
						
					}
					
					if(grepl("limit reached without convergence", mw3)) {
						
						m3 <- myCatch(glmmTMB(my.fm2,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				if(!is.null(m3)) {
					
					if(any(names(m3)=="warning")) {
						
						w3 <- m3$warning} else{ w3 <- NA }
					
					coefs.out2 <- data.frame(summary(m3)$coefficients$cond,
							AIC=rep(AIC(m3),4),
							Warning=c(rep(w3,4)))
					
					colnames(coefs.out2) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out2) <- NULL
					
					out2 <- data.frame(
							REALM=lm2.df[1,"REALM"],
							ECOREGION=lm2.df[1,"ECOREGION"],
							SPECIES=lm2.df[1,"SPECIES"],
							SPECIES.ORIG=lm2.df[1,"SPECIES.ORIG"],
              FEEDING.TYPE=lm2.df[1,"FEEDING.TYPE"],
							TEMP.OPTIM=lm2.df[1,"TEMP.OPTIM"],
							THERMAL.GUILD=lm2.df[1,"THERMAL.GUILD"],
							FlagAnalysis="Abund_Environmental",
							NYEARS=lm2.df[1,"n.yrs"],
							Effect=c(
									"Ref.Int",
									"Ref.Slope",
									"Int.Dev",
									"Slope.Dev"
							),
							coefs.out2
					)
					
				} else if (is.null(m3)) { out2 <- NULL }
				
				# temporal trends
				if(autocor=="no_autocor") {
	
					my.fm3 <- formula(ABUND ~ YEAR*CONTR + (1|group))
	
				} else if(autocor=="ou") {
					
					my.fm3 <- formula(ABUND ~ YEAR*CONTR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm3 <- formula(ABUND ~ YEAR*CONTR + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				# NOTE: abundance values are a constant through time for ref.var (the value estimated
				# at the optimal abundance site in the first year). Thus, the reference slope is zero
				# and the difference in slopes between CONTR levels is actually a test of the hypothesis
				# that temporal trends in optimal fish abundance deviate from zero.
				
				# Since YEAR is centered on the last sampling date, the difference in intercepts
				# between levels of CONTR provides a contrast of fish abundance estimated at the
				# last optimal site (last peak abundance for obs.var) vs the initial abundance
				# at ref.var - i.e. whether shifting latitude through time has determined an increase,
				# a decrease or no change in fish abundance at the final optimal site compared to the
				# abundance esitmated at the first peak site.
				m4 <- myCatch(glmmTMB(my.fm3,
								data=lm2.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family,
								na.action="na.omit"))
				
				if(any(names(m4)=="warning")) {
					
					mw4 <- m4$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw4)|
							grepl("false convergence", mw4)) {
						
						m4 <- myCatch(glmmTMB(my.fm3,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optimizer=optim,
												optArgs=list(method="BFGS"))))
						
						if(any(names(m4)=="warning")) {
							mw4 <- m4$warning
						}
						
					}
					
					if(grepl("limit reached without convergence", mw4)) {
						
						m4 <- myCatch(glmmTMB(my.fm3,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				if(!is.null(m4)) {
					
					if(any(names(m4)=="warning")) {
						
						w4 <- m4$warning} else{ w4 <- NA }
					
					coefs.out3 <- data.frame(summary(m4)$coefficients$cond,
							AIC=rep(AIC(m4),4),
							Warning=c(rep(w4,4)))
					
					colnames(coefs.out3) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out3) <- NULL
					
					out3 <- data.frame(
							REALM=lm2.df[1,"REALM"],
							ECOREGION=lm2.df[1,"ECOREGION"],
							SPECIES=lm2.df[1,"SPECIES"],
							SPECIES.ORIG=lm2.df[1,"SPECIES.ORIG"],
              FEEDING.TYPE=lm2.df[1,"FEEDING.TYPE"],
							TEMP.OPTIM=lm2.df[1,"TEMP.OPTIM"],
							THERMAL.GUILD=lm2.df[1,"THERMAL.GUILD"],
							FlagAnalysis="Abund_Trend",
							NYEARS=lm2.df[1,"n.yrs"],
							Effect=c(
									"Ref.Int",
									"Ref.Slope",
									"Int.Dev",
									"Slope.Dev"
							),
							coefs.out3
					)
					
				} else if(is.null(m4)) { out3 <- NULL }
				
				if(!is.null(out2)&!is.null(out3)) {
					
					out <- rbind(out2,out3)
					
				} else if(!is.null(out2)&is.null(out3)) {
					
					out <- out3
					
				} else if(is.null(out3)&!is.null(out2)) {
					
					out <- out2
					
				}
				
			}
			
	  }
		
		model.out
		
	}  # else { out <- NULL }
	
}

if(!is.null(shifted.trend.analysis.leadtime)) {
	
	mod.parms <- paste(var.type, obs.env.var, autocor, trans, mod.family, sep="_")
	shifted.trend.analysis.leadtime <- shifted.trend.analysis.leadtime |>
			mutate(mod.parms=mod.parms) |>
			relocate(mod.parms, .before="Warning")
	
	outputName=paste(mod.parms, ".RData", sep="")
	outputPath=file.path(file.path(dir_data, "Optim_leadtime_trend"), outputName)
	assign(paste(mod.parms), value=shifted.trend.analysis.leadtime, pos=1, inherits=T)
	save(list=paste(mod.parms), file=outputPath)
	
}













