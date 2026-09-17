#### ---- Analysis of the thermal profiles of fish communities inside and outside MPAs ---- ####

# load libraries
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
# generate and set output directory
if(!file.exists(file.path(dir_data, "Optim_trend")))
	dir.create(file.path(dir_data, "Optim_trend"))

#load(input_file("sp.optim.abund.shift.RData"))
load(input_file("sp.optim.density.shift.RData"))

# Bassian splitted

load(input_file("fish.traits.dat.RData"))

# Trend of the modskurt optimum-location shift data (HPC).
shift.df <- sp.optim.density.shift |>
		rename(ECOREGION=ORIG.ECOREGION)

ftr <- fish.traits.dat |>
		select(SPECIES, FeedingType, sst_mean_TM) 

shifted.temp.df <- shift.df |>
		rename(LAT=SHIFTED.LAT, LON=SHIFTED.LON) |>
		left_join(ftr, by=c("SPECIES.ORIG"="SPECIES")) |>
		group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
		mutate(
				ndays=ifelse(leap_year(YEAR)=="TRUE", 366, 365),
				TEMP.OPTIM=mean(sst_mean_TM),
				THERMAL.GUILD=case_when(
						TEMP.OPTIM<23 ~ "Temperate",
						TEMP.OPTIM>=23 ~ "Tropical"
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
		relocate(SPECIES.ORIG, .after=SPECIES)

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
resp.env <- as.character(args[2])
autocor <- as.character(args[3])
trans <- as.character(args[4])
#mod.family <- as.character(args[5])
mod.family <- gsub("\r", "", as.character(args[5]))

pred.env <- resp.env
min.yrs <- 4

sp.id <- unique(shifted.temp.df$SPECIES) 

shifted.trend.analysis <- foreach(i=1:length(sp.id), .combine=rbind) %dopar% {
	
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
				# select years when the species is present in a given site
				# transform YEAR into decimal values and use function numFactor in glmmTMB t#o 
				# generate factor yrs with irregularly time spaced level#s
				# also generate a fake group (1) to fit a model with temporal autorcorrelation 
				if(trans=="Standardize") {
					
					# standardize only predictor (YEAR)
					lm1.df <- pop.df |>
							filter(abund>0) |>
							mutate(
									RESP.ENV=!!sym(resp.env),
									YEAR=(YEAR-YEAR[1]) + 1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									RESP.ENV=as.vector(standardize(!!sym(resp.env))),
									# standardize YEAR
									YEAR=as.vector(standardize(YEAR)),
									group=1)
					
				} else if(trans=="NULL") {
					
					lm1.df <- pop.df |>
							filter(abund>0) |>
							mutate(
									RESP.ENV=!!sym(resp.env),
									YEAR=(YEAR-YEAR[1])+1,
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1)
					
				}
				
				if(autocor=="no_autocor") {
					
					my.fm1 <- formula(RESP.ENV ~ YEAR + (1|group))
					
				} else if(autocor=="ou") {
					
					my.fm1 <- formula(RESP.ENV ~ YEAR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm1 <- formula(RESP.ENV ~ YEAR + ar1(as.factor(YEAR) + 0|group))
					
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
				
				if(!is.null(m1)) {
					
					if(any(names(m1)=="warning")) {
						w1 <- m1$warning
					} else{ w1 <- NA }
					
					coefs.out1.m1 <- data.frame(summary(m1)$coefficients$cond,
						AIC=rep(AIC(m1),2),
						Warning=rep(w1,2))
					colnames(coefs.out1.m1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out1.m1) <- NULL
					
				} else { 
					
					coefs.out1.m1 <- data.frame(t(rep(NA,4)),
						AIC=rep(NA,2),
						Warning=rep(NA,2))
					colnames(coefs.out1.m1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out1.m1) <- NULL
                                       
				}
				
				out1 <- data.frame(
						REALM=lm1.df[1,"REALM"],
						ECOREGION=lm1.df[1,"ECOREGION"],
						#MPA=lm1.df[1,"MPA"],
						SPECIES=lm1.df[1,"SPECIES"],
						SPECIES.ORIG=lm1.df[1,"SPECIES.ORIG"],
						FEEDING.TYPE=lm1.df[1,"FEEDING.TYPE"],
						TEMP.OPTIM=lm1.df[1, "TEMP.OPTIM"], 
						THERMAL.GUILD=lm1.df[1,"THERMAL.GUILD"],                          
						FlagAnalysis="Trend_Environmental",
						Nyears=lm1.df[1,"n.yrs"],
						Effect=c(
								"Env.Int",
								"Env.Trend"),
						coefs.out1.m1
				)
				
				if(!is.null(out1)) { out <- out1 }
				
			} else if(var.type=="Abund") {
				
				# temporal trend in abundance, including species absences
				# and regression of abundance vs. environmental data
				# transform YEAR into decimal values and use function numFactor in glmmTMB to 
				# generate factor yrs with irregularly time spaced levels
				# also generate a fake group (1) to fit a model with temporal autorcorrelation 
				if(trans=="Standardize") {
					
					# standardize only predictor (PRED.THERM)
					lm2.df <- pop.df |>
							mutate(
									#RESP.ABUND=as.vector(standardize(!!sym(resp.biol))),
									RESP.ABUND=abund,
									PRED.THERM=as.vector(standardize(!!sym(pred.env))),
									YEAR=(YEAR-YEAR[1])+1,		
									n.yrs=n(),
									#pred.env_l=lag(!!sym(pred.env)),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1
							)
					
				} else if (trans=="NULL") {
					
					lm2.df <- pop.df |>
							mutate(
									RESP.ABUND=abund,
									PRED.THERM=!!sym(pred.env),
									YEAR=(YEAR-YEAR[1])+1,		
									n.yrs=n(),
									#pred.env_l=lag(!!sym(pred.env)),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1
							)
					
				}
				
				if(autocor=="no_autocor") {
					
					my.fm2 <- formula(RESP.ABUND ~ YEAR + (1|group))
					
				} else if(autocor=="ou") {
					
					my.fm2 <- formula(RESP.ABUND ~ YEAR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm2 <- formula(RESP.ABUND ~ YEAR + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				m2 <- myCatch(glmmTMB(my.fm2,
								data=lm2.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family,
								na.action="na.omit"))
				
				if(any(names(m2)=="warning")) {
					
					mw2 <- m2$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw2)|
							grepl("false convergence", mw2)) {
						
						m2 <- myCatch(glmmTMB(my.fm2,
										data=lm2.df,
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
						
						m2 <- myCatch(glmmTMB(my.fm2,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				if(!is.null(m2)) {
					
					if(any(names(m2)=="warning")) {
						w2 <- m2$warning
					} else{ w2 <- NA }
					
					coefs.out2.m2 <- data.frame(summary(m2)$coefficients$cond,
						AIC=rep(AIC(m2),2),
						Warning=rep(w2,2))
					colnames(coefs.out2.m2) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out2.m2) <- NULL
					
				} else  { 
					
					coefs.out2.m2 <- data.frame(t(rep(NA,4)),
						AIC=rep(NA,2),						
						Warning=rep(NA,2))
					colnames(coefs.out2.m2) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out2.m2) <- NULL
				
        }
				
				out2 <- data.frame(
						REALM=lm2.df[1,"REALM"],
						ECOREGION=lm2.df[1,"ECOREGION"],
						#MPA=lm2.df[1,"MPA"],
						SPECIES=lm2.df[1,"SPECIES"],
						SPECIES.ORIG=lm2.df[1,"SPECIES.ORIG"],
						FEEDING.TYPE=lm2.df[1,"FEEDING.TYPE"],
						TEMP.OPTIM=lm2.df[1, "TEMP.OPTIM"],
						THERMAL.GUILD=lm2.df[1,"THERMAL.GUILD"],            
						FlagAnalysis="Trend_Abundance",
						Nyears=lm2.df[1,"n.yrs"],
						Effect=c(
								"Abund.Int",
								"Abund.Trend"),
						coefs.out2.m2
				)						
				
				if(autocor=="no_autocor") {
					
					my.fm3 <- formula(RESP.ABUND ~ PRED.THERM + (1|group))
					
				} else if(autocor=="ou") {
					
					my.fm3 <- formula(RESP.ABUND ~ PRED.THERM + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm3 <- formula(RESP.ABUND ~ PRED.THERM + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				m3 <- myCatch(glmmTMB(my.fm3,
								data=lm2.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family,
								na.action="na.omit"))
				
				if(any(names(m3)=="warning")) {
					
					mw3 <- m3$warning
					
					if(grepl("non-positive-definite Hessian matrix", mw3)|
							grepl("false convergence", mw3)) {
						
						m3 <- myCatch(glmmTMB(my.fm3,
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
						
						m3<- myCatch(glmmTMB(my.fm3,
										data=lm2.df,
										ziformula = ~ 0, dispformula = ~ 1,
										family=mod.family,
										control=glmmTMBControl(
												optCtrl=list(iter.max=1e3,eval.max=1e3))))	
					}
					
				}
				
				if(!is.null(m3)) {
					
					if(any(names(m3)=="warning")) {
						w3 <- m3$warning
					} else{ w3 <- NA }
					
					coefs.out3.m3 <- data.frame(summary(m3)$coefficients$cond,
						AIC=rep(AIC(m3),2),						
						Warning=rep(w3,2))
					colnames(coefs.out3.m3) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out3.m3) <- NULL
					
				} else  { 
					
					coefs.out3.m3 <- data.frame(t(rep(NA,4)),
						AIC=rep(NA,2),
						Warning=rep(NA,2))
					colnames(coefs.out3.m3) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out3.m3) <- NULL
                                  
				}
				
				out3 <- data.frame(
						REALM=lm2.df[1,"REALM"],
						ECOREGION=lm2.df[1,"ECOREGION"],
						#MPA=lm2.df[1,"MPA"],
						SPECIES=lm2.df[1,"SPECIES"],
						SPECIES.ORIG=lm2.df[1,"SPECIES.ORIG"],
						FEEDING.TYPE=lm2.df[1,"FEEDING.TYPE"],
						TEMP.OPTIM=lm2.df[1, "TEMP.OPTIM"], 
						THERMAL.GUILD=lm2.df[1,"THERMAL.GUILD"],                                                            
						FlagAnalysis="Abundance_Environmental",
						Nyears=lm2.df[1,"n.yrs"],
						Effect=c(
								"AbundEnv.Int",
								"AbundEnv.Trend"),
						coefs.out3.m3
				
				)
				
				if(!is.null(out2)&!is.null(out3)) {
					
					out <- rbind(out2,out3)
					
				} else if(!is.null(out2)&is.null(out3)) {
					
					out <- out3
					
				} else if(is.null(out3)&!is.null(out2)) {
					
					out <- out2
					
				}
				
			} else if (var.type%in%c("Lat", "Lon")) {
				
				# latitudinal temporal trends
				if(trans=="Standardize") {
					
					# standardize env predictor (PRED.THERM) and YEAR
					lm3.df <- pop.df |>
							mutate(
									#RESP.LAT=abs(LAT),
									RESP=case_when(
										var.type=="Lat" ~ abs(LAT),
										var.type=="Lon" ~ abs(LON)),                  
                  YEAR=(YEAR-YEAR[1])+1,		
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									YEAR=as.vector(standardize(YEAR)),
									group=1
							)
					
				} else if (trans=="NULL") {
					
					lm3.df <- pop.df |>
							mutate(
									#RESP.LAT=abs(LAT),
									RESP=case_when(
										var.type=="Lat" ~ abs(LAT),
										var.type=="Lon" ~ abs(LON)),
                  YEAR=(YEAR-YEAR[1])+1,		
									n.yrs=n(),
									times.dec=as.numeric(paste(YEAR, 1, sep=".")),
									cum.times=lead(cumsum(lag(times.dec, default=0)),
											default=sum(times.dec)),
									yrs=glmmTMB::numFactor(cum.times),
									group=1
							)
					
				}
				
				if(autocor=="no_autocor") {
					
					my.fm4 <- formula(RESP ~ YEAR + (1|group))
					
				} else if(autocor=="ou") {
					
					my.fm4 <- formula(RESP ~ YEAR + ou(yrs + 0|group))
					
				} else if(autocor=="ar1") {
					
					my.fm4 <- formula(RESP ~ YEAR + ar1(as.factor(YEAR) + 0|group))
					
				}
				
				m4 <- myCatch(glmmTMB(my.fm4,
								data=lm3.df,
								ziformula = ~ 0, dispformula = ~ 1,
								family=mod.family,
								na.action="na.omit"))
				
				if(!is.null(m4)) {
					
					if(any(names(m4)=="warning")) {
						w4 <- m4$warning
					} else{ w4 <- NA }
					
					coefs.out4.m4 <- data.frame(summary(m4)$coefficients$cond,
						AIC=rep(AIC(m4),2),
						Warning=rep(w4,2))
					colnames(coefs.out4.m4) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out4.m4) <- NULL
					
				} else  { 
					
					coefs.out4.m4 <- data.frame(t(rep(NA,4)),
						AIC=rep(NA,2),
						Warning=c(rep(NA,2)))
					colnames(coefs.out4.m4) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
					rownames(coefs.out4.m4) <- NULL
                                       
				}
				
				out4 <- data.frame(
						REALM=lm3.df[1,"REALM"],
						ECOREGION=lm3.df[1,"ECOREGION"],
						#MPA=lm3.df[1,"MPA"],
						SPECIES=lm3.df[1,"SPECIES"],
						SPECIES.ORIG=lm3.df[1,"SPECIES.ORIG"],
						FEEDING.TYPE=lm3.df[1,"FEEDING.TYPE"],
						TEMP.OPTIM=lm3.df[1, "TEMP.OPTIM"], 
						THERMAL.GUILD=lm3.df[1,"THERMAL.GUILD"],                                                            
						FlagAnalysis=ifelse(var.type=="Lat", "Trend_Lat", "Trend_Lon"),
						Nyears=lm3.df[1,"n.yrs"],
						Effect=c(paste(var.type, "Int", sep="."), paste(var.type, "Trend", sep=".")),
						coefs.out4.m4
				)
				
				if(!is.null(out4)) { out <- out4 }
			}
			
		}
		
		model.out
		
	}  # else { out <- NULL }
	
}

if(!is.null(shifted.trend.analysis)) {
	
	mod.parms <- paste(var.type, resp.env, autocor, trans, mod.family, sep="_")
	shifted.trend.analysis$mod.parms <- mod.parms
	shifted.trend.analysis <- shifted.trend.analysis |>
			relocate(mod.parms, .before="Warning")
	
	outputName=paste(mod.parms, ".RData", sep="")
	outputPath=file.path(file.path(dir_data, "Optim_trend"), outputName)
	assign(paste(mod.parms), value=shifted.trend.analysis, pos=1, inherits=T)
	save(list=paste(mod.parms), file=outputPath)
	
}














