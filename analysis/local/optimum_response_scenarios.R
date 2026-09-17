# Derive scenarios of response to warming from shifts in optimal location. The input
# file is a data frme originating from temporal models of shift in optimal locaiton
# (e.g., from function fit_optimum_scenario_ecoregion_leadtime). This function is
# called in script sp.optimloc.scenarios
########################################################################################

# multiple, complex, scenarios (scenarios and responses)


optimum_response_scenarios <- function(df, ...) {
	
	p.value <- as.vector(unlist(df[,"P.value"]))
	estimate <- as.vector(unlist(df[,"Estimate"]))
	
	int.ref <-  estimate[1]
	slope.ref <-  estimate[2]
	int.dev <-  estimate[3]
	slope.dev <-  estimate[4]
	realized.int <- int.ref + int.dev
	realized.slope <- slope.ref + slope.dev
	
	p.int.ref <- p.value[1]
	p.slope.ref <- p.value[2]
	p.int.dev <- p.value[3]
	p.slope.dev <- p.value[4]
	
	#flag <- as.vector(unlist(df$FlagAnalysis[1]))
	#model <- as.vector(unlist(df$mod.parms[1]))
	
	df <- df |>
			mutate(
					Scenario=NA,
					Response=NA,
					Scenario=case_when(							
							(slope.ref>0&slope.dev>0&p.slope.dev<0.05)| # Increasing - magnificaiton, slope
									(slope.ref>0&p.slope.ref<0.05&int.dev>0&p.int.dev<0.05&p.slope.dev>=0.05)| # Increasing - magnificaiton, elevation
									(slope.ref>0&slope.dev<0&realized.slope>0&p.slope.dev<0.05)| # Increasing - mitigation, slope
									(slope.ref>0&p.slope.ref<0.05&int.dev<0&p.int.dev<0.05&p.slope.dev>=0.05)| # Increasing - mitigation, elevation
									(slope.ref<0&p.slope.ref<0.05&slope.dev>0&p.slope.dev<0.05&realized.slope>0) ~ "Increasing", # Increasing - opposition
							p.int.dev>=0.05&p.slope.dev>=0.05 ~ "Matching",
							(slope.ref<0&slope.dev<0&p.slope.dev<0.05)| # Buffering - strong, slope
									(slope.ref<0&p.slope.ref<0.05&int.dev<0&p.int.dev<0.05&p.slope.dev>=0.05)| # Buffering - strong, elevation
									(slope.ref<0&slope.dev>0&realized.slope<0&p.slope.dev<0.05)| # Buffering - moderate, slope
									(slope.ref<0&p.slope.ref<0.05&int.dev>0&p.int.dev<0.05&p.slope.dev>=0.05)| # Buffering - moderate, elevation
									(slope.ref>0&p.slope.ref<0.05&slope.dev<0&p.slope.dev<0.05&realized.slope<0)	~ "Buffering",	# Buffering, offsetting
							# Level-shift only: slopes non-significant but int.dev significant
							p.slope.ref>=0.05 & p.slope.dev>=0.05 & p.int.dev<0.05 ~ "Increasing",
							TRUE ~ NA_character_),   # scalar default: all case_when inputs size 1 (no dplyr 1.2 deprecation)
					Response=case_when(
							Scenario=="Increasing"&((slope.ref>0&slope.dev>0&p.slope.dev<0.05)| # Increasing - magnificaiton, slope
										(slope.ref>0&p.slope.ref<0.05&int.dev>0&p.int.dev<0.05&p.slope.dev>=0.05)) ~ "Magnification", # Increasing - magnification, elevation
							Scenario=="Increasing"&((slope.ref>0&slope.dev<0&realized.slope>0&p.slope.dev<0.05)| # Increasing - mitigation, slope
										(slope.ref>0&p.slope.ref<0.05&int.dev<0&p.int.dev<0.05&p.slope.dev>=0.05)) ~ "Mitigation", # Increasing - mitigation, elevation
							Scenario=="Increasing"&(slope.ref<0&p.slope.ref<0.05&slope.dev>0&p.slope.dev<0.05&realized.slope>0) ~ "Opposition",
							# Level-shift magnification / mitigation (non-sig slopes, sig int.dev)
							Scenario=="Increasing"&p.slope.ref>=0.05&p.slope.dev>=0.05&int.dev>0 ~ "Magnification",
							Scenario=="Increasing"&p.slope.ref>=0.05&p.slope.dev>=0.05&int.dev<0 ~ "Mitigation",
							Scenario=="Matching"&p.slope.ref>=0.05 ~ "Matching-Steady",
							Scenario=="Matching"&slope.ref>0&p.slope.ref<0.05 ~ "Matching-Warming",
							Scenario=="Matching"&slope.ref<0&p.slope.ref<0.05 ~ "Matching-Cooling",
							Scenario=="Buffering"&((slope.ref<0&slope.dev<0&p.slope.dev<0.05)| # Buffering - strong, slope
										(slope.ref<0&p.slope.ref<0.05&int.dev<0&p.int.dev<0.05&p.slope.dev>=0.05)) ~ "Acceleration", # Buffering - acceleration, elevation
							Scenario=="Buffering"&((slope.ref<0&slope.dev>0&realized.slope<0&p.slope.dev<0.05)| # Buffering - moderate, slope
										(slope.ref<0&p.slope.ref<0.05&int.dev>0&p.int.dev<0.05&p.slope.dev>=0.05)) ~ "Deceleration", # Buffering - deceleration, elevation
							Scenario=="Buffering"&(slope.ref>0&p.slope.ref<0.05&slope.dev<0&p.slope.dev<0.05&realized.slope<0) ~ "Annihilation", # Buffering, annihilation
							TRUE ~ NA_character_)
			)
	
	df
	
}

# simple scenarios from observed temporal trends
simple_response_scenarios <- function(df, ...) {
	
	p.value <- as.vector(unlist(df[,"P.value"]))
	estimate <- as.vector(unlist(df[,"Estimate"]))
	
	df <- df |>
			mutate(
					Scenario=NA,
					Response=NA,
					Scenario=case_when(
							estimate[2]>0&p.value[2]<0.05 ~ "Warming",
							p.value[2]>0.05 ~ "No trend",
							estimate[2]<0&p.value[2]<0.05 ~ "Cooling",
							TRUE ~ NA_character_)
			)
	
	df
	
}


