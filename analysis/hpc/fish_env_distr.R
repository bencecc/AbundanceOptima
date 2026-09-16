# Function to compute the quantile distribution of environmental variables from
# reef fish geographic distribution (e.g., the species thermal index, oxygen
# and salinity, etc). The function scores generic fish names and resolves
# aggregate estimates up to the genus or family lavel. 
#####################################################################################

fish_env_distr <- function(sp.name, env.layer, supp.occ.df, out.name, ...) {
	
	# Step 1: try robis::occurrence directly on taxon name
	occ.tmp <- robis::occurrence(sp.name)
	
	if(nrow(occ.tmp)==0) {
		tmp.name <- common_to_sci(sp.name)
		
		if(nrow(tmp.name)==1) {
			occ.tmp <- robis::occurrence(tmp.name$Species)
		}
		
		else if(nrow(tmp.name)>1) {
			# if robis::occurrence returns many genera, one could use occ_mult_sp
			# to select the most frequent genus 
			#occ.tmp <- occ_mult_sp(sp=tmp.name$Species)
			occ.tmp1 <- lapply(tmp.name$Species, function(x) robis::occurrence(x))
			occ.tmp <- bind_rows(occ.tmp1)
		}
		
		else {
			
			# Step 2a: validates genera or (sub)families separating text within brackets or by other characters
			chr.split <- str_extract_all(sp.name, boundary("word"))[[1]]
			chr.rm <- c("Unidentified", "unidentified", "fish", "larvae", "larval",
					"sp", "sp.", "spp", "spp.")
			
			if(any(chr.split%in%chr.rm)) {
				chr.split <- chr.split[-which(chr.split%in%chr.rm)]
			}
			
			
			# Step 2b: validates individual genera or (sub)family names
			if(length(chr.split)>0) {
				
				# recheck if the first and second character of chr.split (to capture genus or (sub)family)
				# is a common name; limit to the first two characters to minimize the risk for the algorithm
				# to interpret specific names as the common name of another species.
				tmp.name <- common_to_sci(chr.split[1])
				
				if((nrow(tmp.name)==0)&&(length(chr.split)>1)) {
					tmp.name <- common_to_sci(chr.split[2])
				}
				# if robis::occurrence returns many genera, one could use occ_mult_sp
				# to select the most frequent genus 
#				if(nrow(tmp.name)>0) {
#					occ.tmp <- occ_mult_sp(tmp.name$Species)
#				}
				if(nrow(tmp.name)==1) {
					occ.tmp <- robis::occurrence(tmp.name$Species)
				}
				
				if(nrow(tmp.name)>1) {
					# if robis::occurrence returns many genera, one could use occ_mult_sp
					# to select the most frequent genus
					#occ.tmp <- occ_mult_sp(sp=tmp.name$Species)
					occ.tmp1 <- lapply(tmp.name$Species, function(x) robis::occurrence(x))
					occ.tmp <- bind_rows(occ.tmp1)
				}				
				
				#check for genera or families within any element of chr.split
				if(nrow(tmp.name)==0) {
					tmp.name <- fishbase |> filter(Genus%in%chr.split)
				}
				
				if(nrow(tmp.name)==0) {
					tmp.name <- fishbase |> filter(Family%in%chr.split|
									SubFamily%in%chr.split)
				}				
				
			}
			
			# if tmp.name includes many species within a genus of (sub)family,
			# select those included in the original taxon name
			if((nrow(tmp.name)>0)&&((!exists("occ.tmp"))||nrow(occ.tmp)==0)) {
				
				# filter by other elements (2nd onward) of chr.split if they are species names 
				if(length(chr.split)>1) {
					tmp.name1 <- tmp.name |> filter(Species%in%chr.split[2:length(chr.split)])
					
					if(nrow(tmp.name1)>0) {
						tmp.name <- tmp.name1				
					}
				}
				
				tmp.sp <- tmp.name |> group_by(Genus, Species) |>
						mutate(genspec=paste(Genus, Species, sep=" ")) |>
						ungroup() |> distinct(genspec)
				
				if(nrow(tmp.sp)==1) {
					occ.tmp <- robis::occurrence(tmp.sp$genspec)
				}
				
				if(nrow(tmp.sp>1)) {
					occ.tmp1 <- lapply(tmp.sp$genspec, function(x) robis::occurrence(x))
					occ.tmp <- bind_rows(occ.tmp1)
				}
			}
			
			# if taxon is not yet found, try using the generic (and slow) worms_validate function
			else {
				tmp.name <- worms_validate(sp.name)
				
				if((nrow(tmp.name)==0)&(exists("chr.split"))) {
					tmp.name <- worms_validate(chr.split[1])
					
					if((nrow(tmp.name)==0)&(length(chr.split)>1)) {
						tmp.name1 <- worms_validate(chr.split[2])
						
						if(nrow(tmp.name1)>0) {
							tmp.name <- tmp.name1				
						}
					}	
				}
				
				if((nrow(tmp.name)>0)&&(ncol(tmp.name)==1)) {
					occ.tmp <- robis::occurrence(tmp.name)
				}
				
				if((nrow(tmp.name)>0)&&(ncol(tmp.name)>1)) {
					# check species
					tmp.sp <- tmp.name |> filter(Species%in%chr.split) |> 
							group_by(Genus, Species) |>
							mutate(genspec=paste(Genus, Species, sep=" ")) |>
							ungroup() |> distinct(genspec)
					occ.tmp <- robis::occurrence(tmp.sp$genspec)
					
					# check genus
					if(nrow(tmp.sp)==0) {
						tmp.sp <- tmp.name |> filter(Genus%in%chr.split) |> 
								group_by(Genus, Species) |>
								mutate(genspec=paste(Genus, Species, sep=" ")) |>
								ungroup() |> distinct(genspec)
						occ.tmp <- robis::occurrence(tmp.sp$genspec)
					}
					
					# check family
					if(nrow(tmp.sp)==0) {
						
						tmp.sp <- tmp.name |> filter(Family%in%chr.split) |> 
								group_by(Genus, Species) |>
								mutate(genspec=paste(Genus, Species, sep=" ")) |>
								ungroup() |> distinct(genspec)
						occ.tmp <- robis::occurrence(tmp.sp$genspec)
						
					}
					
				}
				
			}
			
		}
		
	}	
		
	if(nrow(occ.tmp)>0) {
		
		# combine OBIS occurrences with those in fish.occ.df
		occ.df <- rbind(
				occ.tmp |>
				rename(LAT=decimalLatitude, LON=decimalLongitude,
						YEAR=date_year, SPECIES=scientificName) |>
				select(LON,LAT,YEAR,SPECIES),
				fish.occ.df |> filter(SPECIES%in%sp.name)
						) |>
						group_by(YEAR,SPECIES) |>
						distinct(LON,LAT) |>
						ungroup()			
		
	} 

  else {
		
		occ.df <- supp.occ.df |> filter(SPECIES%in%sp.name) |>
				group_by(YEAR,SPECIES) |>
				distinct(LON,LAT) |>
				ungroup()		
	}
	
	# Species oxygen midpoint and quartile distribution
	if(nrow(occ.df)>0) {
		
		env.layer <- unwrap(env.layer)
		crs(env.layer) <- "EPSG:4326" #"+proj=longlat +datum=WGS84"
		time.vec <- terra::time(env.layer)
		timestamp = as.Date(format(time.vec, "%Y-%m-%d"))
		years_all <- as.vector(format(time.vec, "%Y"))
		names(env.layer) <- years_all
		
		sf.dat <- st_as_sf(occ.df, coords=c("LON","LAT"), crs = 4326)
		sf.dat <- st_transform(sf.dat, crs = "epsg:4326")
		#focal_coords <- st_coordinates(sf.dat)
		dat.wide <- terra::extract(env.layer, sf.dat)
		dat.long <- dat.wide |>
				pivot_longer((any_of("ID")+1):last_col(),
						names_to='YEAR', values_to='env')		
		occ.df$ID <- dat.wide$ID
		clim.tmp <- occ.df |>
				left_join(dat.long, by=c("ID","YEAR")) |>
				drop_na(env)
		
    # select min/max latitude for point in range (remove q1 and q99 first)
  	q1 <- quantile(clim.tmp$LAT, probs=0.01, na.rm=T)
	  q99 <- quantile(clim.tmp$LAT, probs=0.99, na.rm=T)
			
		lat.ord <- clim.tmp |> arrange(desc(LAT)) |>
			filter(LAT>q1&LAT<q99)

		env_tmp <- clim.tmp |>
				summarise(
							min.lat = min(lat.ord$LAT),
							max.lat = max(lat.ord$LAT),
							mean.lat = mean(lat.ord$LAT),
							q5lat = quantile(lat.ord$LAT, probs=0.05, na.rm=TRUE),
							q10lat = quantile(lat.ord$LAT, probs=0.10, na.rm=TRUE),
							q25lat = quantile(lat.ord$LAT, probs=0.25, na.rm=TRUE),
							q50lat = quantile(lat.ord$LAT, probs=0.5, na.rm=TRUE),
							q70lat = quantile(lat.ord$LAT, probs=0.70, na.rm=TRUE),
              q75lat = quantile(lat.ord$LAT, probs=0.75, na.rm=TRUE),
							q90lat = quantile(lat.ord$LAT, probs=0.90, na.rm=TRUE),
							q95lat = quantile(lat.ord$LAT, probs=0.95, na.rm=TRUE),
							env_mean = mean(env, na.rm=T), 
							env_q2.5 = quantile(env, probs = 0.025,na.rm=T),
							env_q5 = quantile(env, probs = 0.05,na.rm=T),
							env_q10 = quantile(env, probs = 0.10,na.rm=T),
							env_q50 = quantile(env, probs = 0.50,na.rm=T),
							env_q70 = quantile(env, probs = 0.70,na.rm=T),
						  env_q90 = quantile(env, probs = 0.90,na.rm=T),
						  env_q95 = quantile(env, probs = 0.95,na.rm=T),
						  env_q97.5 = quantile(env, probs = 0.975,na.rm=T),
						  env_mean_TM  = (env_q2.5 + env_q97.5)/2,
						  count   = n())
		
		env_q <- cbind(SPECIES=sp.name, env_tmp)
		
	} 
  
  else {
	  env_q <- data.frame(
				SPECIES=sp.name,
				min.lat=NA,
				max.lat=NA,
				mean.lat=NA,
				q5lat=NA,
				q10lat=NA,
				q25lat=NA,
				q50lat=NA,
				q70lat=NA,
	      q75lat=NA,
				q90lat=NA,
				q95lat=NA,
			  env_mean=NA,
			  env_q2.5=NA,
			  env_q5=NA,
			  env_q10=NA,
			  env_q50=NA,
			  env_q70=NA,
			  env_q90=NA,
			  env_q95=NA,
			  env_q97.5=NA,
			  env_mean_TM=NA,
			  count=NA
		)
	}
	
	env_q <- env_q |>
			rename_with(.fn = function(.x){str_replace(.x,"env", out.name)})
	
	return(env_q)

}


