# This function calculates the exposure to environmental stressful conditions at
# shifted sites determined by function modskurt_analysis, which identifies the optimal
# position along latitudinal gradients. Like sp_optimloc_ecoreg_shift, but it calculates
# environmental statistics at time lags to enable comparison of actual fish exposure
# to stressful conditions to the exposure that would have occurred if peaks in
# abundance: a) did not move in space and occurred repeatedly through time at the
# site where the first peak was observed or b) remained at the peak site of the previous
# sampling year. Comparisons quantify exposure to thermal extremes at the first peak
# site for all the subsqeuent leading years (a) and at the preceeding sampling year
# for all the other sites. The output provides differences of extreme temperatures
# (above 0.975 percentile threshold for each species) and days above threshold with
# respect to previous year and site of peak abundace (diff.temp, diff.days), to
# first year in the time series (diff.firstyear.temp, diff.firstyear.days) and
# cumulative differences across all sampling years in a time series
# (cum.temp, cum.days).

# load libraries
require(dplyr)
require(tidyr)
require(stringr)
require(sf)
require(terra)
require(ncdf4)
require(ncdf4.helpers)

require(foreach, quietly=T)
require(doMC)
registerDoMC(cores=5)
source("config.R")
#require(parallel)
#cl <- makeCluster(5)
#doParallel::registerDoParallel(cl)

# load data
# indices to subset ncdf file by slecting only sites and years as per input data

load(input_file("reef_fish_sti_glorys.RData"))

load(input_file("optim.abund.relocated.RData"))
#load(input_file("optim.density.relocated.RData"))

# indices to subset ncdf file by selecting only sites and years as per input data
load(input_file("idx.optim.abund.relocated.RData"))
#load(input_file("idx.optim.density.relocated.RData"))

# Bassian region splitted by LON

idx.df <- idx.optim.abund.relocated # does not include year 1992, which must be removed from range.df set below to ensure the two data frames match

# prefix name of environmental variable in species index file
env.name <- "thetao"
# desired name for output files
out.name <- "optim"
# desired abbreviated name of environmental variable to be combined with meand, median and percentiles
out.env.name <- "temp"

# generate and set output directory
if(!file.exists(file.path(dir_data, paste0("sp_", out.name, "_shift_leadtime"))))
	dir.create(file.path(dir_data, paste0("sp_", out.name, "_shift_leadtime")))

dir <- file.path(dir_data, paste0("sp_", out.name, "_shift_leadtime"))

range.df <- optim.abund.relocated |> mutate(abund=H.LAT) |> 
		# |> filter(YEAR!=1992) # remove year 1992 for which there is no temp data
		left_join(idx.df, by=c("ECOREGION","SPECIES","SPECIES.ORIG","LAT","LON","YEAR")) |>
		arrange(ECOREGION,SPECIES)

# generate id for each species in each ECOREGION
ngroups <- range.df %>% arrange(ECOREGION, LAT, LON, YEAR) %>%
		group_by(ECOREGION, SPECIES) %>% dplyr::summarise(tlength=n(), .groups="drop")
group_id <- rep(1:nrow(ngroups), times=ngroups$tlength)
range.df <- range.df |> mutate(GROUP_ID = group_id) |>
		relocate(GROUP_ID, .before=SAMPLED_AREA)
#idx.df$GROUP_ID <- group_id

# replace env.name prefix with "env" 
env.df <- reef_fish_sti_glorys |>
		rename_with(.fn = function(.x){str_replace(.x, env.name, "env")})

# combine data on peak location with specie STI thresholds
work.df <-  range.df |>
		left_join(env.df, by=c("SPECIES.ORIG"="SPECIES"))

sf.dat <- st_as_sf(work.df, coords=c("LON","LAT"))
st_crs(sf.dat) <- "+proj=longlat +datum=WGS84"
#sf.dat <- st_shift_longitude(sf.dat) # not needed with glorys data where longtude is already in the range -180,180
coords <- st_coordinates(sf.dat)
work.df$LON <- coords[,1]
work.df$LAT <- coords[,2]

r <- terra::rast(glorys_daily_mean_tif)
time.vec <- terra::time(r)
timestamp = as.Date(format(time.vec, "%Y-%m-%d"))
years_all <- data.frame(year=as.vector(format(time.vec, "%Y")))

r <- wrap(r)

args <- commandArgs(trailingOnly = TRUE)
lw <- as.numeric(args[1])
up <- as.numeric(args[2])
cat(lw,up)

# utility function to extract indices of environmental exposure
# from grouped data of species environmental affinities 
var_idx_lead <- function(data, ...) {
	
	env.mean <- unique(unlist(as.vector(data$env_mean_TM)))	
#	env.q2.5 <- unique(unlist(as.vector(data$env_q2.5)))
#	env.q5 <- unique(unlist(as.vector(data$env_q5)))
#	env.q10 <- unique(unlist(as.vector(data$env_q10)))
	env.q50 <- unique(unlist(as.vector(data$env_q50)))
	env.q70 <- unique(unlist(as.vector(data$env_q70)))
	env.q90 <- unique(unlist(as.vector(data$env_q90)))
	env.q95 <- unique(unlist(as.vector(data$env_q95)))
	env.q97.5 <- unique(unlist(as.vector(data$env_q97.5)))
	
	data.frame(
			SHIFTED.LAT=data[1,"LAT"],
			SHIFTED.LON=data[1,"LON"],
			SAMPLED_AREA=data[1,"SAMPLED_AREA"],
			ECOREGION=data[1,"ECOREGION"],
			YEAR=data[1,"YEAR"],
			SPECIES.ORIG=data[1,"SPECIES.ORIG"],
			abund=data[1,"abund"],
#			shifted.cum.env=sum(data[,"env"], na.rm=T),
			shifted.mean.env=mean(data[,"env"], na.rm=T),
			shifted.max.env=max(data[,"env"], na.rm=T),
			shifted.sd.env=sd(data[,"env"], na.rm=T),
			shifted.cumenv.above.mean=sum(data[which(data$env>env.mean), "env"], na.rm=T),
			shifted.cumenv.above.q97.5=sum(data[which(data$env>env.q97.5), "env"], na.rm=T),
			shifted.cumenv.above.q95=sum(data[which(data$env>env.q95), "env"], na.rm=T),
			shifted.cumenv.above.q90=sum(data[which(data$env>env.q90), "env"], na.rm=T),
			shifted.cumenv.above.q70=sum(data[which(data$env>env.q70), "env"], na.rm=T),
			shifted.cumenv.above.q50=sum(data[which(data$env>env.q50), "env"], na.rm=T),
#			shifted.cumenv.below.mean=sum(data[which(data$env<env.mean), "env"], na.rm=T),
#			shifted.cumenv.below.q50=sum(data[which(data$env<env.q50), "env"], na.rm=T),
#			shifted.cumenv.below.q10=sum(data[which(data$env<env.q10), "env"], na.rm=T),
#			shifted.cumenv.below.q5=sum(data[which(data$env<env.q5), "env"], na.rm=T),
#			shifted.cumenv.below.q2.5=sum(data[which(data$env<env.q2.5), "env"], na.rm=T),
			shifted.days.above.mean=length(data[which(data$env>env.mean), "env"]),
			shifted.days.above.q97.5=length(data[which(data$env>env.q97.5), "env"]),
			shifted.days.above.q95=length(data[which(data$env>env.q95), "env"]),
			shifted.days.above.q90=length(data[which(data$env>env.q90), "env"]),
			shifted.days.above.q70=length(data[which(data$env>env.q70), "env"]),
			shifted.days.above.q50=length(data[which(data$env>env.q50), "env"])
#			shifted.days.below.mean=length(data[which(data$env<env.mean), "env"]),
#			shifted.days.below.q50=length(data[which(data$env<env.q50), "env"]),
#			shifted.days.below.q10=length(data[which(data$env<env.q10), "env"]),
#			shifted.days.below.q5=length(data[which(data$env<env.q5), "env"]),
#			shifted.days.below.q2.5=length(data[which(data$env<env.q2.5), "env"])
	)
	
}

group.id <- unique(work.df$GROUP_ID)

env.indices <- foreach(i=lw:up, .combine="rbind",
				.packages = c("dplyr","terra","foreach","tidyr"),
				#.export = c("r"),
				.inorder = TRUE) %do% {
			
			cat('Doing iter ', i, ' of ', length(group.id), '\n', sep = '')  
			
			r <- unwrap(r)
			
      # here, env_mean etc. are fish species STI thresholds
			df.sub <- work.df |> filter(GROUP_ID%in%group.id[i]) |>
					select(ECOREGION, LAT, LON, SAMPLED_AREA, YEAR, SPECIES, SPECIES.ORIG, abund, env_mean,
							env_q2.5, env_q5, env_q10, env_q50, env_q70, env_q90,
							env_q95, env_q97.5, env_mean_TM)
			
			if(nrow(df.sub)>1) {
				
				idx.group <- work.df |> filter(GROUP_ID%in%group.id[i])
				
				sp.env.stat <- foreach(idx = 1:nrow(idx.group), .combine="rbind") %dopar% {
					
					p <- sf::st_point(as.vector(unlist(df.sub[idx,c("LON","LAT")]))) |> 
							sf::st_sfc(crs = "epsg:4326") |> 
							terra::vect()
					
					lyr <- NULL
					
					for (t in idx:nrow(idx.group)) {
						lyr <- c(lyr, idx.group[t,"time.start"]:
										(idx.group[t,"time.start"]+idx.group[t,"time.step"]-1))
					}
					
					r.sub <- terra::subset(r, subset=lyr)
					dat <- terra::extract(r.sub, p, ID=FALSE) |> as.double()
					
					# if there are more than 10% NAs, increase the spatial extent
					# and extract daily temperature data from increments of 1*res
					# in all directions
					if(length(which(is.na(dat)))>=length(dat)*0.1) {
						
						lon <- as.vector(unlist(df.sub[idx,c("LON")]))
						lat <- as.vector(unlist(df.sub[idx,c("LAT")]))
						res <- res(r)[1]
						m1 <- matrix(data=c(lon-res,lon,lon+res, lat-res,lat,lat+res), ncol=2, nrow=3)
						p1 <- sf::st_multipoint(m1) |> sf::st_sfc(crs = "epsg:4326") |> 
								terra::vect()
						dat.tmp <- terra::extract(r.sub, p1, ID=FALSE)
						
						# if all data are NA (all sites are on land) increase the spatial extent
						# further to 2*res in all directions
						if(all(is.na(dat.tmp))) {
							
							m2 <- matrix(data=c(lon-2*res, lon-res,lon,lon+res, lon+2*res,
											lat-2*res,lat-res,lat,lat+res, lat+2*res), ncol=2, nrow=5)
							p2 <- sf::st_multipoint(m2) |> sf::st_sfc(crs = "epsg:4326") |> 
									terra::vect()
							dat.tmp <- terra::extract(r.sub, p2, ID=FALSE)	
							
						}
						
						# take daily temperature averages
						dat <- apply(dat.tmp, 2, function(x) mean(x, na.rm=T))
						
					} 
          
          if(length(which(is.na(dat)))<length(dat)*0.1) {
						
						yrs  <- NULL
						for(ny in idx:nrow(idx.group)) {
							
							yrs  <- c(yrs ,
									years_all$year[idx.group[ny,"time.start"]:
													(idx.group[ny,"time.start"]+idx.group[ny,"time.step"]-1)]
							)
						}
						
						env.ann.tmp <- data.frame(LAT=rep(df.sub[idx,c("LAT")], length(dat)),
								LON=rep(df.sub[idx,c("LON")], length(dat)),		
								LEAD_YEAR=as.numeric(yrs), env=dat)
						env.ann <- df.sub[idx,] |>
								right_join(env.ann.tmp, by=c("LAT","LON"))
						
						env.stats <- env.ann |>
								group_by(LEAD_YEAR, SPECIES) |>
								group_modify(~ var_idx_lead(as.data.frame(.))) |>
								ungroup() |>
								relocate(c(LEAD_YEAR, SPECIES, SPECIES.ORIG), .before=abund) |> 
								relocate(ECOREGION, .before=(SHIFTED.LAT))
						
						env.stats
						
					}
					
				}

				# NOTE: YEAR indicates the optimal peak sites and thus different
				# YEARS are also different SITES. Whilst YEAR points to the site
				# sampled that year, LEAD_YEAR includes env values for the subsequent
				# yeare at that site. When YEAR = LEAD_YEAR we have the env value observed
				# at a given site in that year, while the following LEAD_YEAR rows
				# report env values at the site corresponding to YEAR for the subsequent
 				# years. We can then compare the env vaule at which any fish species
				# has been exposed at the observed peak site in a given YEAR to that
				# the fish would have experienced in the same year if it remained at the
				# previous peak site (df2) or at the first peak site (df3). This
				# allwos testing the couterfactual.
				
				# cumulative across years and sites
				df1a <- sp.env.stat |> filter(YEAR==LEAD_YEAR) |>
						reframe(
								ambient.env.mean=mean(shifted.mean.env),
								ambient.env.max=max(shifted.mean.env),
								ambient.env.sd=mean(shifted.sd.env),
								cum.env.mean=sum(shifted.cumenv.above.mean),
								cum.env.50=sum(shifted.cumenv.above.q50),
								cum.env.70=sum(shifted.cumenv.above.q70),
								cum.env.90=sum(shifted.cumenv.above.q90),
								cum.env.95=sum(shifted.cumenv.above.q95),
								cum.env.97.5=sum(shifted.cumenv.above.q97.5),
								cum.days.mean=sum(shifted.days.above.mean),
								cum.days.50=sum(shifted.days.above.q50),
								cum.days.70=sum(shifted.days.above.q70),
								cum.days.90=sum(shifted.days.above.q90),
								cum.days.95=sum(shifted.days.above.q95),
								cum.days.97.5=sum(shifted.days.above.q97.5))	
				# cumulative across years in the first peak site
				df1b <- sp.env.stat |> filter(YEAR==min(YEAR)) |>
						reframe(
								ambient.env.mean=mean(shifted.mean.env),
								ambient.env.max=max(shifted.mean.env),
								ambient.env.sd=mean(shifted.sd.env),
								cum.env.mean=sum(shifted.cumenv.above.mean),
								cum.env.50=sum(shifted.cumenv.above.q50),
								cum.env.70=sum(shifted.cumenv.above.q70),
								cum.env.90=sum(shifted.cumenv.above.q90),
								cum.env.95=sum(shifted.cumenv.above.q95),
								cum.env.97.5=sum(shifted.cumenv.above.q97.5),
								cum.days.mean=sum(shifted.days.above.mean),
								cum.days.50=sum(shifted.days.above.q50),
								cum.days.70=sum(shifted.days.above.q70),
								cum.days.90=sum(shifted.days.above.q90),
								cum.days.95=sum(shifted.days.above.q95),
								cum.days.97.5=sum(shifted.days.above.q97.5))
				
				# get the mean env back
				df1 <- df1a-df1b
				
				# differences with respect to previous year and site
				df2 <- sp.env.stat |>
						group_by(YEAR) |> slice_head(n=2) |>
						#group_by(LEAD_YEAR) |>
						ungroup() |>
						mutate(shift.abund=lag(abund),
								shift.ambient.env.mean=lag(shifted.mean.env), # simple mean, not cumulative as below
								shift.ambient.env.max=lag(shifted.max.env),
								shift.ambient.env.sd=lag(shifted.sd.env),
								shift.env.mean=lag(shifted.cumenv.above.mean),
								shift.env.50=lag(shifted.cumenv.above.q50),
								shift.env.70=lag(shifted.cumenv.above.q70),
								shift.env.90=lag(shifted.cumenv.above.q90),
								shift.env.95=lag(shifted.cumenv.above.q95),
								shift.env.97.5=lag(shifted.cumenv.above.q97.5),
								shift.days.mean=lag(shifted.days.above.mean),
								shift.days.50=lag(shifted.days.above.q50),
								shift.days.70=lag(shifted.days.above.q70),
								shift.days.90=lag(shifted.days.above.q90),
								shift.days.95=lag(shifted.days.above.q95),
								shift.days.97.5=lag(shifted.days.above.q97.5)) |>
						filter(YEAR==LEAD_YEAR) |>
						mutate(
								diff.PreviousSite.abund=abund-shift.abund,
								diff.PreviousSite.ambient.env.mean=shifted.mean.env-shift.ambient.env.mean,
								diff.PreviousSite.ambient.env.max=shifted.max.env-shift.ambient.env.max,
								diff.PreviousSite.ambient.env.sd=shifted.sd.env-shift.ambient.env.sd,
								diff.PreviousSite.env.mean=shifted.cumenv.above.mean-shift.env.mean,
								diff.PreviousSite.env.50=shifted.cumenv.above.q50-shift.env.50,
								diff.PreviousSite.env.70=shifted.cumenv.above.q70-shift.env.70,
								diff.PreviousSite.env.90=shifted.cumenv.above.q90-shift.env.90,
								diff.PreviousSite.env.95=shifted.cumenv.above.q95-shift.env.95,
								diff.PreviousSite.env.97.5=shifted.cumenv.above.q97.5-shift.env.97.5,
								diff.PreviousSite.days.mean=shifted.days.above.mean-shift.days.mean,
								diff.PreviousSite.days.50=shifted.days.above.q50-shift.days.50,
								diff.PreviousSite.days.70=shifted.days.above.q70-shift.days.70,
								diff.PreviousSite.days.90=shifted.days.above.q90-shift.days.90,
								diff.PreviousSite.days.95=shifted.days.above.q95-shift.days.95,
								diff.PreviousSite.days.97.5=shifted.days.above.q97.5-shift.days.97.5) |>
						#relocate(c(diff.env, diff.days), .before="abund") |>
						# remove first raw where variables are NA
						#drop_na(diff.env.mean) |>
						slice(-1L) |>
						select(ECOREGION, SHIFTED.LAT, SHIFTED.LON, SAMPLED_AREA, YEAR, SPECIES,SPECIES.ORIG,
								diff.PreviousSite.abund,
								diff.PreviousSite.ambient.env.mean,
								diff.PreviousSite.ambient.env.max,
								diff.PreviousSite.ambient.env.sd,
								diff.PreviousSite.env.mean,
								diff.PreviousSite.env.50,
								diff.PreviousSite.env.70,
								diff.PreviousSite.env.90,
								diff.PreviousSite.env.95,
								diff.PreviousSite.env.97.5,
								diff.PreviousSite.days.mean,
								diff.PreviousSite.days.50,
								diff.PreviousSite.days.70,
								diff.PreviousSite.days.90,
								diff.PreviousSite.days.95,
								diff.PreviousSite.days.97.5)
				
				# differences with respect to first year in time series
				df3 <- sp.env.stat |>
						mutate(SEL=case_when(
										YEAR>min(YEAR)&LEAD_YEAR>YEAR ~ "NO",
										TRUE ~ as.character("YES")
								)) |>
						relocate(SEL, .before="abund") |>
						filter(SEL=="YES") |>
						arrange(LEAD_YEAR) |>
						filter(LEAD_YEAR>min(LEAD_YEAR)) |>
						mutate(
								shift.abund=lag(abund),
								shift.ambient.env.mean=lag(shifted.mean.env), # simple mean, not cumulative as below
								shift.ambient.env.max=lag(shifted.max.env),
								shift.ambient.env.sd=lag(shifted.sd.env),
								shift.env.mean=lag(shifted.cumenv.above.mean),
								shift.env.50=lag(shifted.cumenv.above.q50),
								shift.env.70=lag(shifted.cumenv.above.q70),
								shift.env.90=lag(shifted.cumenv.above.q90),
								shift.env.95=lag(shifted.cumenv.above.q95),
								shift.env.97.5=lag(shifted.cumenv.above.q97.5),
								shift.days.mean=lag(shifted.days.above.mean),
								shift.days.50=lag(shifted.days.above.q50),
								shift.days.70=lag(shifted.days.above.q70),
								shift.days.90=lag(shifted.days.above.q90),
								shift.days.95=lag(shifted.days.above.q95),
								shift.days.97.5=lag(shifted.days.above.q97.5)) |>
						filter(YEAR==LEAD_YEAR) |>
						mutate(diff.FirstSite.abund=abund-shift.abund,
								diff.FirstSite.ambient.env.mean=shifted.mean.env-shift.ambient.env.mean,
								diff.FirstSite.ambient.env.max=shifted.max.env-shift.ambient.env.max,
								diff.FirstSite.ambient.env.sd=shifted.sd.env-shift.ambient.env.sd,
								diff.FirstSite.env.mean=shifted.cumenv.above.mean-shift.env.mean,
								diff.FirstSite.env.50=shifted.cumenv.above.q50-shift.env.50,
								diff.FirstSite.env.70=shifted.cumenv.above.q70-shift.env.70,
								diff.FirstSite.env.90=shifted.cumenv.above.q90-shift.env.90,
								diff.FirstSite.env.95=shifted.cumenv.above.q95-shift.env.95,
								diff.FirstSite.env.97.5=shifted.cumenv.above.q97.5-shift.env.97.5,
								diff.FirstSite.days.mean=shifted.days.above.mean-shift.days.mean,
								diff.FirstSite.days.50=shifted.days.above.q50-shift.days.50,
								diff.FirstSite.days.70=shifted.days.above.q70-shift.days.70,
								diff.FirstSite.days.90=shifted.days.above.q90-shift.days.90,
								diff.FirstSite.days.95=shifted.days.above.q95-shift.days.95,
								diff.FirstSite.days.97.5=shifted.days.above.q97.5-shift.days.97.5) |>
						#relocate(c(diff.FirstSite.env, diff.FirstSite.days), .before="abund") |>
						select(ECOREGION, SHIFTED.LAT, SHIFTED.LON, SAMPLED_AREA, YEAR, SPECIES, SPECIES.ORIG,
								diff.FirstSite.abund,
								diff.FirstSite.ambient.env.mean,
								diff.FirstSite.ambient.env.max,
								diff.FirstSite.ambient.env.sd,
								diff.FirstSite.env.mean,
								diff.FirstSite.env.50,
								diff.FirstSite.env.70,
								diff.FirstSite.env.90,
								diff.FirstSite.env.95,
								diff.FirstSite.env.97.5,
								diff.FirstSite.days.mean,
								diff.FirstSite.days.50,
								diff.FirstSite.days.70,
								diff.FirstSite.days.90,
								diff.FirstSite.days.95,
								diff.FirstSite.days.97.5)
				
				out <- cbind(df2 |> left_join(df3,
								by=c("ECOREGION", "SHIFTED.LAT", "SHIFTED.LON", "SAMPLED_AREA", "YEAR", "SPECIES", "SPECIES.ORIG")),
						df1)
				
		}
	
}

#stopCluster(cl)

env.indices <- env.indices |>
		rename_with(.fn = function(.x){str_replace(.x, "env", out.env.name)})

dim(env.indices)

assign(paste(out.name, lw, '_', up, sep=""), value=env.indices, pos=1, inherits=T)
outputName=paste(out.name, lw, '_', up, ".RData",sep="")
outputPath=file.path(dir, outputName)
save(list=paste(out.name, lw, '_', up, sep=""), file=outputPath)


