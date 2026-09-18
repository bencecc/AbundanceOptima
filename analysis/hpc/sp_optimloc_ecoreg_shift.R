# Calculates the thermal exposure REALISED at the yearly peak-abundance (optimum)
# sites of each population: the optima located by modskurt_analysis and relocated
# onto reef habitat (optim.<kind>.relocated), identified by their coordinates
# (SHIFTED.LAT / SHIFTED.LON), not by survey SITE_ID. For every population-year,
# daily GLORYS subsurface temperature is extracted at that site for that year
# (the raster layers indexed by idx.optim.<kind>.relocated: time.start / time.step);
# if more than 10% of the days are missing (site on land at raster resolution) the
# extraction is widened to the surrounding 1-cell, then 2-cell, neighbourhood and
# averaged. From the daily series it computes: the mean, maximum, SD and sum of
# temperature over the year (shifted.mean/max/sd/cum.temp), and the cumulative
# degree-days and number of days ABOVE each species-specific threshold — the
# species' thermal-index mean and its 50th, 70th, 90th, 95th and 97.5th percentiles
# (from reef_fish_sti_glorys) — and BELOW its mean and 50th, 10th, 5th and 2.5th
# percentiles (shifted.cumtemp.above/below.*, shifted.days.above/below.*). Peak
# abundance (abund) and the ecoregion the population belongs to (ORIG.ECOREGION)
# are carried along. The two counterfactual exposures (first / previous peak site)
# are NOT computed here but in sp_optimloc_ecoreg_shift_lead_time.R.
#
# Args: lw up = the row range of population-years to process (one chunk per
# HPC task; ranges listed in id.optim.<kind>.relocated.txt, launched by
# sp.optim.shift.sh). Data kind (abund | density) = which pair of load() lines is
# active below.
#
# Inputs (paths from config.R): optim.<kind>.relocated.RData,
# idx.optim.<kind>.relocated.RData, reef_fish_sti_glorys.RData, glorys_daily_mean_tif.
#
# Output: one file per chunk, optim<lw>_<up>.RData in sp_optim_shift/ (rename the
# folder to sp_optim_<kind>_shift), reassembled by cluster_summaries.R into
# sp.optim.<kind>.shift — the realised-exposure panel used by every downstream
# trend, counterfactual and figure script.

# load libraries
require(dplyr)
require(tidyr)
require(stringr)
require(sf)
require(terra)

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
# indices to subset ncdf file by slecting only sites and years as per input data
load(input_file("idx.optim.abund.relocated.RData"))
#load(input_file("idx.optim.density.relocated.RData"))

idx.df <- idx.optim.abund.relocated # does not include year 1992, which must be removed from range.df set below to ensure the two data frames match

# prefix name of environmental variable in species index file
env.name <- "thetao"
# desired name for output files
out.name <- "optim"
# desired abbreviated name of environmental variable to be combined with meand, median and percentiles
out.env.name <- "temp"

# generate and set output directory
if(!file.exists(file.path(dir_data, paste0("sp_", out.name, "_shift"))))
	dir.create(file.path(dir_data, paste0("sp_", out.name, "_shift")))

dir <- file.path(dir_data, paste0("sp_", out.name, "_shift"))

range.df <- optim.abund.relocated |> mutate(abund=H.LAT) |> 
		# |> filter(YEAR!=1992) # remove year 1992 for which there is no temp data
		left_join(idx.df, by=c("ECOREGION","SPECIES","SPECIES.ORIG","LAT","LON","YEAR"))

# replace env.name prefix with "env" 
env.df <- reef_fish_sti_glorys |>
		rename_with(.fn = function(.x){str_replace(.x, env.name, "env")})

work.df <-  range.df |>
		left_join(env.df, by=c("SPECIES.ORIG"="SPECIES"))

sf.dat <- st_as_sf(work.df, coords=c("LON","LAT"))
st_crs(sf.dat) <- "+proj=longlat +datum=WGS84"
#sf.dat <- st_shift_longitude(sf.dat) # not needed with glorys data where longitude is already in the range -180,180
coords <- st_coordinates(sf.dat)
work.df$LON <- coords[,1]
work.df$LAT <- coords[,2]

r <- terra::rast(glorys_daily_mean_tif)

#timeInfo(r)
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
var_idx <- function(data, env.df) {
	
	env.mean <- unlist(as.vector(data$env_mean_TM))
	env.q2.5 <- unlist(as.vector(data$env_q2.5))
	env.q5 <- unlist(as.vector(data$env_q5))
	env.q10 <- unlist(as.vector(data$env_q10))
	env.q50 <- unlist(as.vector(data$env_q50))
	env.q70 <- unlist(as.vector(data$env_q70))
	env.q90 <- unlist(as.vector(data$env_q90))
	env.q95 <- unlist(as.vector(data$env_q95))
	env.q97.5 <- unlist(as.vector(data$env_q97.5))
	
	data.frame(
			ORIG.ECOREGION=data$ECOREGION,
			SHIFTED.LAT=data$LAT,
			SHIFTED.LON=data$LON,
			SAMPLED_AREA=data$SAMPLED_AREA,
      		YEAR=data$YEAR,
			SPECIES.ORIG=data$SPECIES.ORIG,
			abund=data$abund,
			shifted.cum.env=sum(env.df[,"env"], na.rm=T),
			shifted.mean.env=mean(env.df[,"env"], na.rm=T),
			shifted.max.env=max(env.df[,"env"], na.rm=T),
			shifted.sd.env=sd(env.df[,"env"], na.rm=T),
			shifted.cumenv.above.mean=sum(env.df[which(env.df$env>env.mean), "env"], na.rm=T),
			shifted.cumenv.above.q97.5=sum(env.df[which(env.df$env>env.q97.5), "env"], na.rm=T),
			shifted.cumenv.above.q95=sum(env.df[which(env.df$env>env.q95), "env"], na.rm=T),
			shifted.cumenv.above.q90=sum(env.df[which(env.df$env>env.q90), "env"], na.rm=T),
			shifted.cumenv.above.q70=sum(env.df[which(env.df$env>env.q70), "env"], na.rm=T),
			shifted.cumenv.above.q50=sum(env.df[which(env.df$env>env.q50), "env"], na.rm=T),
			shifted.cumenv.below.mean=sum(env.df[which(env.df$env<env.mean), "env"], na.rm=T),
			shifted.cumenv.below.q50=sum(env.df[which(env.df$env<env.q50), "env"], na.rm=T),
			shifted.cumenv.below.q10=sum(env.df[which(env.df$env<env.q10), "env"], na.rm=T),
			shifted.cumenv.below.q5=sum(env.df[which(env.df$env<env.q5), "env"], na.rm=T),
			shifted.cumenv.below.q2.5=sum(env.df[which(env.df$env<env.q2.5), "env"], na.rm=T),
			shifted.days.above.mean=length(env.df[which(env.df$env>env.mean), "env"]),
			shifted.days.above.q97.5=length(env.df[which(env.df$env>env.q97.5), "env"]),
			shifted.days.above.q95=length(env.df[which(env.df$env>env.q95), "env"]),
			shifted.days.above.q90=length(env.df[which(env.df$env>env.q90), "env"]),
			shifted.days.above.q70=length(env.df[which(env.df$env>env.q70), "env"]),
			shifted.days.above.q50=length(env.df[which(env.df$env>env.q50), "env"]),
			shifted.days.below.mean=length(env.df[which(env.df$env<env.mean), "env"]),
			shifted.days.below.q50=length(env.df[which(env.df$env<env.q50), "env"]),
			shifted.days.below.q10=length(env.df[which(env.df$env<env.q10), "env"]),
			shifted.days.below.q5=length(env.df[which(env.df$env<env.q5), "env"]),
			shifted.days.below.q2.5=length(env.df[which(env.df$env<env.q2.5), "env"])
	)
	
}

id <- 1:nrow(work.df)

system.time(env.indices <- foreach(i=lw:up, .combine="rbind", .packages = c("dplyr","terra","foreach"),
				#.export = c("r"),
				.inorder = TRUE) %dopar% {
	
  #	cat('Doing iter ', i, ' of ', length(site.id), '\n', sep = '')  
	
	df.sub <- work.df %>% slice(id[i]) |>
							select(ECOREGION, LAT, LON, SAMPLED_AREA, YEAR, SPECIES, SPECIES.ORIG, abund, env_mean,
									env_q2.5, env_q5, env_q10, env_q50, env_q70, env_q90,
									env_q95, env_q97.5, env_mean_TM) 
	
	p <- sf::st_point(as.vector(unlist(df.sub[1,c("LON","LAT")]))) |> 
			sf::st_sfc(crs = "epsg:4326") |> 
			terra::vect()
      
	lyr <-  tryCatch(work.df[i,"time.start"]:
					      (work.df[i,"time.start"]+work.df[i,"time.step"]-1),
  				error=function(e) NULL)

  if(!is.null(lyr)) {

  	r <- unwrap(r)
  	r.sub <- subset(r, subset=lyr)
  	dat <- terra::extract(r.sub, p, ID=FALSE) |> as.double()
	
  	# if there are more than 10% NAs, increase the spatial extent
   	# and extract daily temperature data from increments of 1*res
	  # in all directions
	  if(length(which(is.na(dat)))>length(dat)*0.1) {
		
		  lon <- as.vector(unlist(df.sub[1,c("LON")]))
		  lat <- as.vector(unlist(df.sub[1,c("LAT")]))
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
		
	  	match.yr <- unique(as.vector(unlist(df.sub$YEAR)))
		  env.ann <- data.frame(YEAR=years_all$year[work.df[i,"time.start"]:
				  				(work.df[i,"time.start"]+work.df[i,"time.step"]-1)],
				  env=dat)
		
		  yrs <- unique(as.vector(unlist(df.sub$YEAR)))
		
	  	sp.env.stats <- foreach(j=1:length(yrs), .combine=rbind) %do% {
		  	
			  env.yr <- env.ann |> filter(YEAR%in%yrs[j])
			
			  sp.yr <- df.sub |> filter(YEAR%in%yrs[j]) |> 
					  distinct(SPECIES, .keep_all=T) |>
					  select (ECOREGION, LAT, LON, SAMPLED_AREA, YEAR, SPECIES, SPECIES.ORIG, abund, env_mean,
						  	env_q2.5, env_q5, env_q10, env_q50, env_q70, env_q90,
						  	env_q95, env_q97.5, env_mean_TM) |>
					  group_by(SPECIES) |>
					  group_modify(~ var_idx(., env.df=env.yr)) |>
					  ungroup() |>
					  relocate(c(SPECIES, SPECIES.ORIG), .before=abund)		
	  	}  
	
	    sp.env.stats
	
	  }
   
  }

})

#stopCluster(cl)

env.indices <- env.indices |>
		rename_with(.fn = function(.x){str_replace(.x, "env", out.env.name)})

dim(env.indices)

assign(paste(out.name, lw, '_', up, sep=""), value=env.indices, pos=1, inherits=T)
outputName=paste(out.name, lw, '_', up, ".RData",sep="")
outputPath=file.path(dir, outputName)
save(list=paste(out.name, lw, '_', up, sep=""), file=outputPath)


