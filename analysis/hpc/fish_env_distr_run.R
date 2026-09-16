# Wrap function to compute the quantile distribution of environmental variables from reef fish
# geographic distribution (e.g., the species thermal index, oxygen and salinity, etc).
# Function to compute fish species exposure index to an environmental variable on cluster. It uses
# worms_validate. Called by fish_env_distr.sh. Uses fish_env_distr.R
#### ----------------------------------------------------------------------------------------------- ####
require(sf)
require(terra)
require(dplyr)
require(tidyr)
require(purrr)
require(robis)
require(dbplyr)
require(rfishbase)
require(stringr)
require(worrms)

require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=1)

#library(parallel)
#cl <- makeCluster(10)
#doParallel::registerDoParallel(cl)

#### ---------------------------------------------------------------- LOAD DATA AND SET VARIABLES --------------------------------------------------------- ####

# loaf functions
source("/home/lisandro/workspace/MPA_Timeseries/fish_env_distr.R")
source("/home/lisandro/workspace/MPA_Timeseries/worms_validate.R")

# load fish names and occurrene data
load("/home/lisandro/Lavori/MPA_timeseries/fish_names.RData")
load("/home/lisandro/Lavori/MPA_timeseries/fish.occ.df.RData")

# set name of environmental variable
env.name <- "thetao"

# generate and set output directory
if(!file.exists(paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, "_results", sep="")))
	dir.create(paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, "_results", sep=""))

dir <- paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, "_results", sep="")

# select environmental data (consistent with env.name)
#env_path <- "/home/lisandro/Lavori/MPA_timeseries/EnvData/glorys_salinity_annual_1993_2022_5_10_metres_mean.nc"
#env_path <- "/home/lisandro/Lavori/MPA_timeseries/EnvData/glorys_oxygen_annual_1993_2022_5_10_metres_mean.nc"
env_path <- "/home/lisandro/Lavori/MPA_timeseries/EnvData/glorys_annual_1993_2021_5_10_metres_mean.nc"

# mutate names to some taxa
fish_names <- fish_names |>
		# fish_names on server has three species names modified (rows below)
		# to retreive occurrence from robis;original names have been restored
		# in final fish_sti_res file (see bottom of script) 
		mutate(SPECIES=case_when(
						SPECIES%in%"baitfish, unidentified" ~ "baitfish",
						SPECIES%in%"Pleuronectiformes N.id." ~ "Pleuronectiformes",
						SPECIES%in%"Soleidae n. Id." ~ "Soleidae",
						TRUE ~ as.character(SPECIES)
				))

dat.tmp <- fish_names

#### ---------------------------------------------------------------- START ANALYSIS --------------------------------------------------------- ####

env.rast <- wrap(terra::rast(env_path))

args <- commandArgs(trailingOnly = TRUE)
lw <- as.numeric(args[1])
up <- as.numeric(args[2])

taxa.names <- unique(dat.tmp$SPECIES)[lw:up]

env.tmp <- foreach(i = 1:length(taxa.names), .packages=c("sf","terra","dplyr","tidyr","purrr","robis","rfishbase","stringr","worrms")) %dopar% {
	
	cat("iter = ", i, "of ", length(taxa.names), "\n")
	
	sp.occ <- fish.occ.df |> filter(SPECIES%in%taxa.names[i])
	
	fish.occ.tmp <- try(fish_env_distr(sp.name=taxa.names[i], env.layer=env.rast,
					supp.occ.df=sp.occ, out.name=env.name), silent=T)
	
	if (inherits(fish.occ.tmp, "try-error")) {
		fish.occ <- NULL
	}
	
	else {
		fish.occ <- fish.occ.tmp
	}
	
	fish.occ
}

env.chunk <- bind_rows(purrr::compact(env.tmp))

#stopCluster(cl)

assign(paste(env.name, lw, '_', up, sep=""), value=env.chunk, pos=1, inherits=T)
outputName=paste(env.name, lw, '_', up, ".RData",sep="")
outputPath=file.path(dir, outputName)
save(list=paste(env.name, lw, '_', up, sep=""), file=outputPath)

