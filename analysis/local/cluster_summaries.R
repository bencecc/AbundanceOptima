# cluster_summaries.R — reassemble the HPC per-chunk outputs into the consolidated
# .RData objects used downstream. Single R session (run after the HPC jobs).
# Covers: site/ecoregion temperature (stage 2 in the analysis workflow, see
# Extended Data FIg. 1), modskurt optima (stage 5), unimodal filter (6),
# sea-path relocation (7), spatial-shift + long-term reference panels (8),
# trends (9). Abundance vs density is selected per block (DATA.KIND, or by
# choosing the input folder).
##################################################################################
require(dplyr)
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=5)
source("config.R")

# COMBINE PER-SPECIES MODSKURT FITS -> modskurt.optim.<kind> -------------------
# Each per-species file (from modskurt_analysis.R) ALREADY carries the cluster split:
#   ECOREGION    = base name (e.g. "Bassian")
#   SPECIES      = "<sp>_<suffix>" for split ecoregions (Bassian -> _W/_E/_NW;
#                  Hawaii -> _SE/_NW), else "<sp>"
#   SPECIES.ORIG = clean species name (for trait/STI/sp.df joins)
# So consolidation (FULL run) is just an rbind of every file in the folder — no tag
# function, no separate Bassian folders, no merge. (suffix = directional label from
# split_plan$ecoreg.)

# ABUND -------------------
setwd(file.path(dir_data, "ModskurtOptimEcoreg_Abund"))
file_list <- as.list(list.files())
modskurt.optim.abund <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	d <- mget(load(file_list[[i]]))[[1]]
	if(!"SPECIES.ORIG" %in% names(d)) d$SPECIES.ORIG <- d$SPECIES   # back-compat for any pre-split files
	d
} |> arrange(ECOREGION, SPECIES, YEAR)

#save(modskurt.optim.abund, file="modskurt.optim.abund.RData")

# DENSITY -------------------
setwd(file.path(dir_data, "ModskurtOptimEcoreg_Density"))
file_list <- as.list(list.files())
modskurt.optim.density <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	d <- mget(load(file_list[[i]]))[[1]]
	if(!"SPECIES.ORIG" %in% names(d)) d$SPECIES.ORIG <- d$SPECIES
	d
} |> arrange(ECOREGION, SPECIES, YEAR)

#save(modskurt.optim.density, file="modskurt.optim.density.RData")


# Apply unimodal filter

DATA.KIND <- "abund"
stopifnot(DATA.KIND %in% c("abund", "density"))

modskurt_dir <- dir_data
load(file.path(modskurt_dir, "sp.df.RData"))                              # sampling sites (density guidance)
load(file.path(modskurt_dir, paste0("modskurt.optim.", DATA.KIND, ".RData")))

# -- Drop weak unimodal-fit populations, then RE-SAVE the filtered optima ---
# Relocation is a separate CPU/memory-heavy step that RELOADS the filtered optima, so the QC filter
# is applied HERE and saved as modskurt.optim.<kind>.unimodal.RData (the original
# modskurt.optim.<kind>.RData is kept intact).
# Filter the ORIGINAL object (keeps m.LAT/m.LON, so the relocation reload+rename
# still works) — NOT a renamed copy. Requires unimodal_fit_check.R to have produced
# the decision CSV. sp_df = sp.df gives the correct per-species ecoregion-id map;
# suffixed split species (SPECIES != SPECIES.ORIG, _W/_E/_NW/_SE) are handled by the
# helper's split branch (id 1). Defaults: all deviation classes, >=50% threshold.
source(file.path(PROJ, "analysis/local/unimodal_filter_helper.R"))
optim_obj    <- paste0("modskurt.optim.", DATA.KIND)             # loaded (unfiltered) object
unimodal_obj <- paste0(optim_obj, ".unimodal")                   # filtered object: NAME == file base
assign(unimodal_obj, unimodal_filter_helper(
  get(optim_obj), sp_df = sp.df,
  decision_csv = file.path(modskurt_dir,
    if (DATA.KIND == "abund") "unimodal_fit_check_pop_decision.csv"
                         else "unimodal_fit_check_pop_decision_density.csv")))
save(list = unimodal_obj, file = file.path(modskurt_dir, paste0(unimodal_obj, ".RData")))

#### RELOCATE MISPLACED SITES ####
require(dplyr)
require(tidyr)
require(sf)
require(terra)
require(RANN)                   # fast NN (density surface in relocation)
require(foreach, quietly=T)
require(doMC, quietly=T)
require(gdistance, quietly=T)   # least-cost sea path (sea-path relocation)
require(raster, quietly=T)      # accCost / extract over the transition
# Cores drive the PARALLEL INNER point loop inside replace_to_bathy_seapath
# (ecoregions are looped serially). Set registerDoMC() cores to suit your machine.
registerDoMC(cores=64)

# -- RELOCATION MODEL: SEA-PATH (adopted 2026-06-05) -----------------------
# This relocation block is CPU- and memory-heavy (run it on a large-memory
# multicore machine). replace_to_bathy_seapath = least-cost sea path,
# land-avoiding, latitude-preserving (with a built-in straight-line fallback).
# Memory model: GEBCO pre-aggregated ONCE to a ~1 km grid
# (cached tif), ecoregions SERIAL, per-point accCost parallelised, transition
# forked read-only (COW).
source(file.path(PROJ, "analysis/local/replace_to_bathy_seapath.R"))
gebco_rast <- rast(gebco_file)  # .nc/.tif/.grd
terra::crs(gebco_rast) <- "epsg:4326"

# Load the QC-FILTERED optima saved by the filter block above (this RELOADS the freshly
# filtered modskurt.optim.<kind>.unimodal.RData) and build the relocation input.
# For a standalone run, set DATA.KIND + modskurt_dir first.
load(file.path(modskurt_dir, paste0("modskurt.optim.", DATA.KIND, ".unimodal.RData")))
optim.df <- get(paste0("modskurt.optim.", DATA.KIND, ".unimodal")) |>
  rename(LAT = m.LAT, LON = m.LON)
ecoreg <- unique(optim.df$ECOREGION)

# -- Pre-aggregate GEBCO ONCE to the ~1 km working grid (native x2) ---------
# Cached as a tif so each foreach worker crops a SMALL grid instead of
# re-aggregating native GEBCO every cycle (that per-worker native aggregation
# was the memory blow-up). Computed only the first time; ~1 km is far finer
# than the ~9 km GLORYS pixel, so exposure is unaffected.
gebco_work_tif <- file.path(modskurt_dir, "GEBCO_work_x2.tif")
if (!file.exists(gebco_work_tif)) {
  message("Aggregating GEBCO x2 to the ~1 km working grid (one-time)...")
  terra::writeRaster(terra::aggregate(gebco_rast, fact = 2, fun = "mean", na.rm = TRUE),
                     gebco_work_tif, overwrite = TRUE)
}

# Ecoregions SERIAL (outer); points PARALLEL inside replace_to_bathy_seapath
# (its per-point accCost runs %dopar% over the 64 registered cores). The
# gdistance sea transition is built ONCE per ecoregion and forked read-only.
optim.shift.ecoreg <- do.call(rbind, lapply(seq_along(ecoreg), function(i) {
	# Memory-light: open the SMALL pre-aggregated grid (lazy crop), not native GEBCO.
	gebco_w <- terra::rast(gebco_work_tif); terra::crs(gebco_w) <- "epsg:4326"
	optim.ecoreg <- optim.df |> filter(ECOREGION%in%ecoreg[i]) |>
			mutate(DATA.TYPE="Initial modeled peak sites") |>
			relocate(DATA.TYPE, .before="ECOREGION")
	sp.df.ecoreg <- sp.df |> filter(ECOREGION%in%optim.ecoreg$ECOREGION&SPECIES%in%optim.ecoreg$SPECIES.ORIG) |>
			dplyr::select("SITE_ID","ECOREGION","SPECIES","YEAR","LAT","LON") |>
			mutate(DATA.TYPE="Original sampled sites") |>
			relocate(DATA.TYPE, .before="ECOREGION")

	t0 <- Sys.time()
	check.coords <- replace_to_bathy_seapath(
    	modeled_df        = optim.ecoreg,    # peaks to relocate (LON/LAT)
    	sample_points     = sp.df.ecoreg,    # density guidance = species' sampling sites
    	bathy_rast        = gebco_w,         # pre-aggregated ~1 km working grid (NOT native)
    	target_depth      = 10,              # relocate to ~10 m isobath
    	tol_m             = 1,               # band half-width = target_depth +/- 1 m
    	target_range      = c(0, 30),        # kept unchanged if depth already in 0-30 m
    	sea_is_negative   = TRUE,            # GEBCO: sea depths negative
    	depth_layer       = 1,               # bathy layer
    	aoi_pad_deg       = 0.75,            # crop pad (auto-enlarged to cover max_dist_km)
    	agg_fact_sea      = 1,               # NO further aggregation (grid already ~1 km)
    	density_weight    = 0.5,             # nearest vs sample-density tie-break
    	k_neighbors       = 10,              # neighbours for the density surface
    	lat_band_deg      = 1/12,            # "same temperature pixel" = GLORYS lat-res (1/12 deg)
    	max_dist_km       = 100              # least-cost sea-path search buffer (km)
  	)
	cat(sprintf("[%d/%d] %-40s %5d pts  %4.1f min\n", i, length(ecoreg), ecoreg[i],
	            nrow(optim.ecoreg), as.numeric(difftime(Sys.time(), t0, units="mins"))))

	check.coords$points_df |>
			mutate(DATA.TYPE="Final modeled peak sites") |>
			drop_na(LAT) |>
			relocate(c(LON,LAT), .before=H.LAT)
}))

setwd(dir_data)

# Drop the incomplete first year (1992). The unimodal-fit QC filter was applied to
# optim.df above (unimodal_filter_helper), so the relocated set is already clean.
optim.relocated <- optim.shift.ecoreg |> filter(YEAR != 1992)

# -- Save with the kind-specific name (optim.abund.relocated / optim.density.relocated)
out_name <- paste0("optim.", DATA.KIND, ".relocated")
assign(out_name, optim.relocated)
save(list = out_name, file = paste0(out_name, ".RData"))
cat("Saved", paste0(out_name, ".RData"), "\n")

# LONGTERM REFERENCE: using temperature timeseries from 1993 to the sampling year for each population
# Run the counterfactual against long-term temperature trends

#setwd(file.path(dir_data, "sp_optim_abund_firstsite_longterm_ref"))
#setwd(file.path(dir_data, "sp_optim_abund_previoussite_longterm_ref"))
#setwd(file.path(dir_data, "sp_optim_density_firstsite_longterm_ref"))
setwd(file.path(dir_data, "sp_optim_density_previoussite_longterm_ref"))

file_list <- as.list(list.files())

shift.tmp <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

sp.optimloc.shift.tmp <- shift.tmp

# rename modified fish names to those originally included in master.fish.dat for matching
# sp.optimloc.shift.tmp <- shift.tmp |>
# 		mutate(SPECIES=case_when(
# 						SPECIES%in%"baitfish" ~ "baitfish, unidentified",
# 						SPECIES%in%"Pleuronectiformes" ~ "Pleuronectiformes N.id.",
# 						SPECIES%in%"Soleidae" ~ "Soleidae n. Id.",
# 						SPECIES%in%"Mugilidae n.id." ~ "Mugilidae",
# 						TRUE ~ as.character(SPECIES)
# 				),
# 				SAMPLED_AREA=as.numeric(SAMPLED_AREA))

setwd(dir_data)

# sp.optim.abund.firstsite.longterm <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.abund.firstsite.longterm, file=file.path(dir_data, "sp.optim.abund.firstsite.longterm.RData"))

# sp.optim.abund.previoussite.longterm <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.abund.previoussite.longterm, file=file.path(dir_data, "sp.optim.abund.previoussite.longterm.RData"))
	
# sp.optim.density.firstsite.longterm <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.density.firstsite.longterm, file=file.path(dir_data, "sp.optim.density.firstsite.longterm.RData"))

sp.optim.density.previoussite.longterm <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.density.previoussite.longterm, file=file.path(dir_data, "sp.optim.density.previoussite.longterm.RData"))

#### ---- Effects of SPATIAL SHIFTS in optimum location ---- ####

#setwd(file.path(dir_data, "sp_optim_abund_shift"))
#setwd(file.path(dir_data, "sp_optim_density_shift"))
#setwd(file.path(dir_data, "sp_optim_abund_shift_leadtime"))
setwd(file.path(dir_data, "sp_optim_density_shift_leadtime"))


file_list <- as.list(list.files())

shift.tmp <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

sp.optimloc.shift.tmp <- shift.tmp

# # rename modified fish names to those originally included in master.fish.dat for matching
# sp.optimloc.shift.tmp <- shift.tmp |>
# 		mutate(SPECIES=case_when(
# 						SPECIES%in%"baitfish" ~ "baitfish, unidentified",
# 						SPECIES%in%"Pleuronectiformes" ~ "Pleuronectiformes N.id.",
# 						SPECIES%in%"Soleidae" ~ "Soleidae n. Id.",
# 						SPECIES%in%"Mugilidae n.id." ~ "Mugilidae",
# 						TRUE ~ as.character(SPECIES)
# 				))



setwd(dir_data)

# load("sp.optim.abund.shift.RData")
# load("sp.optim.abund.shift.leadtime.RData")
# load("sp.optim.density.shift.RData")

# load("optim.abund.relocated.RData")
# load("optim.density.relocated.RData")

# load("sp.optim.abund.shift.leadtime.RData")
# load("sp.optim.density.shift.leadtime.RData")

#### Save files --------------------------------------------------------------------------------------
# sp.optim.abund.shift <- sp.optimloc.shift.tmp |> arrange(ORIG.ECOREGION, SPECIES, YEAR)
# save(sp.optim.abund.shift, file=file.path(dir_data, "sp.optim.abund.shift.RData"))

# sp.optim.density.shift <- sp.optimloc.shift.tmp |> arrange(ORIG.ECOREGION, SPECIES, YEAR)
# save(sp.optim.density.shift, file=file.path(dir_data, "sp.optim.density.shift.RData"))

# sp.optim.abund.shift.leadtime <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.abund.shift.leadtime, file=file.path(dir_data, "sp.optim.abund.shift.leadtime.RData"))

# sp.optim.density.shift.leadtime <- sp.optimloc.shift.tmp |> arrange(ECOREGION, SPECIES, YEAR)
# save(sp.optim.density.shift.leadtime, file=file.path(dir_data, "sp.optim.density.shift.leadtime.RData"))
	
####

# LONGTERM TRENDS ---------------
# setwd(file.path(dir_data, "Optim_abund_trend_firstsite_longterm"))
# setwd(file.path(dir_data, "Optim_abund_trend_previousSite_longterm"))

# setwd(file.path(dir_data, "Optim_density_trend_firstsite_longterm"))
setwd(file.path(dir_data, "Optim_density_trend_previousSite_longterm"))

file_list <- as.list(list.files())

trend.tmp <- foreach(i=1:length(file_list), .combine="rbind") %dopar% {
	
	mget(load(file_list[[i]]))[[1]]
	#out <- mget(load(file_list[[i]]))[[1]]
	#if(is.data.frame(out)) {return(out)}
}

optimloc.trend.tmp <- trend.tmp |>
		arrange(REALM,SPECIES,FEEDING.TYPE,TEMP.OPTIM,THERMAL.GUILD,
				mod.parms, FlagAnalysis) # |> drop_na(SE)

setwd(dir_data)

# optim.abund.firstsite.trend.longterm <- optimloc.trend.tmp
# save(optim.abund.firstsite.trend.longterm, file="optim.abund.firstsite.trend.longterm.RData")

# optim.abund.previoussite.trend.longterm <- optimloc.trend.tmp
# save(optim.abund.previoussite.trend.longterm, file="optim.abund.previoussite.trend.longterm.RData")

# optim.density.firstsite.trend.longterm <- optimloc.trend.tmp
# save(optim.density.firstsite.trend.longterm, file="optim.density.firstsite.trend.longterm.RData")

optim.density.previoussite.trend.longterm <- optimloc.trend.tmp
# save(optim.density.previoussite.trend.longterm, file="optim.density.previoussite.trend.longterm.RData")

# LEADTIME LONGTERM TRENDS ---------------
# setwd(file.path(dir_data, "Optim_abund_leadtime_trend_firstsite_longterm"))
# setwd(file.path(dir_data, "Optim_abund_leadtime_trend_previoussite_longterm"))

# setwd(file.path(dir_data, "Optim_density_leadtime_trend_firstsite_longterm"))
setwd(file.path(dir_data, "Optim_density_leadtime_trend_previoussite_longterm"))


file_list <- as.list(list.files())

trend.tmp <- foreach(i=1:length(file_list), .combine="rbind") %dopar% {
	
	mget(load(file_list[[i]]))[[1]]
	#out <- mget(load(file_list[[i]]))[[1]]
	#if(is.data.frame(out)) {return(out)}
}

optimloc.trend.tmp <- trend.tmp |>
		arrange(REALM,SPECIES,FEEDING.TYPE,TEMP.OPTIM,THERMAL.GUILD,
				mod.parms, FlagAnalysis) # |> drop_na(SE)

setwd(dir_data)

# optim.abund.leadtime.firstsite.trend.longterm <- optimloc.trend.tmp
# save(optim.abund.leadtime.firstsite.trend.longterm, file="optim.abund.leadtime.firstsite.trend.longterm.RData")

# optim.abund.leadtime.previoussite.trend.longterm <- optimloc.trend.tmp
# save(optim.abund.leadtime.previoussite.trend.longterm, file="optim.abund.leadtime.previoussite.trend.longterm.RData")

# optim.density.leadtime.firstsite.trend.longterm <- optimloc.trend.tmp
# save(optim.density.leadtime.firstsite.trend.longterm, file="optim.density.leadtime.firstsite.trend.longterm.RData")

optim.density.leadtime.previoussite.trend.longterm <- optimloc.trend.tmp
# save(optim.density.leadtime.previoussite.trend.longterm, file="optim.density.leadtime.previoussite.trend.longterm.RData")


#### ---- SHORT-TERM TRENDS in optimal location ---- ####

#setwd(file.path(dir_data, "Optim_abund_trend"))
#setwd(file.path(dir_data, "Optim_density_trend"))
#setwd(file.path(dir_data, "Optim_abund_leadtime_trend"))
#setwd(file.path(dir_data, "Optim_density_leadtime_trend"))

### optimloc previous site ###
#setwd(file.path(dir_data, "Optim_abund_leadtime_prevsite_trend"))
setwd(file.path(dir_data, "Optim_density_leadtime_prevsite_trend"))

file_list <- as.list(list.files())

trend.tmp <- foreach(i=1:length(file_list), .combine="rbind") %dopar% {
	
	mget(load(file_list[[i]]))[[1]]
	#out <- mget(load(file_list[[i]]))[[1]]
	#if(is.data.frame(out)) {return(out)}
}

optimloc.trend.tmp <- trend.tmp |>
		arrange(REALM,SPECIES,FEEDING.TYPE,TEMP.OPTIM,THERMAL.GUILD,
				mod.parms, FlagAnalysis) # |> drop_na(SE)

setwd(dir_data)

#optim.abund.trend <- optimloc.trend.tmp
#save(optim.abund.trend, file="optim.abund.trend.RData")
#optim.density.trend <- optimloc.trend.tmp
#save(optim.density.trend, file="optim.density.trend.RData")
#optim.abund.leadtime.firstsite.trend <- optimloc.trend.tmp
#save(optim.abund.leadtime.firstsite.trend, file="optim.abund.leadtime.firstsite.trend.RData")
#optim.density.leadtime.firstsite.trend <- optimloc.trend.tmp
#save(optim.density.leadtime.firstsite.trend, file="optim.density.leadtime.firstsite.trend.RData")

# optim.abund.leadtime.trend <- optimloc.trend.tmp
# save(optim.abund.leadtime.trend, file="optim.abund.leadtime.trend.RData")
# optim.density.leadtime.trend <- optimloc.trend.tmp
# save(optim.density.leadtime.trend, file="optim.density.leadtime.trend.RData")

### Save optimloc previous site ###
# optim.abund.leadtime.prevsite.trend <- optimloc.trend.tmp
# save(optim.abund.leadtime.prevsite.trend, file="optim.abund.leadtime.prevsite.trend.RData")
optim.density.leadtime.prevsite.trend <- optimloc.trend.tmp
# save(optim.density.leadtime.prevsite.trend, file="optim.density.leadtime.prevsite.trend.RData")


#### All sites temperatures to calculate warming trends at shfted sites; ####
#### thetao values at all shifted sites from1 1993 to the last peak year ####
setwd(file.path(dir_data, "thetao_all_sites"))

file_list1 <- as.list(list.files())

thetao.tmp <- foreach(i=1:length(file_list1), .combine='rbind') %dopar% {
	mget(load(file_list1[[i]]))[[1]]
}

temperature_all_sites <- thetao.tmp

# -- Apply the Bassian/Hawaii ecoregion split to the site rows --------------
# The thetao_all_sites files carry the CLEAN species name, but downstream warming
# trends are joined to the split-suffixed populations (…_W/_E/_NW). Suffix each
# split-ecoregion site by its W/E/NW group (nearest split_plan$sites site, matched
# WITHIN the ecoregion) and rebuild SPECIES; non-split ecoregions are unchanged.
# Without this, Bassian/Hawaii populations get NA warming in Fig 4.
load(input_file("split_plan.RData"))
suffix_split_sites <- function(df, plan) {
  df$.suffix <- NA_character_
  for (E in unique(plan$ECOREGION)) {
    idx <- which(df$ECOREGION == E)
    if (!length(idx)) next
    pl <- plan[plan$ECOREGION == E, ]
    nn <- RANN::nn2(pl[, c("LON","LAT")],
                    df[idx, c("SHIFTED_LON","SHIFTED_LAT")], k = 1)$nn.idx[, 1]
    df$.suffix[idx] <- pl$suffix[nn]
  }
  df$SPECIES <- ifelse(is.na(df$.suffix), df$SPECIES, paste0(df$SPECIES, "_", df$.suffix))
  df$.suffix <- NULL
  df
}
temperature_all_sites <- suffix_split_sites(temperature_all_sites, split_plan$sites)

save(temperature_all_sites, file=file.path(dir_data, "temperature_all_sites.RData"))

# -------------------------------------------------------------------------------- #

#### All ecoregions temperatures to calculate warming trends for ecoregions; ####
#### thetao values at all shifted sites from1 1993 to 2021 ####
setwd(file.path(dir_data, "thetao_all_ecoregions_coastal"))

file_list1 <- as.list(list.files())

thetao.tmp <- foreach(i=1:length(file_list1), .combine='rbind') %dopar% {
	mget(load(file_list1[[i]]))[[1]]
}

temperature_all_ecoregions_coastal <- thetao.tmp
save(temperature_all_ecoregions_coastal, file=file.path(dir_data, "temperature_all_ecoregions_coastal.RData"))

setwd(file.path(dir_data, "thetao_all_ecoregions_all"))

file_list1 <- as.list(list.files())

thetao.tmp <- foreach(i=1:length(file_list1), .combine='rbind') %dopar% {
	mget(load(file_list1[[i]]))[[1]]
}

temperature_all_ecoregions_all <- thetao.tmp
save(temperature_all_ecoregions_all, file=file.path(dir_data, "temperature_all_ecoregions_all.RData"))

