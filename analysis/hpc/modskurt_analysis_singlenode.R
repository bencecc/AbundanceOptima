#### ---- modskurt optima within ecoregions: SINGLE-NODE version (one machine, many cores) ---- ####
#### ---- Same analysis as modskurt_analysis.R, but instead of one species per HPC task (sp.id  ---- ####
#### ---- from the command line) it loops over ALL species in parallel on one node: species are ---- ####
#### ---- distributed over the cores (doMC/fork) and the years within a species run sequentially ---- ####
#### ---- (the reverse of the HPC version, where a task = one species and years run in parallel). ---- ####
#### ---- Outputs are identical: one <species>.RData per species in out_optim (+ plot data).    ---- ####
#### ---- Usage: Rscript modskurt_analysis_singlenode.R [n_cores] [species-id file]              ---- ####
#### ----   n_cores          default 120 (one worker per species; each fit is single-threaded)   ---- ####
#### ----   species-id file  optional, one id per line (e.g. spID.txt / spID_pooledEAus.txt);   ---- ####
#### ----                    default = all species. Already-finished species are skipped, so    ---- ####
#### ----                    an interrupted run can simply be restarted.                        ---- ####
require(dplyr)
require(cmdstanr)
require(modskurt1)

require(foreach, quietly=T)
require(doMC, quietly=T)
source("config.R")
args    <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 120L
registerDoMC(cores = n_cores)

# load data -- two alternatives:
# 1) MAIN ANALYSIS, one population per species x ecoregion (with the split_plan sub-regions):
load(input_file("sp.df.RData"))
# 2) SENSITIVITY "pooled East Australia": species with sites in >= 2 adjacent ecoregions of
#    the Torres Strait -> Bassian (eastern group) coast, relabelled to ONE ecoregion so the
#    optimum is fit over the whole 29-degree gradient (co-author query on the ecoregion split).
#    Uncomment the two lines below instead of 1); keep run_ecoregions <- NULL (a Bassian/Hawaii
#    subset would empty the pooled data); the split_plan block below then matches nothing and
#    leaves SPECIES unchanged. ALSO switch the output folders (out_optim / out_plot below) to
#    the _pooledEAus alternatives, or the main outputs are overwritten, and point modskurt1.sh
#    at spID_pooledEAus.txt (1..624) in its last line ("done < ...") -- do not overwrite spID.txt.
#load(input_file("sp.df.pooled.EAus.RData"))
#sp.df <- sp.df.pooled.EAus
#load(input_file("species.id.RData"))

# ---- Apply the ecoregion sub-regions from split_plan --------------------------------
# split_plan (from cluster_ecoregions.R) flags which ecoregions are split and maps each
# site to its group. For every split ecoregion (Bassian + Hawaii) the sub-group suffix
# (directional: _W/_E/_NW) goes on SPECIES ("X" -> "X_W"), ECOREGION keeps its base name,
# and the clean name is stored in SPECIES.ORIG. Each suffixed SPECIES is then fit as its
# own population, so the output already has ECOREGION = "Bassian", SPECIES = "X_W",
# SPECIES.ORIG = "X" — no downstream relabelling needed.
load(input_file("split_plan.RData"))

# Select parameters
resp.var <- 'abund'
# save data for later plotting?
save.plot.data <- TRUE 
# plot the results as they are generated?
plot <- FALSE
# set the minumum number of years in a timeseries to proceed
nyrs <- 4

# RUN SCOPE (keeps the original "all or a subset of ecoregions" flexibility):
#   run_ecoregions <- NULL          -> ALL ecoregions (full re-run; splits auto-applied).
#   run_ecoregions <- c("Bassian")  -> only these BASE ecoregions (re-fit just what split).
run_ecoregions <- NULL # NULL = ALL ecoregions (full run); c("Bassian","Hawaii") = re-fit only the split ecoregions

if (!is.null(run_ecoregions)) sp.df <- sp.df |> filter(ECOREGION %in% run_ecoregions)

split_ecos <- unique(split_plan$ecoreg$ECOREGION[split_plan$ecoreg$n_groups > 1])
ckey <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
plan_sites <- split_plan$sites |> mutate(coord_key = ckey(LON, LAT)) |>
		dplyr::select(ECOREGION, coord_key, group, suffix)     # suffix = directional label (_W/_E/_NW)
sp.df <- sp.df |>
		mutate(coord_key = ckey(LON, LAT)) |>
		left_join(plan_sites, by = c("ECOREGION", "coord_key")) |>
		mutate(SPECIES.ORIG = SPECIES,                                  # clean name for trait/STI joins
		       SPECIES = ifelse(ECOREGION %in% split_ecos & !is.na(group),
						paste0(SPECIES, "_", suffix), SPECIES)) |>      # e.g. "X_W"; ECOREGION unchanged
		dplyr::select(-coord_key, -group, -suffix)

# Output folders: SINGLE fixed names (as in modskurt_analysis_old). Every species'
# files go into ONE optim folder + ONE plot folder; the ecoregion is encoded in the
# filename (_ecoreg_<i>_), not as a folder. RENAME these manually after each run
# (e.g. -> ModskurtOptimEcoreg_Abund / _Density) so abund and density don't collide,
# since both runs write to the same fixed folder.
modskurt_dir <- dir_data
out_optim <- file.path(modskurt_dir, "ModskurtOptimEcoreg")          # 1) main analysis
out_plot  <- file.path(modskurt_dir, "ModskurtOptimEcoregPlot")
#out_optim <- file.path(modskurt_dir, "ModskurtOptimEcoreg_pooledEAus")     # 2) pooled sensitivity
#out_plot  <- file.path(modskurt_dir, "ModskurtOptimEcoregPlot_pooledEAus")
                   

species.id <- sp.df |> distinct(SPECIES) 

# generate and set output directory (per sub-region; see out_optim/out_plot above)
if(!dir.exists(out_optim)) dir.create(out_optim, recursive = TRUE)
if(!dir.exists(out_plot))  dir.create(out_plot,  recursive = TRUE)

stan_path <- system.file("stan", "modskurt1.stan", package = "modskurt1")
mod <- cmdstanr::cmdstan_model(stan_path)
#mod <- cmdstan_model(stan_path, compile = FALSE)
#mod$compile(dir = "/tmp")

if(resp.var=="density") {
	
	sp.df <- sp.df |>
			mutate(density=abund/SAMPLED_AREA)

}

# species to run: all, or the ids listed in the optional file (one per line); skip the ones
# whose output file already exists (restartable)
sp.ids <- if (length(args) >= 2) as.integer(readLines(args[2])) else seq_len(nrow(species.id))
done   <- file.exists(file.path(out_optim, paste0(species.id$SPECIES[sp.ids], ".RData")))
sp.ids <- sp.ids[!done]
cat("species to run:", length(sp.ids), "(", sum(done), "already done ) on", n_cores, "cores
")

# ---- OUTER PARALLEL LOOP over species (each worker writes its own files; errors in one
# species are recorded and do not stop the others) ----------------------------------------
run_log <- foreach(sp.id = sp.ids, .errorhandling = "pass") %dopar% {

sp.name <- species.id$SPECIES[sp.id]
cat(sp.id, sp.name, "
")

test.dat <- sp.df |> filter(SPECIES%in%species.id$SPECIES[sp.id])

# Ecoregions
ecoreg <- unique(test.dat$ECOREGION)

sp.optim.res <- foreach(i=1:length(ecoreg), .combine="rbind") %do% {
	
	sp.ecoreg <- test.dat |> filter(ECOREGION%in%ecoreg[i])
	
	yr.df <- sp.ecoreg |> group_by(YEAR) |>
			summarise(n=n(), .groups="drop")
	
	unique.sp <- unique(sp.ecoreg$SPECIES)
	
	latlon.range <- sp.df |> filter(SPECIES%in%unique.sp) |>
			reframe(max.lat=max(LAT), min.lat=min(LAT),
					max.lon=max(LON), min.lon=min(LON))
	max.resp <- sp.df |> filter(SPECIES%in%unique.sp) |>
			reframe(max.resp=max(!!sym(resp.var)))
	
	if(!is.null(yr.df)) {
		
		mod.res <- foreach(j=1:nrow(yr.df), .combine="rbind") %do% {   # sequential: the parallelism is over species
			
			cat('Doing YEAR ', j, ' of ', nrow(yr.df),
					' for SPECIES ', i, ' of ', length(sp.id), '\n', sep = '')
			
			# proceed if there are at least n (>1 or >4) years
			yrs.check <- yr.df[j,"n"]
			
			if(yrs.check$n>nyrs) {
				
				yr <- yr.df[j, "YEAR"]
				sub.df <- sp.ecoreg |> filter(YEAR%in%yr$YEAR)
				
				# check for non-zeros abundances
				perc.not.zeros <- nrow(sub.df |> filter(!!sym(resp.var)>0))/nrow(sub.df)					
				
				# proceed if there are at least 5 sites and
				# the percentage of non-zero values is > 0
				if(nrow(sub.df)>4&perc.not.zeros>0) {						
					
					sub.test.dat <- as.data.frame(sub.df)
										
					mod.lat <- modskurt_flow(df=sub.test.dat, resp.var=resp.var,
							pred.var='LAT', mod=mod, latlon.range=latlon.range,
							max.resp=max.resp, plot.res=plot, custom.dist=NULL) 
					
					mod.lon <- modskurt_flow(df=sub.test.dat, resp.var=resp.var,
							pred.var='LON', mod=mod, latlon.range=latlon.range,
							max.resp=max.resp, plot.res=plot, custom.dist=NULL) 
					
					if(!is.null(mod.lat)&!is.null(mod.lon)) {
						
						mod.comb <- mod.lat[[1]] |>
								left_join(mod.lon[[1]][1,c("SPECIES","H.LON","m.LON","SELECTED_MODEL_LON")],
										by="SPECIES") |>
								relocate(c("H.LON","m.LON"), .after=m.LAT) |>
								relocate("SELECTED_MODEL_LON", .after=SELECTED_MODEL_LAT)
						
						if(save.plot.data==TRUE) {
							
							# save LAT summary data for plotting
							lat.df <- mod.lat[[2]]
							plot.lat.name <- file.path(out_plot,
									paste(gsub(" ", ".", sp.name), "_ecoreg_", i, "_", yr$YEAR, "_LAT.RData", sep=""))
							save(lat.df, file=plot.lat.name)
							
							# save LON summary data for plotting
							lon.df <- mod.lon[[2]]
							plot.lon.name <- file.path(out_plot,
									paste(gsub(" ", ".", sp.name), "_ecoreg_", i, "_", yr$YEAR, "_LON.RData", sep=""))
							save(lon.df, file=plot.lon.name)
							
						}
						
						# return model output
						return(mod.comb)
						
					} else if (!is.null(mod.lat)&is.null(mod.lon)) {
						
						if(save.plot.data==TRUE) {
							
							# save LAT summary data for plotting
							lat.df <- mod.lat[[2]]
							
							plot.lat.name <- file.path(out_plot,
									paste(gsub(" ", ".", sp.name), "_ecoreg_", i, "_", yr$YEAR, "_LAT.RData", sep=""))
							save(lat.df, file=plot.lat.name)
							
						}
						# return model output
						return(mod.lat)
						
					} else if (is.null(mod.lat)&!is.null(mod.lon)) {
						
						if(save.plot.data==TRUE) {
							
							# save LON summary data for plotting
							lon.df <- mod.lon[[2]]
							plot.lon.name <- file.path(out_plot,
									paste(gsub(" ", ".", sp.name), "_ecoreg_", i, "_", yr$YEAR, "_LON.RData", sep=""))
							save(lon.df, file=plot.lon.name)
							
						}
						
						# return model output
						return(mod.lon)
						
					} else (NULL)
					
				}
				
			} 
			
		}
		
	}	
	
	if(!is.null(mod.res)) mod.res |> tidyr::drop_na(m.LAT)
}	

if(!is.null(sp.optim.res)) {

	# carry the clean species name (constant for this sp.id) so the output has
	# ECOREGION=base, SPECIES=suffixed, SPECIES.ORIG=clean
	sp.optim.res$SPECIES.ORIG <- test.dat$SPECIES.ORIG[1]

	outputName=paste(sp.name, ".RData", sep="")
	outputPath=file.path(out_optim, outputName)
	assign(paste(sp.name), value=sp.optim.res)
	save(list=paste(sp.name), file=outputPath)

}

sp.name   # value returned to run_log
}   # ---- end of the outer species loop ----

# report species whose worker raised an error (run_log holds the error object for them)
errs <- sapply(run_log, function(x) inherits(x, "error"))
cat("finished:", sum(!errs), "species ok,", sum(errs), "with errors
")
if (any(errs)) for (k in which(errs)) cat("  species id", sp.ids[k], ":", conditionMessage(run_log[[k]]), "
")

