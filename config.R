# ============================================================================
# config.R — single source of every path used by the AbundanceOptima analysis.
# Source this at the top of any script:   source("config.R")
# It DEFINES path variables; it does not run anything. Scripts then read files
# via these variables (dir_data, env_dir, gebco_file, input_file(), ...) so no
# absolute path is ever hard-coded in a script.
# ============================================================================

# ---- EDIT THIS ONE LINE to point at your project root --------------------
PROJ <- path.expand("~/abundance_optima")
# For development on this machine, uncomment the override block at the BOTTOM
# of this file instead (points PROJ + data/grids at their real locations).
# --------------------------------------------------------------------------

dir_data     <- file.path(PROJ, "data")      # inputs: provided .RData (Zenodo); also the working data dir
dir_packages <- file.path(PROJ, "packages")  # modskurt1 + CockR sources
dir_out      <- file.path(PROJ, "output")    # freshly generated .RData (optional; see input_file)
dir_results  <- file.path(PROJ, "results")   # FINAL figures & tables (each script's own subfolder)
dir_figures  <- file.path(PROJ, "figures")   # figure scripts + helpers (plot_sites.R)
for (d in c(dir_out, dir_results)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- Large user-supplied grids (NOT redistributed — see README) ----------
# Point env_dir (and the files) at wherever you stored the big GLORYS/GEBCO/MEOW files.
env_dir        <- dir_data   # GLORYS rasters live here (edit if you keep them elsewhere)
gebco_file     <- file.path(dir_data, "GEBCO_2024.nc")
meow_shapefile <- file.path(dir_data, "Marine_Ecoregions_Of_the_World__MEOW_.shp")   # + .shx/.dbf/.prj
glorys_daily_nc        <- file.path(env_dir, "glorys_daily_1993_2021_5_10_metres_mean.nc")            # raw daily (input to make_tiled)
glorys_daily_mean_tif  <- file.path(env_dir, "raster_glorys_daily_1993_2021_5_10_metres_mean.tif")    # daily mean raster
glorys_daily_tiled_tif <- file.path(env_dir, "raster_glorys_daily_1993_2021_5_10_metres_mean_TILED.tif")  # tiled (make_tiled_from_nc.R)
glorys_annual_mean_nc  <- file.path(env_dir, "glorys_annual_1993_2021_5_10_metres_mean.nc")           # annual-mean layers

# Per-stage SCRATCH / OUTPUT folders (HPC chunk in/out dirs, consolidated saves)
# live under dir_data — a script builds them as file.path(dir_data, "<name>"),
# e.g. file.path(dir_data, "thetao_all_sites"). No separate variable needed.

# Resolve an input .RData by name: prefer a freshly generated copy in output/,
# else the provided (Zenodo / working) copy in data/.
input_file <- function(name) {
  p_out <- file.path(dir_out, name)
  if (file.exists(p_out)) p_out else file.path(dir_data, name)
}

if (!dir.exists(dir_data))
  stop("dir_data not found: ", dir_data, "\nEdit PROJ in config.R — see README.")

# ============================================================================
# DEV OVERRIDE (this project's machines) — uncomment to point at the real data.
# ----------------------------------------------------------------------------
# .od <- c("C:/Users/LBenedettiCecchi/OneDrive - University of Pisa",
#          "C:/Users/lisan/OneDrive - University of Pisa")
# .od <- .od[file.exists(.od)][1]
# PROJ           <- "G:/My Drive/workspace/BiodGlob/Fish_Thermal/AbundanceOptima"
# dir_data       <- "G:/My Drive/Lavori/MPA_timeseries/Modskurt"
# dir_packages   <- file.path(PROJ, "packages"); dir_figures <- file.path(PROJ, "figures")
# dir_out        <- file.path(PROJ, "output");   dir_results <- file.path(PROJ, "results")
# env_dir        <- file.path(.od, "Lavori/MPA_timeseries/EnvData")
# gebco_file     <- file.path(.od, "Lavori/Shorelines/GEBCO/GEBCO_2024.nc")
# meow_shapefile <- file.path(.od, "Lavori/MEOWs/Marine_Ecoregions_Of_the_World__MEOW_.shp")
# glorys_daily_nc        <- file.path(env_dir, "glorys_daily_1993_2021_5_10_metres_mean.nc")
# glorys_daily_mean_tif  <- file.path(env_dir, "raster_glorys_daily_1993_2021_5_10_metres_mean.tif")
# glorys_daily_tiled_tif <- file.path(env_dir, "raster_glorys_daily_1993_2021_5_10_metres_mean_TILED.tif")
# glorys_annual_mean_nc  <- file.path(env_dir, "glorys_annual_1993_2021_5_10_metres_mean.nc")
# ============================================================================
