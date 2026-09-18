##############################################################################
# all_ecoregions_temperature.R
#
# Extract ECOREGION-level warming rates (1993-2021) from the GLORYS raster.
#
# Design:
#   - Clip the GLORYS raster to each focal MEOW ecoregion polygon
#   - For every grid cell inside the ecoregion:
#       * extract daily subsurface temperature 1993-2021
#       * compute annual means
#       * fit OLS slope (°C / yr) of annual mean temp vs YEAR
#   - Aggregate grid-cell slopes to ecoregion-level warming rate (mean + SE)
#
# Input:
#   - raster_glorys_daily_1993_2021_5_10_metres_mean.tif   (GLORYS daily raster)
#   - Marine_Ecoregions_Of_the_World__MEOW_.shp            (MEOW polygons)
#   - focal_ecoregions.RData                               (character vector of
#                                                           focal ecoregion names,
#                                                           or taken from sp.focal)
#
# Output:
#   - per-chunk .RData written into thetao_all_ecoregions_<coastal|all>/, reassembled
#     into temperature_all_ecoregions_{coastal,all}.RData by cluster_summaries.R.
#     data.frame with columns:
#       ECOREGION, warming_rate, warming_se, n_cells, n_years_mean
#
# HPC (doMC). Chunked by ecoregion index via CLI args.
# Usage:
#   Rscript all_ecoregions_temperature.R [lw] [up]
##############################################################################

require(dplyr)

require(tidyr)
require(sf)
require(terra)
require(foreach, quietly = TRUE)
require(doMC)
source("config.R")
registerDoMC(cores = 5)

# ---- 1. Paths ---------------------------------------------------------
# prefix name of environmental variable in species index file
env.name <- "thetao"
#  name for output files
out.name <- "_all_ecoregions"

BASE_DIR   <- dir_data
MEOW_PATH  <- meow_shapefile
GLORYS_TIF <- glorys_daily_mean_tif
GEBCO_NC   <- gebco_file

# ---- Spatial domain flag: "coastal" (shelf only) or "all" (whole ecoregion) ----
# Can be overridden on the command line as the 3rd argument:
#   Rscript all_ecoregions_temperature.R [lw] [up] [domain]
# e.g.:  Rscript ... 1 50 coastal      →  apply -30 - 0 m bathymetry mask
#        Rscript ... 1 50 all          →  use all cells inside the polygon
DOMAIN <- "all"    # c("coastal", "all")

# Depth range — defined only when DOMAIN == "coastal"
if (DOMAIN == "coastal") {
  DEPTH_MIN <- -30     # deepest allowed (shallow reef zone)
  DEPTH_MAX <- 0       # shallowest allowed (coastline; excludes land = positive)
}

# Output directory (suffix by domain so the two runs don't overwrite each other)
out.dir.suffix <- paste0(out.name, "_", DOMAIN)
if(!file.exists(file.path(dir_data, paste0(env.name, out.dir.suffix))))
	dir.create(file.path(dir_data, paste0(env.name, out.dir.suffix)))
OUT_DIR <- file.path(dir_data, paste0(env.name, out.dir.suffix))

# ---- 2. Identify focal ecoregions   ---------------------------------------------
# Use the ecoregions present in the modskurt fish dataset
load(file.path(BASE_DIR, "sp.optim.abund.shift.RData"))
focal_ecos <- sort(unique(sp.optim.abund.shift$ORIG.ECOREGION))
cat("Focal ecoregions:", length(focal_ecos), "\n")

# ---- 3. Load GLORYS raster and identify annual layer groups ---------------------
r <- terra::rast(GLORYS_TIF)
time.vec  <- terra::time(r)
years_all <- as.numeric(format(time.vec, "%Y"))

# Keep layers in 1993-2021
keep_lyrs <- which(years_all >= 1993 & years_all <= 2021)
r <- terra::subset(r, keep_lyrs)
years_all <- years_all[keep_lyrs]
year_layers <- split(seq_along(years_all), years_all)
available_years <- as.numeric(names(year_layers))
cat("GLORYS years:", min(available_years), "-", max(available_years), "\n")
cat("GLORYS CRS:", terra::crs(r, describe = TRUE)$name, "\n")
cat("GLORYS extent:", paste(as.vector(terra::ext(r)), collapse = ", "), "\n")

# ---- 3b. GEBCO bathymetry (only announced if coastal domain) ----
# Each worker opens GEBCO lazily inside the foreach loop, only when needed.
if (DOMAIN == "coastal") {
  cat("GEBCO path:", GEBCO_NC, "(opened per-worker, coastal mask active)\n")
}

# Flag longitude convention: GLORYS may be in 0-360 or -180 to 180
glorys_ext <- as.vector(terra::ext(r))
glorys_0_360 <- glorys_ext[2] > 180   # xmax > 180 → 0-360 convention

# ---- 4. Load MEOW polygons and reconcile CRS/longitude with raster -----------------
meow_raw <- sf::read_sf(MEOW_PATH)
cat("MEOW CRS:", sf::st_crs(meow_raw)$Name, "\n")

# Reproject MEOW to GLORYS CRS (both should be WGS84 but make sure)
meow_raw <- sf::st_transform(meow_raw, terra::crs(r))

# Wrap longitudes if raster uses 0-360 convention (shift dateline)
if (glorys_0_360) {
  cat("Shifting MEOW longitudes to 0-360 convention\n")
  meow_raw <- sf::st_wrap_dateline(meow_raw,
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  # st_shift_longitude moves features so lon is in [0, 360]
  meow_raw <- sf::st_shift_longitude(meow_raw)
}

# Make geometries valid (MEOW has some complex multipolygons)
meow_raw <- sf::st_make_valid(meow_raw)

# Filter to focal ecoregions
meow <- meow_raw |> filter(ECOREGION %in% focal_ecos)
cat("MEOW polygons matched:", nrow(meow), "\n")

missing <- setdiff(focal_ecos, meow$ECOREGION)
if (length(missing) > 0) {
  cat("WARNING — ecoregions missing from MEOW shapefile:\n")
  print(missing)
}

# Report each polygon's extent vs raster extent (diagnostic)
for (k in seq_len(min(3, nrow(meow)))) {
  pe <- as.vector(sf::st_bbox(meow[k, ]))
  cat(sprintf("Sample poly %s bbox: %s\n",
    meow$ECOREGION[k], paste(round(pe, 2), collapse = ", ")))
}

# Wrap raster for parallel workers
r_wrap <- terra::wrap(r)

# -- 5. CLI args for cluster chunking and domain flag -------------
#   arg 1 : lower ecoregion index
#   arg 2 : upper ecoregion index
#   arg 3 : domain ("coastal" or "all") — optional, overrides the default
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  lw <- as.numeric(args[1])
  up <- as.numeric(args[2])
} else {
  lw <- 1
  up <- nrow(meow)
}
# arg[3] is the container name
if (length(args) >= 4) {
  DOMAIN <- tolower(args[4])
  stopifnot(DOMAIN %in% c("coastal", "all"))
  # Rebuild output directory with overridden domain
  out.dir.suffix <- paste0(out.name, "_", DOMAIN)
  OUT_DIR <- file.path(dir_data, paste0(env.name, out.dir.suffix))
  if (!file.exists(OUT_DIR)) dir.create(OUT_DIR)
  # Ensure depth range exists if CLI switched us into coastal mode
  if (DOMAIN == "coastal" && !exists("DEPTH_MIN")) {
    DEPTH_MIN <- -30
    DEPTH_MAX <- 0
  }
}
cat("Processing ecoregions", lw, "to", up, "| domain =", DOMAIN, "\n")
if (DOMAIN == "coastal") {
  cat("Applying bathymetry mask: depth", DEPTH_MIN, "to", DEPTH_MAX, "m\n")
} else {
  cat("No bathymetry mask — using all cells inside the ecoregion polygon\n")
}

# -- 6. Helper: compute per-cell annual means inside a polygon -------------
# If domain == "coastal", apply the bathymetry depth mask; otherwise skip it
# and use all cells inside the polygon.
cell_annual_means <- function(r_local, bathy_local = NULL, poly,
                              years_all, year_layers,
                              depth_min = NULL, depth_max = NULL,
                              domain = DOMAIN) {

  # Canonical CRS for this worker — use the GLORYS WKT as ground truth
  ref_crs <- terra::crs(r_local)

  # Convert polygon to SpatVector and force CRS to match raster
  poly_v <- terra::vect(poly)
  if (!terra::same.crs(poly_v, r_local)) {
    poly_v <- terra::project(poly_v, ref_crs)
  }
  terra::crs(poly_v) <- ref_crs   # harmonise WKT string exactly

  # Crop + mask the raster to the ecoregion polygon.
  # tryCatch handles the rare case where a polygon lies entirely outside
  # the raster extent (e.g., polar ecoregion below GLORYS lat range).
  r.poly <- tryCatch({
      terra::crop(r_local, poly_v, snap = "out") |>
        terra::mask(poly_v)
    }, error = function(e) {
      cat("  → crop/mask failed for this polygon:", conditionMessage(e), "\n")
      NULL
    })
  if (is.null(r.poly)) return(NULL)
  terra::crs(r.poly) <- ref_crs

  # -- Bathymetry mask (only if domain == "coastal") ------------------------
  if (domain == "coastal") {
    # Reproject GEBCO to the reference CRS if needed, then crop + resample.
    if (!terra::same.crs(bathy_local, r_local)) {
      bathy_local <- terra::project(bathy_local, ref_crs)
    }
    bathy.poly <- tryCatch({
        terra::crop(bathy_local, poly_v, snap = "out") |>
          terra::mask(poly_v) |>
          terra::resample(terra::subset(r.poly, 1), method = "average")
      }, error = function(e) {
        cat("  → bathymetry crop/resample failed:", conditionMessage(e), "\n")
        NULL
      })
    if (is.null(bathy.poly)) return(NULL)
    terra::crs(bathy.poly) <- ref_crs

    coastal_mask <- bathy.poly >= depth_min & bathy.poly <= depth_max
    terra::crs(coastal_mask) <- ref_crs
    r.poly <- terra::mask(r.poly, coastal_mask, maskvalue = FALSE)
  }

  # Grid cell coordinates (non-NA cells only)
  vals_ref <- terra::values(terra::subset(r.poly, 1))
  cell_ok  <- which(!is.na(vals_ref))
  if (length(cell_ok) == 0) {
    msg <- if (domain == "coastal")
      paste0("  → no coastal cells (", depth_min, " to ", depth_max, " m)")
    else
      "  → no valid cells in polygon"
    cat(msg, "\n")
    return(NULL)
  }

  xy <- terra::xyFromCell(r.poly, cell_ok)

  # Daily value matrix: rows = cells, cols = days
  day_vals <- terra::extract(r.poly, xy)      # data.frame [ncells × ndays]

  # Compute annual means per cell
  ann <- vapply(seq_along(year_layers), function(k) {
    lyr <- year_layers[[k]]
    rowMeans(day_vals[, lyr, drop = FALSE], na.rm = TRUE)
  }, numeric(length(cell_ok)))
  colnames(ann) <- names(year_layers)

  data.frame(
    cell_id = cell_ok,
    lon     = xy[, 1],
    lat     = xy[, 2],
    ann,                                      # one column per year
    check.names = FALSE)
}

# ---- 7. Main loop: per ecoregion → per-cell slopes → ecoregion mean -------------
results <- foreach(i = lw:up, .combine = "rbind",
                   .packages = c("terra", "sf", "dplyr"),
                   .errorhandling = "pass",
                   .inorder = TRUE) %dopar% {

  tryCatch({

    eco_name <- meow$ECOREGION[i]
    cat("Ecoregion", i, ":", eco_name, "\n")

    poly <- meow[i, ]
    r_local <- terra::unwrap(r_wrap)

    # Only open GEBCO when needed (coastal domain)
    bathy_local <- NULL
    if (DOMAIN == "coastal") {
      bathy_local <- terra::rast(GEBCO_NC)
      terra::crs(bathy_local) <- "epsg:4326"
    }

    if (DOMAIN == "coastal") {
      ann.df <- cell_annual_means(r_local, bathy_local, poly,
                                   years_all, year_layers,
                                   DEPTH_MIN, DEPTH_MAX,
                                   domain = DOMAIN)
    } else {
      ann.df <- cell_annual_means(r_local, poly = poly,
                                   years_all = years_all,
                                   year_layers = year_layers,
                                   domain = DOMAIN)
    }
    if (is.null(ann.df)) return(NULL)

  # Long format: cell × year → fit OLS per cell
  ann.long <- ann.df |>
    tidyr::pivot_longer(cols = -c(cell_id, lon, lat),
                        names_to = "YEAR", values_to = "temp") |>
    mutate(YEAR = as.numeric(YEAR)) |>
    filter(!is.nan(temp) & !is.na(temp))

  cell_slopes <- ann.long |>
    group_by(cell_id, lon, lat) |>
    filter(n() >= 10) |>
    summarise(
      n_years    = n(),
      cell_slope = coef(lm(temp ~ YEAR))[2],
      cell_se    = summary(lm(temp ~ YEAR))$coefficients[2, 2],
      .groups    = "drop")

    if (nrow(cell_slopes) == 0) return(NULL)

    # Aggregate to ecoregion
    out <- data.frame(
      ECOREGION    = eco_name,
      warming_rate = mean(cell_slopes$cell_slope, na.rm = TRUE),
      warming_se   = sd(cell_slopes$cell_slope, na.rm = TRUE) / sqrt(nrow(cell_slopes)),
      n_cells      = nrow(cell_slopes),
      n_years_mean = mean(cell_slopes$n_years))

    out

  }, error = function(e) {
    cat("Ecoregion", i, "failed:", conditionMessage(e), "\n")
    NULL
  })
}

# ---- 8. Save combined output ---------------------------------------------

assign(paste(env.name, lw, '_', up, sep=""), value=results, pos=1, inherits=T)
outputName=paste(env.name, lw, '_', up, ".RData",sep="")
outputPath=file.path(OUT_DIR, outputName)
save(list=paste(env.name, lw, '_', up, sep=""), file=outputPath)
