##############################################################################
# extract_standardized_warming.R
#
# Extract annual mean subsurface temperature at ALL optimal peak sites
# visited by each population, continuously from 1993 to the last fish
# sampling year.
#
# Design:
#   - For each population (ECOREGION × SPECIES), identify all unique sites
#     (LAT, LON) where peaks were observed, in chronological order.
#   - For EACH site, extract daily GLORYS temperature for every year from
#     1993 to the population's last sampling year, then compute annual means.
#   - Flag each site-year with:
#       SITE_ORDER  : chronological rank of this site (1 = first peak site)
#       FISH_SAMPLED: TRUE if fish were actually surveyed at this site in
#                     this year (from the original shift data)
#
# This enables flexible counterfactual analyses:
#   - Filter SITE_ORDER == 1 → first-site counterfactual (current approach)
#   - For each site, use SITE_ORDER == (current - 1) → preceding-site
#     counterfactual (tracks movement more closely)
#   - Compare trends across counterfactual choices as sensitivity analysis
#
# Input:
#   - optim.abund.relocated.RData       (peak locations per year)
#   - idx.optim.abund.relocated.RData   (raster layer indices)
#   - reef_fish_sti_glorys.RData        (species thermal indices)
#   - raster_glorys_daily_1993_2021_5_10_metres_mean.tif (GLORYS raster)
#
# Output:
#   - standardized_warming_allsites_[lw]_[up].RData
#     Columns: ECOREGION, SPECIES, SITE_LAT, SITE_LON, SITE_ORDER,
#              YEAR, annual_mean_temp, FISH_SAMPLED
#
# Designed for Linux cluster (doMC). Adjust cores and paths as needed.
# Usage: Rscript extract_standardized_warming.R [lw] [up]
##############################################################################

require(tidyverse)
require(sf)
require(terra)

require(foreach, quietly = TRUE)
require(doMC)
registerDoMC(cores = 20)

# ── 1. Generate and set output directory ─────────────────────────────────────────────────────────────

# prefix name of environmental variable in species index file
env.name <- "thetao"
# desired name for output files
out.name <- "_all_sites"

if(!file.exists(paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, out.name, sep="")))
	dir.create(paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, out.name, sep=""))

dir <- paste("/home/lisandro/Lavori/MPA_timeseries/Modskurt/", env.name, out.name, sep="")
  
# ── 2. Load data ─────────────────────────────────────────────────────────────
load("/home/lisandro/Lavori/MPA_timeseries/Modskurt/optim.abund.relocated.RData")
load("/home/lisandro/Lavori/MPA_timeseries/Modskurt/idx.optim.abund.relocated.RData")

idx.df <- idx.optim.abund.relocated

# ── 3. Build population-level site catalogue ─────────────────────────────────
range.df <- optim.abund.relocated |>
  mutate(abund = H.LAT) |>
  left_join(idx.df, by = c("ECOREGION", "SPECIES", "LAT", "LON", "YEAR")) |>
  arrange(ECOREGION, SPECIES, YEAR)

# For each population: unique sites in order of first appearance,
# plus the years fish were actually sampled at each site
pop_sites <- range.df |>
  group_by(ECOREGION, SPECIES) |>
  mutate(
    n_fish_years = n(),
    last_year    = max(YEAR)
  ) |>
  filter(n_fish_years >= 3) |>
  ungroup()

# Site catalogue: unique (LAT, LON) per population, ordered by first visit
# SAMPLED_YEAR records the fish survey year(s) at each site (comma-separated
# if the peak returned to the same site in multiple years)
site_catalogue <- pop_sites |>
  group_by(ECOREGION, SPECIES, LAT, LON) |>
  reframe(
    first_visit  = min(YEAR),
    SAMPLED_YEAR = paste(sort(unique(YEAR)), collapse = ",")
  ) |>
  arrange(ECOREGION, SPECIES, first_visit) |>
  group_by(ECOREGION, SPECIES) |>
  mutate(SITE_ORDER = row_number()) |>
  ungroup()

# Fish-sampled flags: which (site × year) combinations have actual fish data
fish_sampled_flags <- pop_sites |>
  distinct(ECOREGION, SPECIES, LAT, LON, YEAR) |>
  mutate(FISH_SAMPLED = TRUE)

# Population-level summary (for iteration)
pop_summary <- pop_sites |>
  group_by(ECOREGION, SPECIES) |>
  reframe(
    last_year    = max(YEAR),
    n_fish_years = first(n_fish_years)
  ) 
  

#cat("Populations:", nrow(pop_summary), "\n")
#cat("Total unique sites across all populations:", nrow(site_catalogue), "\n")
#cat("Mean sites per population:",
#    round(mean(table(paste(site_catalogue$ECOREGION, site_catalogue$SPECIES))), 1), "\n")

# ── 4. Load GLORYS raster ────────────────────────────────────────────────────
r <- terra::rast("/home/lisandro/Lavori/MPA_timeseries/EnvData/raster_glorys_daily_1993_2021_5_10_metres_mean.tif")
time.vec  <- terra::time(r)
years_all <- as.numeric(format(time.vec, "%Y"))

# Lookup: for each year, which raster layers correspond to it
year_layers <- split(seq_along(years_all), years_all)
available_years <- as.numeric(names(year_layers))

#cat("GLORYS years available:", min(available_years), "-", max(available_years), "\n")

r <- wrap(r)

# ── 5. Command-line args for cluster chunking ────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  lw <- as.numeric(args[1])
  up <- as.numeric(args[2])
} else {
  lw <- 1
  up <- nrow(pop_summary)
}
cat("Processing populations", lw, "to", up, "\n")

# ── 6. Helper: extract annual mean temp at a point, with NA expansion ────────
extract_annual_temp <- function(r_local, lon, lat, lyrs, years_all, year_layers) {

  p <- sf::st_point(c(lon, lat)) |>
    sf::st_sfc(crs = "epsg:4326") |>
    terra::vect()

  r.sub <- terra::subset(r_local, subset = lyrs)
  dat <- terra::extract(r.sub, p, ID = FALSE) |> as.double()

  # If >10% NAs, expand spatial extent (1×res in all directions)
  if (sum(is.na(dat)) >= length(dat) * 0.1) {
    res_r <- res(r_local)[1]
    m1 <- matrix(c(lon - res_r, lon, lon + res_r,
                    lat - res_r, lat, lat + res_r),
                 ncol = 2, nrow = 3)
    p1 <- sf::st_multipoint(m1) |>
      sf::st_sfc(crs = "epsg:4326") |>
      terra::vect()
    dat.tmp <- terra::extract(r.sub, p1, ID = FALSE)

    # If still all NA, expand to 2×res
    if (all(is.na(dat.tmp))) {
      m2 <- matrix(c(lon - 2*res_r, lon - res_r, lon, lon + res_r, lon + 2*res_r,
                      lat - 2*res_r, lat - res_r, lat, lat + res_r, lat + 2*res_r),
                   ncol = 2, nrow = 5)
      p2 <- sf::st_multipoint(m2) |>
        sf::st_sfc(crs = "epsg:4326") |>
        terra::vect()
      dat.tmp <- terra::extract(r.sub, p2, ID = FALSE)
    }

    dat <- apply(dat.tmp, 2, function(x) mean(x, na.rm = TRUE))
  }

  # Too many NAs even after expansion
  if (sum(is.na(dat)) >= length(dat) * 0.1) return(NULL)

  # Compute annual means
  day_years <- years_all[lyrs]
  ann_means <- tapply(dat, day_years, mean, na.rm = TRUE)

  data.frame(
    YEAR             = as.numeric(names(ann_means)),
    annual_mean_temp = as.numeric(ann_means)
  )
}

# ── 7. Main extraction loop: all sites per population ────────────────────────

results <- foreach(i = lw:up, .combine = "rbind",
                   .packages = c("terra", "sf", "dplyr"),
                   .inorder = TRUE) %dopar% {

  pop <- pop_summary[i, ]
  cat("Pop", i, ":", pop$SPECIES, "@", pop$ECOREGION, "\n")

  r_local <- unwrap(r)

  # All sites for this population
  sites_i <- site_catalogue |>
    filter(ECOREGION == pop$ECOREGION, SPECIES == pop$SPECIES) |>
    arrange(SITE_ORDER)

  if (nrow(sites_i) == 0) return(NULL)

  # Target years: 1993 to last fish sampling year (continuous)
  target_years <- seq(1993, pop$last_year)
  target_years <- target_years[target_years %in% available_years]

  if (length(target_years) < 3) return(NULL)

  # Raster layers for all target years
  lyrs <- unlist(year_layers[as.character(target_years)])

  # Extract at each site
  site_results <- lapply(seq_len(nrow(sites_i)), function(s) {

    site <- sites_i[s, ]

    ann <- extract_annual_temp(r_local, site$LON, site$LAT,
                               lyrs, years_all, year_layers)
    if (is.null(ann)) return(NULL)

    ann$ECOREGION    <- pop$ECOREGION
    ann$SPECIES      <- pop$SPECIES
    ann$SHIFTED_LAT     <- site$LAT
    ann$SHIFTED_LON     <- site$LON
    ann$SITE_ORDER   <- site$SITE_ORDER
    ann$SAMPLED_YEAR <- site$SAMPLED_YEAR

    ann
  })

  do.call(rbind, site_results)
}

#cat("\nExtraction complete. Rows:", nrow(results), "\n")

# ── 7. Add fish-sampled flags ────────────────────────────────────────────────
results <- results |>
  left_join(
    fish_sampled_flags,
    by = c("ECOREGION", "SPECIES", "SHIFTED_LAT" = "LAT",
           "SHIFTED_LON" = "LON", "YEAR")
  ) |>
  mutate(FISH_SAMPLED = ifelse(is.na(FISH_SAMPLED), FALSE, TRUE))

# Reorder columns
results <- results |>
  select(ECOREGION, SPECIES, SHIFTED_LAT, SHIFTED_LON, SITE_ORDER,
         SAMPLED_YEAR, YEAR, annual_mean_temp, FISH_SAMPLED)

# ── 8. Save ──────────────────────────────────────────────────────────────────
dim(results)

assign(paste(env.name, lw, '_', up, sep=""), value=results, pos=1, inherits=T)
outputName=paste(env.name, lw, '_', up, ".RData",sep="")
outputPath=file.path(dir, outputName)
save(list=paste(env.name, lw, '_', up, sep=""), file=outputPath)