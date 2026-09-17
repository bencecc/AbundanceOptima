# Build per-population long-term reference exposure panels at either the
# FIRST optimum site (constant per population) or the PREVIOUS optimum site
# (second-to-last sampling year's optimum, also constant per population).
# For each year 1993 → last sampling year of the population, daily
# temperatures are extracted at the chosen ref site and the same annual
# exposure metrics used by sp_optimloc_ecoreg_shift_lead_time.R are
# recomputed (shifted.cumtemp.above.{mean,q50,q70,q90,q95,q97.5},
# shifted.days.above.{...}, shifted.{mean,max,sd}.temp).
#
# Args:
#   ref.mode  firstsite  | previoussite
#   lw        lower index of populations to process (1-based)
#   up        upper index of populations to process
#
# Output (one file per ref.mode):
#   sp.optim.abund.firstsite.longterm.RData      — list element name: sp.optim.abund.firstsite.longterm
#   sp.optim.abund.previoussite.longterm.RData   — list element name: sp.optim.abund.previoussite.longterm
# Each is a long panel (one row per ECOREGION × SPECIES × YEAR) with the
# same column names as sp.optim.abund.shift, so the downstream HPC and
# plotting scripts can be pointed at it via a `load()` swap.

require(dplyr)

require(tidyr)

require(stringr)
require(sf)
require(terra)
require(foreach, quietly = TRUE)
require(doMC,    quietly = TRUE)
registerDoMC(cores = 15)
source("config.R")
# 1. Load -----------------------------------------------------------------
load(input_file("optim.abund.relocated.RData"))
load(input_file("idx.optim.abund.relocated.RData"))
#load(input_file("optim.density.relocated.RData"))
#load(input_file("idx.optim.density.relocated.RData"))

load(input_file("reef_fish_sti_glorys.RData"))

input.df <- optim.abund.relocated
idx.df       <- idx.optim.abund.relocated
env.name     <- "thetao"
out.env.name <- "temp"

if (!file.exists(file.path(dir_data, "sp_optim_longterm_ref")))
  dir.create(file.path(dir_data, "sp_optim_longterm_ref"))

# 2. Parse args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ref.mode <- gsub("\r", "", as.character(args[1]))   # firstsite | previoussite
lw       <- as.numeric(args[2])
up       <- as.numeric(args[3])
stopifnot(ref.mode %in% c("firstsite", "previoussite"))
cat("ref.mode =", ref.mode, " | lw =", lw, " | up =", up, "\n")

# 3. Identify the reference site per population --------------------------
# firstsite     = first sampled YEAR per (ECOREGION, SPECIES)
# previoussite  = second-to-last sampled YEAR per (ECOREGION, SPECIES)
range.df <- input.df |>
  mutate(abund = H.LAT,
         SAMPLED_AREA = as.character(SAMPLED_AREA)) |>   # force one type
  filter(YEAR != 1992) |>
  left_join(idx.df, by = c("ECOREGION", "SPECIES", "SPECIES.ORIG", "LAT", "LON", "YEAR")) |>
  arrange(ECOREGION, SPECIES, YEAR)

# nth(x, -2) returns NA of the column's own type for groups with n() < 2,
# so no ifelse wrappers are needed.
ref_sites <- range.df |>
  group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
  arrange(YEAR, .by_group = TRUE) |>
  summarise(
    first_year = first(YEAR),
    last_year  = last(YEAR),
    first_LAT  = first(LAT),
    first_LON  = first(LON),
    first_AREA = first(SAMPLED_AREA),
    prev_year  = nth(YEAR,  -2),
    prev_LAT   = nth(LAT,   -2),
    prev_LON   = nth(LON,   -2),
    prev_AREA  = nth(SAMPLED_AREA, -2),
    .groups = "drop") |>
  mutate(
    ref_LAT       = if (ref.mode == "firstsite") first_LAT  else prev_LAT,
    ref_LON       = if (ref.mode == "firstsite") first_LON  else prev_LON,
    ref_AREA      = if (ref.mode == "firstsite") first_AREA else prev_AREA,
    ref_anchor_yr = if (ref.mode == "firstsite") first_year else prev_year) |>
  filter(!is.na(ref_LAT), !is.na(ref_LON), !is.na(last_year))

# Join species thermal index thresholds
env.df <- reef_fish_sti_glorys |>
  rename_with(.fn = function(.x) str_replace(.x, env.name, "env"))

ref_sites <- ref_sites |>
  left_join(env.df, by=c("SPECIES.ORIG"="SPECIES"))

# Build the (pop � evaluation year) grid 1993 at last_year per population
eval_grid <- ref_sites |>
  rowwise() |>
  mutate(years = list(1993:last_year)) |>
  unnest(years) |>
  rename(EVAL_YEAR = years) |>
  ungroup()

# 4. Load the GLORYS daily raster ----------------------------------------
r <- terra::rast(file.path(
  glorys_daily_mean_tif))
time.vec  <- terra::time(r)
years_all <- as.integer(format(time.vec, "%Y"))
r         <- terra::wrap(r)

# 5. Per-year metric calculator ------------------------------------------
var_idx_lt <- function(daily_env, env_thr) {
  daily_env <- daily_env[!is.na(daily_env)]
  if (length(daily_env) == 0) return(rep(NA_real_, 15))
  m <- env_thr
  c(shifted.mean.temp = mean(daily_env),
    shifted.max.temp  = max(daily_env),
    shifted.sd.temp   = sd(daily_env),
    shifted.cumtemp.above.mean   = sum(daily_env[daily_env > m["env_mean"]],   na.rm = TRUE),
    shifted.cumtemp.above.q50    = sum(daily_env[daily_env > m["env_q50"]],    na.rm = TRUE),
    shifted.cumtemp.above.q70    = sum(daily_env[daily_env > m["env_q70"]],    na.rm = TRUE),
    shifted.cumtemp.above.q90    = sum(daily_env[daily_env > m["env_q90"]],    na.rm = TRUE),
    shifted.cumtemp.above.q95    = sum(daily_env[daily_env > m["env_q95"]],    na.rm = TRUE),
    shifted.cumtemp.above.q97.5  = sum(daily_env[daily_env > m["env_q97.5"]],  na.rm = TRUE),
    shifted.days.above.mean      = sum(daily_env > m["env_mean"]),
    shifted.days.above.q50       = sum(daily_env > m["env_q50"]),
    shifted.days.above.q70       = sum(daily_env > m["env_q70"]),
    shifted.days.above.q90       = sum(daily_env > m["env_q90"]),
    shifted.days.above.q95       = sum(daily_env > m["env_q95"]),
    shifted.days.above.q97.5     = sum(daily_env > m["env_q97.5"]))
}

# 6. Loop over populations ----------------------------------------------
pops <- ref_sites |>
  mutate(GROUP_ID = row_number()) |>
  select(GROUP_ID, ECOREGION, SPECIES, SPECIES.ORIG, ref_LAT, ref_LON, ref_AREA,
         last_year,
         env_mean, env_q50, env_q70, env_q90, env_q95, env_q97.5)

batch_pops <- pops |> filter(GROUP_ID >= lw, GROUP_ID <= up)

cat("Populations in this batch:", nrow(batch_pops), "\n")

out_panel <- foreach(i = seq_len(nrow(batch_pops)), .combine = bind_rows,
                     .packages = c("terra", "sf", "dplyr")) %dopar% {

  rr <- terra::unwrap(r)
  pop <- batch_pops[i, ]
  cat("[", i, "/", nrow(batch_pops), "] ",
      pop$ECOREGION, " | ", pop$SPECIES, "\n", sep = "")

  p <- sf::st_point(c(pop$ref_LON, pop$ref_LAT)) |>
       sf::st_sfc(crs = "epsg:4326") |> terra::vect()

  env_thr <- c(env_mean = pop$env_mean,  env_q50 = pop$env_q50,
               env_q70  = pop$env_q70,   env_q90 = pop$env_q90,
               env_q95  = pop$env_q95,   env_q97.5 = pop$env_q97.5)

  years_to_eval <- 1993:pop$last_year

  per_year <- lapply(years_to_eval, function(yy) {
    lyr_idx  <- which(years_all == yy)
    if (length(lyr_idx) == 0) return(NULL)
    r.sub    <- terra::subset(rr, subset = lyr_idx)
    daily    <- as.double(terra::extract(r.sub, p, ID = FALSE))

    # Fallback to neighbour cells if too many NAs at the point
    if (sum(is.na(daily)) >= length(daily) * 0.1) {
      res_xy <- res(rr)[1]
      m1 <- matrix(c(pop$ref_LON - res_xy, pop$ref_LON, pop$ref_LON + res_xy,
                     pop$ref_LAT - res_xy, pop$ref_LAT, pop$ref_LAT + res_xy),
                   ncol = 2, nrow = 3)
      p1 <- sf::st_multipoint(m1) |> sf::st_sfc(crs = "epsg:4326") |> terra::vect()
      dat.tmp <- terra::extract(r.sub, p1, ID = FALSE)
      if (all(is.na(dat.tmp))) return(NULL)
      daily <- apply(dat.tmp, 2, function(x) mean(x, na.rm = TRUE))
    }

    stats <- var_idx_lt(daily, env_thr)

    tibble(
      ECOREGION    = pop$ECOREGION,
      SPECIES      = pop$SPECIES,
      SPECIES.ORIG = pop$SPECIES.ORIG,
      YEAR         = yy,
      SHIFTED.LAT  = pop$ref_LAT,
      SHIFTED.LON  = pop$ref_LON,
      SAMPLED_AREA = pop$ref_AREA,
      !!!as.list(stats))
  })

  bind_rows(per_year)
}

# 7. Save ----------------------------------------------------------------
out_name <- paste0("sp.optim.", ref.mode, ".longterm")
assign(out_name, out_panel)

outfile <- file.path(file.path(dir_data, "sp_optim_longterm_ref"),
                     paste0(out_name, "_", lw, "_", up, ".RData"))
save(list = out_name, file = outfile)
cat("Saved:", outfile, "\n")
