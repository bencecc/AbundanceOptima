# =============================================================================
# bestpath_optimum_search.r
# GLOBAL-OPTIMUM "best reachable lowest-exposure path" per population — an exhaustive,
# provably-optimal dynamic-programming search.
#
# Per population: DYNAMIC PROGRAMMING (DP) over the year x reef-cell combinations.
#   Candidate positions = 0-30 m reef cells (reef_xy); a transition between consecutive years is
#   allowed only when the LEAST-COST-BY-SEA distance <= d_max (tr_sea: LAND impassable, deep water
#   PASSABLE for transit). Year 1 is pinned to the anchor. The DP returns the trajectory that minimises
#   TOTAL exposure over the sampled years.
#   Tie-break: among equal-cost predecessors take the NEAREST (staying wins).
#   Records distance travelled: OBSERVED (realized peak shifts) vs BEST PATH (obs/best_dist_km).
# Deviation on the SAME ruler as Slope.Dev via a JOINT leadtime model (reference + observed +
# best, stacked; ENV.VAR ~ YEAR*CONTR + autocor). One fit gives best_dev (best vs ref), obs_dev (obs
# vs ref), AND the best-vs-obs contrast with the CORRECT covariance (obs & best share the reference
# and data -> NOT independent). One-sided tests: best < ref (p_best_lt_ref), best < obs (p_best_lt_obs).
#
# CONSTANT-EXPOSURE guard: if a best path's exposure is CONSTANT across years (typically best_expo = 0
# every year — a cell that never exceeds the STI), the joint deviation fit is singular (-> NaN), but the
# path's slope is 0 by construction, so best_dev = 0 - ref_slope = -ref_slope analytically (= 0 ONLY if
# the counterfactual is itself flat; NEGATIVE / buffering when it warms). We set this directly and flag it
# (the `refugium` column) instead of returning NaN. NB a best path that merely STAYS at the first site is
# NOT constant (its exposure varies year to year like the reference), so it fits normally and is not
# flagged — the flag captures only the constant-exposure case.
#
# Flexible over DATA.KIND x ref.mode. PROTOTYPE by default (PROTO_N populations); set RUN_FULL.
# =============================================================================
suppressMessages({
  require(dplyr,      quietly = TRUE)
  require(tidyr,      quietly = TRUE)
  require(foreach,    quietly = TRUE)
  require(doMC,       quietly = TRUE)
  require(terra,      quietly = TRUE)
  require(raster,     quietly = TRUE)   # gdistance works on Raster*
  require(gdistance,  quietly = TRUE)   # transition / costDistance (least-cost by sea)
  require(glmmTMB,    quietly = TRUE)
  require(datawizard, quietly = TRUE)   # standardize()
  require(geosphere,  quietly = TRUE)   # distGeo (reachable-buffer prefilter)
  require(RANN,       quietly = TRUE)   # nearest-neighbour (reachable-buffer prefilter)
})
# ----------- SERVER RUN CONFIG — set, then source the whole file -----------
DATA.KIND <- "abund"        # "abund" | "density"
ref.mode  <- "firstsite"    # "firstsite" | "previoussite"
NCORES    <- 40L            # forked workers (DP needs an N x N least-cost matrix per pop -> heavier RAM/CPU;
                            # start modest, watch `free -g`, raise if comfortable)
# PROTOTYPE controls: run a small subset first to inspect behaviour + TIME it (the DP is the slow part;
# check the printed EST full-run before setting RUN_FULL <- TRUE).
RUN_FULL  <- FALSE          # TRUE = all populations; FALSE = a PROTO_N spread (prototype)
PROTO_N   <- 12L            # populations for the prototype spread (ignored if PROTO_POPS set)
PROTO_POPS<- NULL           # optional explicit c("ECOREGION|SPECIES", ...) to target specific pops
SEARCH_RES_DEG <- NA_real_  # reef search grid (deg); NA = use the GLORYS pixel (~1/12 deg)
# ---------------------------------------------------------------------------------------------
stopifnot(DATA.KIND %in% c("abund", "density"),
          ref.mode  %in% c("firstsite", "previoussite"))
registerDoMC(cores = NCORES)
terra::terraOptions(memmax = 2)

source("config.R")
modskurt_dir <- dir_data
daily_tif <- glorys_daily_tiled_tif
if (!file.exists(daily_tif))
  stop("Tiled daily raster not found:\n  ", daily_tif, "\nBuild it once from the NetCDF, then re-run.")
bathy_nc  <- gebco_file
DENS_W  <- 0.5
BAND    <- c(0, 30)     # coastal habitat band (m depth, positive down)
AGG_DEG <- 0.05         # ~5 km movement grid (least-cost sea surface)
min.yrs <- 4

shift_obj_name <- paste0("sp.optim.", DATA.KIND, ".shift")

# ---- joint leadtime-model deviation — reproduces optimloc_ecoregion_leadtime_trend_longterm.R.
# Stacks the reference trajectory (CONTR "ref.var", full series) and a candidate path (CONTR
# "obs.var", sampled years), fits ENV.VAR ~ YEAR*CONTR + autocor, and returns the YEAR:CONTR
# interaction = the path's deviation from the reference — the SAME quantity as the stored Slope.Dev.
# STANDARDIZED units when spec$trans=="Standardize" (the caller back-transforms with
# sd_env_stacked / sd_year, exactly as Fig 4 does for Slope.Dev).
fit_leadtime_dev <- function(ref_env, ref_year, path_env, path_year, spec) {
  d <- bind_rows(
        data.frame(ENV.VAR = ref_env,  YEAR = ref_year,  CONTR = "ref.var"),
        data.frame(ENV.VAR = path_env, YEAR = path_year, CONTR = "obs.var")) |>
    tidyr::drop_na(ENV.VAR) |> arrange(CONTR, YEAR)
  if (sum(d$CONTR == "obs.var") < min.yrs || sum(d$CONTR == "ref.var") < min.yrs)
    return(c(dev = NA_real_, se = NA_real_))
  d <- d |> mutate(YEAR_c    = (YEAR - min(YEAR)) + 1,
                   times.dec = as.numeric(paste(YEAR_c, 1, sep = ".")),
                   cum.times = lead(cumsum(lag(times.dec, default = 0)), default = sum(times.dec)),
                   yrs   = glmmTMB::numFactor(cum.times),
                   group = 1,
                   CONTR = factor(CONTR, levels = c("ref.var", "obs.var")))
  if (spec$trans == "Standardize")
    d <- d |> mutate(YEAR = as.vector(standardize(YEAR_c)), ENV.VAR = as.vector(standardize(ENV.VAR)))
  else d$YEAR <- d$YEAR_c
  fm <- switch(spec$autocor,
    no_autocor = ENV.VAR ~ YEAR * CONTR + (1 | group),
    ou         = ENV.VAR ~ YEAR * CONTR + ou(yrs + 0 | group),
    ar1        = ENV.VAR ~ YEAR * CONTR + ar1(as.factor(YEAR) + 0 | group),
    ENV.VAR ~ YEAR * CONTR + (1 | group))
  m <- tryCatch(suppressWarnings(glmmTMB(fm, data = d, ziformula = ~ 0,
        dispformula = ~ 1, family = spec$family)), error = function(e) NULL)
  if (is.null(m))
    m <- tryCatch(suppressWarnings(glmmTMB(fm, data = d, ziformula = ~ 0,
          dispformula = ~ 1, family = spec$family,
          control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))),
          error = function(e) NULL)
  if (is.null(m)) return(c(dev = NA_real_, se = NA_real_))
  co <- tryCatch(summary(m)$coefficients$cond["YEAR:CONTRobs.var", c("Estimate", "Std. Error")],
                 error = function(e) c(NA_real_, NA_real_))
  c(dev = unname(co[1]), se = unname(co[2]))
}

# ---- Joint leadtime model: reference (baseline) + observed + best, stacked. ONE fit gives
# best_dev (best vs ref) and obs_dev (obs vs ref) as the two YEAR:CONTR interactions, plus the
# best-vs-obs contrast (diff = best - obs) with its SE from the model COVARIANCE — correct, since obs
# and best share the reference and data (no independence assumption). Standardized units when
# spec$trans=="Standardize" (caller back-transforms). Returns NA list on failure.
fit_leadtime_dev3 <- function(ref_env, ref_year, obs_env, obs_year, best_env, best_year, spec) {
  na3 <- list(best_dev = NA_real_, se_best = NA_real_, obs_dev = NA_real_, se_obs = NA_real_,
              diff = NA_real_, se_diff = NA_real_)
  d <- bind_rows(
        data.frame(ENV.VAR = ref_env,  YEAR = ref_year,  CONTR = "ref"),
        data.frame(ENV.VAR = obs_env,  YEAR = obs_year,  CONTR = "obs"),
        data.frame(ENV.VAR = best_env, YEAR = best_year, CONTR = "best")) |>
    tidyr::drop_na(ENV.VAR) |> arrange(CONTR, YEAR)
  if (any(tapply(d$ENV.VAR, d$CONTR, length) < min.yrs)) return(na3)
  d <- d |> mutate(YEAR_c    = (YEAR - min(YEAR)) + 1,
                   times.dec = as.numeric(paste(YEAR_c, 1, sep = ".")),
                   cum.times = lead(cumsum(lag(times.dec, default = 0)), default = sum(times.dec)),
                   yrs   = glmmTMB::numFactor(cum.times), group = 1,
                   CONTR = factor(CONTR, levels = c("ref", "obs", "best")))   # ref = baseline
  if (spec$trans == "Standardize")
    d <- d |> mutate(YEAR = as.vector(standardize(YEAR_c)), ENV.VAR = as.vector(standardize(ENV.VAR)))
  else d$YEAR <- d$YEAR_c
  fm <- switch(spec$autocor,
    no_autocor = ENV.VAR ~ YEAR * CONTR + (1 | group),
    ou         = ENV.VAR ~ YEAR * CONTR + ou(yrs + 0 | group),
    ar1        = ENV.VAR ~ YEAR * CONTR + ar1(as.factor(YEAR) + 0 | group),
    ENV.VAR ~ YEAR * CONTR + (1 | group))
  m <- tryCatch(suppressWarnings(glmmTMB(fm, data = d, ziformula = ~ 0, dispformula = ~ 1,
        family = spec$family)), error = function(e) NULL)
  if (is.null(m))
    m <- tryCatch(suppressWarnings(glmmTMB(fm, data = d, ziformula = ~ 0, dispformula = ~ 1,
          family = spec$family, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))),
          error = function(e) NULL)
  if (is.null(m)) return(na3)
  co <- tryCatch(summary(m)$coefficients$cond, error = function(e) NULL)
  Vc <- tryCatch(vcov(m)$cond, error = function(e) NULL)
  ob <- "YEAR:CONTRobs"; be <- "YEAR:CONTRbest"
  if (is.null(co) || is.null(Vc) || !all(c(ob, be) %in% rownames(co)) || !all(c(ob, be) %in% rownames(Vc)))
    return(na3)
  b_obs <- co[ob, "Estimate"];   b_best <- co[be, "Estimate"]
  se_ob <- co[ob, "Std. Error"]; se_be  <- co[be, "Std. Error"]
  vdiff <- Vc[be, be] + Vc[ob, ob] - 2 * Vc[be, ob]              # Var(best - obs) with the correlation
  list(best_dev = b_best, se_best = se_be, obs_dev = b_obs, se_obs = se_ob,
       diff = b_best - b_obs, se_diff = sqrt(max(vdiff, 0)))
}
parse_spec <- function(mp) list(
  autocor = if (grepl("_ou_", mp)) "ou" else if (grepl("_ar1_", mp)) "ar1" else "no_autocor",
  trans   = if (grepl("_Standardize_", mp)) "Standardize" else "NULL",
  family  = if (grepl("_lognormal", mp)) "lognormal" else "gaussian")

# ---- data (same objects as the null script) --------------------------------
load(file.path(modskurt_dir, paste0(shift_obj_name, ".RData")))
load(file.path(modskurt_dir, "reef_fish_sti_glorys.RData"))
sp.focal <- get(shift_obj_name) |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON, ECOREGION = ORIG.ECOREGION)
sti_thr <- reef_fish_sti_glorys |>
  dplyr::select(SPECIES, thr = thetao_mean_TM) |> distinct()

# ---- observed reference + deviation, from the leadtime longterm trend file ------------------------
# optim.<kind>.leadtime.<ref>.trend.longterm already stores, per candidate model and population,
# BOTH the reference-trajectory slope (Effect == "Ref.Slope") and the observed deviation from it
# (Effect == "Slope.Dev", identical to Fig 4's trends_df$Slope.Dev). Per population we take the
# AIC-best cumtemp.above.mean model and read Ref.Slope + Slope.Dev (+SE) from that ONE model, so
# reference and deviation are mutually consistent (obs slope = Ref.Slope + Slope.Dev). No scenario
# data, no re-fit. Standardize models are back-transformed to per-year units with the population's
# YEAR sd (as in the reference fit and the null).
lt_obj <- sprintf("optim.%s.leadtime.%s.trend.longterm", DATA.KIND, ref.mode)
load(file.path(modskurt_dir, paste0(lt_obj, ".RData")))
lt_df <- get(lt_obj) |>
  filter(FlagAnalysis == "Trend_Environmental",
         grepl("shifted.cumtemp.above.mean_", mod.parms, fixed = TRUE),
         Effect %in% c("Ref.Slope", "Slope.Dev"), is.finite(AIC))
lt_best <- lt_df |> group_by(ECOREGION, SPECIES) |>
  slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup() |>
  dplyr::select(ECOREGION, SPECIES, sel.mod.parms = mod.parms)
obs_ref_tab <- lt_df |> inner_join(lt_best, by = c("ECOREGION", "SPECIES")) |>
  filter(mod.parms == sel.mod.parms) |>
  dplyr::select(ECOREGION, SPECIES, mod.parms, Effect, Estimate, SE) |>
  tidyr::pivot_wider(names_from = Effect, values_from = c(Estimate, SE)) |>
  rename(ref.slope = `Estimate_Ref.Slope`, obs.dev = `Estimate_Slope.Dev`,
         ref.se = `SE_Ref.Slope`, obs.se = `SE_Slope.Dev`)

ref_panel_obj <- paste0("sp.optim.", DATA.KIND, ".", ref.mode, ".longterm")
load(file.path(modskurt_dir, paste0(ref_panel_obj, ".RData")))
ref_panel <- get(ref_panel_obj)
if ("ORIG.ECOREGION" %in% names(ref_panel)) ref_panel <- ref_panel |> rename(ECOREGION = ORIG.ECOREGION)
# Back-transform for Standardize models — EXACTLY as Fig 4. The leadtime fit standardizes BOTH
# ENV.VAR and YEAR, so raw slope = Estimate * sd_env_stacked / sd_year. sd_year from the reference
# panel; sd_env_stacked = sd of the pooled reference + realized-observed exposure (long-term stack).
resp.env <- "shifted.cumtemp.above.mean"
sd_year_lu <- ref_panel |> group_by(ECOREGION, SPECIES) |>
  summarise(sd_year = sd((YEAR - min(YEAR)) + 1), .groups = "drop")
sd_env_lu <- bind_rows(
    ref_panel |> dplyr::select(ECOREGION, SPECIES, v = all_of(resp.env)),
    sp.focal  |> dplyr::select(ECOREGION, SPECIES, v = all_of(resp.env))) |>
  filter(!is.na(v)) |> group_by(ECOREGION, SPECIES) |>
  summarise(sd_env_stacked = sd(v), .groups = "drop")
obs_ref_tab <- obs_ref_tab |>
  left_join(sd_year_lu, by = c("ECOREGION", "SPECIES")) |>
  left_join(sd_env_lu,  by = c("ECOREGION", "SPECIES")) |>
  mutate(is_std = grepl("_Standardize_", mod.parms),
         bt     = if_else(is_std & is.finite(sd_year) & sd_year > 0, sd_env_stacked / sd_year, 1),
         ref.slope = ref.slope * bt, obs.dev = obs.dev * bt,
         ref.se    = ref.se    * bt, obs.se  = obs.se  * bt) |>
  dplyr::select(ECOREGION, SPECIES, ref.slope, obs.dev, ref.se, obs.se, sd_year, ref.mod.parms = mod.parms)

pops <- sp.focal |> distinct(ECOREGION, SPECIES) |> arrange(ECOREGION, SPECIES) |>
  mutate(GROUP_ID = row_number())
# --- prototype subset (deterministic spread, or an explicit list) -----------
if (!RUN_FULL) {
  if (!is.null(PROTO_POPS)) {
    key <- paste(pops$ECOREGION, pops$SPECIES, sep = "|")
    batch_pops <- pops[key %in% PROTO_POPS, ]
  } else {
    batch_pops <- pops[unique(round(seq(1, nrow(pops), length.out = PROTO_N))), ]
  }
} else batch_pops <- pops
cat("bestpath | DATA.KIND =", DATA.KIND, "| ref.mode =", ref.mode,
    "| pops =", nrow(batch_pops), if (RUN_FULL) "(FULL)" else "(PROTO)", "| cores =", NCORES, "\n")

daily_full <- rast(daily_tif)
yr_layer   <- as.numeric(format(time(daily_full), "%Y"))
if (!all(is.finite(as.numeric(time(daily_full)))))      # time-axis safety (as in the null)
  yr_layer <- as.numeric(format(time(rast(glorys_daily_nc)), "%Y"))
glorys_lat_res <- terra::res(daily_full)[2]
stopifnot(any(is.finite(yr_layer)), is.finite(glorys_lat_res), glorys_lat_res > 0)
yr_layer_set <- unique(yr_layer[is.finite(yr_layer)])
search_res <- if (is.finite(SEARCH_RES_DEG)) SEARCH_RES_DEG else glorys_lat_res

reason_row <- function(eco, sp, why)
  data.frame(ECOREGION = eco, SPECIES = sp, resolvable = NA, d_max_km = NA_real_, n_years = NA_integer_,
             ref = NA_real_, best_se = NA_real_,
             best_dev = NA_real_, obs_dev = NA_real_, obs_dev_joint = NA_real_,
             obs_expo = NA_real_, best_expo = NA_real_,
             obs_dist_km = NA_real_, best_dist_km = NA_real_,
             obs_dist_mean_km = NA_real_, best_dist_mean_km = NA_real_, best_net_disp_km = NA_real_,
             p_best_lt_ref = NA_real_, p_best_lt_obs = NA_real_,
             best_poleward = NA, refugium = NA, secs = NA_real_, drop_reason = why, stringsAsFactors = FALSE)

# ── population loop ────────────────────────────────────────────────────────
out_panel <- foreach(i = seq_len(nrow(batch_pops)), .combine = bind_rows,
                     .errorhandling = "remove") %dopar% {
  t0 <- proc.time()[3]
  this_eco <- batch_pops$ECOREGION[i]; this_sp <- batch_pops$SPECIES[i]
  tryCatch({
    pop_df <- sp.focal |> filter(ECOREGION == this_eco, SPECIES == this_sp) |> arrange(YEAR)
    n_all  <- nrow(pop_df); pop_df <- pop_df |> filter(YEAR %in% yr_layer_set)
    if (nrow(pop_df) < min.yrs)
      return(reason_row(this_eco, this_sp, sprintf("few_years: %d/%d", n_all, nrow(pop_df))))
    this_sp_orig <- pop_df$SPECIES.ORIG[1]
    thr <- sti_thr |> filter(SPECIES == this_sp_orig) |> pull(thr)
    if (length(thr) != 1 || is.na(thr)) return(reason_row(this_eco, this_sp, "no_STI"))
    or_row <- obs_ref_tab |> filter(ECOREGION == this_eco, SPECIES == this_sp)
    if (nrow(or_row) != 1) return(reason_row(this_eco, this_sp, "no_ref"))
    ref_slope  <- or_row$ref.slope; obs_dev <- or_row$obs.dev; obs_se_use <- or_row$obs.se
    spec <- parse_spec(or_row$ref.mod.parms)
    obs_expo <- mean(pop_df$shifted.cumtemp.above.mean, na.rm = TRUE)   # observed mean annual exposure

    yrs <- pop_df$YEAR; n_years <- length(yrs)
    anchor_lon <- pop_df$LON[1]; anchor_lat <- pop_df$LAT[1]
    pad <- 3
    ext_reg <- ext(min(pop_df$LON) - pad, max(pop_df$LON) + pad,
                   min(pop_df$LAT) - pad, max(pop_df$LAT) + pad)

    # --- least-cost-by-sea surface + d_max (largest realized annual step, by sea)
    bathy_w <- terra::rast(bathy_nc); terra::crs(bathy_w) <- "epsg:4326"
    bz    <- crop(bathy_w, ext_reg); depth <- -bz
    fact  <- max(1, round(AGG_DEG / mean(res(bz))))
    sea   <- aggregate(classify(depth, rbind(c(-Inf, 0, NA), c(0, Inf, 1))), fact = fact,
                       fun = function(x) if (any(!is.na(x))) 1 else NA_real_)
    sea_r <- raster::raster(sea)
    tr_sea <- tryCatch(gdistance::geoCorrection(
      gdistance::transition(sea_r, transitionFunction = function(x) mean(x), directions = 8), type = "c"),
      error = function(e) conditionMessage(e))
    if (is.null(tr_sea) || is.character(tr_sea)) return(reason_row(this_eco, this_sp, "no_transition"))
    sea_xy <- raster::rasterToPoints(sea_r)[, 1:2, drop = FALSE]
    if (nrow(sea_xy) < 5) return(reason_row(this_eco, this_sp, "few_sea_cells"))
    snap_sea <- function(lon, lat) sea_xy[which.min((sea_xy[,1]-lon)^2 + (sea_xy[,2]-lat)^2), ]
    real_sea <- t(apply(cbind(pop_df$LON, pop_df$LAT), 1, function(p) snap_sea(p[1], p[2])))
    step_d <- tryCatch(diag(gdistance::costDistance(tr_sea,
                real_sea[-nrow(real_sea), , drop = FALSE], real_sea[-1, , drop = FALSE])),
                error = function(e) NA_real_)
    d_max <- suppressWarnings(max(step_d[is.finite(step_d)], na.rm = TRUE))

    resolvable <- is.finite(d_max) && d_max > 0
    # daily crop for the sampled years (reused for every candidate cell)
    daily_w   <- terra::rast(daily_tif)
    need_lays <- which(yr_layer %in% yrs)
    daily_tf  <- tempfile(fileext = ".tif")
    daily_reg <- terra::crop(daily_w[[need_lays]], ext_reg, filename = daily_tf, overwrite = TRUE)
    yr_reg    <- yr_layer[need_lays]
    drop_id <- function(x) { x <- as.data.frame(x); x[["ID"]] <- NULL; as.matrix(x) }
    extract_daily <- function(rstack, pts) {
      pts <- matrix(pts, ncol = 2); v <- drop_id(terra::extract(rstack, pts))
      bad <- which(rowSums(!is.na(v)) == 0)
      if (length(bad) > 0) { dd <- mean(terra::res(rstack))
        for (b in bad) { nb <- rbind(c(pts[b,1]-dd, pts[b,2]), c(pts[b,1]+dd, pts[b,2]),
                                     c(pts[b,1], pts[b,2]-dd), c(pts[b,1], pts[b,2]+dd))
          v[b, ] <- colMeans(drop_id(terra::extract(rstack, nb)), na.rm = TRUE) } }
      v
    }
    expo_of <- function(v) if (all(is.na(v))) NA_real_ else sum(v[v > thr], na.rm = TRUE)
    expo_by_year <- function(vmat)                    # vmat: n_pts x n_layers -> n_pts x n_years
      vapply(yrs, function(y) apply(vmat[, yr_reg == y, drop = FALSE], 1, expo_of), numeric(nrow(vmat)))

    # --- reef candidate cells (0-30 m band), at the search resolution, within reach
    reef_bin <- classify(depth, rbind(c(-Inf, BAND[1], NA), c(BAND[1], BAND[2], 1), c(BAND[2], Inf, NA)))
    rf <- max(1, round(search_res / mean(res(depth))))
    reef_ag  <- aggregate(reef_bin, fact = rf, fun = function(x) if (any(!is.na(x))) 1 else NA_real_)
    reef_xy  <- as.matrix(as.data.frame(reef_ag, xy = TRUE)[, 1:2])
    if (!nrow(reef_xy)) return(reason_row(this_eco, this_sp, "no_reef_cells"))
    reach_buf <- max(d_max, 1) * (n_years + 1) * 1.3   # cumulative straight-line reach (m), generous
    keep <- geosphere::distGeo(c(anchor_lon, anchor_lat), reef_xy) <= reach_buf
    reef_xy <- reef_xy[keep, , drop = FALSE]
    if (!nrow(reef_xy)) return(reason_row(this_eco, this_sp, "no_reef_in_reach"))

    # exposure of every candidate reef cell, per year (extract ONCE, reuse in the walk)
    reef_expo <- expo_by_year(extract_daily(daily_reg, reef_xy))     # n_reef x n_years
    start_idx <- which.min(geosphere::distGeo(c(anchor_lon, anchor_lat), reef_xy))

    if (!resolvable) {                                # no reachable alternative -> best = stay
      best_idx_path <- rep(start_idx, n_years)
      best_dist_total <- 0; best_dist_mean <- 0; best_net_disp <- 0
    } else {
      # ---- GLOBAL OPTIMUM by DYNAMIC PROGRAMMING over the (year x reef-cell) trellis -----------------
      # Seascape: positions = 0-30 m reef cells (reef_xy); a transition is allowed
      # only when the LEAST-COST-BY-SEA distance <= d_max (tr_sea: land impassable, deep water passable
      # for transit). Minimise TOTAL exposure over the sampled years; year 1 pinned to the anchor.
      # Tie-break: among equal-cost predecessors take the NEAREST (staying wins, self-distance = 0).
      Dm <- as.matrix(gdistance::costDistance(tr_sea, reef_xy))   # N x N least-cost by sea (m)
      reachable <- is.finite(Dm) & Dm <= d_max
      N  <- nrow(reef_xy)
      LAMBDA <- 1e-9                                            # tiny distance penalty: breaks exposure ties
      # toward the SHORTEST movement (staying, self-dist 0, wins) WITHOUT changing the exposure optimum.
      Vc <- rep(Inf, N); Vc[start_idx] <- { e1 <- reef_expo[start_idx, 1]; if (is.finite(e1)) e1 else 0 }
      back <- matrix(NA_integer_, nrow = N, ncol = n_years)
      for (t in 2:n_years) {                                    # vectorised forward pass (C-level, not R loops)
        et <- reef_expo[, t]
        C  <- matrix(Vc, N, N) + LAMBDA * Dm                    # C[i,j] = cost-to-reach-i + tiny*dist(i->j)
        C[!reachable] <- Inf
        pj <- suppressWarnings(max.col(-t(C), ties.method = "first"))   # arg-min predecessor per destination j
        cj <- C[cbind(pj, seq_len(N))]                          # its (min) value
        Vn <- cj + et
        ok_j <- is.finite(et) & is.finite(cj)
        Vn[!ok_j] <- Inf
        back[ok_j, t] <- pj[ok_j]
        Vc <- Vn
      }
      endj <- if (all(!is.finite(Vc))) start_idx else which.min(Vc)      # min-total-exposure endpoint
      best_idx_path <- integer(n_years); best_idx_path[n_years] <- endj
      for (t in n_years:2) best_idx_path[t - 1] <- back[best_idx_path[t], t]
      best_idx_path[1] <- start_idx                                      # year 1 fixed at the anchor
      if (anyNA(best_idx_path)) best_idx_path <- rep(start_idx, n_years) # safety on a broken backtrack
      st <- if (n_years >= 2) Dm[cbind(best_idx_path[-n_years], best_idx_path[-1])] else numeric(0)
      st <- st[is.finite(st)]
      best_dist_total <- if (length(st)) sum(st)  else 0
      best_dist_mean  <- if (length(st)) mean(st) else 0
      best_net_disp   <- if (n_years >= 2) { d0 <- Dm[best_idx_path[1], best_idx_path[n_years]]
                                             if (is.finite(d0)) d0 else 0 } else 0
    }
    best_expo_vec <- reef_expo[cbind(best_idx_path, seq_len(n_years))]
    best_expo     <- mean(best_expo_vec, na.rm = TRUE)         # best-path mean annual exposure
    best_lat      <- reef_xy[best_idx_path, 2]
    best_polew    <- tryCatch(unname(coef(lm(abs(best_lat) ~ yrs))[2]) > 0, error = function(e) NA)
    # OBSERVED movement distance (realized peak shifts) for the best-vs-real check.
    obs_steps      <- step_d[is.finite(step_d)]
    obs_dist_total <- if (length(obs_steps)) sum(obs_steps)  else NA_real_
    obs_dist_mean  <- if (length(obs_steps)) mean(obs_steps) else NA_real_

    rm(tr_sea, sea_r, daily_reg); unlink(daily_tf); gc(FALSE)
    # --- deviation + tests via a JOINT leadtime model (reference + observed + best): one fit
    # gives best_dev (best vs ref), the joint obs_dev (cross-check vs stored Slope.Dev), and the
    # best-vs-obs contrast with the CORRECT covariance. Back-transform Standardize models with
    # sd(pooled env)/sd_year, as Fig 4 does for Slope.Dev.
    ref_series <- ref_panel |> filter(ECOREGION == this_eco, SPECIES == this_sp) |> arrange(YEAR)
    ref_env <- ref_series$shifted.cumtemp.above.mean; ref_yr <- ref_series$YEAR
    obs_env <- pop_df$shifted.cumtemp.above.mean
    is_std  <- spec$trans == "Standardize"; sd_year_pop <- or_row$sd_year
    bt <- function(x, envs) { if (!is_std) return(x); sde <- sd(unlist(envs), na.rm = TRUE)
      if (is.finite(sd_year_pop) && sd_year_pop > 0) x * sde / sd_year_pop else NA_real_ }
    # CONSTANT-EXPOSURE guard (see header): a CONSTANT best-path exposure (best_expo = 0 every year — a
    # cell never above the STI) makes the joint fit singular. The flat path has slope 0, so analytically
    # best_dev = -ref_slope (0 only when the counterfactual is flat; NEGATIVE/buffering when it warms);
    # best is then EXACT, so best-vs-obs uses obs's own SE. Flag `refugium`, never NaN, never a blanket 0.
    bv_fin   <- best_expo_vec[is.finite(best_expo_vec)]
    refugium <- length(bv_fin) < min.yrs || length(unique(bv_fin)) < 2
    if (refugium) {
      best_dev <- -ref_slope; best_se <- or_row$ref.se; obs_dev_joint <- NA_real_
      z1 <- best_dev / best_se
      se_d <- if (is.finite(or_row$obs.se) && or_row$obs.se > 0) or_row$obs.se else best_se
      z2 <- (best_dev - obs_dev) / se_d
    } else {
      j3   <- fit_leadtime_dev3(ref_env, ref_yr, obs_env, yrs, best_expo_vec, yrs, spec)
      envs <- list(ref_env, obs_env, best_expo_vec)
      best_dev      <- bt(j3$best_dev, envs); best_se <- bt(j3$se_best, envs)
      obs_dev_joint <- bt(j3$obs_dev,  envs)                 # cross-check: should ~ stored Slope.Dev
      z1 <- if (is.finite(best_se) && best_se > 0)   best_dev / best_se else NA_real_
      z2 <- if (is.finite(j3$se_diff) && j3$se_diff > 0) j3$diff / j3$se_diff else NA_real_
    }
    data.frame(ECOREGION = this_eco, SPECIES = this_sp, resolvable = resolvable,
               d_max_km = d_max / 1000, n_years = n_years, ref = ref_slope,
               best_se = best_se, best_dev = best_dev, obs_dev = obs_dev,
               obs_dev_joint = obs_dev_joint,                         # 3-arm cross-check vs stored Slope.Dev
               obs_expo = obs_expo, best_expo = best_expo,
               obs_dist_km = obs_dist_total / 1000, best_dist_km = best_dist_total / 1000,
               obs_dist_mean_km = obs_dist_mean / 1000, best_dist_mean_km = best_dist_mean / 1000,
               best_net_disp_km = best_net_disp / 1000,
               p_best_lt_ref = unname(pnorm(z1)), p_best_lt_obs = unname(pnorm(z2)),
               best_poleward = best_polew, refugium = refugium, secs = proc.time()[3] - t0,
               drop_reason = NA_character_, stringsAsFactors = FALSE)
  }, error = function(e) reason_row(this_eco, this_sp, paste("ERROR:", conditionMessage(e))))
}

ok <- out_panel |> filter(is.na(drop_reason))
cat(sprintf("\nDONE: %d/%d populations produced a result (%.0f%%). Median time/pop: %.1f s.\n",
            nrow(ok), nrow(out_panel), 100 * nrow(ok) / max(1, nrow(out_panel)), median(ok$secs, na.rm = TRUE)))
if (nrow(out_panel) > nrow(ok)) { cat("drop reasons:\n"); print(table(sub(":.*", "", out_panel$drop_reason[!is.na(out_panel$drop_reason)]))) }
cat(sprintf("EST. full run (%d pops, %d cores): %.0f min\n",
            nrow(pops), NCORES, nrow(pops) * median(ok$secs, na.rm = TRUE) / NCORES / 60))
print(ok |> dplyr::select(ECOREGION, SPECIES, n_years, d_max_km, ref,
                          best_dev, obs_dev, obs_expo, best_expo, obs_dist_km, best_dist_km,
                          p_best_lt_ref, p_best_lt_obs, best_poleward, refugium, secs) |>
      mutate(across(where(is.numeric), \(x) round(x, 4))) |> as.data.frame())
cat(sprintf("[distance] best-path vs observed movement: mean best_dist %.0f km vs obs_dist %.0f km (median over pops)\n",
            median(ok$best_dist_km, na.rm = TRUE), median(ok$obs_dist_km, na.rm = TRUE)))

out_name <- paste0("bestpath.opt.", DATA.KIND, ".", ref.mode, if (!RUN_FULL) ".proto" else "")
assign(out_name, out_panel)
save(list = out_name, file = file.path(dir_data, paste0(out_name, ".RData")))
cat("Wrote", file.path(dir_data, paste0(out_name, ".RData")), "\n")
