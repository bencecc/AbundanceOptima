# ==========================================================================
# cluster_ecoregions.R  —  STEP 1 of the population-splitting workflow
#
# For each ecoregion, split its sampling sites into sub-regions ONLY where a
# single modskurt fit would otherwise violate the analysis's core requirement of
# a CONTINUOUS LATITUDINAL GRADIENT. Two barrier types are detected, recursively.
# The split criterion is purely a BARRIER (can a population span the gap?), NOT
# temperature — the thermal gradient is the signal the analysis tracks.
# Temperature is reported (meanT/dT) for inspection only.
#
#   (a) LATITUDINAL GAP  — a wide unsampled latitude band (>= LAT_GAP_MIN_DEG)
#       separating two DIFFERENT land masses (e.g. Bass Strait: mainland-Victoria
#       vs Tasmania): -> split N / S. Gate = gap + groups_diff_landmass (barrier),
#       NO thermal requirement. (A gap on the SAME continuous coast is just a
#       sampling gap on a continuous gradient -> not split).
#   (b) LONGITUDINAL DOUBLING — two longitude clusters at OVERLAPPING latitudes
#       with LAND between them (e.g. E/W Tasmania).
#       Pooling them lands the peak ON LAND and relocation would swap coasts
#       (a spurious LONGITUDINAL shift):  -> split E / W. Land barrier alone.
#
# Justification: Martins et al. 2024 (PNAS; same group as modskurt) split
# North America E/W to preserve continuous latitudinal gradients, justified by a
# topographic barrier + climatic distinctness. This is the marine analogue.
#
# NOT used: by-sea connectivity (lets fish swim around -> misses E/W Tasmania);
# community clustering (can split on non-climatic grounds); Hawaii exclusion
# (Hawaii passed the unimodal QC -> its marginals are unimodal; its off-diagonal
# raw peak is a separability artifact handled by relocation, latitude valid).
#
# Output = a SPLIT PLAN.
# ==========================================================================

suppressMessages({
  require(dplyr); require(tidyr); require(terra)
  require(foreach); require(doMC)
})
source("config.R")
registerDoMC(cores = 64)   # adjust to your core count (Windows: silently sequential)
terra::terraOptions(memmax = 8)   # cap terra at 8 GB/op so big-RAM nodes STREAM rasters
                                  #   to disk instead of loading the whole GEBCO into memory

# ── Tunables (the TEST BATTERY tunes these) ───────────────────────────────
# DT_MIN          (REMOVED as a split criterion — temperature no longer gates any
#                 split; meanT/dT are reported for inspection only. See header.)
MIN_SITES       <- 5      # a kept group needs >= this many unique sites
MIN_YEARS       <- 4      # per-species/group fitting guard (= nyrs in modskurt_analysis)
LAT_GAP_MIN_DEG <- 1.0    # unsampled lat gap >= this (~110 km) to consider a lat split.
                          #   >0.75 because reef-fish LARVAE bridge ~80 km gaps (adults
                          #   don't), so a smaller gap is not a true population break; the
                          #   land-mass barrier test is the real safeguard anyway.
LAT_GAP_FAC     <- 4      #   ...and >= this multiple of the median lat spacing
LON_GAP_MIN_DEG <- 0.30   # a longitudinal doubling needs a lon gap >= this between clusters
LAT_OVERLAP_MIN <- 0.30   #   ...the two lon clusters must share >= this frac of their lat span
LAND_FRAC_MIN   <- 0.30   #   ...and >= this frac of the lon gap (over the overlap) must be LAND
LARGE_LAND_KM2  <- 3000   # a "large" land mass for the lat-split disconnection test
NEARBY_LAND_KM  <- 200    # a group near a SMALL island <= this far from a large land mass is
                          #   assigned to that mainland gradient (keeps GBR cays / Cape Howe
                          #   Bass-Strait islands whole); FAR isolated islands (Hawaii NW chain)
                          #   or large islands (Tasmania) keep their own identity -> split
CONN_RES_DEG    <- 1/30    # ~3.7 km GEBCO grid for the land-between check
AOI_PAD_DEG     <- 0.5    # small crop pad (only local land/sea needed now)

# ── Paths (local work machine / server) ────────────────────────────────────
modskurt_dir <- dir_data
gebco_nc     <- gebco_file
glorys_path  <- glorys_daily_mean_tif   # temperature is report-only here (not a split criterion)

# NULL = all 50 -> split_plan.RData ; a subset -> split_plan_TEST.RData
# test battery: c("Bassian","Hawaii","Ningaloo","Western Sumatra","Tweed-Moreton")
TEST_ECOREGIONS <- NULL

# ── CONFIRMED biogeographic barriers (final gate) ──────────────────────────
# The detectors above PROPOSE candidate splits; the saved plan keeps only the
# candidates that correspond to a recognised barrier. This is a deliberate gate,
# not a threshold, because no single geometric/bathymetric cut separates the real
# splits from the spurious ones — the gap-corridor depths rank them BACKWARDS:
#   Bassian (Bass Strait, ~ -75 m shelf)  -> SPLIT  (cool-temperate Tasmania vs
#       warm-temperate mainland; a recognised faunal-province break, soft-sediment
#       strait the rocky-reef fauna does not cross)
#   Hawaii  (NW chain, ~ -4000 m abyss)   -> SPLIT  (main islands isolated from the
#       NW Hawaiian chain by deep open ocean)
#   Red Sea (~ -850 m axial trough but groups at the SAME longitude, 0.21 deg apart)
#       -> NOT split: a continuous N-S coastal reef gradient with only a sampling
#       gap; the deep axis runs BESIDE the gradient, not across it; the mid-sea
#       offset is relocation-handled (latitude/exposure preserved).
#   Sunda  (Java vs Sumatra, ~1100 km but ~ -15..-50 m continuous Sunda Shelf)
#       -> NOT split: one shallow-shelf fauna, larval-connected across the shelf;
#       the pooled peak lands on shelf habitat, not unsuitable water.
# Set to NULL to save every detected split (no confirmation gate).
CONFIRMED_BARRIERS <- c("Bassian", "Hawaii")

# ── Inputs ─────────────────────────────────────────────────────────────────
load(file.path(modskurt_dir, "sp.df.RData"))

# CRITICAL (server memory): pre-aggregate the FULL GEBCO to the coarse connectivity
# grid ONCE, cached to disk — exactly like relocate_seapath's GEBCO_work_x2.tif.
# terra::aggregate(filename=...) streams in blocks, so this is memory-bounded.
# Parallel workers then crop this SMALL cached raster, NEVER the large native nc
# (importing the native into many workers is what blows up memory).
gebco_work <- file.path(modskurt_dir, "GEBCO_cluster_work.tif")
if (!file.exists(gebco_work)) {
  message("Pre-aggregating GEBCO to the ~", round(CONN_RES_DEG * 111, 1),
          " km grid (one-time, cached)...")
  g    <- terra::rast(gebco_nc); terra::crs(g) <- "epsg:4326"
  fact <- max(1, round(CONN_RES_DEG / mean(terra::res(g))))
  # fun="max" is LAND-PRESERVING (GEBCO land = positive): a coarse cell is land if
  # ANY native cell is land, so small islands/atolls (Hawaii NW chain) survive the
  # aggregation. "mean" would average them below sea level and the landmass test
  # would misfire (resolution-sensitive). cluster only needs the land/sea split.
  terra::aggregate(g, fact = fact, fun = "max", na.rm = TRUE,
                   filename = gebco_work, overwrite = TRUE)
  rm(g); gc()
}

# Extract GLORYS mean T for every unique site ONCE here (serial; points only, cheap)
# and attach to sp.df, so the parallel workers NEVER open GLORYS — removing it as a
# per-worker memory cost (64 workers each opening a daily .nc was a second RAM hog).
# GLORYS long-term mean T per site. extract() reads ONLY the cells under the sites
# (never the whole grid), and rowMeans collapses the layers — so this is fast on the
# single-layer .tif (local) OR the ~29-layer ANNUAL .nc (server). Do NOT point this at
# the ~10k-layer DAILY file: summarising that is what made the step crawl.
usites_all <- sp.df |> distinct(LON, LAT)
# Temperature is now REPORT-ONLY (no longer a split criterion), so GLORYS is
# optional: if the raster is unavailable (e.g. a OneDrive online-only stub), run
# with T = NA — the barrier-based splits are unaffected.
glorys <- tryCatch({ g <- terra::rast(glorys_path); terra::crs(g) <- "epsg:4326"; g },
                   error = function(e) NULL)
if (!is.null(glorys)) {
  gv <- terra::extract(glorys, as.matrix(usites_all[, c("LON","LAT")]))          # nearest cell, points only
  usites_all$T <- if (ncol(gv) > 2) rowMeans(gv[, -1], na.rm = TRUE) else gv[, 2]  # mean over annual layers
  rm(glorys); gc()
} else {
  message("GLORYS raster unavailable -> proceeding with T = NA (report-only; splits unaffected).")
  usites_all$T <- NA_real_
}
.tk <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
sp.df <- sp.df |> mutate(.tk = .tk(LON, LAT)) |>
  left_join(usites_all |> mutate(.tk = .tk(LON, LAT)) |> dplyr::select(.tk, T), by = ".tk") |>
  dplyr::select(-.tk)

ecoreg_all <- sort(unique(sp.df$ECOREGION))
ecoreg <- if (is.null(TEST_ECOREGIONS)) ecoreg_all else intersect(ecoreg_all, TEST_ECOREGIONS)
cat("Running on:", paste(ecoreg, collapse = " | "), "\n\n")

# ── Detectors ──────────────────────────────────────────────────────────────
# (a) largest unsampled latitudinal gap -> candidate N/S cut latitude (or NA)
find_lat_gap <- function(s) {
  lat <- sort(unique(round(s$LAT, 3)))
  if (length(lat) < 2) return(NA_real_)
  d <- diff(lat); i <- which.max(d); gap <- d[i]
  if (gap >= LAT_GAP_MIN_DEG && gap >= LAT_GAP_FAC * median(d)) (lat[i] + lat[i + 1]) / 2 else NA_real_
}

# (b) longitudinal doubling: the sites split into two longitude clusters that
#     OVERLAP in latitude (parallel coasts) with LAND between them (coastline
#     "doubles back" around a land mass, e.g. E/W Tasmania). -> candidate E/W cut
#     longitude (or NA). Pooled across latitude (robust to patchy per-band
#     sampling); the lat-overlap test rejects a coastline that merely bends.
find_lon_doubling <- function(s, gz) {
  lon <- sort(unique(round(s$LON, 3)))
  if (length(lon) < 2) return(NA_real_)
  d <- diff(lon); i <- which.max(d)
  if (d[i] < LON_GAP_MIN_DEG) return(NA_real_)
  divide <- (lon[i] + lon[i + 1]) / 2
  L <- s[s$LON <= divide, ]; R <- s[s$LON > divide, ]
  if (nrow(L) < MIN_SITES || nrow(R) < MIN_SITES) return(NA_real_)
  loL <- range(L$LAT); loR <- range(R$LAT)
  ov  <- min(loL[2], loR[2]) - max(loL[1], loR[1])                      # latitude overlap of the two coasts
  if (ov <= 0 || ov / (max(loL[2], loR[2]) - min(loL[1], loR[1])) < LAT_OVERLAP_MIN) return(NA_real_)
  lat_seq <- seq(max(loL[1], loR[1]), min(loL[2], loR[2]), length.out = 20)
  lon_seq <- seq(lon[i], lon[i + 1], length.out = 10)
  dep <- terra::extract(gz, as.matrix(expand.grid(LON = lon_seq, LAT = lat_seq)))[, 1]
  if (mean(dep >= 0, na.rm = TRUE) >= LAND_FRAC_MIN) divide else NA_real_   # land between -> doubling
}

# TRUE open-water disconnection test for a latitudinal split: are the two groups
# on DIFFERENT land masses (mainland Australia vs Tasmania) rather than the same
# continuous coast (GBR north/south, California — both the mainland)?
# Each group is assigned to its nearest LARGE land mass if one is within
# NEARBY_LAND_KM (so a group sitting by small offshore islands — GBR cays, Cape
# Howe's Bass-Strait islets — is folded into the adjacent mainland gradient);
# otherwise it keeps its own (small/far) island identity (Hawaii NW chain). The
# split fires only if the two groups' assigned land masses differ.
groups_diff_landmass <- function(A, B, gz) {
  comp <- terra::patches(terra::ifel(gz >= 0, 1, NA), directions = 8)  # label land masses
  lc <- terra::as.data.frame(comp, xy = TRUE, na.rm = TRUE)
  if (nrow(lc) == 0) return(TRUE)                      # no land at all -> open ocean -> disconnected
  names(lc)[3] <- "pid"
  cell_km2  <- prod(terra::res(comp)) * 111.32^2 * cos(mean(lc$y) * pi / 180)
  large_ids <- as.integer(names(which(table(lc$pid) * cell_km2 >= LARGE_LAND_KM2)))
  assign_land <- function(lon, lat) {
    if (length(large_ids)) {                            # nearest LARGE land mass within NEARBY_LAND_KM?
      L <- lc[lc$pid %in% large_ids, ]
      i <- which.min((L$x - lon)^2 + (L$y - lat)^2)
      dkm <- sqrt(((L$x[i] - lon) * cos(lat * pi / 180))^2 + (L$y[i] - lat)^2) * 111.32
      if (dkm <= NEARBY_LAND_KM) return(L$pid[i])
    }
    lc$pid[which.min((lc$x - lon)^2 + (lc$y - lat)^2)]  # else own (small/far) island
  }
  assign_land(mean(A$LON), mean(A$LAT)) != assign_land(mean(B$LON), mean(B$LAT))
}

# recursive split: returns a list of site data.frames (one per kept group).
# Two gates, applied in order:
#  - LATITUDINAL gap: split only if ALSO thermally distinct (ΔT >= DT_MIN) — a lat
#    gap without a thermal step is just a sampling gap on a continuous gradient.
#  - LONGITUDINAL doubling: split on the LAND BARRIER alone (Δlon gap + land
#    between), NO temperature requirement — a land mass between two coasts is a
#    hard barrier (fish can't cross), the peak would otherwise land ON LAND and
#    relocation would swap coasts, and the coasts are distinct oceanographic /
#    warming regimes (e.g. EAC-warm east vs cool west Tasmania) even when their
#    mean temperature matches. Justification follows Martins et al. (2024).
# Lat is tried first so a land mass is separated latitudinally before any E/W cut
# (keeps the mainland whole instead of cutting it E/W).
# `tag` carries a directional label for the leaf group: a LON doubling stamps the
# two children "W"/"E" (their distinguishing axis); a LAT split leaves the tag
# empty (those groups are labelled by compass-from-ecoregion-centroid later, so
# the northern mainland reads e.g. "NW" rather than just "N"). Result for Bassian:
# W-Tas -> _W, E-Tas -> _E, mainland -> _NW.
split_group <- function(s, gz, tag = "") {
  ok_size <- function(A, B) nrow(A) >= MIN_SITES && nrow(B) >= MIN_SITES
  lc <- find_lat_gap(s)
  if (!is.na(lc)) { A <- s[s$LAT <= lc, ]; B <- s[s$LAT > lc, ]
    if (ok_size(A, B) && groups_diff_landmass(A, B, gz))  # spatial gap + barrier; NO thermal gate
      return(c(split_group(A, gz, tag), split_group(B, gz, tag))) }
  nc <- find_lon_doubling(s, gz)
  if (!is.na(nc)) { A <- s[s$LON <= nc, ]; B <- s[s$LON > nc, ]
    if (ok_size(A, B))                                   # land barrier — no ΔT gate
      return(c(split_group(A, gz, "W"), split_group(B, gz, "E"))) }
  s$gtag <- tag
  list(s)
}

# ── Per-ecoregion processing (one self-contained chunk) ────────────────────
# Rasters are passed in already-opened so the caller controls fork-safety.
process_ecoregion <- function(eco, gebco) {
  usite <- sp.df |> filter(ECOREGION == eco) |> distinct(LON, LAT, T) |> as.data.frame()  # T pre-attached

  # crop the (already coarse, pre-aggregated) GEBCO to this ecoregion's AOI chunk
  ext_pad <- terra::ext(min(usite$LON)-AOI_PAD_DEG, max(usite$LON)+AOI_PAD_DEG,
                        min(usite$LAT)-AOI_PAD_DEG, max(usite$LAT)+AOI_PAD_DEG)
  gz <- terra::crop(gebco, ext_pad)

  groups <- split_group(usite, gz)
  for (k in seq_along(groups)) groups[[k]]$group <- k
  usite_g <- bind_rows(groups); K <- length(groups)

  # ── Directional suffix per group ─────────────────────────────────────────
  # LON-doubling children keep their "W"/"E" tag (their distinguishing axis);
  # LAT-split / unsplit groups get an 8-point compass label from the ecoregion
  # centroid (so the northern mainland reads "NW", the Tasmania coasts "W"/"E").
  # K == 1 (no split) -> empty suffix (SPECIES unchanged downstream).
  eco_cLON <- mean(usite$LON); eco_cLAT <- mean(usite$LAT)
  sp_lon <- diff(range(usite$LON)); sp_lat <- diff(range(usite$LAT))
  compass <- function(clon, clat) {
    ns <- if (clat - eco_cLAT >  0.15 * sp_lat) "N" else if (clat - eco_cLAT < -0.15 * sp_lat) "S" else ""
    ew <- if (clon - eco_cLON >  0.15 * sp_lon) "E" else if (clon - eco_cLON < -0.15 * sp_lon) "W" else ""
    lab <- paste0(ns, ew); if (nzchar(lab)) lab else "C"
  }
  gtag_k <- vapply(seq_len(K), function(k) groups[[k]]$gtag[1], character(1))
  cl     <- usite_g |> group_by(group) |> summarise(clon = mean(LON), clat = mean(LAT), .groups = "drop")
  suf <- vapply(seq_len(K), function(k)
    if (nzchar(gtag_k[k])) gtag_k[k] else compass(cl$clon[cl$group == k], cl$clat[cl$group == k]),
    character(1))
  if (K == 1) suf <- ""
  if (any(duplicated(suf))) suf <- make.unique(suf, sep = "")
  suffix_map <- tibble::tibble(group = seq_len(K), suffix = suf)
  usite_g <- usite_g |> left_join(suffix_map, by = "group") |> dplyr::select(-gtag)

  key <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
  usite_g$.k <- key(usite_g$LON, usite_g$LAT)
  s <- sp.df |> filter(ECOREGION == eco) |> mutate(.k = key(LON, LAT)) |>   # sp.df already carries T
       left_join(usite_g |> dplyr::select(.k, group) |> distinct(.k, .keep_all = TRUE), by = ".k")

  spy <- s |> group_by(group, YEAR) |> summarise(ns = n_distinct(.k), .groups = "drop") |>
    group_by(group) |> summarise(med_site_yr = median(ns), min_site_yr = min(ns),
                                 yrs_ge5 = sum(ns >= 5), .groups = "drop")
  gsum <- s |> group_by(group) |>
    summarise(n_sites = n_distinct(.k), n_years = n_distinct(YEAR),
              meanT = round(mean(T, na.rm = TRUE), 2),
              cLON = round(mean(LON), 2), cLAT = round(mean(LAT), 2), .groups = "drop") |>
    left_join(spy, by = "group") |>
    left_join(suffix_map, by = "group") |>
    mutate(fit = ifelse(yrs_ge5 >= MIN_YEARS, "OK",
                 ifelse(yrs_ge5 >= MIN_YEARS - 1, "MARGINAL", "NO")))
  sp_occ <- s |> group_by(SPECIES, group, YEAR) |>
    summarise(ns = n_distinct(.k), .groups = "drop") |> filter(ns >= 5) |>
    group_by(SPECIES, group) |> summarise(n_years = n_distinct(YEAR), .groups = "drop") |>
    filter(n_years >= MIN_YEARS) |>
    group_by(SPECIES) |>
    summarise(n_groups_occupied = n_distinct(group),
              groups = paste(sort(unique(group)), collapse = ","), .groups = "drop") |>
    mutate(ECOREGION = eco, spans_multiple = n_groups_occupied > 1)
  dT <- if (nrow(gsum) > 1) round(diff(range(gsum$meanT, na.rm = TRUE)), 2) else 0

  # Per-group viability in the one-liner: <suffix>:<n_sites>s/<yrs_ge5>y[<fit>]
  #   n_sites = unique sites; yrs_ge5 = years with >=5 sites (modskurt needs
  #   >= MIN_YEARS of these); fit = OK / MARGINAL / NO. A "NO"/"MARGINAL" group
  #   means the split left too few site-years to fit reliably -> reconsider.
  grp_lbl <- paste(sprintf("%s:%ds/%dy[%s]",
                           ifelse(nzchar(gsum$suffix), gsum$suffix, paste0("g", gsum$group)),
                           gsum$n_sites, gsum$yrs_ge5, gsum$fit), collapse = " ")
  report <- paste0(
    sprintf("[%-28s] %d group(s) | %3d uniq sites | refit %d/%d sp | %s\n",
            eco, K, nrow(usite), sum(sp_occ$spans_multiple), nrow(sp_occ), grp_lbl),
    paste(capture.output(print(as.data.frame(gsum), row.names = FALSE)), collapse = "\n"), "\n")
  list(sites = usite_g |> mutate(ECOREGION = eco) |> dplyr::select(-.k),
       eco_summary = gsum |> mutate(ECOREGION = eco, n_groups = K, dT = dT),
       species = sp_occ, report = report)
}

# ── Run: ecoregions in PARALLEL. Each worker RE-OPENS GEBCO/GLORYS from disk and
#    crops its own AOI window — terra's C++ raster pointers do NOT survive forking,
#    so they cannot be shared from the parent (this is why a plain `for` was used
#    before). doMC forks the rest of the env (helpers, sp.df, params) copy-on-write.
results <- foreach(eco = ecoreg, .packages = c("dplyr","tidyr","terra")) %dopar% {
  gebco <- terra::rast(gebco_work); terra::crs(gebco) <- "epsg:4326"   # only the SMALL cached tif
  tryCatch(process_ecoregion(eco, gebco),
           error = function(e) list(sites = NULL, eco_summary = NULL, species = NULL,
                                    report = sprintf("[%-28s] ERROR: %s\n", eco, conditionMessage(e))))
}
for (r in results) cat(r$report)

# The saved plan holds ONLY the ecoregions that actually split (n_groups > 1) —
# the 46 single-group ecoregions need no action and are omitted (the stdout report
# above still lists all 50 that were evaluated). modskurt_analysis applies a suffix
# only to ecoregions present here; everything else is fit as-is.
all_sites   <- bind_rows(lapply(results, `[[`, "sites"))
all_ecoreg  <- bind_rows(lapply(results, `[[`, "eco_summary"))
all_species <- bind_rows(lapply(results, `[[`, "species"))
split_ecos  <- all_ecoreg |> filter(n_groups > 1) |> distinct(ECOREGION) |> pull(ECOREGION)
# Confirmation gate: keep only recognised barriers (see CONFIRMED_BARRIERS above).
if (!is.null(CONFIRMED_BARRIERS)) {
  dropped    <- setdiff(split_ecos, CONFIRMED_BARRIERS)
  split_ecos <- intersect(split_ecos, CONFIRMED_BARRIERS)
  if (length(dropped))
    cat("Detected but NOT confirmed (excluded from plan):",
        paste(dropped, collapse = " | "), "\n")
}

split_plan <- list(
  sites   = all_sites   |> filter(ECOREGION %in% split_ecos),
  ecoreg  = all_ecoreg  |> filter(ECOREGION %in% split_ecos),
  species = all_species |> filter(ECOREGION %in% split_ecos),   # spans_multiple flags which to refit
  params  = list(MIN_SITES=MIN_SITES, MIN_YEARS=MIN_YEARS,
                 LAT_GAP_MIN_DEG=LAT_GAP_MIN_DEG, LAT_GAP_FAC=LAT_GAP_FAC,
                 LON_GAP_MIN_DEG=LON_GAP_MIN_DEG, LAT_OVERLAP_MIN=LAT_OVERLAP_MIN,
                 LAND_FRAC_MIN=LAND_FRAC_MIN, LARGE_LAND_KM2=LARGE_LAND_KM2, NEARBY_LAND_KM=NEARBY_LAND_KM)
)
out <- file.path(modskurt_dir, if (is.null(TEST_ECOREGIONS)) "split_plan.RData" else "split_plan_TEST.RData")
save(split_plan, file = out)
cat("\nSaved", out, "\nSplit ecoregions:",
    paste(split_plan$ecoreg |> distinct(ECOREGION, n_groups) |> filter(n_groups > 1) |> pull(ECOREGION),
          collapse = " | "), "\n")
