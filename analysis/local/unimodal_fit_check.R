# ==========================================================================
# unimodal_fit_check.R
#
# Quality check on modskurt unimodal fits per population × year × dimension.
#   DIM = the gradient the curve was fitted along, LAT or LON: modskurt_analysis.R
#   fits abundance separately against latitude and against longitude for every
#   population x year, saved as ..._<year>_LAT.RData and ..._<year>_LON.RData.
# Inspects the FITTED CURVE (mu_mean) from each .RData file (p.df[[2]] in
# cluster_summaries.R convention) and classifies its shape, regardless of
# the number of underlying data points: sample size is gated upstream, since
# modskurt_analysis.R fits a year only with at least 5 surveyed sites and a
# population  only with > 4 years, so the QC judges curve SHAPE only.
#
# Classes:
#   clear_unimodal : 1 internal local maximum, peak away from boundaries,
#                    prominence (max-min)/max above threshold
#   bimodal        : >= 2 internal local maxima
#   edge_peak      : global max within 5% of either x-boundary
#   flat           : single max but prominence below threshold (no real peak)
#
# Per-population aggregate (ECOREGION × SPECIES × DIM) reports % years in
# each class. Populations with consistent deviation from clear_unimodal
# (>= 50% of years) are flagged as candidates for sensitivity-test exclusion.
#
# Outputs (saved to dir_data, from config.R):
#   unimodal_fit_check_per_year.csv  — one row per file
#   unimodal_fit_check_summary.csv   — one row per (SPECIES × ECOREGION × DIM)
# ==========================================================================

require(dplyr)
require(stringr)
require(foreach)
require(doMC)
source("config.R")

# ---- KIND switch: set ONCE; drives the scan folder AND the output suffix ----
KIND     <- "abund"                                  # "abund" | "density"
stopifnot(KIND %in% c("abund", "density"))
KIND_DIR <- if (KIND == "density") "Density" else "Abund"
sfx      <- if (KIND == "density") "_density" else ""   # output filename suffix

mod_root <- dir_data
src_dir  <- file.path(mod_root, paste0("ModskurtOptimEcoregPlot_", KIND_DIR))
out_dir  <- mod_root

setwd(src_dir)
file_list <- list.files(pattern = "\\.RData$")
cat("N files to scan:", length(file_list), "\n")

# Parallel backend:
registerDoMC(cores = max(1, parallel::detectCores() - 1))

# --- Tunable thresholds --------------------------------------------------
edge_thr  <- 0.05   # peak within 5% of either x-boundary => edge_peak
flat_thr  <- 0.20   # (max - min)/max < 20% of peak => flat
consis_thr <- 50    # >=50% deviating years => consistent_deviation flag

# --- Peak classifier on the fitted curve ---------------------------------
classify_fit <- function(x, mu) {
  if (length(mu) < 3 || any(!is.finite(mu))) return(NULL)
  ord <- order(x)
  x <- x[ord]; mu <- mu[ord]
  span <- diff(range(x))
  if (span <= 0 || max(mu) <= 0) return(NULL)

  # Internal local maxima: mu[i] > mu[i-1] & mu[i] > mu[i+1] for i = 2..n-1
  n <- length(mu)
  is_internal_max <- c(FALSE,
                       (mu[2:(n-1)] > mu[1:(n-2)]) & (mu[2:(n-1)] > mu[3:n]),
                       FALSE)
  n_internal_max <- sum(is_internal_max)

  pk_idx <- which.max(mu)
  pk_x   <- x[pk_idx]
  pk_rel <- (pk_x - x[1]) / span
  is_edge <- pk_rel < edge_thr | pk_rel > 1 - edge_thr

  prom <- (max(mu) - min(mu)) / max(mu)

  cls <- if (n_internal_max >= 2)  "bimodal"
         else if (is_edge)         "edge_peak"
         else if (prom < flat_thr) "flat"
         else                      "clear_unimodal"

  list(class = cls, n_modes = n_internal_max,
       peak_x = pk_x, peak_rel = pk_rel, prominence = prom)
}

# --- Per-file scan ------------------------------------------------------
fit_diag <- foreach(i = seq_along(file_list),
                    .combine = "rbind") %dopar% {
  fname <- file_list[[i]]
  m <- stringr::str_match(fname,
        "^(.+)_ecoreg_(\\d+)_(\\d{4})_(LAT|LON)\\.RData$")
  if (is.na(m[1, 1])) return(NULL)
  sp  <- m[1, 2]
  eco <- as.integer(m[1, 3])
  yr  <- as.integer(m[1, 4])
  dim <- m[1, 5]

  obj <- tryCatch(mget(load(fname))[[1]], error = function(e) NULL)
  if (is.null(obj) || length(obj) < 2) return(NULL)
  fit <- obj[[2]]
  if (!all(c(dim, "mu_mean") %in% names(fit))) return(NULL)

  res <- classify_fit(fit[[dim]], fit$mu_mean)
  if (is.null(res)) return(NULL)

  data.frame(SPECIES = sp, ECOREGION_ID = eco, YEAR = yr, DIM = dim,
             class = res$class, n_modes = res$n_modes,
             peak_x = res$peak_x, peak_rel = res$peak_rel,
             prominence = res$prominence,
             stringsAsFactors = FALSE)
}

cat("Files classified:", nrow(fit_diag), "\n\n")
cat("Class breakdown by dimension:\n")
print(table(fit_diag$class, fit_diag$DIM))

# --- Per-population (× DIM) aggregate -----------------------------------
pop_summary <- fit_diag |>
  group_by(SPECIES, ECOREGION_ID, DIM) |>
  summarise(
    n_years            = n(),
    pct_unimodal       = mean(class == "clear_unimodal") * 100,
    pct_bimodal        = mean(class == "bimodal")        * 100,
    pct_edge_peak      = mean(class == "edge_peak")      * 100,
    pct_flat           = mean(class == "flat")           * 100,
    pct_deviating      = 100 - pct_unimodal,
    consistent_deviation = pct_deviating >= consis_thr,
    .groups = "drop") |>
  arrange(desc(pct_deviating))

# --- Per-population (collapsed across DIM): flagged if either LAT or LON
# is consistently deviating, since the modskurt fit is run separately on
# each axis and either failure compromises the peak location.
pop_flag <- pop_summary |>
  group_by(SPECIES, ECOREGION_ID) |>
  summarise(any_dim_consistent_dev = any(consistent_deviation),
            worst_pct_deviating    = max(pct_deviating),
            n_years_min            = min(n_years),
            .groups = "drop") |>
  arrange(desc(worst_pct_deviating))

# --- Per-population, pooling LAT + LON instances (the inclusion unit)
# Each fit (year × DIM) is one instance. A population is KEPT if >=50%
# of its instances across both axes are clear_unimodal; otherwise REMOVE.
pop_keep <- fit_diag |>
  group_by(SPECIES, ECOREGION_ID) |>
  summarise(
    n_instances    = n(),                                   # years × dims
    pct_unimodal   = mean(class == "clear_unimodal") * 100,
    pct_bimodal    = mean(class == "bimodal")        * 100,
    pct_edge_peak  = mean(class == "edge_peak")      * 100,
    pct_flat       = mean(class == "flat")           * 100,
    decision       = ifelse(pct_unimodal >= consis_thr, "KEEP", "REMOVE"),
    .groups = "drop") |>
  arrange(decision, desc(pct_unimodal - 100))               # REMOVEs first, worst first

# --- Save (filenames carry the KIND-driven suffix: "" for abund, "_density") ---
write.csv(fit_diag,
          file.path(out_dir, paste0("unimodal_fit_check_per_year", sfx, ".csv")),
          row.names = FALSE)
write.csv(pop_summary,
          file.path(out_dir, paste0("unimodal_fit_check_summary", sfx, ".csv")),
          row.names = FALSE)
write.csv(pop_flag,
          file.path(out_dir, paste0("unimodal_fit_check_pop_flag", sfx, ".csv")),
          row.names = FALSE)
write.csv(pop_keep,
          file.path(out_dir, paste0("unimodal_fit_check_pop_decision", sfx, ".csv")),
          row.names = FALSE)

#### -----------------------------------------------------------

cat("\nPopulations × DIM summarised:", nrow(pop_summary), "\n")
cat("Populations (× DIM) with consistent deviation (>=", consis_thr,
    "% deviating years):", sum(pop_summary$consistent_deviation), "\n")
cat("Populations flagged in EITHER LAT or LON:",
    sum(pop_flag$any_dim_consistent_dev), "\n")
cat("\nOutputs written to:", out_dir, "\n")


#### FILTER THE OPTIMA ---------------------------------------------------- ####
# Drop populations with dubious (non-unimodal) modskurt peaks and save the filtered
# optima as modskurt.optim.<kind>.unimodal.RData — the ORIGINAL modskurt.optim.<kind>
# .RData is kept intact (NOT overwritten). Object name == file base (.unimodal),
# matching cluster_summaries.R so the relocation reload+rename still works.
# Driven by the SAME KIND switch set at the top (KIND, mod_root, sfx). The
# id->name mapping is sp.df-based (the only order that matches the decision CSV)
# and Bassian E/W split-aware. EXCLUDE_CLASSES deliberately omits edge_peak.
EXCLUDE_CLASSES <- c("bimodal", "flat")          # edge_peak intentionally kept
source(file.path(PROJ, "analysis/local/unimodal_filter_helper.R"))

obj_name     <- paste0("modskurt.optim.", KIND)            # original (unfiltered) object
unimodal_obj <- paste0(obj_name, ".unimodal")             # filtered object: NAME == file base
load(file.path(mod_root, paste0(obj_name, ".RData")))
optim.df <- get(obj_name)
load(file.path(mod_root, "sp.df.RData"))         # ecoregion ORDER for the id map
dec_csv  <- file.path(mod_root, paste0("unimodal_fit_check_pop_decision", sfx, ".csv"))

optim.filt <- unimodal_filter_helper(optim.df, classes = EXCLUDE_CLASSES,
                                     decision_csv = dec_csv, sp_df = sp.df)
assign(unimodal_obj, optim.filt)
save(list = unimodal_obj, file = file.path(mod_root, paste0(unimodal_obj, ".RData")))

# ==========================================================================
# OPTIONAL DROP DIAGNOSTIC
# ==========================================================================
# Identifies files on disk that the main scan dropped (returned NULL on) and
# attributes each to a drop reason. Useful for understanding why N classified
# << N on disk (e.g., density yielded 8,401 classified vs 55,658 on disk).
#
# Set RUN_DROP_DIAGNOSTIC <- FALSE to skip.
# Output: unimodal_fit_check_drop_diagnostic_<density|abund>.csv
# ==========================================================================
RUN_DROP_DIAGNOSTIC <- TRUE
DROP_SAMPLE_N       <- 2000   # sample size; set to Inf to scan all dropped

if (RUN_DROP_DIAGNOSTIC) {
  data_tag <- KIND

  # All file keys on disk
  disk_keys <- stringr::str_match(file_list,
      "^(.+)_ecoreg_(\\d+)_(\\d{4})_(LAT|LON)\\.RData$")
  disk_df <- data.frame(SPECIES = disk_keys[, 2],
                        ECOREGION_ID = as.integer(disk_keys[, 3]),
                        YEAR = as.integer(disk_keys[, 4]),
                        DIM = disk_keys[, 5],
                        file = file_list,
                        stringsAsFactors = FALSE)

  # Anti-join against classified fits to get dropped files
  dropped_df <- disk_df |>
    dplyr::anti_join(fit_diag,
                     by = c("SPECIES", "ECOREGION_ID", "YEAR", "DIM"))
  cat("\n>>> Drop diagnostic\n")
  cat(">>> on disk    :", nrow(disk_df), "\n")
  cat(">>> classified :", nrow(fit_diag), "\n")
  cat(">>> dropped    :", nrow(dropped_df), "\n")

  # Sample dropped files and re-classify with reason tracking
  n_samp <- min(DROP_SAMPLE_N, nrow(dropped_df))
  set.seed(1)
  samp_idx <- sample(seq_len(nrow(dropped_df)), n_samp)
  samp_df  <- dropped_df[samp_idx, ]
  cat(">>> Re-scanning", n_samp, "sampled dropped files for reason ...\n")

  drop_reasons <- foreach(i = seq_len(n_samp), .combine = c) %dopar% {
    fname <- samp_df$file[i]; dm <- samp_df$DIM[i]
    if (is.na(fname))      return("name_regex_fail")
    obj <- tryCatch(mget(load(fname))[[1]], error = function(e) NULL)
    if (is.null(obj))      return("load_fail")
    if (length(obj) < 2)   return("obj_len_lt_2")
    fit <- obj[[2]]
    if (!(dm %in% names(fit)))        return(paste0("missing_col_", dm))
    if (!("mu_mean" %in% names(fit))) return("missing_col_mu_mean")
    mu <- fit$mu_mean; x <- fit[[dm]]
    if (length(mu) < 3)               return("mu_len_lt_3")
    if (any(!is.finite(mu)))          return("mu_nonfinite")
    if (diff(range(x)) <= 0)          return("x_zero_span")
    if (max(mu) <= 0)                 return("max_mu_le_0")
    "classify_other"   # would have classified — rare; suggests staleness
  }

  cat(">>> Drop-reason tally (sample of", n_samp, "):\n")
  print(sort(table(drop_reasons), decreasing = TRUE))

  samp_df$drop_reason <- drop_reasons
  diag_out <- file.path(out_dir,
                paste0("unimodal_fit_check_drop_diagnostic_", data_tag, ".csv"))
  write.csv(samp_df, diag_out, row.names = FALSE)
  cat(">>> Diagnostic written to:", diag_out, "\n")
}
