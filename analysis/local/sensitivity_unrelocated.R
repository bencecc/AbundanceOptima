# ============================================================================
# sensitivity_unrelocated.R
#
# SENSITIVITY ANALYSIS — do the main conclusions hold for populations whose
# abundance optimum was essentially NOT relocated?
#
# Subset:  populations with MAX relocation distance == 0 km across ALL their
#          years (never moved onto the reef band).  n = 429 of 2,880.
#          [switch the SUBSET rule below to use status=='ok' (n=319) or a
#           distance cutoff if desired.]
#
# Re-tallies/refits the four core analyses on the subset vs the full dataset:
#   1. Directional shifts  — significant LATITUDINAL trend (Fig 2)          [filter + re-tally]
#   2. P(poleward shift) ~ previous-year exposure (Fig 3)                    [REFIT: GLMM pools populations]
#   3. Counterfactual response — Scenario Buffering/Matching/Increasing (Fig 4a) [filter + re-tally]
#   4. Best mitigating-path outcome (Fig 4b)                                 [filter + re-tally]
#
# Prints FULL vs SUBSET side by side for each. Nothing is saved.
# Data kind = abundance, reference = firstsite (the primary panels).
# ============================================================================

.libPaths(c("C:/R_libs/win-library/4.5", .libPaths()))
suppressMessages({require(dplyr); require(glmmTMB)})
source("config.R")

lo    <- function(f) get(load(input_file(f))[1])
strip <- function(x) sub("_(W|E|NW|SE)$", "", x)   # not needed (keys already orig ECO + suffixed SPECIES)
DIST_CUTOFF <- 0        # km; a population is "unrelocated" if max(dist_km) <= this (0 -> 429 pops)

# ---- 0. Define the unrelocated subset + per-population relocation distance -
rel <- lo("optim.abund.relocated.RData") |>
  transmute(ECOREGION, SPECIES, YEAR, dist_km = ifelse(status == "ok", 0, dist_km))
popdist <- rel |> group_by(ECOREGION, SPECIES) |>
  summarise(maxd = max(dist_km, na.rm = TRUE), meand = mean(dist_km, na.rm = TRUE), .groups = "drop")
sub_pops <- popdist |> filter(maxd <= DIST_CUTOFF) |> distinct(ECOREGION, SPECIES)
cat(sprintf("Unrelocated subset (max move <= %g km): %d populations\n\n",
            DIST_CUTOFF, nrow(sub_pops)))
in_sub <- function(df) semi_join(df, sub_pops, by = c("ECOREGION", "SPECIES"))

# sp.optim.abund.shift — realized optima + exposure (used to build transitions)
sh <- lo("sp.optim.abund.shift.RData") |>
  rename(ECOREGION = ORIG.ECOREGION, LAT = SHIFTED.LAT, LON = SHIFTED.LON)

# ===========================================================================
# 1. DIRECTIONAL SHIFTS — significant latitudinal trend
# ===========================================================================
lat <- lo("optim.abund.trend.RData") |>
  filter(FlagAnalysis == "Trend_Lat", Effect == "Lat.Trend") |>
  group_by(ECOREGION, SPECIES) |> slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup() |>
  filter(!is.na(Estimate)) |>
  mutate(dir = case_when(
    P.value < 0.05 & Estimate > 0 ~ "Sig_poleward",
    P.value < 0.05 & Estimate < 0 ~ "Sig_equatorward",
    TRUE ~ "NS"))

lat_summary <- function(d, lab) {
  n <- nrow(d)
  cat(sprintf("%-8s n=%4d | sig %5.1f%% (poleward %4.1f%% / equatorward %4.1f%%) | NS %4.1f%%\n",
      lab, n,
      100 * mean(d$dir != "NS"),
      100 * mean(d$dir == "Sig_poleward"),
      100 * mean(d$dir == "Sig_equatorward"),
      100 * mean(d$dir == "NS")))
}
cat("\n=== 1. DIRECTIONAL (LATITUDINAL) SHIFTS ===\n")
lat_summary(lat,          "FULL")
lat_summary(in_sub(lat),  "SUBSET")

# ===========================================================================
# 2. P(POLEWARD SHIFT) ~ previous-year exposure   (REFIT on the subset)
# ===========================================================================
lt <- lo("sp.optim.abund.shift.leadtime.RData")
if ("ORIG.ECOREGION" %in% names(lt)) lt <- lt |> rename(ECOREGION = ORIG.ECOREGION)
lt <- lt |> rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON, diff_prev = diff.PreviousSite.temp.mean)

dat <- sh |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON, realized_expo = shifted.cumtemp.above.mean) |>
  left_join(lt, by = c("ECOREGION", "SPECIES", "YEAR", "LAT", "LON"))

trans <- dat |> arrange(ECOREGION, SPECIES, YEAR) |>
  group_by(ECOREGION, SPECIES) |>
  mutate(abs_lat = abs(LAT), prev_abs_lat = lag(abs_lat), prev_year = lag(YEAR),
         yrs = YEAR - prev_year, prev_exposure = realized_expo - diff_prev,
         poleward_km = (abs_lat - prev_abs_lat) * 111) |> ungroup() |>
  filter(!is.na(prev_abs_lat), !is.na(prev_exposure), yrs == 1) |>
  mutate(pop = paste(ECOREGION, SPECIES, sep = "::"),
         poleward_yn = as.numeric(poleward_km > 0))

# GLOBAL model P(poleward) ~ poly2(prev exposure) + (1|pop), binomial. The strict
# 429-pop subset is underpowered, so we report a POWER GRADIENT across increasingly
# relaxed relocation thresholds to show the effect recovers with sample size (i.e.
# it is not an artifact of relocation).
fit_pole <- function(keys, lab) {
  d <- trans |> semi_join(keys, by = c("ECOREGION", "SPECIES"))
  d$prev_exposure_z <- as.numeric(scale(d$prev_exposure))
  m  <- try(glmmTMB(poleward_yn ~ poly(prev_exposure_z, 2) + (1 | pop), data = d, family = binomial), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("%-22s MODEL FAILED\n", lab)); return(invisible()) }
  m0 <- glmmTMB(poleward_yn ~ (1 | pop), data = d, family = binomial)
  p  <- anova(m0, m)$`Pr(>Chisq)`[2]
  q  <- summary(m)$coefficients$cond["poly(prev_exposure_z, 2)2", 1]
  cat(sprintf("%-22s pops=%4d  transitions=%5d  LRT P=%-9.3g quad=%+.2f (%s)\n",
      lab, dplyr::n_distinct(d$pop), nrow(d), p, q, ifelse(q < 0, "concave-down", "up")))
}
cat("\n=== 2. P(POLEWARD SHIFT) ~ PREVIOUS-YEAR EXPOSURE — power gradient ===\n")
fit_pole(popdist,                          "Full (2880)")
fit_pole(filter(popdist, maxd  == 0),      "Never displaced (429)")
fit_pole(filter(popdist, meand <= 1),      "Mean <=1 km (650)")
fit_pole(filter(popdist, meand <= 2),      "Mean <=2 km (922)")

# ===========================================================================
# 3. COUNTERFACTUAL RESPONSE — Mitigation vs Matching + Magnification
# ===========================================================================
# Two headline RESPONSE categories, matching the main text: strict Mitigation, and
# Matching + Magnification (31 + 41 = 72% in the full analysis). The remaining ~4%
# (Acceleration + Annihilation) fall in neither and are footnoted in Extended Data Table 2.
MM_RESP <- c("Matching-Steady", "Matching-Warming", "Matching-Cooling", "Magnification")
scen <- lo("scenario.abund.leadtime.firstsite.longterm.RData") |>
  filter(FlagAnalysis == "Trend_Environmental",
         grepl("cumtemp.above.mean", mod.parms, fixed = TRUE), Effect == "Slope.Dev") |>
  group_by(ECOREGION, SPECIES) |> slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup() |>
  filter(!is.na(Response))

cf_summary <- function(d, lab) {
  cat(sprintf("%-8s n=%4d | mitigate (Mitigation) %4.1f%% | matching+magnification %4.1f%% | residual (accel+annih) %4.1f%%\n",
      lab, nrow(d),
      100 * mean(d$Response == "Mitigation"),
      100 * mean(d$Response %in% MM_RESP),
      100 * mean(!(d$Response %in% c("Mitigation", MM_RESP)))))
}
cat("\n\n=== 3. COUNTERFACTUAL RESPONSE (mitigate vs matched/increased) ===\n")
cf_summary(scen,         "FULL")
cf_summary(in_sub(scen), "SUBSET")

# ===========================================================================
# 4. BEST MITIGATING-PATH OUTCOME (Fig 4b)
# ===========================================================================
bp <- lo("bestpath.opt.abund.firstsite.RData") |> filter(is.na(drop_reason), resolvable,
                        is.finite(p_best_lt_ref), is.finite(p_best_lt_obs)) |>
  mutate(outcome = case_when(
    !(p_best_lt_ref < 0.05) ~ "No available path",
    (p_best_lt_obs < 0.05)  ~ "Available, not followed",
    TRUE                    ~ "Available and followed"))

path_summary <- function(d, lab) {
  n <- nrow(d)
  cat(sprintf("%-8s n=%4d | no-path %4.1f%% | not-followed %4.1f%% | followed %4.1f%%\n",
      lab, n,
      100 * mean(d$outcome == "No available path"),
      100 * mean(d$outcome == "Available, not followed"),
      100 * mean(d$outcome == "Available and followed")))
}
cat("\n=== 4. BEST MITIGATING-PATH OUTCOME ===\n")
cat("(file:", bpf, ")\n")
path_summary(bp,         "FULL")
path_summary(in_sub(bp), "SUBSET")
