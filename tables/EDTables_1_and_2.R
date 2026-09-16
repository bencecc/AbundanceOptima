# ==========================================================================
# Extended Data Tables 1 & 2  (reef-fish abundance-optima paper)
#
#   TABLE 1 — Reef-fish sampling, peak-abundance shifts and analysis coverage
#             by marine ecoregion (one row per ecoregion entering >=1 analysis,
#             plus a whole-dataset total; + the taxa count for the caption).
#   TABLE 2 — Sensitivity of the four main results to optimum relocation
#             (never-displaced subset vs the full dataset).
#
#   Both tables share the same inputs and the same AIC-best model selection.
#   Data kind = abundance, reference = firstsite.
#
#   CONVENTIONS
#   * Table 1 total and Table 2 base are 2,876 populations / 30 ecoregions.
#   * Table 2 row 3 (counterfactual) uses the MAIN-TEXT split — Mitigation vs
#     matching + magnification (24% / 72% full; 15% / 82% subset) — with the ~4%
#     residual (acceleration + annihilation) footnoted.
# ==========================================================================

.libPaths(c("C:/R_libs/win-library/4.5", .libPaths()))
suppressMessages({library(dplyr); library(tidyr); library(glmmTMB)})
setwd("~/Lavori/MPA_timeseries/Modskurt")

lo    <- function(f) get(load(f)[1])
comma <- function(n) formatC(n, big.mark = ",", format = "d")
# P-value formatted as in Table 2: "4 x 10^-21" for tiny P, else 3 decimals.
fmt_p <- function(p) if (p < 1e-3) {
  e <- floor(log10(p)); m <- round(p / 10^e); sprintf("%d x 10^%d", m, e)
} else sprintf("%.3f", p)

# ==========================================================================
# SHARED DATA + DERIVED OBJECTS (used by both tables)
# ==========================================================================
load("sp.df.RData")                                  # raw survey (sites, years)
load("sp.optim.abund.shift.RData")                   # analysed optima
load("optim.abund.trend.RData")                      # latitudinal + env trends
scenario <- lo("scenario.abund.leadtime.firstsite.longterm.RData")   # counterfactual
bpf <- c("bestpath_optimum/bestpath.opt.abund.firstsite.RData",
         "bestpath.opt.abund.firstsite.RData")
bestpath <- lo(bpf[file.exists(bpf)][1])             # best-path

# realised optima, two views: `shift` keeps SHIFTED.* (Table 1 counts);
# `sh` renames to LAT/LON (Table 2 transitions).
shift <- sp.optim.abund.shift |> rename(ECOREGION = ORIG.ECOREGION)
sh    <- sp.optim.abund.shift |> rename(ECOREGION = ORIG.ECOREGION, LAT = SHIFTED.LAT, LON = SHIFTED.LON)

# AIC-best latitudinal trend per population (+ direction / significance)
lat <- optim.abund.trend |>
  filter(FlagAnalysis == "Trend_Lat", Effect == "Lat.Trend") |>
  group_by(ECOREGION, SPECIES) |> slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup() |>
  filter(!is.na(Estimate)) |>
  mutate(dir = case_when(P.value < 0.05 & Estimate > 0 ~ "Sig_poleward",
                         P.value < 0.05 & Estimate < 0 ~ "Sig_equatorward",
                         TRUE ~ "NS"))

# AIC-best counterfactual response per population (cumtemp.above.mean Slope.Dev)
scen <- scenario |>
  filter(FlagAnalysis == "Trend_Environmental",
         grepl("cumtemp.above.mean", mod.parms, fixed = TRUE), Effect == "Slope.Dev") |>
  group_by(ECOREGION, SPECIES) |> slice_min(AIC, n = 1, with_ties = FALSE) |> ungroup() |>
  filter(!is.na(Response))

# ==========================================================================
# ======================  PART A — EXTENDED DATA TABLE 1  ==================
# ==========================================================================
# ---- A1. Survey effort — Sites, Time series --------------------------------
sites_eco <- sp.df |> group_by(ECOREGION) |>            # ALL sites in the ecoregion
  summarise(Sites = n_distinct(SITE_ID), .groups = "drop")
# Time series = distinct site x species among analysed pops. Split ecoregions
# (Bassian, Hawaii) suffix SPECIES, so join the raw survey on SPECIES.ORIG.
pops_key <- shift |> distinct(ECOREGION, SPECIES, SPECIES.ORIG)
ts_eco <- sp.df |>
  inner_join(pops_key, by = c("ECOREGION", "SPECIES" = "SPECIES.ORIG"),
             relationship = "many-to-many") |>
  group_by(ECOREGION) |>
  summarise(TimeSeries = n_distinct(paste(SITE_ID, SPECIES)), .groups = "drop")

# ---- A2. Populations, Years, significant latitudinal shifts ----------------
popyr <- shift |> group_by(ECOREGION) |>
  summarise(Populations = n_distinct(SPECIES),
            yr0 = min(YEAR), yr1 = max(YEAR), .groups = "drop")
lat_eco <- lat |> group_by(ECOREGION) |>
  summarise(estimable = n(), LatShifts = sum(P.value < 0.05, na.rm = TRUE), .groups = "drop")

# ---- A3. Analysis coverage  L / P / C / B ----------------------------------
# L (Fig 2) & C (Fig 4a) use the >=5-population per-ecoregion display threshold;
# P (Fig 3) & B (Fig 4b) include any ecoregion with >=1 contributing population.
L_eco <- lat_eco |> filter(estimable >= 5) |> pull(ECOREGION)
P_eco <- shift |> arrange(ECOREGION, SPECIES, YEAR) |>
  group_by(ECOREGION, SPECIES) |> mutate(dy = YEAR - lag(YEAR)) |> ungroup() |>
  filter(dy == 1) |> distinct(ECOREGION) |> pull(ECOREGION)
C_eco <- scen |> count(ECOREGION) |> filter(n >= 5) |> pull(ECOREGION)
B_eco <- bestpath |>
  filter(is.na(drop_reason), resolvable, is.finite(p_best_lt_ref), is.finite(p_best_lt_obs)) |>
  distinct(ECOREGION) |> pull(ECOREGION)
cov_code <- function(e) {
  f <- c(if (e %in% L_eco) "L", if (e %in% P_eco) "P",
         if (e %in% C_eco) "C", if (e %in% B_eco) "B")
  if (length(f) == 4) "All" else paste(f, collapse = ", ")
}

# ---- A4. Assemble ----------------------------------------------------------
full <- popyr |>
  left_join(sites_eco, by = "ECOREGION") |>
  left_join(ts_eco,    by = "ECOREGION") |>
  left_join(lat_eco,   by = "ECOREGION") |>
  rowwise() |> mutate(Analyses = cov_code(ECOREGION)) |> ungroup() |>
  # Lat shifts is NA when the ecoregion is not in the latitudinal (L) analysis
  # (< 5 estimable trends); an ecoregion IN L with no significant shift shows 0.
  transmute(Ecoregion = ECOREGION, Sites, `Time series` = TimeSeries, Populations,
            `Lat shifts` = ifelse(ECOREGION %in% L_eco, LatShifts, NA_integer_),
            Years = paste0(yr0, "-", yr1), Analyses)
# Displayed rows drop Mascarene (enters no analysis) and are ordered by latitude
# (N->S, mean latitude of the optima), matching the manuscript. Total = displayed rows.
eco_lat <- shift |> group_by(ECOREGION) |>
  summarise(mlat = mean(SHIFTED.LAT, na.rm = TRUE), .groups = "drop")
tab <- full |> filter(Analyses != "") |>
  left_join(eco_lat, by = c("Ecoregion" = "ECOREGION")) |>
  arrange(desc(mlat)) |> select(-mlat)
total_row <- tibble(
  Ecoregion = "Total", Sites = sum(tab$Sites), `Time series` = sum(tab$`Time series`),
  Populations = sum(tab$Populations), `Lat shifts` = sum(tab$`Lat shifts`, na.rm = TRUE),
  Years = paste0(min(shift$YEAR), "-", max(shift$YEAR)), Analyses = "")
n_taxa <- shift |> filter(ECOREGION %in% tab$Ecoregion) |> distinct(SPECIES.ORIG) |> nrow()

cat("==================== EXTENDED DATA TABLE 1 (body) ====================\n")
print(as.data.frame(bind_rows(tab, total_row)), row.names = FALSE)
cat(sprintf("\nCaption taxa: %d taxa across %d populations (manuscript: 1,205 species + 52 higher-level; curated).\n",
            n_taxa, total_row$Populations))
cat("Ecoregions entering no analysis (excluded from table):",
    paste(setdiff(full$Ecoregion, tab$Ecoregion), collapse = ", "), "\n\n")

# ==========================================================================
# ======================  PART B — EXTENDED DATA TABLE 2  ==================
# ==========================================================================
DIST_CUTOFF <- 0        # km; "unrelocated" = max(dist_km) <= this

# ---- B0. Unrelocated subset (Mascarene excluded -> base 2,876, subset 428) --
rel <- lo("optim.abund.relocated.RData") |>
  transmute(ECOREGION, SPECIES, YEAR, dist_km = ifelse(status == "ok", 0, dist_km))
popdist <- rel |> group_by(ECOREGION, SPECIES) |>
  summarise(maxd = max(dist_km, na.rm = TRUE), meand = mean(dist_km, na.rm = TRUE), .groups = "drop") |>
  filter(!grepl("Mascaren", ECOREGION, ignore.case = TRUE))
sub_pops <- popdist |> filter(maxd <= DIST_CUTOFF) |> distinct(ECOREGION, SPECIES)
in_sub   <- function(df) semi_join(df, sub_pops, by = c("ECOREGION", "SPECIES"))
cat(sprintf("Table 2 — unrelocated subset (max move <= %g km): %d of %d populations\n\n",
            DIST_CUTOFF, nrow(sub_pops), nrow(popdist)))

# ---- ROW 1. Directional (latitudinal) shifts (Fig 2) -----------------------
lat_cell <- function(d) sprintf("%.1f%% (%.1f%%; %.1f%%); n = %s",
  100 * mean(d$dir != "NS"), 100 * mean(d$dir == "Sig_poleward"),
  100 * mean(d$dir == "Sig_equatorward"), comma(nrow(d)))

# ---- ROW 2. P(poleward shift) ~ previous-year exposure (Fig 3) -------------
# Sensitivity cell fitted on the minimally relocated subset (mean <= 2 km).
lt <- lo("sp.optim.abund.shift.leadtime.RData")
if ("ORIG.ECOREGION" %in% names(lt)) lt <- lt |> rename(ECOREGION = ORIG.ECOREGION)
lt <- lt |> rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON, diff_prev = diff.PreviousSite.temp.mean)
dat <- sh |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON, realized_expo = shifted.cumtemp.above.mean) |>
  left_join(lt, by = c("ECOREGION", "SPECIES", "YEAR", "LAT", "LON"))
trans <- dat |> arrange(ECOREGION, SPECIES, YEAR) |> group_by(ECOREGION, SPECIES) |>
  mutate(abs_lat = abs(LAT), prev_abs_lat = lag(abs_lat), prev_year = lag(YEAR),
         yrs = YEAR - prev_year, prev_exposure = realized_expo - diff_prev,
         poleward_km = (abs_lat - prev_abs_lat) * 111) |> ungroup() |>
  filter(!is.na(prev_abs_lat), !is.na(prev_exposure), yrs == 1) |>
  mutate(pop = paste(ECOREGION, SPECIES, sep = "::"), poleward_yn = as.numeric(poleward_km > 0))
pole_cell <- function(keys) {
  d <- trans |> semi_join(keys, by = c("ECOREGION", "SPECIES"))
  d$z <- as.numeric(scale(d$prev_exposure))
  m  <- glmmTMB(poleward_yn ~ poly(z, 2) + (1 | pop), data = d, family = binomial)
  m0 <- glmmTMB(poleward_yn ~ (1 | pop), data = d, family = binomial)
  list(cell = sprintf("Significant, concave-down (P = %s)", fmt_p(anova(m0, m)$`Pr(>Chisq)`[2])),
       trans = nrow(d))
}

# ---- ROW 3. Counterfactual response (Fig 4a) — Mitigation vs matching+magnif -
MM_RESP <- c("Matching-Steady", "Matching-Warming", "Matching-Cooling", "Magnification")
cf_cell <- function(d) sprintf("%.0f%%; %.0f%%; n = %s",
  100 * mean(d$Response == "Mitigation"), 100 * mean(d$Response %in% MM_RESP), comma(nrow(d)))

# ---- ROW 4. Best mitigating-path outcome (Fig 4b) --------------------------
bp <- bestpath |>
  filter(is.na(drop_reason), resolvable, is.finite(p_best_lt_ref), is.finite(p_best_lt_obs)) |>
  mutate(outcome = case_when(!(p_best_lt_ref < 0.05) ~ "No available path",
                             (p_best_lt_obs < 0.05)  ~ "Available, not followed",
                             TRUE                    ~ "Available and followed"))
path_cell <- function(d) sprintf("%.0f%%; %.0f%%; %.0f%%; n = %s",
  100 * mean(d$outcome == "No available path"),
  100 * mean(d$outcome == "Available, not followed"),
  100 * mean(d$outcome == "Available and followed"), comma(nrow(d)))

# ---- Emit — ED Table 2 cells (Full | Unrelocated) --------------------------
p_full <- pole_cell(popdist)
p_sub  <- pole_cell(filter(popdist, meand <= 2))
n_sub922 <- nrow(filter(popdist, meand <= 2))   # subset SIZE (incl. pops w/o transitions)
row <- function(lab, full, sub) cat(sprintf("%-28s FULL: %s\n%-28s UNREL: %s\n\n", lab, full, "", sub))

cat("==================== EXTENDED DATA TABLE 2 (cells) ====================\n\n")
row("1. Directional (Fig 2)",     lat_cell(lat),   lat_cell(in_sub(lat)))
row("2. P(poleward) (Fig 3)",     p_full$cell,     p_sub$cell)
cat(sprintf("     (Fig 3 subset = minimally relocated, mean <= 2 km: n = %d populations, %d transitions)\n\n",
            n_sub922, p_sub$trans))
row("3. Counterfactual (Fig 4a)", cf_cell(scen),   cf_cell(in_sub(scen)))
row("4. Best-path (Fig 4b)",      path_cell(bp),   path_cell(in_sub(bp)))
cat("Row 3 = Mitigation ; matching + magnification (main-text split); residual ~4% (accel + annih) footnoted.\n")

# ---- optional Word export for Table 1 (uncomment; set your own folder) ------
# require(flextable); require(officer)
# print(flextable(bind_rows(tab, total_row)),
#       target = "~/Lavori/MPA_timeseries/Modskurt/Tables/ExtendedDataTable1.docx")
