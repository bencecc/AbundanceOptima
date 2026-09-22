# ==========================================================================
# Fig3.R
#
# Within-population sequential test of H2: does thermal exposure at the
# previous peak site (year t-1) trigger a poleward shift from t-1 to t?
#
# Predictor : prev_exposure_z = shifted.cumtemp.above.mean at year t-1
#             back-derived as shifted.cumtemp.above.mean[t] - diff.PreviousSite.temp.mean[t]
# Response  : poleward step shift = (|LAT_t| - |LAT_{t-1}|) * 111  (km)
# Model     : glmmTMB, random intercept + slope by population
# Output    : Fig3 (observed points + marginal predicted trend with 95% CI)
#
# PRODUCES (main + Extended Data, selected by the switches at the top):
#   DATA.KIND = "abund", LAG = 1  -> Fig 3 (main).
#   ED Fig 3 (sensitivity) — re-run once for each of the three panels:
#     (a) DATA.KIND = "abund",   LAG = 2   (abundance, lag-2)
#     (b) DATA.KIND = "density", LAG = 1   (density,   lag-1)
#     (c) DATA.KIND = "density", LAG = 2   (density,   lag-2)
#     (d) abundance, lag-1, EXCLUDING transitions whose previous-year optimum is an
#         edge peak — the SENSITIVITY block after the main figure (any switch setting).
#   The ggsave() call is commented — uncomment and set the output folder.
#   Console prints report the N quoted in the Fig 3 / ED Fig 3 legend.
# ==========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
require(ggh4x)
library(glmmTMB)
require(DHARMa)
library(ggeffects)
library(CockR)
require(cluster)   # daisy()
require(ape)       # pcoa()
require(vegan)     # envfit()

# paths via config.R (run with the working directory set to the project root)
source("config.R")

# ---- INPUT SWITCH ------------------------------------------------------------
# DATA.KIND : "abund" | "density". prev_exposure = the exposure realized at the
# previous peak site that year (realized exposure - diff.PreviousSite.temp.mean).
# (Inputs are ALREADY unimodal-fit-QC filtered upstream in the analysis scripts.)
DATA.KIND  <- "abund"     # "abund" | "density"
LAG        <- 1L          # 1 | 2 — exposure lag (years) used as predictor + step length
stopifnot(DATA.KIND %in% c("abund", "density"), LAG %in% c(1L, 2L))

load(input_file("reef_fish_sti_glorys.RData"))   # THERMAL.GUILD + RANGE classification
load(input_file("fish.traits.dat.RData"))        # (abandoned trait supplement only)

load_obj <- function(f) { e <- new.env(); get(load(input_file(f), envir = e)[1], envir = e) }
to_eco <- function(d) {
  if ("ORIG.ECOREGION" %in% names(d)) dplyr::rename(d, ECOREGION = ORIG.ECOREGION) else d
}

# Realized-optimum panel (peak LAT/LON + realized exposure) — ALWAYS the shift file.
sh <- load_obj(paste0("sp.optim.", DATA.KIND, ".shift.RData")) |> to_eco() |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, SPECIES.ORIG, YEAR, LAT, LON,
         realized_expo = shifted.cumtemp.above.mean)

# prev_exposure = exposure realized at the previous site that year
# (realized exposure - diff.PreviousSite.temp.mean).
lt <- load_obj(paste0("sp.optim.", DATA.KIND, ".shift.leadtime.RData")) |> to_eco() |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON,
         diff_prev = diff.PreviousSite.temp.mean)
dat <- sh |>
  left_join(lt, by = c("ECOREGION", "SPECIES", "YEAR", "LAT", "LON")) |>
  mutate(prev_exposure = realized_expo - diff_prev)
cat(sprintf("[%s] rows: %d | prev_exposure NAs: %d\n",
            DATA.KIND, nrow(dat), sum(is.na(dat$prev_exposure))))

# ============================================================================
# (1) Build step-level transitions (consecutive years) ---------------------
# prev_exposure is precomputed above (mode-specific); split species keep their
# suffix on SPECIES and carry SPECIES.ORIG for the STI/trait joins below.
# ----------------------------------------------------------------------------
transitions <- dat |>
  arrange(ECOREGION, SPECIES, YEAR) |>
  group_by(ECOREGION, SPECIES) |>
  mutate(
    abs_lat        = abs(LAT),
    prev_abs_lat   = lag(abs_lat, LAG),
    prev_year      = lag(YEAR, LAG),
    yrs_since_prev = YEAR - prev_year,
    # predictor = exposure LAG years earlier. LAG 1: exposure at the previous site
    # (realized - diff.PreviousSite); LAG 2: realized exposure two years prior.
    prev_exposure  = if (LAG == 1L) realized_expo - diff_prev else lag(realized_expo, LAG),
    poleward_km    = (abs_lat - prev_abs_lat) * 111
  ) |>
  ungroup() |>
  filter(!is.na(prev_abs_lat), !is.na(prev_exposure),
         yrs_since_prev == LAG) |>                      # exactly LAG years apart
  mutate(population_id   = paste(ECOREGION, SPECIES, sep = "::"),
         prev_exposure_z = as.numeric(scale(prev_exposure)),
         yrs             = glmmTMB::numFactor(YEAR))

# --- Add THERMAL.GUILD and RANGE for faceting ----------------------------
sti <- reef_fish_sti_glorys |>
  select(SPECIES, q25lat, q75lat, thetao_mean_TM) |>
  rename(TEMP.OPTIM = thetao_mean_TM)

pop.traits <- transitions |>
  group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
  reframe(MEAN.LAT = mean(LAT, na.rm = TRUE)) |>
  # Split ecoregions suffix SPECIES (Bassian _W/_E/_NW, Hawaii _SE/_NW); reef_fish_sti_glorys
  # is keyed by the CLEAN name, so join on SPECIES.ORIG (else split populations get NA
  # guild/range and are dropped below).
  left_join(sti, by = c("SPECIES.ORIG" = "SPECIES")) |>
  mutate(
    THERMAL.GUILD = factor(ifelse(TEMP.OPTIM >= 23, "Tropical", "Temperate"),
                           levels = c("Tropical", "Temperate")),
    RANGE = case_when(
      is.na(q25lat) | is.na(q75lat) ~ NA_character_,
      q25lat > 0 & q75lat > 0 ~ case_when(
        MEAN.LAT < q25lat ~ "Warm edge",
        MEAN.LAT > q75lat ~ "Cold edge", TRUE ~ "Mid range"),
      q25lat < 0 & q75lat < 0 ~ case_when(
        MEAN.LAT > q75lat ~ "Warm edge",
        MEAN.LAT < q25lat ~ "Cold edge", TRUE ~ "Mid range"),
      q25lat < 0 & q75lat > 0 & abs(q25lat) <  abs(q75lat) ~ case_when(
        MEAN.LAT < q25lat ~ "Warm edge",
        MEAN.LAT > q75lat ~ "Cold edge", TRUE ~ "Mid range"),
      q25lat < 0 & q75lat > 0 & abs(q25lat) >= abs(q75lat) ~ case_when(
        MEAN.LAT > q75lat ~ "Warm edge",
        MEAN.LAT < q25lat ~ "Cold edge", TRUE ~ "Mid range"),
      TRUE ~ NA_character_),
    RANGE = factor(RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))) |>
  select(ECOREGION, SPECIES, THERMAL.GUILD, RANGE)

transitions <- transitions |>
  left_join(pop.traits, by = c("ECOREGION", "SPECIES")) |>
  filter(!is.na(THERMAL.GUILD), !is.na(RANGE))

# Fig 3 legend N: populations with >=1 consecutive-year transition AND an assigned
# thermal-guild + range-position category (expect 2,509 for abund, LAG 1).
cat(sprintf("Fig 3 legend N = %d populations (>=1 transition + assigned guild & range)\n",
            n_distinct(transitions$population_id)))

# Back-transform constants for the plot
pe_mean <- mean(transitions$prev_exposure, na.rm = TRUE)
pe_sd   <- sd(transitions$prev_exposure,   na.rm = TRUE)

cat("N transitions:", nrow(transitions),
    " | N populations:", n_distinct(transitions$population_id), "\n")
cat("\nprev_exposure summary:\n");  print(summary(transitions$prev_exposure))
cat("poleward_km   summary:\n");    print(summary(transitions$poleward_km))
cat("Transitions per population (quantiles):\n")
print(quantile(table(transitions$population_id), c(0, .25, .5, .75, 1)))

# ============================================================================
# (2) Linear binary model: P(poleward shift) ~ prev exposure
# ----------------------------------------------------------------------------
transitions$poleward_yn <- as.numeric(transitions$poleward_km > 0)

m_h2_bin <- glmmTMB(poleward_yn ~ prev_exposure_z * THERMAL.GUILD * RANGE +
            (1 | population_id) +
            ar1(as.factor(YEAR) + 0 | population_id),
            data = transitions, family = binomial(), REML = FALSE)
cat("\nBinary H2 model (P(poleward shift)):\n"); print(summary(m_h2_bin))
car::Anova(m_h2_bin)

# sim_bin <- simulateResiduals(m_h2_bin, n = 250); plot(sim_bin)

# ============================================================================
# (3) Quadratic binary model — captures non-linear log-odds curves per facet
# ----------------------------------------------------------------------------
# transitions$poleward_yn <- as.numeric(transitions$poleward_km > 0)

binned <- transitions |>
  group_by(THERMAL.GUILD, RANGE) |>
  #mutate(bin = ntile(prev_exposure, 10)) |>
  mutate(bin = cut(prev_exposure,                                                                         
                   breaks = unique(quantile(prev_exposure,
                                            probs = seq(0, 1, length.out = 11),
                                            na.rm = TRUE)),
                   include.lowest = TRUE, labels = FALSE)) |>
  group_by(THERMAL.GUILD, RANGE, bin) |>
  summarise(prev_exposure = mean(prev_exposure),
            p_obs = mean(poleward_yn),
            n     = n(),
            se    = sqrt(p_obs * (1 - p_obs) / n),
            .groups = "drop") |>
  mutate(lo = pmax(p_obs - se, 0),
         hi = pmin(p_obs + se, 1))

facet_ranges <- binned |>
  group_by(THERMAL.GUILD, RANGE) |>
  summarise(x_min = min(prev_exposure, na.rm = TRUE),
            x_max = max(prev_exposure, na.rm = TRUE),
            .groups = "drop")

m_h2_bin_q <- glmmTMB(poleward_yn ~ poly(prev_exposure_z, 2) * THERMAL.GUILD * RANGE +
            (1 | population_id) +
            ar1(as.factor(YEAR) + 0 | population_id),
            data = transitions, family = binomial(), REML = FALSE)
cat("\nQuadratic binary H2 model:\n"); print(summary(m_h2_bin_q))
car::Anova(m_h2_bin_q)

# Likelihood-ratio test vs linear
cat("\nLR test (linear vs quadratic):\n")
print(anova(m_h2_bin, m_h2_bin_q))

sim_bin_q <- simulateResiduals(m_h2_bin_q, n = 250)
plot(sim_bin_q)

# Predictions on the quadratic surface, per Guild × Range
pred_bin_q <- ggeffects::ggpredict(m_h2_bin_q,
            terms = c("prev_exposure_z [all]", "THERMAL.GUILD", "RANGE")) |>
  as.data.frame() |>
  rename(prev_exposure_z = x, fit = predicted,
         lwr = conf.low, upr = conf.high,
         THERMAL.GUILD = group, RANGE = facet) |>
  mutate(prev_exposure = prev_exposure_z * pe_sd + pe_mean,
         THERMAL.GUILD = factor(THERMAL.GUILD, levels = levels(transitions$THERMAL.GUILD)),
         RANGE         = factor(RANGE, levels = levels(transitions$RANGE)))

# Trim quadratic curve to per-facet binned-points range
pred_bin_q <- pred_bin_q |>
  inner_join(facet_ranges, by = c("THERMAL.GUILD", "RANGE")) |>
  filter(prev_exposure >= x_min, prev_exposure <= x_max)

# Warm-to-cold range palette from CockR (matches Fig 2 / Fig 3 conventions).
thermal   <- taster_palettes_continuous()[["Thermal_continuous"]]
maitai    <- taster_palettes_continuous()[["Mai Tai_continuous"]]
blueangel <- taster_palettes_continuous()[["Blue Angel_continuous"]]
range_pal <- c("Warm edge" = thermal[200],
               "Mid range" = maitai[70],
               "Cold edge" = blueangel[120])

plot.h2.bin.q <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey70") +
  geom_ribbon(data = pred_bin_q,
              aes(x = prev_exposure, ymin = lwr, ymax = upr, fill = RANGE),
              alpha = 0.3) +
  geom_line(data = pred_bin_q,
            aes(x = prev_exposure, y = fit, colour = RANGE),
            linewidth = 0.9) +
  geom_errorbar(data = binned,
                aes(x = prev_exposure, ymin = lo, ymax = hi),
                colour="grey40", width = 0, linewidth = 0.3) +
  geom_point(data = binned,
             aes(x = prev_exposure, y = p_obs), colour="grey40", size = 1.6) +
  scale_colour_manual(values = range_pal, guide = "none") +
  scale_fill_manual(values = range_pal, guide = "none") +
  facet_grid2(THERMAL.GUILD ~ RANGE, scales = "free_x", independent = "x") +
  scale_y_continuous(limits = c(0.15, 0.75),
                     breaks = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
  labs(x = expression("Cumulative temperature exposure ("*degree*"C" %*% "day)"),
       y = "Probability of poleward shift") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

plot(plot.h2.bin.q)
outdir <- file.path(dir_results, "Fig3_poleward")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
save.path <- file.path(outdir,
                paste0("poleward_shift_binary_quadratic_", DATA.KIND,
                       if (LAG == 2L) "_lag2" else "", ".pdf"))
#ggsave(save.path, plot.h2.bin.q, width = 6, height = 5)

# ============================================================================
# SENSITIVITY — ED Fig 3d: abundance, lag 1, EXCLUDING transitions whose previous-year
# latitude optimum was classed "edge_peak" by the unimodal QC (peak at the boundary of the
# surveyed sites, from where no further poleward move can be recorded). Co-author query:
# is the decline of P(poleward) at high exposure an artefact of optima pinned at the edge?
# Self-contained (uses load_obj, to_eco, sti, range_pal, outdir from above); independent of the
# DATA.KIND / LAG switches, so it runs after any main-figure setting.
# ----------------------------------------------------------------------------
sh_d <- load_obj("sp.optim.abund.shift.RData") |> to_eco() |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, SPECIES.ORIG, YEAR, LAT, LON, realized_expo = shifted.cumtemp.above.mean)
lt_d <- load_obj("sp.optim.abund.shift.leadtime.RData") |> to_eco() |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  select(ECOREGION, SPECIES, YEAR, LAT, LON, diff_prev = diff.PreviousSite.temp.mean)
tr_d <- sh_d |> left_join(lt_d, by = c("ECOREGION", "SPECIES", "YEAR", "LAT", "LON")) |>
  arrange(ECOREGION, SPECIES, YEAR) |> group_by(ECOREGION, SPECIES) |>
  mutate(abs_lat = abs(LAT), prev_abs_lat = lag(abs_lat), prev_year = lag(YEAR),
         prev_exposure = realized_expo - diff_prev, poleward_yn = as.numeric(abs_lat > prev_abs_lat)) |>
  ungroup() |> filter(!is.na(prev_abs_lat), !is.na(prev_exposure), YEAR - prev_year == 1) |>
  mutate(population_id = paste(ECOREGION, SPECIES, sep = "::"))

# previous-year LAT class from the QC (SPECIES with dots; ECOREGION_ID = per-species
# first-appearance index in sp.df; split species -> id 1), then drop edge peaks
load(input_file("sp.df.RData"))
qc_d <- read.csv(input_file("unimodal_fit_check_per_year.csv")) |> filter(DIM == "LAT") |>
  select(SPECIES_dot = SPECIES, ECOREGION_ID, prev_year = YEAR, prev_class = class)
eco_lookup <- bind_rows(
    sp.df |> distinct(SPECIES, ECOREGION) |> group_by(SPECIES) |> mutate(ECOREGION_ID = row_number()) |> ungroup(),
    tr_d |> filter(SPECIES != SPECIES.ORIG) |> distinct(SPECIES, ECOREGION) |> mutate(ECOREGION_ID = 1L)) |>
  mutate(SPECIES_dot = gsub(" ", ".", SPECIES))
qc_d <- qc_d |> inner_join(eco_lookup, by = c("SPECIES_dot", "ECOREGION_ID")) |> select(ECOREGION, SPECIES, prev_year, prev_class)
n0 <- nrow(tr_d)
tr_d <- tr_d |> left_join(qc_d, by = c("ECOREGION", "SPECIES", "prev_year")) |> filter(!(prev_class %in% "edge_peak"))
cat(sprintf("ED Fig 3d: dropped %d of %d transitions (previous-year optimum = edge peak)\n", n0 - nrow(tr_d), n0))

# guild x range per population (same rule as above)
traits_d <- tr_d |> group_by(ECOREGION, SPECIES, SPECIES.ORIG) |> reframe(MEAN.LAT = mean(LAT)) |>
  left_join(sti, by = c("SPECIES.ORIG" = "SPECIES")) |>
  mutate(THERMAL.GUILD = factor(ifelse(TEMP.OPTIM >= 23, "Tropical", "Temperate"), levels = c("Tropical", "Temperate")),
         RANGE = case_when(
           is.na(q25lat) | is.na(q75lat) ~ NA_character_,
           q25lat > 0 & q75lat > 0 ~ case_when(MEAN.LAT < q25lat ~ "Warm edge", MEAN.LAT > q75lat ~ "Cold edge", TRUE ~ "Mid range"),
           q25lat < 0 & q75lat < 0 ~ case_when(MEAN.LAT > q75lat ~ "Warm edge", MEAN.LAT < q25lat ~ "Cold edge", TRUE ~ "Mid range"),
           abs(q25lat) < abs(q75lat) ~ case_when(MEAN.LAT < q25lat ~ "Warm edge", MEAN.LAT > q75lat ~ "Cold edge", TRUE ~ "Mid range"),
           TRUE ~ case_when(MEAN.LAT > q75lat ~ "Warm edge", MEAN.LAT < q25lat ~ "Cold edge", TRUE ~ "Mid range")),
         RANGE = factor(RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))) |>
  select(ECOREGION, SPECIES, THERMAL.GUILD, RANGE)
tr_d <- tr_d |> left_join(traits_d, by = c("ECOREGION", "SPECIES")) |> filter(!is.na(THERMAL.GUILD), !is.na(RANGE))
cat(sprintf("ED Fig 3d: %d transitions, %d populations\n", nrow(tr_d), n_distinct(tr_d$population_id)))

# same quadratic binomial GLMM as Fig 3, binned points, curves trimmed to the binned range
pe_mean_d <- mean(tr_d$prev_exposure); pe_sd_d <- sd(tr_d$prev_exposure)
tr_d$prev_exposure_z <- (tr_d$prev_exposure - pe_mean_d) / pe_sd_d
m_d <- glmmTMB(poleward_yn ~ poly(prev_exposure_z, 2) * THERMAL.GUILD * RANGE +
                 (1 | population_id) + ar1(as.factor(YEAR) + 0 | population_id),
               data = tr_d, family = binomial(), REML = FALSE)
binned_d <- tr_d |> group_by(THERMAL.GUILD, RANGE) |>
  mutate(bin = cut(prev_exposure, breaks = unique(quantile(prev_exposure, probs = seq(0, 1, length.out = 11))),
                   include.lowest = TRUE, labels = FALSE)) |>
  group_by(THERMAL.GUILD, RANGE, bin) |>
  summarise(prev_exposure = mean(prev_exposure), p_obs = mean(poleward_yn), n = n(),
            se = sqrt(p_obs * (1 - p_obs) / n), .groups = "drop") |>
  mutate(lo = pmax(p_obs - se, 0), hi = pmin(p_obs + se, 1))
pred_d <- tr_d |> group_by(THERMAL.GUILD, RANGE) |>
  reframe(prev_exposure = seq(min(prev_exposure), max(prev_exposure), length.out = 80)) |>
  mutate(prev_exposure_z = (prev_exposure - pe_mean_d) / pe_sd_d)
X_d <- model.matrix(delete.response(terms(m_d)), pred_d)
eta_d <- as.vector(X_d %*% fixef(m_d)$cond); se_d <- sqrt(rowSums((X_d %*% vcov(m_d)$cond) * X_d))
pred_d <- pred_d |> mutate(fit = plogis(eta_d), lwr = plogis(eta_d - 1.96 * se_d), upr = plogis(eta_d + 1.96 * se_d)) |>
  inner_join(binned_d |> group_by(THERMAL.GUILD, RANGE) |> summarise(x_min = min(prev_exposure), x_max = max(prev_exposure), .groups = "drop"),
             by = c("THERMAL.GUILD", "RANGE")) |>
  filter(prev_exposure >= x_min, prev_exposure <= x_max)

plot.h2.bin.q.noedge <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey70") +
  geom_ribbon(data = pred_d, aes(x = prev_exposure, ymin = lwr, ymax = upr, fill = RANGE), alpha = 0.3) +
  geom_line(data = pred_d, aes(x = prev_exposure, y = fit, colour = RANGE), linewidth = 0.9) +
  geom_pointrange(data = binned_d, aes(x = prev_exposure, y = p_obs, ymin = lo, ymax = hi),
                  colour = "grey30", size = 0.3, linewidth = 0.4) +
  scale_colour_manual(values = range_pal, guide = "none") +
  scale_fill_manual(values = range_pal, guide = "none") +
  facet_grid2(THERMAL.GUILD ~ RANGE, scales = "free_x", independent = "x") +
  labs(x = "Previous-year cumulative exposure (degree-days)", y = "P(poleward shift)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
plot(plot.h2.bin.q.noedge)
save.path.d <- file.path(outdir, "poleward_shift_binary_quadratic_abund_noedge.pdf")
#ggsave(save.path.d, plot.h2.bin.q.noedge, width = 6, height = 5)






