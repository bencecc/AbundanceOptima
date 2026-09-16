# ==========================================================================
# Fig2.r
#
# Set DATA_KIND ("abund" | "density") ONCE at the top of section 1: it drives
# BOTH the input objects loaded and the output-folder/file suffix ("abundance" /
# "density"), so the whole figure (roses, inset, map, table) is produced per kind.
#
# Figure 2: Ecoregion warming (polygon fill) × movement roses (insets)
#
# Background: MEOW ecoregion polygons filled by SST warming rate
# Overlay:    Directional roses at ecoregion centroids showing
#             movement direction (bearing) + poleward/equatorward colour
#             (blue = poleward/cooling, red = equatorward/warming)
#             + significance layering (solid = sig, faded = all);
#             overlay is done manually in Adobe Illustrator.
#
# PRODUCES (main + Extended Data, selected by DATA_KIND at the top of section 1):
#   DATA_KIND = "abund"   -> Fig 2 (main): warming map + directional roses + inset.
#   DATA_KIND = "density" -> ED Fig 2: density sensitivity — the same map + roses
#                            + inset built on density optima.
#   Re-run once per kind; outputs carry the "_abundance" / "_density" suffix.
#   Console prints report every N quoted in the Fig 2 / ED Fig 2 legend and Results.
#
# Requires: sf, rnaturalearth, rnaturalearthdata, dplyr, ggplot2,
#           ggforce (for inset roses), RColorBrewer
#
# NOTE: You need your MEOW ecoregion shapefile. Adjust the path below.
# ==========================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(CockR)

# ==========================================================================
# 1. LOAD DATA
# ==========================================================================
setwd("~/Lavori/MPA_timeseries/Modskurt")

# ---- DATA-KIND SWITCH: set ONCE here (abundance vs density) --------------------------
# Drives BOTH the input objects loaded below AND the output-folder/file suffix,
# so every downstream section (roses, inset, map, table) uses the chosen data.
DATA_KIND   <- "abund"                                          # "abund" | "density"
stopifnot(DATA_KIND %in% c("abund", "density"))
kind_suffix <- if (DATA_KIND == "density") "density" else "abundance"   # output suffix

# Movement + trend data for the chosen kind
load(sprintf("sp.optim.%s.shift.RData", DATA_KIND))
load(sprintf("optim.%s.trend.RData",   DATA_KIND))

# load species traits and temperature affinity (kind-independent)
load("fish.traits.dat.RData")
load("reef_fish_sti_glorys.RData")

# Set up data frames for analysis (both derived from DATA_KIND)
trend.focal <- get(sprintf("optim.%s.trend", DATA_KIND)) |>
    select(-c(THERMAL.GUILD,TEMP.OPTIM,FEEDING.TYPE)) # thermal properties recomputed from thetao_mean_TM below

sp.focal <- get(sprintf("sp.optim.%s.shift", DATA_KIND)) |>
  rename(LAT=SHIFTED.LAT, LON=SHIFTED.LON, ECOREGION=ORIG.ECOREGION)

# ---- Fig 2 legend N: modskurt-retained input set -----------------------------
# The Mascarene Islands ecoregion (4 populations) yields no estimable latitudinal
# trend and enters no downstream analysis, so the retained set excludes it:
# expect 2,876 populations / 30 ecoregions (the raw shift file holds 2,880 / 31).
masc_eco <- grepl("Mascaren", sp.focal$ECOREGION, ignore.case = TRUE)
cat(sprintf("[Fig 2 legend] input set (retained, ex-Mascarene): %d populations / %d ecoregions\n",
            n_distinct(paste(sp.focal$ECOREGION, sp.focal$SPECIES)[!masc_eco]),
            n_distinct(sp.focal$ECOREGION[!masc_eco])))

# ---- MEOW shapefile ------------------------------------------------------------------
meow <- sf::read_sf('C:/Users/LBenedettiCecchi/OneDrive - University of Pisa/Lavori/MEOWs/Marine_Ecoregions_Of_the_World__MEOW_.shp') 

# ==========================================================================
# 2. ECOREGION-LEVEL WARMING RATE (for polygon fill)
# ==========================================================================
# Warming rates come from the HPC script all_ecoregions_temperature.R,
# which fits a per-cell OLS slope of annual mean SST (1993-2021) from GLORYS and
# aggregates to ecoregion. Two outputs are available:
#   * all_ecoregion_temperatures_coastal.RData  — 0 to -30 m cells only
#   * all_ecoregion_temperatures_all.RData      — all cells inside the polygon

load("temperature_all_ecoregions_coastal.RData")
load("temperature_all_ecoregions_all.RData")

# ==========================================================================
# 3. ECOREGION CENTROIDS
# ==========================================================================

eco_centroids <- sp.focal |>
  group_by(ECOREGION) |>
  summarise(
    cen_lon = mean(LON, na.rm = TRUE),
    cen_lat = mean(LAT, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================================
# 4. POPULATION BEARINGS + EXPOSURE
# ==========================================================================

# Compute SD of centered YEAR per population to back-transform standardized slopes.
# In model: YEAR = (YEAR - YEAR[1]) + 1, then standardize() divides by sd.
sd_year_lookup <- sp.focal |>
  group_by(ECOREGION, SPECIES) |>
  reframe(sd_year  = sd((YEAR - min(YEAR)) + 1),
          mabs_lat = mean(abs(LAT), na.rm = TRUE),   # for lognormal back-transform (response = |lat|)
          mabs_lon = mean(abs(LON), na.rm = TRUE))    # for lognormal back-transform (response = |lon|)

# Lat and Lon slopes per population (for bearing)
lat_trends <- trend.focal |> filter(FlagAnalysis=="Trend_Lat"&Effect=="Lat.Trend") |>
    group_by(ECOREGION,SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    filter(!is.na(Estimate)) |>
    ungroup() |>
    left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
    mutate(
      is_std = grepl("Standardize", mod.parms),
      is_log = grepl("lognormal",  mod.parms),
      Estimate = ifelse(is_std, Estimate / sd_year, Estimate),
      SE       = ifelse(is_std, SE / sd_year, SE),
      # lognormal fits model log(|lat|): slope is d log(|lat|)/dyr. Convert to
      # degrees/yr at the series mean: d|lat|/dyr = beta_log * mean(|lat|).
      Estimate = ifelse(is_log, Estimate * mabs_lat, Estimate),
      SE       = ifelse(is_log, SE * mabs_lat, SE)) |>
    rename(lat_slope = Estimate, lat_se = SE, lat_p = P.value)

lon_trends <- trend.focal |> filter(FlagAnalysis=="Trend_Lon"&Effect=="Lon.Trend") |>
   group_by(ECOREGION,SPECIES) |>
   slice_min(AIC, n = 1, with_ties = FALSE) |>
   filter(!is.na(Estimate)) |>
   ungroup() |>
   left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
   mutate(
     is_std = grepl("Standardize", mod.parms),
     is_log = grepl("lognormal",  mod.parms),
     Estimate = ifelse(is_std, Estimate / sd_year, Estimate),
     SE       = ifelse(is_std, SE / sd_year, SE),
     # lognormal fits model log(|lon|): convert slope to degrees/yr at the mean.
     Estimate = ifelse(is_log, Estimate * mabs_lon, Estimate),
     SE       = ifelse(is_log, SE * mabs_lon, SE)) |>
   rename(lon_slope = Estimate, lon_se = SE, lon_p = P.value) |>
   select(ECOREGION,SPECIES,lon_slope,lon_se,lon_p)


mean.focal <- sp.focal |> group_by(ECOREGION, SPECIES) |>
  reframe(mean_lat = mean(LAT, na.rm = TRUE),
          mean_lon = mean(LON, na.rm = TRUE)) |>   # for signed E/W recovery from the |lon| slope
  left_join(lon_trends, by = c("ECOREGION", "SPECIES"))

# 0-meridian straddlers (Western Med): abs(LON) folds E/W and sign(mean_lon) is
# unstable near 0 deg, so |lon| x sign mis-signs their bearing. Recompute their
# signed longitude slope with the pipeline's estimator (glmmTMB, gaussian,
# autocorr none/OU/AR1, AIC-best; lognormal dropped as it is invalid on signed
# values). Every other population is within one hemisphere, where
# |lon| x sign(mean_lon) already equals the signed slope, so they are untouched.
require(glmmTMB)
fit_signed_lon <- function(d) {  # same estimator as the pipeline, on signed LON
  d <- d |> arrange(YEAR) |> mutate(Yc = (YEAR - YEAR[1]) + 1, grp = 1,
                                    yrs = glmmTMB::numFactor(Yc))
  cand <- list(
    none = try(glmmTMB(LON ~ Yc + (1 | grp), data = d), silent = TRUE),
    ou   = try(glmmTMB(LON ~ Yc + ou(yrs + 0 | grp), data = d), silent = TRUE),
    ar1  = try(glmmTMB(LON ~ Yc + ar1(as.factor(Yc) + 0 | grp), data = d), silent = TRUE))
  aic <- sapply(cand, function(m) if (inherits(m, "try-error") || !is.finite(logLik(m))) Inf else AIC(m))
  co  <- summary(cand[[which.min(aic)]])$coefficients$cond
  tibble(lon_slope_signed = co["Yc", 1], lon_p_signed = co["Yc", 4])
}
straddlers <- sp.focal |>
  group_by(ECOREGION, SPECIES) |>
  filter(any(LON > 0) & any(LON < 0) & max(abs(LON)) < 30) |>
  group_modify(~ fit_signed_lon(.x)) |>
  ungroup()

# Merge with lat trends
pop_bear <- lat_trends |>
  inner_join(mean.focal, by = c("SPECIES", "ECOREGION")) |>
  left_join(straddlers, by = c("ECOREGION", "SPECIES")) |>
  mutate(
    cos_lat  = cos(abs(mean_lat) * pi / 180),
    # Model uses abs(LAT): positive slope = poleward (away from equator).
    # For bearing we need geographic direction, so flip sign for S hemisphere.
    dy_km    = lat_slope * 111 * sign(mean_lat),
    # lon model uses |lon| with sign(mean_lon) recovery, EXCEPT 0-deg straddlers,
    # whose directly-fitted signed slope is used instead (lon_slope_signed non-NA).
    dx_km    = ifelse(!is.na(lon_slope_signed), lon_slope_signed * 111 * cos_lat,
                      lon_slope * 111 * cos_lat * sign(mean_lon)),
    lon_p    = ifelse(!is.na(lon_p_signed), lon_p_signed, lon_p),
    bearing  = (atan2(dx_km, dy_km) * 180 / pi) %% 360,
    speed_km = sqrt(dy_km^2 + dx_km^2),
    # Significance: either lat or lon P < 0.05
    lat_sig  = !is.na(lat_p) & lat_p < 0.05,
    lon_sig  = !is.na(lon_p) & lon_p < 0.05,
    sig      = lat_sig | lon_sig
  )
  
# Mean exposure per population
pop_exposure <- sp.focal |>
  group_by(SPECIES, ECOREGION) |>
  reframe(exposure = mean(shifted.cumtemp.above.mean, na.rm = TRUE))

pop <- pop_bear |>
  left_join(pop_exposure, by = c("SPECIES", "ECOREGION"))

cat("Total populations:", nrow(pop), "\n")
cat("Lat sig:", sum(pop$lat_sig, na.rm = TRUE), "\n")
cat("Lon sig:", sum(pop$lon_sig, na.rm = TRUE), "\n")
cat("Either sig:", sum(pop$sig, na.rm = TRUE), "\n")

# ---- Numbers quoted in the Results text (LATITUDINAL shifts) -----------------
# Direction = sign of the latitudinal slope (lat_slope > 0 = poleward, i.e. away
# from the equator, since the model fits abs(LAT)). Velocity = |lat_slope| in
# km/decade (deg/yr * 111 * 10), averaged over the SIGNIFICANT shifts only, with
# SE of the mean. Expected (abund): 54% / 46% split, ~20% sig (11% / 9%),
# poleward 58 +/- 5 (n=273), equatorward 72 +/- 7 (n=216) km/dec.
lat_pop <- pop |>
  filter(!is.na(lat_slope)) |>
  mutate(dir        = ifelse(lat_slope > 0, "poleward", "equatorward"),
         vel_km_dec = abs(lat_slope) * 111 * 10)
.se <- function(x) sd(x) / sqrt(length(x))
cat(sprintf("\n[Results text] latitudinal shifts, n = %d populations\n", nrow(lat_pop)))
cat(sprintf("  direction : poleward %.0f%% (n=%d) | equatorward %.0f%% (n=%d)\n",
            100 * mean(lat_pop$dir == "poleward"),  sum(lat_pop$dir == "poleward"),
            100 * mean(lat_pop$dir == "equatorward"), sum(lat_pop$dir == "equatorward")))
cat(sprintf("  significant: %.0f%% total | poleward %.0f%% (n=%d) | equatorward %.0f%% (n=%d)\n",
            100 * mean(lat_pop$lat_sig),
            100 * mean(lat_pop$lat_sig & lat_pop$dir == "poleward"),
            sum(lat_pop$lat_sig & lat_pop$dir == "poleward"),
            100 * mean(lat_pop$lat_sig & lat_pop$dir == "equatorward"),
            sum(lat_pop$lat_sig & lat_pop$dir == "equatorward")))
for (dd in c("poleward", "equatorward")) {
  v <- lat_pop$vel_km_dec[lat_pop$lat_sig & lat_pop$dir == dd]
  cat(sprintf("  velocity (sig %-11s): %.0f +/- %.0f km/dec (n=%d)\n",
              dd, mean(v), .se(v), length(v)))
}

eco_summary <- pop |>
  group_by(ECOREGION) |>
  reframe(
    n_total   = n(),
    n_sig     = sum(sig, na.rm = TRUE),
    n_lat_sig = sum(lat_sig, na.rm = TRUE),
    n_lon_sig = sum(lon_sig, na.rm = TRUE)
  )

eco_info <- eco_centroids |>
  inner_join(eco_summary, by = "ECOREGION") |>
  filter(n_total >= 5)

# Order by ABUNDANCE latitude so density uses the SAME ordering as abundance
# (density's own peaks can flip near-tied ecoregions, e.g. Hawaii/Southern China).
ref_lat <- if (DATA_KIND == "abund") {
  eco_centroids |> select(ECOREGION, ref_lat = cen_lat)
} else {
  get(load("sp.optim.abund.shift.RData")[1]) |>
    group_by(ORIG.ECOREGION) |>
    summarise(ref_lat = mean(SHIFTED.LAT, na.rm = TRUE), .groups = "drop") |>
    rename(ECOREGION = ORIG.ECOREGION)
}
eco_info <- eco_info |> left_join(ref_lat, by = "ECOREGION") |>
  arrange(desc(ref_lat)) |> mutate(eco_id = sprintf("%02d", row_number()))

cat("\nEcoregions:", nrow(eco_info), "\n")
# Fig 2 legend N (roses): populations with an estimable latitudinal trend that
# fall in the ecoregions retained for display (n_total >= 5). Expect 2,511 / 28.
cat(sprintf("[Fig 2 legend] roses: %d populations across %d ecoregions (estimable lat trend, n>=5)\n",
            sum(pop$ECOREGION %in% eco_info$ECOREGION), nrow(eco_info)))

# ==========================================================================
# 5. COMPUTE ROSE DATA PER ECOREGION
# ==========================================================================

nbins <- 8
bin_width <- 360 / nbins
bin_edges <- seq(0, 360, length.out = nbins + 1)

bin_bearing <- function(bearing, edges) {
  # Returns bin index 1..nbins
  b <- findInterval(bearing, edges, rightmost.closed = TRUE)
  ifelse(b > length(edges) - 1, length(edges) - 1, b)
}

# Compute binned data for each ecoregion
rose_data <- pop |>
  mutate(bin = bin_bearing(bearing, bin_edges)) |>
  group_by(ECOREGION, bin) |>
  reframe(
    count_all = n(),
    count_sig = sum(sig, na.rm = TRUE)) |>
  mutate(bearing_mid = (bin - 0.5) * bin_width)

# Normalise counts within each ecoregion
rose_data <- rose_data %>%
  group_by(ECOREGION) %>%
  mutate(
    max_count = max(count_all),
    norm_all  = count_all / max_count,
    norm_sig  = count_sig / max_count
  ) %>%
  ungroup()


# ==========================================================================
# 6. DIRECTIONAL COLOUR SCALE (poleward = blue, equatorward = red)
# ==========================================================================

# Poleward score: +1 = straight poleward, -1 = straight equatorward
poleward_score <- function(bearing_deg, cen_lat) {
  score <- cos(bearing_deg * pi / 180)  # +1 at 0° (N), -1 at 180° (S)
  if (cen_lat < 0) score <- -score      # flip for southern hemisphere
  score
}

# 4 discrete direction colours derived from CockR continuous palettes
cockr_cont <- taster_palettes_continuous()
blueangel  <- cockr_cont[["Blue Angel_continuous"]]
aperol     <- cockr_cont[["Aperol Spritz_continuous"]]
dir_cols <- c(
  equatorward         = aperol[103],     # deep red-orange
  equatorward_lateral = aperol[31],      # pink
  poleward_lateral    = blueangel[43],   # light blue
  poleward            = blueangel[156]   # deep blue
)

direction_to_colour <- function(bearing_deg, cen_lat) {
  score <- poleward_score(bearing_deg, cen_lat)
  dplyr::case_when(
    score >= 0.5  ~ dir_cols[["poleward"]],
    score >= 0    ~ dir_cols[["poleward_lateral"]],
    score >= -0.5 ~ dir_cols[["equatorward_lateral"]],
    TRUE          ~ dir_cols[["equatorward"]]
  )
}

# ==========================================================================
# 7. GENERATE INDIVIDUAL ROSE PNGs + PDFs (base R graphics)
# ==========================================================================

outdir <- file.path("~/Lavori/MPA_timeseries/Modskurt/Figs", paste0("Fig2_roses_", kind_suffix))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Helper: draw one rose on the current device
draw_rose <- function(rd, eco_name, cen_lat, n_total, n_sig, half_width) {

  par(mar = c(0, 0, 0, 0))

  plot.new()
  plot.window(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), asp = 1)

  # Background disc: pure white for contrast against dark map
  angles <- seq(0, 2 * pi, length.out = 200)
  polygon(1.05 * cos(angles), 1.05 * sin(angles),
          col = "white", border = "grey50", lwd = 1.5)

  # Concentric reference circles (drawn before petals so petals overlay)
  for (r in c(0.25, 0.5, 0.75, 1.0)) {
    lines(r * cos(angles), r * sin(angles),
          col = "grey80", lwd = 0.6, lty = 1)
  }

  # Petals: convert bearing (0=N, clockwise) to math angle (0=E, counter-clockwise)
  for (j in seq_len(nrow(rd))) {
    bearing_rad <- rd$bearing_mid[j] * pi / 180
    theta_mid <- pi / 2 - bearing_rad

    theta_lo <- theta_mid - half_width
    theta_hi <- theta_mid + half_width

    petal_col <- direction_to_colour(rd$bearing_mid[j], cen_lat)

    # --- All shifts (faded = non-significant) ---
    r_all <- rd$norm_all[j]
    if (r_all > 0) {
      arc <- seq(theta_lo, theta_hi, length.out = 20)
      x_poly <- c(0, r_all * cos(arc), 0)
      y_poly <- c(0, r_all * sin(arc), 0)
      polygon(x_poly, y_poly, col = adjustcolor(petal_col, alpha.f = 0.45),
              border = adjustcolor("grey40", 0.5), lwd = 0.8)
    }

    # --- Significant shifts (solid) ---
    r_sig <- rd$norm_sig[j]
    if (r_sig > 0) {
      arc <- seq(theta_lo, theta_hi, length.out = 20)
      x_poly <- c(0, r_sig * cos(arc), 0)
      y_poly <- c(0, r_sig * sin(arc), 0)
      polygon(x_poly, y_poly, col = adjustcolor(petal_col, alpha.f = 0.95),
              border = "black", lwd = 2)
    }
  }

}

for (i in seq_len(nrow(eco_info))) {

  eco_name <- eco_info$ECOREGION[i]
  eco_id   <- eco_info$eco_id[i]
  cen_lat  <- eco_info$cen_lat[i]
  n_total  <- eco_info$n_total[i]
  n_sig    <- eco_info$n_sig[i]

  # Get binned data for this ecoregion
  rd <- rose_data %>% filter(ECOREGION == eco_name)

  if (nrow(rd) == 0) next

  # File names
  fname_base <- paste0(eco_id, "_", gsub("[/ ]", "_", eco_name))
  half_w <- (bin_width / 2) * pi / 180

  # --- PNG (transparent background) ---
  png(file.path(outdir, paste0(fname_base, ".png")),
      width = 400, height = 400, bg = "transparent")
  draw_rose(rd, eco_name, cen_lat, n_total, n_sig, half_w)
  dev.off()

  # --- PDF (transparent background) ---
  cairo_pdf(file.path(outdir, paste0(fname_base, ".pdf")),
            width = 3, height = 3, bg = "transparent")
  draw_rose(rd, eco_name, cen_lat, n_total, n_sig, half_w)
  dev.off()
}

cat("\nSaved", nrow(eco_info), "roses (PNG + PDF) to", outdir, "/\n")

# ==========================================================================
# 8. STANDALONE COLOUR LEGEND (for map overlay)
# ==========================================================================

draw_legend <- function() {
  par(mar = c(1, 1, 1, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # --- 4 discrete direction swatches ---
  swatch_w <- 0.08
  swatch_h <- 0.15
  labels <- c("Equatorward", "Equatorward-\nlateral", "Poleward-\nlateral", "Poleward")
  cols <- unname(dir_cols)
  x_pos <- c(0.10, 0.30, 0.55, 0.75)

  for (k in seq_along(cols)) {
    rect(x_pos[k], 0.60, x_pos[k] + swatch_w, 0.60 + swatch_h,
         col = cols[k], border = "grey30", lwd = 1.2)
    text(x_pos[k] + swatch_w / 2, 0.50, labels[k], cex = 1.2, adj = 0.5, font = 1)
  }

  # --- Significance key ---
  # Solid swatch
  rect(0.20, 0.18, 0.28, 0.30, col = adjustcolor(dir_cols[["poleward"]], alpha.f = 0.95),
       border = "grey30", lwd = 1)
  text(0.30, 0.24, "Significant shift", cex = 1.2, adj = 0, font = 1)

  # Faded swatch
  rect(0.55, 0.18, 0.63, 0.30, col = adjustcolor(dir_cols[["poleward"]], alpha.f = 0.45),
       border = adjustcolor("grey50", 0.4), lwd = 1)
  text(0.65, 0.24, "Non-significant", cex = 1.2, adj = 0, font = 1)
}

# Save as PNG
png(file.path(outdir, "legend_direction.png"),
    width = 800, height = 250, bg = "transparent")
draw_legend()
dev.off()

# Save as PDF
cairo_pdf(file.path(outdir, "legend_direction.pdf"),
          width = 8, height = 2.5, bg = "transparent")
draw_legend()
dev.off()

cat("Saved legend to", outdir, "/legend_direction.{png,pdf}\n")

# ==========================================================================
# 9. MAP RESULTS
# ==========================================================================

# World map
world_map <- rnaturalearth::ne_countries(scale = 10, returnclass = c("sf"))

ftr <- fish.traits.dat |>
		select(SPECIES, FeedingType)
        
point.range <- reef_fish_sti_glorys |>
		arrange(SPECIES) |>
		select(SPECIES, min.lat, max.lat, q25lat, q75lat,thetao_mean_TM) |>
        rename(TEMP.OPTIM=thetao_mean_TM)

pop.traits <- sp.focal |>
		group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
		reframe(MEAN.LAT=mean(LAT), MEAN.LON=mean(LON)) |>
		ungroup() |>
		select(ECOREGION,MEAN.LAT,MEAN.LON,SPECIES,SPECIES.ORIG) |>
        # split ecoregions suffix SPECIES; traits/STI are keyed by the clean name,
        # so join on SPECIES.ORIG (else Bassian/Hawaii get NA guild+range and drop)
        left_join(ftr, by=c("SPECIES.ORIG"="SPECIES")) |>
        left_join(point.range, by=c("SPECIES.ORIG"="SPECIES")) |>
        select(-SPECIES.ORIG) |>
        mutate(
        THERMAL.GUILD = ifelse(TEMP.OPTIM >= 23, "Tropical", "Temperate"),
        THERMAL.GUILD=factor(THERMAL.GUILD, levels=c("Tropical","Temperate")),
        RANGE = case_when(
            is.na(q25lat) | is.na(q75lat) ~ NA,
            # Both positive (Northern Hemisphere species)
            q25lat > 0 & q75lat > 0 ~ case_when(
            MEAN.LAT < q25lat ~ "Warm edge",
            MEAN.LAT > q75lat ~ "Cold edge",
            TRUE              ~ "Mid range"
        ),
            # Both negative (Southern Hemisphere species)
            q25lat < 0 & q75lat < 0 ~ case_when(
            MEAN.LAT > q75lat ~ "Warm edge",
            MEAN.LAT < q25lat ~ "Cold edge",
            TRUE              ~ "Mid range"
        ),
            # Trans-equatorial (q25lat < 0, q75lat > 0)
            q25lat < 0 & q75lat > 0 & abs(q25lat) < abs(q75lat) ~ case_when(
            MEAN.LAT < q25lat ~ "Warm edge",
            MEAN.LAT > q75lat ~ "Cold edge",
            TRUE              ~ "Mid range"
        ),
            q25lat < 0 & q75lat > 0 & abs(q25lat) >= abs(q75lat) ~ case_when(
            MEAN.LAT > q75lat ~ "Warm edge",
            MEAN.LAT < q25lat ~ "Cold edge",
            TRUE              ~ "Mid range"
        ),
      TRUE ~ NA),
      RANGE=factor(RANGE, levels=c("Warm edge","Mid range","Cold edge"))
    ) |>
    filter(!is.na(THERMAL.GUILD)&!is.na(RANGE))

temp.exposure.by.time <- sp.focal |>
    left_join(pop.traits, by=c("ECOREGION","SPECIES")) |>
	rename(mean=shifted.cumtemp.above.mean, q50=shifted.cumtemp.above.q50,
            q70=shifted.cumtemp.above.q70, q90=shifted.cumtemp.above.q90,
            q95=shifted.cumtemp.above.q95, q97.5=shifted.cumtemp.above.q97.5) |>
    mutate(EXPOSURE.QUALITATIVE=case_when(q90>0 ~ "Extreme exposure", TRUE ~ "Moderate exposure"),
            EXPOSURE.QUALITATIVE = factor(EXPOSURE.QUALITATIVE, levels=c("Moderate exposure", "Extreme exposure"))) |>        
    pivot_longer(cols=c(mean, q50, q70, q90, q95, q97.5), names_to="EXPOSURE.THRESHOLD", values_to="EXPOSURE") |>
    filter(!is.na(THERMAL.GUILD)&!is.na(RANGE))

temp.exposure.agg <- temp.exposure.by.time |>
    pivot_wider(names_from=EXPOSURE.THRESHOLD, values_from=EXPOSURE) |>
    group_by(ECOREGION,SPECIES,THERMAL.GUILD,RANGE) |>
    reframe(mean=sum(mean, na.rm=TRUE),
            q50=sum(q50, na.rm=TRUE),
            q70=sum(q70, na.rm=TRUE),
            q90=sum(q90, na.rm=TRUE),
            q95=sum(q95, na.rm=TRUE),
            q97.5=sum(q97.5, na.rm=TRUE))

trend.lat.df <- trend.focal |> filter(FlagAnalysis=="Trend_Lat"&Effect=="Lat.Trend") |>
    group_by(REALM,ECOREGION,SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    filter(!is.na(Estimate)) |>
    ungroup() |>
    left_join(pop.traits, by=c("ECOREGION","SPECIES")) |>
    left_join(temp.exposure.agg, by=c("ECOREGION","SPECIES","THERMAL.GUILD","RANGE")) |>
    relocate(c(MEAN.LAT,MEAN.LON), .after=ECOREGION) |>
    filter(!is.na(THERMAL.GUILD)&!is.na(RANGE))     

mean.trendlat.df <- trend.lat.df |>
    group_by(ECOREGION) |>
    reframe(
    expo=mean(mean),
    n = n()
    ) |>
    filter(n >= 5) |>
    mutate(Eco_ID=1:n()) |>
    arrange(ECOREGION)

warming_all <- temperature_all_ecoregions_all |>
    rename(warming_all = warming_rate,
           warming_all_se = warming_se,
           n_all = n_cells) |>
    select(ECOREGION, warming_all, warming_all_se, n_all)

focal.meow <- meow |>
    select(-REALM) |>
    arrange(ECOREGION)|>
		left_join(mean.trendlat.df,
				by=c("ECOREGION")) |>
		filter(ECOREGION %in% mean.trendlat.df$ECOREGION) |>
		st_make_valid() |>
    left_join(warming_all, by = "ECOREGION")

# Diagnostic: check the data before plotting
cat("focal.meow rows:", nrow(focal.meow), "\n")
cat("Geometry types:", paste(unique(st_geometry_type(focal.meow)), collapse = ", "), "\n")
cat("Any empty geometries:", sum(st_is_empty(focal.meow)), "\n")
cat("warming_coastal range:",
    paste(round(range(focal.meow$warming_coastal, na.rm = TRUE), 4), collapse = " - "), "\n")

map.data <- ggplot() +
		geom_sf(data = world_map, fill = "white", col = NA) +
    geom_sf(data = focal.meow, aes(fill = warming_all), col = "grey70",
		        linewidth = 0.2, alpha = 1) +
    geom_sf(data = world_map, fill = "white", col = NA) +
    scale_fill_taster_d("aperol_spritz_blue_angel_mimosa",
        cont = TRUE, n = 8, dir = -1, na.value = NA,
        name = expression("Warming rate ("*degree*"C yr"^-1*")"),
        # For expo (cumulative temperature) variant, use:
        # trans = "log10",
        # breaks = c(25000, 50000, 100000),
        # labels = expression(2.5%*%10^4, 5%*%10^4, 10^5),
        # name = "Cumulative temperature (\u00b0C\u00b7d)",
        guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    theme_void() +
		theme(
				legend.position="bottom",
				legend.background = element_rect(fill = "white"),
				legend.key = element_rect(fill = "black"),
				legend.text = element_text(color = "black"),
				legend.title = element_text(color = "black"),
				panel.background = element_rect(fill = "#c2e7f7", color = NA),
        legend.key.size = unit(1.6, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.key.width = unit(2.5, "cm"),
        legend.margin = margin(5, 10, 5, 10)
        ) +
		coord_sf(expand = FALSE)

plot(map.data)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/Fig2_map_%s.pdf", kind_suffix),
#        map.data, width = 18, height = 9)

unique(focal.meow$ECOREGION)

#### inset ####
dir.df <- trend.lat.df |>
  mutate(
    direction = case_when(
      !is.na(P.value) & P.value < 0.05 & Estimate > 0 ~ "Sig poleward",
      !is.na(P.value) & P.value < 0.05 & Estimate < 0 ~ "Sig equatorward",
      Estimate > 0 ~ "NS poleward",
      Estimate < 0 ~ "NS equatorward",
      TRUE ~ NA_character_),                      # Estimate == 0 / NA
    direction = factor(direction, levels = c("Sig poleward", "NS poleward",
                                             "NS equatorward", "Sig equatorward")),
    THERMAL.GUILD = factor(THERMAL.GUILD, levels = c("Tropical", "Temperate")),
    RANGE = factor(RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))) |>
  filter(!is.na(direction))

# Fig 2 legend N (inset): populations with an assignable realm (guild) AND
# edge-range position. Expect 2,463.
cat(sprintf("[Fig 2 legend] inset: %d populations (realm + range assignable)\n", nrow(dir.df)))

multi.mod <- nnet::multinom(direction ~ THERMAL.GUILD * RANGE, data = dir.df)
print(car::Anova(multi.mod, type = 3))

dir.props <- dir.df |>
    count(THERMAL.GUILD, RANGE, direction) |>
    group_by(THERMAL.GUILD, RANGE) |>
    mutate(pct = n / sum(n) * 100) |>
    ungroup()

# per-bar population totals, placed above each stacked bar
dir.bar.n <- dir.df |>
    count(THERMAL.GUILD, RANGE, name = "n")

# Direction × significance palette — fully CockR-derived. Sig extremes reuse
# the rose colours (dir_cols[["poleward"]] and dir_cols[["equatorward"]]) so
# the inset reads as part of Fig 2. NS categories are neutral greys with a
# subtle cool/warm tint (Blue Lagoon vs Espresso Martini) — visually distinct
# from the rose laterals (light blue + pink) to avoid reader confusion.
cockr_cont  <- taster_palettes_continuous()
blueangel   <- cockr_cont[["Blue Angel_continuous"]]
aperol      <- cockr_cont[["Aperol Spritz_continuous"]]
espresso    <- cockr_cont[["Espresso Martini_continuous"]]

# Both NS categories use Espresso Martini (warm neutral grey) at different
# lightness — no blue/pink hue, no collision with the rose laterals.
dir.pal <- c(
  "Sig poleward"    = blueangel[156],   # = dir_cols[["poleward"]]
  "NS poleward"     = espresso[25],     # pale warm grey
  "NS equatorward"  = espresso[80],     # darker warm grey
  "Sig equatorward" = aperol[103])      # = dir_cols[["equatorward"]]

p.lat.dir <- ggplot(dir.props,
    aes(x = RANGE, y = pct, fill = direction)) +
  geom_col(width = 0.75, colour = NA) +
  geom_text(data = dir.bar.n, aes(x = RANGE, y = 102, label = n),
            inherit.aes = FALSE, size = 2.8, vjust = 0) +
  facet_wrap(~ THERMAL.GUILD) +
  scale_fill_manual(values = dir.pal) +
  coord_cartesian(ylim = c(0, 104), clip = "off") +
  labs(x = "", y = "% of populations", fill = "Direction") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(size = 10),
		axis.title.x = element_text(size = 10))

plot(p.lat.dir)
# ggsave(p.lat.dir, filename = sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/Fig2_inset_%s.pdf", kind_suffix),
  # width = 6, height = 4)

# The density version of this inset (and of every panel) is produced simply by setting
# DATA_KIND <- "density" at the top and re-running — outputs carry the "_density" suffix.
# (The old manual density block that reused the abundance `dir.props` is removed.)
