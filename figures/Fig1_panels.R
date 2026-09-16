# ============================================================================
# Fig 1 — illustrative panels for a single example species (Choerodon fasciatus)
#   a: abundance optimum identified along coordinates (modskurt LAT & LON curves)
#      + a year-to-year trajectory map of the relocated optimum (supporting).
#   b: LATITUDE of the optimum vs YEAR (latitudinal trend through time).
#   c: TEMPERATURE vs YEAR — observed at shifted peak sites vs no-shift counterfactual.
#
# PRODUCES:
#   MAIN: Fig 1 (panels a-c above) — the worked example on one species.
#   EXTENDED DATA: none from this script. ED Fig 1 (the seven-step analytical
#     workflow + the modskurt candidate-curve schematic) is a separate methods
#     figure, not generated here.
# ============================================================================

# ---- paths: run with the working directory set to the project root (where
#      config.R lives). config.R defines dir_data/dir_out/dir_results/dir_figures
#      and input_file(); nothing below hard-codes a path.
source("config.R")

require(dplyr)
require(tidyr)
require(sf)
require(terra)
require(tidyterra)
require(ggplot2)
require(patchwork)
require(CockR)
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturalearthhires)
require(ggOceanMaps)
require(modskurt1)
library(cmdstanr)
library(posterior)

# plot_sites() override (basemap helper for element 3/3 of panel a)
source(file.path(dir_figures, "plot_sites.R"))

# modskurt Stan model
stan_path <- system.file("stan", "modskurt1.stan", package = "modskurt1")
mod <- cmdstanr::cmdstan_model(stan_path)

# ---- inputs (resolved via input_file: output/ if regenerated, else data/) ----
load(input_file("sp.df.RData"))                    # raw survey data (abundance)
load(input_file("optim.abund.relocated.RData"))    # relocated optima + SELECTED_MODEL_LAT/LON

# ---- target species (illustrative) ---------------------------------------
fish.name      <- "Choerodon fasciatus"
target.ecoreg  <- "Central and Southern Great Barrier Reef"
year           <- 2018                 # year fitted for the a/b curves
focal.response <- "count"              # abundance = counts

sp.df <- sp.df |> mutate(focal.var = abund)

# helper: pull a numeric shape parameter (Mr, d, p) out of a model-name string
extract_param <- function(param_name, string) {
  pattern <- paste0(param_name, "([0-9]+\\.?[0-9]*)")
  match   <- regmatches(string, regexec(pattern, string))[[1]]
  if (length(match) > 1) as.numeric(match[2]) else 1  # default = 1
}

# common minimal theme for the abundance curves
theme_curve <- theme(
  legend.position  = "none",
  panel.border     = element_rect(colour = "grey60", fill = NA),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.y = element_text(size = 12, colour = "black"),
  axis.title.x = element_text(size = 12, colour = "black"))

# output folder (its own leaf under results/)
outdir <- file.path(dir_results, "Fig1_abund_select")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# data for the fitted year
df2 <- sp.df |> filter(ECOREGION == target.ecoreg, SPECIES == fish.name, YEAR == year)

# ==========================================================================
# Panel a — element 1/3: modskurt abundance curve along LATITUDE
# ==========================================================================
mod.sel.lat <- (optim.abund.relocated |>
  filter(ECOREGION == target.ecoreg, SPECIES == fish.name, YEAR == year))$SELECTED_MODEL_LAT
parts.lat      <- strsplit(mod.sel.lat, "_")[[1]]
params.str.lat <- parts.lat[1]
model.lat      <- parts.lat[2]
Mr.lat <- extract_param("Mr", params.str.lat)
d.lat  <- extract_param("d",  params.str.lat)
p.lat  <- extract_param("p",  params.str.lat)

spec.lat <- modskurt_spec(
  x = df2$LAT, y = df2$focal.var,
  response_type = focal.response, distribution = model.lat,
  use_r = Mr.lat, use_d = d.lat, use_p = p.lat)
fit.lat <- modskurt_fit(spec.lat, mod = mod, chains = 4)

p.tlat <- modskurt_plot(fit = fit.lat, spec.lat, type = "mean", interval = "none",
  quantiles = NULL, apply_zi = FALSE, show_draws = FALSE, combine = FALSE,
  xlab = "Latitude", plot = FALSE)

p.trend.lat <- ggplot(p.tlat$trend.df, aes(x = Latitude, y = mu_mean)) +
  geom_line(linewidth = 1.2, col = "red") +
  geom_point(data = p.tlat$point.df, aes(y = y), color = "black") +
  scale_y_continuous(name = "Abundance") + theme_curve

print(p.trend.lat)                                          # panel a, element 1 (latitude curve)
#ggsave(p.trend.lat, file = file.path(outdir, "Choerodon_trend.lat.pdf"), width = 4, height = 3, dpi = 300)

# ==========================================================================
# Panel a — element 2/3: modskurt abundance curve along LONGITUDE
# ==========================================================================
mod.sel.lon <- (optim.abund.relocated |>
  filter(ECOREGION == target.ecoreg, SPECIES == fish.name, YEAR == year))$SELECTED_MODEL_LON
parts.lon      <- strsplit(mod.sel.lon, "_")[[1]]
params.str.lon <- parts.lon[1]
model.lon      <- parts.lon[2]
Mr.lon <- extract_param("Mr", params.str.lon)
d.lon  <- extract_param("d",  params.str.lon)
p.lon  <- extract_param("p",  params.str.lon)

spec.lon <- modskurt_spec(
  x = df2$LON, y = df2$focal.var,
  response_type = focal.response, distribution = model.lon,
  use_r = Mr.lon, use_d = d.lon, use_p = p.lon)
fit.lon <- modskurt_fit(spec.lon, mod = mod, chains = 4)

p.tlon <- modskurt_plot(fit = fit.lon, spec.lon, type = "mean", interval = "none",
  quantiles = NULL, apply_zi = FALSE, show_draws = FALSE, combine = FALSE,
  xlab = "Longitude", plot = FALSE)

p.trend.lon <- ggplot(p.tlon$trend.df, aes(x = Longitude, y = mu_mean)) +
  geom_line(linewidth = 1.2, col = "red") +
  geom_point(data = p.tlon$point.df, aes(y = y), color = "black") +
  scale_y_continuous(name = "Abundance") + theme_curve

print(p.trend.lon)                                          # panel a, element 2 (longitude curve)
#ggsave(p.trend.lon, file = file.path(outdir, "Choerodon_trend.lon.pdf"), width = 4, height = 3, dpi = 300)

# ==========================================================================
# Panel a — element 3/3: trajectory of the relocated optimum across years (map)
# ==========================================================================
world <- ne_countries(scale = 10, returnclass = "sf")
world <- st_transform(world, crs = st_crs(4326))

relocated.df <- optim.abund.relocated |>
  filter(ECOREGION == target.ecoreg, SPECIES == fish.name)
yr <- relocated.df |> distinct(YEAR) |>
  filter(YEAR %in% c(2000, 2004, 2006, 2010, 2015, 2018, 2020))
plot.relocated.df <- relocated.df |> filter(YEAR %in% yr$YEAR)

my.pal <- "aperol_spritz_blue_angel_blue_lagoon"

plot.trajectories <- plot_sites(
  plot.dat = plot.relocated.df, plot.map = world, title.col = "ECOREGION",
  basemap = "ggom", point.size = 0, alpha = 0.8, shape = 21, col = "black",
  pad.km = 100, facet = FALSE,
  ggom.bathy = TRUE, ggom.bathy.style = "rcb", ggom.args = list()) +
  geom_point(data = plot.relocated.df,
    aes(x = LON, y = LAT, color = factor(YEAR), fill = factor(YEAR)), size = 5) +
  geom_segment(data = plot.relocated.df,
    aes(x = LON, y = LAT, xend = lead(LON), yend = lead(LAT)),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "grey60") +
  scale_color_taster_d(palette = my.pal, n = 8, dir = -1) +
  scale_fill_taster_d(palette = my.pal, n = 8, dir = -1)

print(plot.trajectories)                                    # panel a, element 3 (trajectory map)
#ggsave(plot.trajectories, file = file.path(outdir, "Choerodon_trajectories.pdf"), width = 6, height = 4, dpi = 300)

# ==========================================================================
# Panels b & c — temporal trends. Need observed vs counterfactual temperature:
#   reference = long-term first-site panel (1993-2021), stacked with observed.
# ==========================================================================
load(input_file("sp.optim.abund.shift.RData"))                # observed peak series
load(input_file("sp.optim.abund.firstsite.longterm.RData"))   # long-term (no-shift) reference

obs.series <- sp.optim.abund.shift |>
  rename(LAT = SHIFTED.LAT, ECOREGION = ORIG.ECOREGION) |>
  filter(ECOREGION == target.ecoreg, SPECIES == fish.name) |>
  transmute(ECOREGION, SPECIES, YEAR, LAT, Counterf = "shifted.temp", Temp = shifted.mean.temp)
ref.series <- sp.optim.abund.firstsite.longterm |>
  rename(LAT = SHIFTED.LAT) |>
  filter(ECOREGION == target.ecoreg, SPECIES == fish.name) |>
  transmute(ECOREGION, SPECIES, YEAR, LAT, Counterf = "initial.temp", Temp = shifted.mean.temp)
shifted.df <- bind_rows(obs.series, ref.series) |> drop_na(Temp)

mycols <- taster_palettes_discrete()[["alexander_blue_lagoon_negroni"]][c(1, 4)]
theme_ts <- theme_bw() + theme(
  panel.border = element_rect(colour = "grey60", fill = NA),
  panel.background = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.title.y = element_text(size = 12, colour = "black"),
  axis.title.x = element_text(size = 12, colour = "black"))

# Panel b — LATITUDE of the optimum vs YEAR (observed shifted series)
# (named p.lat.year to avoid clashing with the modskurt shape parameter p.lat above)
p.lat.year <- ggplot(shifted.df |> filter(YEAR >= 1995, YEAR < 2021, Counterf == "shifted.temp"),
                aes(x = YEAR, y = LAT)) +
  geom_point(shape = 21, color = mycols[2], fill = mycols[2]) +
  stat_smooth(method = "lm", se = TRUE, color = mycols[2], fill = mycols[2], alpha = 0.2) +
  scale_x_continuous(name = NULL, breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(name = "Latitude") + theme_ts

plot(p.lat.year)                                           # panel b on screen
#ggsave(p.lat.year, file = file.path(outdir, "Choerodon_lat.pdf"), width = 5, height = 3, dpi = 300)

# Panel c — TEMPERATURE vs YEAR: observed (shifted) vs no-shift counterfactual
p.cf <- ggplot(shifted.df |> filter(YEAR >= 1995, YEAR < 2021), aes(x = YEAR, y = Temp, group = Counterf)) +
  geom_point(aes(col = Counterf, fill = Counterf), shape = 21) +
  stat_smooth(aes(color = Counterf, fill = Counterf), method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = mycols) + scale_fill_manual(values = mycols) +
  scale_x_continuous(name = NULL, breaks = c(1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(name = "Temperature (°C)") + theme_ts
plot(p.cf)                                                 # panel c on screen
#ggsave(p.cf, file = file.path(outdir, "Choerodon_cf.pdf"), width = 5, height = 3, dpi = 300)
