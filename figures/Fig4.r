# ==========================================================================
# Figure 4 — counterfactual response (4a) + best-path seascape use (4b)
#
# PRODUCES (main + several Extended Data panels, all from one run):
#   MAIN
#     p_rg_all  -> Fig 4a: counterfactual Response by thermal guild x range.
#     p_track   -> Fig 4b: was a better-buffering path reachable, and taken,
#                  by guild x range (from the DP best-path search).
#   EXTENDED DATA — variants/re-aggregations of Fig 4 (same run unless noted; ED numbers per current draft)
#    p_ecoreg        -> ED Fig 4: Fig 4a (counterfactual) by ecoregion.
#    p_ecoreg_track  -> ED Fig 6: Fig 4b (best-path) by ecoregion.
#    SENSITIVITY ANALYSES:
#    re-run p_rg_all / p_track under the alternative SCENARIO settings:
#    p_rg_all -> ED Fig 5: counterfactual sensitivity  (a) abund/previoussite,
#                 (b) density/firstsite, (c) density/previoussite.
#    p_track  -> ED Fig 7: best-path sensitivity (same three settings as above).
#   ALTERNATIVE EXPLANATIONS (ED Fig 8) — a SEPARATE question, not a Fig-4 variant;
#   kept in this script only because panel a reuses the counterfactual Response built here:
#     p_sti_tol       -> ED Fig 8a: STI thermal-tolerance falsification (reuses
#                        trends_df$Response; + STI table).
#     p_hmax/p_slopes -> ED Fig 8b: depth-emigration proxy (self-contained; loads
#                        modskurt.optim.<kind>.unimodal, independent of Fig 4).
#
# SWITCHES: SCENARIO (data kind abund/density x reference first/previous site x
#   short/long term) at the top; USE_OPTIMUM (DP global optimum vs greedy twin)
#   in the Fig 4b block. Every ggsave() is commented — uncomment + set a folder.
# ==========================================================================

require(dplyr)
require(tidyr)
require(forcats)
require(ggplot2)
require(ggpmisc)
require(patchwork)
require(glmmTMB)
require(ggeffects)
library(CockR)

# 1. SCENARIO SELECTION (set before loading so data_kind drives which files load)
# -------------------------------------------------------------------------------
SCENARIO <- "scenario.abund.leadtime.firstsite.longterm"
#SCENARIO <- "scenario.density.leadtime.previoussite.longterm"

# Options:
#   scenario.abund.leadtime.firstsite              (short, firstsite, abundance)
#   scenario.abund.leadtime.prevsite               (short, previoussite, abundance)
#   scenario.abund.leadtime.firstsite.longterm      (long, firstsite, abundance)
#   scenario.abund.leadtime.previoussite.longterm   (long, previoussite, abundance)
#   scenario.density.leadtime.firstsite             (short, firstsite, density)
#   scenario.density.leadtime.previoussite          (short, previoussite, density)
#   scenario.density.leadtime.firstsite.longterm    (long, firstsite, density)
#   scenario.density.leadtime.previoussite.longterm (long, previoussite, density)

# Derive flags + companion-object names from SCENARIO string.
is_longterm <- grepl("longterm$",       SCENARIO)
is_prevsite <- grepl("prev(ious)?site", SCENARIO)
is_density  <- grepl("\\.density\\.",   SCENARIO)
data_kind   <- if (is_density) "density" else "abund"
ref_tag <- if (is_prevsite) {
  if (is_longterm) "previoussite" else "prevsite"
} else "firstsite"

# Output-name key: inherits data kind + reference site + time window from SCENARIO,
# so every ggsave() below is uniquely named per run (mirrors the null-model script).
fig_key <- sprintf("%s_%s_%s", data_kind,
                   if (is_prevsite) "previoussite" else "firstsite",
                   if (is_longterm) "longterm" else "shortterm")
realized_panel_name <- paste0("sp.optim.", data_kind, ".shift")
realized_trend_name <- paste0("optim.",    data_kind, ".trend")
if (is_longterm) {
  ref_panel_name <- paste0("sp.optim.", data_kind, ".", ref_tag, ".longterm")
  ref_trend_name <- paste0("optim.",    data_kind, ".", ref_tag, ".trend.longterm")
  leadtime_name  <- NULL
} else {
  ref_panel_name <- NULL
  ref_trend_name <- NULL
  leadtime_name  <- paste0("sp.optim.", data_kind, ".shift.leadtime")
}

# 2. LOAD DATA -------------------------------------------------------------------------------
setwd("~/Lavori/MPA_timeseries/Modskurt")

# ── Shared files (always needed regardless of data_kind) ------------------
load("reef_fish_sti_glorys.RData")
load("fish.traits.dat.RData")
load("temperature_all_sites.RData")
load("sp.optim.abund.shift.RData")     # unimodal_filter_helper needs it for ECOREGION_ID mapping

# ── Core files — switch on data_kind ------------------------------------------
if (data_kind == "abund") {
  # Realised panels + leadtime (sp.optim.abund.shift already loaded above)
  load("sp.optim.abund.shift.leadtime.RData")
  load("optim.abund.trend.RData")
  # Long-term reference panels
  load("sp.optim.abund.firstsite.longterm.RData")
  load("sp.optim.abund.previoussite.longterm.RData")
  load("optim.abund.firstsite.trend.longterm.RData")
  load("optim.abund.previoussite.trend.longterm.RData")
} else {
  # Realised panels + leadtime (density equivalents)
  load("sp.optim.density.shift.RData")
  load("sp.optim.density.shift.leadtime.RData")
  load("optim.density.trend.RData")
  # Long-term reference panels (density equivalents)
  load("sp.optim.density.firstsite.longterm.RData")
  load("sp.optim.density.previoussite.longterm.RData")
  load("optim.density.firstsite.trend.longterm.RData")
  load("optim.density.previoussite.trend.longterm.RData")
}

# ---- Scenario classification file ---------------------------------------------------------------
load(paste0(SCENARIO, ".RData"))

# 3. SCENARIO ALIASES (derived from SCENARIO + flags above)
# Use rename(any_of(...)) so only columns actually present get renamed —
# longterm panels often already use LAT/LON/ECOREGION instead of SHIFTED.*/ORIG.*.
rename_map <- c(LAT = "SHIFTED.LAT", LON = "SHIFTED.LON", ECOREGION = "ORIG.ECOREGION")

# Realized panel — descriptive stats (pop.traits, eco centroids, ts lengths).
sp.focal.obs <- get(realized_panel_name) |> rename(any_of(rename_map))

# Lat trends — always from the realized trend output (peak latitude doesn't
# exist in the long-term ref panel, which is anchored at first/previous site).
optim.trend.lat <- get(realized_trend_name)

# Scenario object itself
scenario.focal <- get(SCENARIO)

# Env-trend back-transform inputs + leadtime / ref panel
if (is_longterm) {
  sp.focal.env      <- get(ref_panel_name) |> rename(any_of(rename_map))
  sp.focal.leadtime <- sp.focal.env
  optim.trend.env   <- get(ref_trend_name)
} else {
  sp.focal.env      <- sp.focal.obs
  sp.focal.leadtime <- get(leadtime_name) |> rename(any_of(rename_map))
  optim.trend.env   <- get(realized_trend_name)
}

# Diff-from-ref column name (used only in the short-term sd_env_stacked branch;
# absent in long-term scenarios, where sd_env_stacked auto-branches on its
# absence).
contrast.var <- if (is_prevsite) "diff.PreviousSite.temp.mean" else "diff.FirstSite.temp.mean"
  
# Population traits: Thermal guild + Range position
sti <- reef_fish_sti_glorys |>
  select(SPECIES, q25lat, q75lat, thetao_mean_TM) |>
  rename(TEMP.OPTIM = thetao_mean_TM)

pop.traits <- sp.focal.obs |>
  group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
  reframe(
    MEAN.LAT = mean(LAT, na.rm = TRUE),
    SD.cumtemp=sd(shifted.cumtemp.above.mean, na.rm=TRUE),
    CV.cumtemp=SD.cumtemp/mean(shifted.cumtemp.above.mean, na.rm=TRUE)
    ) |>
  # split ecoregions suffix SPECIES; STI is keyed by the clean name, so join on
  # SPECIES.ORIG (else Bassian/Hawaii get NA guild+range and drop from facets)
  left_join(sti, by = c("SPECIES.ORIG" = "SPECIES")) |>
  select(-SPECIES.ORIG) |>
  mutate(
    THERMAL.GUILD = ifelse(TEMP.OPTIM >= 23, "Tropical", "Temperate"),
    THERMAL.GUILD = factor(THERMAL.GUILD, levels = c("Tropical", "Temperate")),
    RANGE = case_when(
      is.na(q25lat) | is.na(q75lat) ~ NA_character_,
      q25lat > 0 & q75lat > 0 ~ case_when(
        MEAN.LAT < q25lat ~ "Warm edge",
        MEAN.LAT > q75lat ~ "Cold edge",
        TRUE              ~ "Mid range"),
      q25lat < 0 & q75lat < 0 ~ case_when(
        MEAN.LAT > q75lat ~ "Warm edge",
        MEAN.LAT < q25lat ~ "Cold edge",
        TRUE              ~ "Mid range"),
      q25lat < 0 & q75lat > 0 & abs(q25lat) < abs(q75lat) ~ case_when(
        MEAN.LAT < q25lat ~ "Warm edge",
        MEAN.LAT > q75lat ~ "Cold edge",
        TRUE              ~ "Mid range"),
      q25lat < 0 & q75lat > 0 & abs(q25lat) >= abs(q75lat) ~ case_when(
        MEAN.LAT > q75lat ~ "Warm edge",
        MEAN.LAT < q25lat ~ "Cold edge",
        TRUE              ~ "Mid range"),
      TRUE ~ NA_character_),
    RANGE = factor(RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))) |>
  filter(!is.na(THERMAL.GUILD) & !is.na(RANGE))

# back-transform estimates from standardized models
sd_year_lookup <- sp.focal.env |>
  group_by(ECOREGION, SPECIES) |>
  reframe(sd_year = sd((YEAR - min(YEAR)) + 1))

# sd of ENV.VAR as used in env_trends fit
# Ratinoale for back-transforming regression coefficnets originated from a standardized regression env over year:.afterFor a regression Y = β X + ε, if both are standardized:
# - Y* = (Y - mean_Y) / sd_Y                                                                                                          
# - X* = (X - mean_X) / sd_X                                                                                                                                
# The standardized slope β* relates to the raw slope β by:                                                                          
# β_raw = β_std × (sd_Y / sd_X)
# Here sd_Y = sd_env, sd_X = sd_year. Hence the back-transform is Estimate × sd_env / sd_year.  
# Set env variable + contrast once; used by both lookups and (implicitly) by
# mod.parms.env / mod.parms.slope below. Change here to switch env variable.
resp.env <- "shifted.cumtemp.above.mean"
# contrast.var is set by the SCENARIO switch at the top of this script

sd_env_obs_lookup <- sp.focal.env |>
  group_by(ECOREGION, SPECIES) |>
  reframe(sd_env_obs = sd(.data[[resp.env]], na.rm = TRUE))

# sd of ENV.VAR as used in counterfactual fit (stacked ref + obs vector).
# Two paths:
#   SHORT-TERM data: leadtime file carries the diff variable `contrast.var` but
#     NOT resp.env, so we join resp.env in from sp.focal and back-derive
#     ref.var = obs - diff (original logic).
#   LONG-TERM data: leadtime panel IS the ref panel and already carries
#     resp.env directly (no diff variable, ref.env.var = obs.env.var per the
#     longterm build script). We stack ref values from sp.focal.leadtime with
#     obs values from sp.focal and take the SD.
if (contrast.var %in% names(sp.focal.leadtime)) {
  # SHORT-TERM: leadtime carries diff variable; bring resp.env in from
  # sp.focal.env (which in short scenarios equals sp.focal.obs) and back-derive
  # ref.var = obs - diff.
  sd_env_stacked_lookup <- sp.focal.leadtime |>
    left_join(
      sp.focal.env |> select(ECOREGION, SPECIES, YEAR, LAT, LON, all_of(resp.env)),
      by = c("ECOREGION", "SPECIES", "YEAR", "LAT", "LON")) |>
    mutate(ref.var = .data[[resp.env]] - .data[[contrast.var]]) |>
    filter(!is.na(ref.var)) |>
    group_by(ECOREGION, SPECIES) |>
    reframe(sd_env_stacked = sd(c(ref.var, .data[[resp.env]]), na.rm = TRUE))
} else {
  # LONG-TERM: leadtime IS the ref panel (1993→last). Stack with the realized
  # obs (sampling-year peaks) from sp.focal.obs.
  ref_vals <- sp.focal.leadtime |>
    select(ECOREGION, SPECIES, v = all_of(resp.env))
  obs_vals <- sp.focal.obs |>
    select(ECOREGION, SPECIES, v = all_of(resp.env))
  sd_env_stacked_lookup <- bind_rows(ref_vals, obs_vals) |>
    filter(!is.na(v)) |>
    group_by(ECOREGION, SPECIES) |>
    reframe(sd_env_stacked = sd(v, na.rm = TRUE))
}

#unique(optim.trend.lat$mod.parms)
#mod.parm.lat <- "Lat_none_ou_NULL_gaussian"  #"Lat_none_no_autocor_NULL_gaussian"
lat_trends <- optim.trend.lat |>
    filter(
        FlagAnalysis%in%c("Trend_Lat")&
        Effect=="Lat.Trend"#&mod.parms%in%c(mod.parm.lat)
    ) |>
    group_by(ECOREGION, SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
    mutate(
      is_std = grepl("Standardize", mod.parms),
      Estimate = ifelse(is_std, Estimate / sd_year, Estimate), # see rationale for back-transformation of regression coefficients at line 69
      SE = ifelse(is_std, SE / sd_year, SE),
      poleward_km = Estimate * 111,
      direction = case_when(
      !is.na(P.value) & P.value < 0.05 & Estimate > 0 ~ "Sig_poleward",
      !is.na(P.value) & P.value < 0.05 & Estimate < 0 ~ "Sig_equatorward",
      TRUE ~ "NS"),REALM,
      direction = factor(direction, levels = c("NS", "Sig_poleward", "Sig_equatorward"))
    ) |>
    rename(Trend.Lat = Estimate, SE.Lat=SE, P.value.Lat = P.value, Lat.mod.parms=mod.parms) |>
    select(-c(FEEDING.TYPE,TEMP.OPTIM,THERMAL.GUILD,FlagAnalysis, Effect, t.value, AIC, Warning))

#unique(optim.trend.env$mod.parms)
#mod.parms.env <- "Env_shifted.cumtemp.above.mean_ou_Standardize_gaussian" #selected for Slope.Dev
#"Env_shifted.cumtemp.above.mean_no_autocor_NULL_gaussian"
#"Env_shifted.cumtemp.above.mean_ou_NULL_gaussian"
# "Env_shifted.sd.temp_ou_Standardize_gaussian"
env_trends <- optim.trend.env |>
    filter(
        FlagAnalysis%in%c("Trend_Environmental")&
        Effect=="Env.Trend"&
        grepl("cumtemp\\.above\\.mean", mod.parms)#&mod.parms%in%c(mod.parms.env)
    ) |>
    group_by(ECOREGION, SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
    left_join(sd_env_obs_lookup, by = c("ECOREGION", "SPECIES")) |>
    mutate(
      is_std = grepl("Standardize", mod.parms),
      Estimate = ifelse(is_std, Estimate * sd_env_obs / sd_year, Estimate), # see rationale for back-transformation of regression coefficients at line 69
      SE       = ifelse(is_std, SE * sd_env_obs / sd_year, SE)
    ) |>
    rename(Trend.Env = Estimate, SE.Env=SE, P.value.Env = P.value, Env.mod.parms=mod.parms) |>
    select(-c(FEEDING.TYPE,TEMP.OPTIM,THERMAL.GUILD,FlagAnalysis, Effect, t.value, AIC, Warning))
    
warming_ecoregion <- temperature_all_sites |>
  group_by(ECOREGION, SPECIES, SITE_ORDER) |>
  reframe(slope = coef(lm(annual_mean_temp ~ YEAR))[2]) |>
  group_by(ECOREGION, SPECIES) |>
  reframe(Trend.Env.Ecoregion.long = mean(slope))

fit_warming_ou <- function(df) {
  df <- df |>
    arrange(YEAR) |>
    mutate(YEAR_c    = (YEAR - min(YEAR)) + 1,
           times.dec = as.numeric(paste(YEAR_c, 1, sep = ".")),
           cum.times = lead(cumsum(lag(times.dec, default = 0)),
                            default = sum(times.dec)),
           yrs       = glmmTMB::numFactor(cum.times),
           group     = 1)
  m <- tryCatch(
    glmmTMB::glmmTMB(annual_mean_temp ~ YEAR_c + ou(yrs + 0 | group),
                     data = df),
    error   = function(e) NULL,
    warning = function(w) NULL)
  if (is.null(m) || !is.finite(glmmTMB::fixef(m)$cond["YEAR_c"]))
    return(unname(coef(lm(annual_mean_temp ~ YEAR_c, data = df))[2]))   # fallback
  unname(glmmTMB::fixef(m)$cond["YEAR_c"])
}

warming_sampled <- temperature_all_sites |>
  filter(FISH_SAMPLED) |>
  group_by(ECOREGION, SPECIES) |>
  filter(n() >= 5) |>
  group_modify(~ data.frame(Trend.Env.Sampled.long = fit_warming_ou(.x))) |>
  ungroup()

env_trend_long <- warming_ecoregion |>
  left_join(warming_sampled, by = c("ECOREGION", "SPECIES"))

#unique(scenario.abund.leadtime.firstsite$mod.parms)
#mod.parms.slope <- "Env_shifted.cumtemp.above.mean_ou_Standardize_gaussian"
#"Env_shifted.cumtemp.above.mean_no_autocor_NULL_gaussian"
#"Env_shifted.cumtemp.above.mean_ar1_NULL_gaussian" 
#"Env_shifted.sd.temp_ou_Standardize_gaussian" 
slope_dev <- scenario.focal |>
    filter(
        FlagAnalysis == "Trend_Environmental"&
        grepl("cumtemp\\.above\\.mean", mod.parms)#&mod.parms%in%c(mod.parms.slope)
    ) |>
    #group_by(ECOREGION, SPECIES, mod.parms) |>
    #Take the Slope.Dev row (row 4) which carries the Scenario/Response info
    #from the optimum_response_scenarios function
    filter(Effect%in%c("Slope.Dev")) |>
    #ungroup() |>
    rename(Slope.Dev=Estimate, SE.Dev=SE, P.Dev=P.value) |>
    group_by(ECOREGION, SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
    left_join(sd_env_stacked_lookup, by = c("ECOREGION", "SPECIES")) |>
    mutate(
      is_std = grepl("Standardize", mod.parms),
      Slope.Dev = ifelse(is_std, Slope.Dev * sd_env_stacked / sd_year, Slope.Dev), # see rationale for back-transformation of regression coefficients at line 69
      SE.Dev    = ifelse(is_std, SE.Dev    * sd_env_stacked / sd_year, SE.Dev)
    ) |>
    select(ECOREGION, SPECIES, Slope.Dev, SE.Dev, P.Dev, Scenario, Response, mod.parms) |>
    rename(SlopeDev.mod.parms = mod.parms)
  
slope_ref <- scenario.focal |>
    filter(
        FlagAnalysis == "Trend_Environmental"&
        grepl("cumtemp\\.above\\.mean", mod.parms)#&mod.parms%in%c(mod.parms.slope)
    ) |>
    #group_by(ECOREGION, SPECIES, mod.parms) |>
    #Take the Slope.Dev row (row 4) which carries the Scenario/Response info
    #from the optimum_response_scenarios function
    filter(Effect%in%c("Ref.Slope")) |>
    #ungroup() |>
    rename(Slope.Ref=Estimate, SE.Ref=SE, P.Ref=P.value) |>
    group_by(ECOREGION, SPECIES) |>
    slice_min(AIC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    left_join(sd_year_lookup, by = c("ECOREGION", "SPECIES")) |>
    left_join(sd_env_stacked_lookup, by = c("ECOREGION", "SPECIES")) |>
    mutate(
      is_std = grepl("Standardize", mod.parms),
      Slope.Ref = ifelse(is_std, Slope.Ref * sd_env_stacked / sd_year, Slope.Ref), # see rationale for back-transformation of regression coefficients at line 69
      SE.Ref    = ifelse(is_std, SE.Ref    * sd_env_stacked / sd_year, SE.Ref)
    ) |>
    select(ECOREGION, SPECIES, Slope.Ref, SE.Ref, P.Ref, mod.parms) |>
    rename(SlopeRef.mod.parms = mod.parms)
  
trends_df <- lat_trends |>
    left_join(slope_dev, by=c("ECOREGION","SPECIES")) |>
    left_join(slope_ref, by=c("ECOREGION", "SPECIES")) |>
    left_join(env_trends |> select(-c("sd_year", "is_std")), by = c("REALM","ECOREGION", "SPECIES","n.yrs")) |>
    left_join(env_trend_long, by = c("ECOREGION", "SPECIES")) |>
    left_join(pop.traits |> select(-c(q25lat, q75lat)), by=c("ECOREGION","SPECIES")) |>
    relocate(c(direction, Trend.Env,Trend.Env.Ecoregion.long, Trend.Env.Sampled.long), .after = Trend.Lat) |>
    relocate(c(Slope.Dev,Slope.Ref), .before="Trend.Lat") |>
    relocate(c(Response, Scenario), .after="P.value.Env") |>
    relocate(RANGE, .after="THERMAL.GUILD") |>
    mutate(
        Slope.Obs=Slope.Ref+Slope.Dev,
        M=1-Slope.Dev/Slope.Ref,
        Response=factor(Response, levels=c(
            "Magnification", "Mitigation", "Opposition",
            "Matching-Steady", "Matching-Warming", "Matching-Cooling",
            "Acceleration", "Deceleration", "Annihilation"
        )),
        Scenario=factor(Scenario, levels=c("Increasing","Matching","Buffering")) 
        ) |>
        relocate(c(Slope.Obs,M), .after="Slope.Ref") |>
        relocate(c(Lat.mod.parms,Env.mod.parms,SlopeDev.mod.parms,SlopeRef.mod.parms), .after="RANGE")

#save(trends_df, file="trends_df.RData")
#load("trends_df.RData")  

# PLOTS
resp_order <- c(
  "Magnification", "Mitigation", "Opposition",
  "Matching-Steady", "Matching-Warming", "Matching-Cooling",
  "Acceleration", "Deceleration", "Annihilation"
)

# Warm-to-cool gradient (orange -> yellow -> blue) from CockR palettes
thermal   <- taster_palettes_continuous()[["Thermal_continuous"]]
maitai    <- taster_palettes_continuous()[["Mai Tai_continuous"]]
blueangel <- taster_palettes_continuous()[["Blue Angel_continuous"]]

soft_orange <- colorRampPalette(c("white", thermal[200]))(16)
soft_yellow <- colorRampPalette(c("white", maitai[70]))(16)
soft_blue   <- colorRampPalette(c("white", blueangel[120]))(16)

response_palette <- c(
  "Magnification"    = soft_orange[16],
  "Mitigation"       = soft_orange[11],
  "Opposition"       = soft_orange[7],
  "Matching-Steady"  = soft_yellow[10],
  "Matching-Warming" = soft_yellow[5],
  "Matching-Cooling" = soft_blue[3],
  "Acceleration"     = soft_blue[7],
  "Deceleration"     = soft_blue[11],
  "Annihilation"     = soft_blue[16]
)

# helper proportions by ecoregion
compute_ecoregion_props <- function(df, resp_order, eco_lat) {
  # Drop populations with no assigned counterfactual Response (NA from the
  # left-join to slope_dev); otherwise they inflate the denominator below and
  # leave white gaps in the bars. Mirrors the !is.na(Response) filter used for
  # the Range x Guild panel.
  df <- df |> filter(!is.na(Response))

  # Count populations per ecoregion (for annotation)
  n_per_eco <- df |>
    group_by(ECOREGION) |>
    reframe(n_pop = n())

  # Proportion of each Response per ecoregion
  props <- df |>
    group_by(ECOREGION, Response) |>
    reframe(count = n()) |>
    group_by(ECOREGION) |>
    mutate(Perc = count / sum(count)) |>
    ungroup() |>
    select(-count)

  # Ensure all Response levels present for every ecoregion
  all_combos <- expand.grid(
    ECOREGION = unique(df$ECOREGION),
    Response  = resp_order,
    stringsAsFactors = FALSE
  )
  props <- all_combos |>
    left_join(props, by = c("ECOREGION", "Response")) |>
    mutate(Perc = replace_na(Perc, 0))

  # Join latitude for sorting and population counts
  props <- props |>
    left_join(eco_lat, by = "ECOREGION") |>
    left_join(n_per_eco, by = "ECOREGION") |>
    mutate(
      Response  = factor(Response, levels = resp_order),
      ECOREGION = fct_reorder(ECOREGION, eco_mean_lat))

  # Numbering: 1 = northernmost, N = southernmost (matches map convention).
  # Factor levels arranged so Eco_ID 1 sits at TOP of a horizontal y-axis →
  # reading top→bottom = N→S.
  eco_ids <- props |>
    distinct(ECOREGION, eco_mean_lat) |>
    arrange(desc(eco_mean_lat)) |>
    mutate(Eco_ID = row_number(),
           Eco_label = sprintf("%d. %s", Eco_ID, as.character(ECOREGION)))
  # Factor levels: southernmost first (bottom of y), northernmost last (top of y)
  props <- props |>
    left_join(eco_ids |> select(ECOREGION, Eco_ID, Eco_label), by = "ECOREGION") |>
    mutate(Eco_label = factor(Eco_label,
      levels = eco_ids$Eco_label[order(eco_ids$eco_mean_lat)]))

  props
}

# Lat centroids
eco_lat <- sp.focal.obs |>
  group_by(ECOREGION) |>
  reframe(eco_mean_lat = mean(LAT, na.rm = TRUE))

# Filter ecoregions with >= 5 populations for robustness
eco_sel <- trends_df |> count(ECOREGION) |> filter(n >= 5)
eco_focal <- trends_df |> filter(!is.na(Response)&ECOREGION %in% eco_sel$ECOREGION)
eco_props <- compute_ecoregion_props(eco_focal, resp_order, eco_lat)

# N-and-time-series labels per ecoregion
ts_length_pop <- sp.focal.obs |>
  group_by(ECOREGION, SPECIES) |>
  summarise(ts_years = n_distinct(YEAR), .groups = "drop")
ts_length_eco <- ts_length_pop |>
  group_by(ECOREGION) |>
  summarise(mean_ts = mean(ts_years), .groups = "drop")
  
n_labels_all <- eco_props |>
  distinct(Eco_label, ECOREGION, n_pop, eco_mean_lat) |>
  left_join(ts_length_eco, by = "ECOREGION") |>
  mutate(label = sprintf("n=%d, %.0fyr", n_pop, mean_ts))

# Scenarios actually OBSERVED (Perc>0 in any ecoregion), in resp_order — so the
# legend lists only observed scenarios, not all 9 possible.
obs_eco <- resp_order[resp_order %in% eco_props$Response[eco_props$Perc > 0]]

# CockR palette shortcuts (from response_palette): blue = buffered/counteract, red = amplify.
cock_blue <- unname(response_palette["Annihilation"])
cock_red  <- unname(response_palette["Magnification"])

# Fig4a by THERMAL_GUILD x RANGE
props_rg <- trends_df |>
    filter(!is.na(Response), !is.na(THERMAL.GUILD), !is.na(RANGE)) |>
    group_by(THERMAL.GUILD, RANGE, Response) |>
    reframe(count = n()) |>
    group_by(THERMAL.GUILD, RANGE) |>
    mutate(Perc = count / sum(count), n_pop = sum(count)) |>
    ungroup() |>
    mutate(Response = factor(Response, levels = resp_order))

all_combos_rg <- expand.grid(
  THERMAL.GUILD = c("Tropical","Temperate"),
  RANGE = c("Warm edge","Mid range","Cold edge"),
  Response = resp_order, stringsAsFactors = FALSE)

props_rg <- all_combos_rg |>
  left_join(props_rg, by = c("THERMAL.GUILD","RANGE","Response")) |>
  mutate(Perc = tidyr::replace_na(Perc, 0),
         THERMAL.GUILD = factor(THERMAL.GUILD, levels = c("Tropical","Temperate")),
         RANGE = factor(RANGE, levels = c("Warm edge","Mid range","Cold edge")),
         Response = factor(Response, levels = resp_order))

n_labels_rg <- props_rg |>
  distinct(THERMAL.GUILD, RANGE, n_pop) |>
  filter(!is.na(n_pop))

# Scenarios actually OBSERVED (Perc>0 in any guild x range cell), in resp_order —
# so the legend lists only observed scenarios, not all 9 possible.
obs_rg <- resp_order[resp_order %in% props_rg$Response[props_rg$Perc > 0]]

p_rg_all <- ggplot(props_rg, aes(x = RANGE, y = Perc, fill = Response)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  geom_text(data = n_labels_rg, inherit.aes = FALSE,
    aes(x = RANGE, y = 1.02, label = paste0(n_pop)),
    size = 5, vjust = 0, colour = "grey30") +
  facet_wrap(~ THERMAL.GUILD) +
  scale_fill_manual(values = response_palette, breaks = obs_rg, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_discrete(
    labels = c("Warm edge" = "Warm\nedge", "Mid range" = "Mid\nrange","Cold edge" = "Cold\nedge")) +
  labs(x = NULL, y = "Proportion of populations") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
		axis.title.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
		axis.text.x = element_text(size = 12, colour = "black"),
    strip.text = element_text(size=12, face = "bold"))

p_rg_all
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/Fig4a_counterfactual_range_guild_%s.pdf", fig_key),
#   p_rg_all, width = 7, height = 4)

# ==========================================================================================================
# (Fig 4b) WAS A BETTER-BUFFERING PATH REACHABLE, AND DID THE POPULATION TAKE IT? by guild x range.
#   INGESTS the DETERMINISTIC best-reachable-cooling-path search (bestpath_greedy_search.r / bestpath_optimum_search.r):
#   bestpath.<kind>.<ref>.RData. Two ONE-SIDED tests (alpha 0.05), relative to the counterfactual
#   (stay at the time-1 site) and to the realised path — both on the Slope.Dev ruler:
#     p_best_lt_ref = best path buffers below the counterfactual (best < ref)
#     p_best_lt_obs = best path buffers below the realised path   (best < obs)
#   Three classes (availability tested FIRST):
#     - No available mitigation path          : NOT(best < ref)                  (no reachable path beats staying)
#     - Available but not followed : best < ref  AND  best < obs      (a cooler path than realised existed)
#     - Available and followed : best < ref  AND  NOT(best < obs) (realised is as cool as the best)
#   Constant-exposure best paths (permanent sub-STI refuges) carry refugium=TRUE and
#   best_dev = -ref_slope (see bestpath header). guild x RANGE from trends_df (as Fig 4a).
# ==========================================================================================================
BEST_KIND <- data_kind                                    # data kind + reference site for the best-path load
BEST_REF  <- if (is_prevsite) "previoussite" else "firstsite"
USE_OPTIMUM <- TRUE                                        # TRUE = DP global optimum (primary); FALSE = greedy twin
bp_dir  <- if (USE_OPTIMUM) "bestpath_optimum" else "bestpath"
bp_stem <- if (USE_OPTIMUM) "bestpath.opt" else "bestpath"
bp_file <- path.expand(sprintf("~/Lavori/MPA_timeseries/Modskurt/%s/%s.%s.%s.RData",
                               bp_dir, bp_stem, BEST_KIND, BEST_REF))
if (!file.exists(bp_file))                                # fallback: flat Modskurt/ (no subfolder)
  bp_file <- path.expand(sprintf("~/Lavori/MPA_timeseries/Modskurt/%s.%s.%s.RData",
                                 bp_stem, BEST_KIND, BEST_REF))
BP_ALPHA <- 0.05
bp <- get(load(bp_file)[1]) |> filter(is.na(drop_reason), resolvable)
cat(sprintf("[bestpath] %s.%s: %d classifiable pops | %d refugia\n",
            BEST_KIND, BEST_REF, sum(is.finite(bp$p_best_lt_ref)), sum(bp$refugium, na.rm = TRUE)))
dv <- bp |>
  left_join(distinct(trends_df, ECOREGION, SPECIES, THERMAL.GUILD, RANGE),
            by = c("ECOREGION", "SPECIES")) |>
  rename(guild = THERMAL.GUILD) |>
  filter(is.finite(p_best_lt_ref), is.finite(p_best_lt_obs), !is.na(guild), !is.na(RANGE)) |>
  mutate(better_ref = p_best_lt_ref < BP_ALPHA,    # a cooler-than-counterfactual path was reachable
         better_obs = p_best_lt_obs < BP_ALPHA,    # a cooler path than realised existed
         outcome = dplyr::case_when(
           !better_ref ~ "No available mitigation path",
           better_obs  ~ "Available but not followed",
           TRUE        ~ "Available and followed"))
cls_lv <- c("Tropical · Warm edge", "Tropical · Mid range", "Tropical · Cold edge",
            "Temperate · Warm edge", "Temperate · Mid range", "Temperate · Cold edge")
dv$class <- factor(paste(dv$guild, dv$RANGE, sep = " · "), levels = cls_lv)
dv$RANGE <- factor(dv$RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))
dv$guild <- factor(dv$guild, levels = c("Tropical", "Temperate"))

diag12 <- dv |> group_by(class) |>
  summarise(n = dplyr::n(),
            pct_nobetter = round(100 * mean(outcome == "No available mitigation path")),
            pct_matched  = round(100 * mean(outcome == "Available and followed")),
            pct_notused  = round(100 * mean(outcome == "Available but not followed")), .groups = "drop")
cat("[Fig 4b] better reachable path & whether taken, by class:\n")
print(as.data.frame(diag12))

lv3  <- c("Available and followed", "Available but not followed", "No available mitigation path")
pal3 <- c("Available and followed" = cock_blue,
          "Available but not followed" = cock_red,
          "No available mitigation path"          = "grey65")
bars <- dv |> count(guild, RANGE, outcome) |> group_by(guild, RANGE) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
bars$outcome <- factor(bars$outcome, levels = lv3)
n_track <- dv |> count(guild, RANGE, name = "n_pop")     # populations per class (bar labels, as Fig 4a)
p_track <- ggplot(bars, aes(RANGE, pct, fill = outcome)) +
  geom_col(width = .7) +
  geom_text(data = n_track, inherit.aes = FALSE,
            aes(x = RANGE, y = 102, label = paste0(n_pop)), size = 5, vjust = 0) +
  facet_wrap(~ guild) +
  scale_fill_manual(values = pal3, name = NULL, labels = lv3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_discrete(
    labels = c("Warm edge" = "Warm\nedge", "Mid range" = "Mid\nrange","Cold edge" = "Cold\nedge")) +
  labs(x = NULL, y = "Proportion of populations", title = "") +
  theme_bw(base_size = 10) +
  scale_x_discrete(
    labels = c("Warm edge" = "Warm\nedge", "Mid range" = "Mid\nrange","Cold edge" = "Cold\nedge")) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 12, colour = "black"),
		axis.title.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
		axis.text.x = element_text(size = 12, colour = "black"),
    strip.text = element_text(size=12, face = "bold"))
    
plot(p_track)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/Fig4b_available_path_use_%s.pdf", fig_key),
#   p_track, width = 7, height = 4)

# ==========================================================================================================
# TEXT STATS — % of best paths that head POLEWARD, by category.
#   best path = the DETERMINISTIC best reachable cooling path (bestpath); poleward if its |lat| trend
#   > 0 (best_poleward). Reported for the two non-barrier categories, by REALM x RANGE and overall —
#   numbers to quote in the text, no figure.
# ==========================================================================================================
pole_df <- dv |>                                             # dv carries best_poleward (from bestpath)
  left_join(distinct(trends_df, ECOREGION, SPECIES, REALM), by = c("ECOREGION", "SPECIES")) |>
  filter(outcome %in% c("Available but not followed", "Available and followed"),
         !is.na(best_poleward))
pole_df$RANGE <- factor(pole_df$RANGE, levels = c("Warm edge", "Mid range", "Cold edge"))
pct_pole <- function(df, ...) df |> group_by(outcome, ...) |>
  summarise(n = dplyr::n(), pct_poleward = round(100 * mean(best_poleward)), .groups = "drop")
cat("\n[Fig 4b] % of best-buffering trajectories that are POLEWARD\n-- by REALM x RANGE --\n")
print(as.data.frame(pct_pole(pole_df, REALM, RANGE)))
cat("-- overall (per category) --\n")
print(as.data.frame(pct_pole(pole_df)))

# ==========================================================================================================
# RESULTS text — percentages quoted under 'Seascape constraints to mitigating shifts' (Fig 4b).
#   Ranges span the three range-edge classes within each realm; "followed" spans all six classes.
#   Expect (abund / firstsite): tropics no-path 60-73%, tropics not-followed 16-28%,
#   temperate no-path 25-37%, temperate not-followed 49-66%, followed 9-15%, poleward 86% / 72%.
# ==========================================================================================================
rng <- function(x) sprintf("%d-%d%%", min(x), max(x))
tro <- diag12 |> filter(grepl("^Tropical",  class))
tem <- diag12 |> filter(grepl("^Temperate", class))
po  <- pct_pole(pole_df)                                     # % poleward per category, overall
cat("\n[Seascape constraints] main-text percentages (Fig 4b):\n")
cat(sprintf("  Tropics   no available path        : %s (warm->cold %s)\n",
            rng(tro$pct_nobetter), paste(tro$pct_nobetter, collapse = "/")))
cat(sprintf("  Tropics   available, not followed  : %s\n", rng(tro$pct_notused)))
cat(sprintf("  Temperate no available path        : %s\n", rng(tem$pct_nobetter)))
cat(sprintf("  Temperate available, not followed  : %s\n", rng(tem$pct_notused)))
cat(sprintf("  Followed (all six classes)         : %s\n", rng(diag12$pct_matched)))
cat(sprintf("  Poleward: available-not-followed %d%%, followed %d%%\n",
            po$pct_poleward[po$outcome == "Available but not followed"],
            po$pct_poleward[po$outcome == "Available and followed"]))

# Fig4a by ecoregion vertical
p_ecoreg <- ggplot(eco_props, aes(x = Perc, y = Eco_label, fill = Response)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  geom_text(data = n_labels_all,
    aes(x = 1.05, y = Eco_label, label = label, fill = NULL),
    size = 2.8, colour = "grey40", hjust = 0.1) +
  scale_fill_manual(values = response_palette, breaks = obs_eco, drop = FALSE,
    name = "Counterfactual\nresponse") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.17)),
    breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent_format()) +
  labs(y = "Ecoregion", x = "Proportion of populations",
    title = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"))

print(p_ecoreg)

# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/FigS_Fig4a_counterfactual_by_ecoregion_%s.pdf", fig_key),
#   p_ecoreg, width = 9, height = 10)

# ==========================================================================================================
# SUPPLEMENTARY — Fig 4b BY ECOREGION (as Fig 4a's ecoregion panel), N->S order.
#   Lets ecoregions with clear physical/oceanographic barriers (e.g. Mediterranean — no sea
#   north of the continent; Bassian/Tasmania — no reef habitat further poleward) be read off:
#   a high "available but not followed" share there is NOT a physical barrier (the null found a
#   reachable path), pointing instead to little need to buffer.
# ==========================================================================================================
eco_lat_b <- sp.focal.obs |> group_by(ECOREGION) |> summarise(eco_mean_lat = mean(LAT, na.rm = TRUE), .groups = "drop")
# Filter ecoregions with >= 5 populations — consistency with the Fig 4a ecoregion supplement (eco_sel, >= 5).
# Applied only to the per-ecoregion figure; the guild x range p_track keeps the full dv.
eco_keep_b <- dv |> count(ECOREGION) |> filter(n >= 5) |> pull(ECOREGION)
dv_eco     <- dv |> filter(ECOREGION %in% eco_keep_b)
eco_ids_b <- eco_lat_b |> filter(ECOREGION %in% dv_eco$ECOREGION) |> arrange(desc(eco_mean_lat)) |>
  mutate(Eco_ID = row_number(), Eco_label = sprintf("%d. %s", Eco_ID, ECOREGION))
eco_track <- dv_eco |> count(ECOREGION, outcome, name = "n") |> group_by(ECOREGION) |>
  mutate(Perc = n / sum(n), n_pop = sum(n)) |> ungroup() |>
  left_join(eco_ids_b, by = "ECOREGION") |>
  mutate(Eco_label = factor(Eco_label, levels = eco_ids_b$Eco_label[order(eco_ids_b$eco_mean_lat)]),
         outcome   = factor(outcome, levels = lv3))
n_lab_b <- eco_track |> distinct(Eco_label, n_pop)
p_ecoreg_track <- ggplot(eco_track, aes(Perc, Eco_label, fill = outcome)) +
  geom_col(width = .8) +
  geom_text(data = n_lab_b, aes(x = 1.02, y = Eco_label, label = paste0("n=", n_pop)),
            hjust = 0, size = 2.4, inherit.aes = FALSE) +
  scale_fill_manual(values = pal3, name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, .12))) +
  labs(x = "Proportion of populations", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8))

plot(p_ecoreg_track)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/FigS_Fig4b_by_ecoregion_%s.pdf", fig_key), p_ecoreg_track, width = 7, height = 8)

# ==========================================================================================================
# SUPPLEMENTARY — STI adaptation / thermal-tolerance hypothesis (falsification)
#   Reuses trends_df built above (same scenario/data/thresholds). Prediction: populations
#   that MATCH or MAGNIFY warming (did not mitigate) are MORE warm-adapted (higher STI)
#   than MITIGATORS. Tested WITHIN realm (guild) — the POOLED contrast is confounded
#   (tropical = high-STI & magnify; temperate = low-STI & mitigate). This produces a
#   SEPARATE supplementary figure, NOT an inset into Fig 4.
# ==========================================================================================================
sti_test <- trends_df |>
  dplyr::filter(!is.na(Response), !is.na(THERMAL.GUILD)) |>
  dplyr::left_join(reef_fish_sti_glorys |> dplyr::transmute(SPECIES, STI = thetao_mean_TM),
                   by = c("SPECIES.ORIG" = "SPECIES")) |>
  dplyr::filter(is.finite(STI)) |>
  dplyr::mutate(grp = dplyr::case_when(
      Response == "Mitigation" ~ "Mitigate",
      Response %in% c("Magnification","Matching-Steady","Matching-Warming","Matching-Cooling") ~ "Match/Magnify",
      TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(grp)) |>
  dplyr::mutate(grp = factor(grp, levels = c("Match/Magnify","Mitigate")),
                THERMAL.GUILD = factor(THERMAL.GUILD, levels = c("Tropical","Temperate")),
                RANGE = factor(RANGE, levels = c("Warm edge","Mid range","Cold edge")))

# within each guild x range class (6) — Wilcoxon PRIMARY (robust to multimodality); Welch t = p_t
sti_grid <- sti_test |> dplyr::group_by(THERMAL.GUILD, RANGE) |>
  dplyr::summarise(n_mm = sum(grp=="Match/Magnify"), n_mit = sum(grp=="Mitigate"),
    med_mm = round(median(STI[grp=="Match/Magnify"]),1), med_mit = round(median(STI[grp=="Mitigate"]),1),
    dSTI = round(median(STI[grp=="Match/Magnify"]) - median(STI[grp=="Mitigate"]), 2),
    p_wilcox = signif(suppressWarnings(wilcox.test(STI[grp=="Match/Magnify"], STI[grp=="Mitigate"])$p.value), 3),
    p_t = signif(t.test(STI[grp=="Match/Magnify"], STI[grp=="Mitigate"])$p.value, 3), .groups = "drop")
cat("\n[STI tolerance] within guild x range (Wilcoxon PRIMARY; Welch t = p_t):\n")
print(as.data.frame(sti_grid), row.names = FALSE)
sti_pool <- function(s, lab){ a <- s$STI[s$grp=="Match/Magnify"]; b <- s$STI[s$grp=="Mitigate"]
  data.frame(panel=lab, dSTI=round(median(a)-median(b),2),
             p_wilcox=signif(suppressWarnings(wilcox.test(a,b)$p.value),3)) }
cat("pooled / by-realm context (confounded):\n")
print(dplyr::bind_rows(sti_pool(sti_test, "All (pooled)"),
  sti_pool(dplyr::filter(sti_test, THERMAL.GUILD=="Tropical"),  "Tropical"),
  sti_pool(dplyr::filter(sti_test, THERMAL.GUILD=="Temperate"), "Temperate")), row.names = FALSE)
cat("facet-controlled lm  STI ~ grp + guild*range  (grpMitigate > 0 => mitigators warmer-adapted, OPPOSITE to tolerance):\n")
print(round(summary(lm(STI ~ grp + THERMAL.GUILD*RANGE, data = sti_test))$coefficients["grpMitigate", , drop=FALSE], 4))
# write.csv(sti_grid, "~/Lavori/MPA_timeseries/Modskurt/Tables/TableS_STI_tolerance.csv", row.names = FALSE)

# supplementary figure: overlapping STI densities per guild x range class + group medians (CockR colours)
col_mit <- unname(taster_palettes_continuous()[["Blue Angel_continuous"]][156])     # Mitigate      (buffer/cool)
col_mm  <- unname(taster_palettes_continuous()[["Aperol Spritz_continuous"]][103])  # Match/Magnify (amplify/warm)
sti_meds <- sti_test |> dplyr::group_by(THERMAL.GUILD, RANGE, grp) |>
  dplyr::summarise(med = median(STI), .groups = "drop")
sti_ann <- sti_grid |>
  dplyr::transmute(THERMAL.GUILD, RANGE, lab = sprintf("ΔSTI = %+.2f °C\np = %.3g", dSTI, p_wilcox))
p_sti_tol <- ggplot(sti_test, aes(STI, fill = grp, colour = grp)) +
  geom_density(alpha = .40, linewidth = .55) +
  geom_vline(data = sti_meds, aes(xintercept = med, colour = grp),
             linetype = "dashed", linewidth = .55, show.legend = FALSE) +
  facet_wrap(THERMAL.GUILD ~ RANGE, ncol = 3, scales = "free") +
  geom_text(data = sti_ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3, size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c("Match/Magnify"=col_mm, "Mitigate"=col_mit),
                    labels = c("Match/Magnify"="Not mitigating", "Mitigate"="Mitigating"), name = NULL) +
  scale_colour_manual(values = c("Match/Magnify"=col_mm, "Mitigate"=col_mit),
                      labels = c("Match/Magnify"="Not mitigating", "Mitigate"="Mitigating"), name = NULL) +
  labs(x = "Species temperature index, STI (°C)", y = "Density",
       title = "") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "bottom")
print(p_sti_tol)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/FigS_STI_tolerance_%s.pdf", fig_key),
#        p_sti_tol, width = 9, height = 5.5, device = cairo_pdf)

# ==========================================================================================================
# SUPPLEMENTARY — depth-migration proxy: trend through time in peak modal abundance.
#   We cannot test depth shifts directly, but if fish descend below the -5 to -10 m survey
#   band over time, the peak modal abundance (Hmax = H.LAT, the modskurt peak height) at the
#   sampled peak sites should DECLINE. A glmmTMB model with population (species x ecoregion)
#   as a random intercept + slope estimates the population-average trend of log(Hmax) through
#   time; no general decline => little evidence for systematic deepening. Individual trends and
#   the % declining significantly (per-population OLS) show whether decliners are a minority.
#   (Same log(H.LAT) ~ year form as the pipeline's sp.lat.shift.analysis.R.)
# ==========================================================================================================
# Hmax = H.LAT (modskurt peak height) lives in the per-year optima, NOT the shift panel.
# NB: modskurt.optim.<kind>.unimodal.RData must be materialised locally to load — if it errors
# with "cannot open compressed file ... Invalid argument", it is online-only on this machine:
# right-click it in Google Drive -> "Available offline" (or run this section on the work computer).
opt_hmax <- get(load(sprintf("modskurt.optim.%s.unimodal.RData", data_kind))[1])
stopifnot("H.LAT" %in% names(opt_hmax))
hd <- opt_hmax |>
  dplyr::transmute(pop_id = interaction(ECOREGION, SPECIES, drop = TRUE), YEAR, Hmax = H.LAT) |>
  filter(is.finite(Hmax), Hmax > 0, is.finite(YEAR)) |>
  group_by(pop_id) |> filter(dplyr::n() >= 4) |> ungroup() |>
  mutate(YEAR_c = YEAR - min(YEAR), logHmax = log(Hmax))
# Random intercept + slope: each population has its OWN baseline AND its OWN estimated trend.
m_hmax <- glmmTMB(logHmax ~ YEAR_c + (YEAR_c | pop_id), data = hd)
cat(sprintf("[Fig S] Hmax trend: n=%d peaks, %d populations\n", nrow(hd), dplyr::n_distinct(hd$pop_id)))
print(summary(m_hmax)$coef$cond)
b_hmax <- summary(m_hmax)$coef$cond["YEAR_c", ]
yr0 <- min(hd$YEAR)

# Population-average marginal trend + 95% CI via ggeffects::ggpredict, on the real log scale.
pa <- as.data.frame(ggpredict(m_hmax, terms = "YEAR_c [all]"))
pa <- transform(pa, YEAR = x + yr0)

# Individual population trends by independent OLS (transparent per-population significance):
# intercept a, slope b, slope p-value; "significant decline" = b < 0 & p < 0.05 (potential deepening).
# also record each population's OWN sampled year span (y0, y1) so its line spans only those years.
popfit <- hd |> group_by(pop_id) |>
  group_modify(~ {
    s <- tryCatch(summary(lm(logHmax ~ YEAR_c, .x))$coefficients, error = function(e) NULL)
    if (is.null(s) || !"YEAR_c" %in% rownames(s))
      data.frame(a = NA_real_, b = NA_real_, p = NA_real_, y0 = NA_real_, y1 = NA_real_)
    else data.frame(a = s["(Intercept)", "Estimate"], b = s["YEAR_c", "Estimate"], p = s["YEAR_c", "Pr(>|t|)"],
                    y0 = min(.x$YEAR_c), y1 = max(.x$YEAR_c)) }) |>
  ungroup() |> filter(is.finite(a), is.finite(b), is.finite(p)) |>
  mutate(trend = dplyr::case_when(b < 0 & p < 0.05 ~ "significant decline",
                                  b > 0 & p < 0.05 ~ "significant increase",
                                  TRUE             ~ "not significant"))
popfit$trend <- factor(popfit$trend, levels = c("not significant", "significant decline", "significant increase"))
pct_signeg <- 100 * mean(popfit$trend == "significant decline")
pct_sigpos <- 100 * mean(popfit$trend == "significant increase")

# individual estimated trend lines over EACH population's actual sampled years (y0..y1),
# at their real heights (own intercept + slope), coloured by decline.
ind <- data.frame(pop_id = rep(popfit$pop_id, each = 2),
                  YEAR_c = as.vector(t(cbind(popfit$y0, popfit$y1))),
                  a = rep(popfit$a, each = 2), b = rep(popfit$b, each = 2),
                  trend  = rep(popfit$trend, each = 2))
ind$YEAR    <- ind$YEAR_c + yr0
ind$logHmax <- ind$a + ind$b * ind$YEAR_c
trend_pal <- c("significant decline" = cock_red, "not significant" = "grey80", "significant increase" = cock_blue)

p_hmax <- ggplot() +
  geom_line(data = ind, aes(YEAR, logHmax, group = pop_id, colour = trend, alpha = trend), linewidth = .3) +
  geom_ribbon(data = pa, aes(YEAR, ymin = conf.low, ymax = conf.high), fill = "grey45", alpha = .5) +
  geom_line(data = pa, aes(YEAR, predicted), colour = "black", linewidth = 1.3) +
  scale_colour_manual(values = trend_pal, name = NULL) +
  scale_alpha_manual(values = c("significant decline" = .5, "not significant" = .10, "significant increase" = .5), guide = "none") +
  labs(x = "Year", y = "log peak modal abundance",
       subtitle = sprintf("population-average %+.1f%% yr⁻¹ (p = %s); significant declines %.0f%%, increases %.0f%%",
                          100 * (exp(b_hmax[1]) - 1), format.pval(b_hmax[4], digits = 2, eps = 1e-3),
                          pct_signeg, pct_sigpos)) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
plot(p_hmax)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/FigS_Hmax_trend_%s.pdf", fig_key), p_hmax, width = 6, height = 4.5)

# Alternative view — distribution of per-population trends; significant decline red, increase blue.
p_slopes <- ggplot(popfit, aes(b, fill = trend)) +
  geom_histogram(bins = 60, colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = trend_pal, name = NULL) +
  coord_cartesian(xlim = quantile(popfit$b, c(.005, .995), na.rm = TRUE)) +
  labs(x = "Per-population trend in log peak modal abundance  (yr⁻¹)", y = "Populations",
       subtitle = sprintf("significantly negative %.0f%%, significantly positive %.0f%% (p < 0.05)",
                          pct_signeg, pct_sigpos)) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
plot(p_slopes)
# ggsave(sprintf("~/Lavori/MPA_timeseries/Modskurt/Figs/FigS_Hmax_slopes_%s.pdf", fig_key), p_slopes, width = 6, height = 4)
