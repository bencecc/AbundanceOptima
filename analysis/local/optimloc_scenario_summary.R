# optimloc_scenario_summary.R — build the counterfactual SCENARIO files used by Fig 4.
# For each population it classifies the response to warming (Mitigation / Magnification /
# Matching / ...) with optimum_response_scenarios(), then saves one scenario.*.RData per
# combination. ALL eight combinations are produced (abundance x density) x (first-site x
# previous-site) x (short-term x long-term) — short-term included for consistency checks
# even where the paper shows only the long-term panels.
#
# Single R session. Inputs = the stage-9 trend outputs (from optimloc_ecoregion_*trend*.R,
# reassembled by cluster_summaries.R). Paths are still the original server paths —
# generalise via config.R when finalising.
###############################################################################

require(dplyr)
require(tidyr)

source("config.R")
# the response-classification function (bundled sibling)
source(file.path(PROJ, "analysis/local/optimum_response_scenarios.R"))

# ---- stage-9 trend inputs (all by ecoregion) --------------------------------
# short-term (leadtime): "firstsite" is the default leadtime trend; prevsite is explicit
load(input_file("optim.abund.leadtime.trend.RData"))
load(input_file("optim.density.leadtime.trend.RData"))
load(input_file("optim.abund.leadtime.prevsite.trend.RData"))
load(input_file("optim.density.leadtime.prevsite.trend.RData"))
# long-term
load(input_file("optim.abund.leadtime.firstsite.trend.longterm.RData"))
load(input_file("optim.abund.leadtime.previoussite.trend.longterm.RData"))
load(input_file("optim.density.leadtime.firstsite.trend.longterm.RData"))
load(input_file("optim.density.leadtime.previoussite.trend.longterm.RData"))

# ---- helper: classify each population, drop unresolved, save ----------------
# Groups by ecoregion x species x analysis flag x model, applies the scenario
# classifier, and writes <out>.RData (object named <out>).
build_scenario <- function(trend_df, out) {
  sc <- trend_df |>
    group_by(ECOREGION, SPECIES, FlagAnalysis, mod.parms) |>
    group_modify(~ optimum_response_scenarios(.)) |>
    ungroup() |>
    drop_na(Response)
  assign(out, sc)
  save(list = out, file = file.path(dir_data, paste0(out, ".RData")))
  cat("saved", out, "-", nrow(sc), "rows\n")
}

# ---- SHORT-TERM scenarios ---------------------------------------------------
build_scenario(optim.abund.leadtime.trend,            "scenario.abund.leadtime.firstsite")
build_scenario(optim.abund.leadtime.prevsite.trend,   "scenario.abund.leadtime.prevsite")
build_scenario(optim.density.leadtime.trend,          "scenario.density.leadtime.firstsite")
build_scenario(optim.density.leadtime.prevsite.trend, "scenario.density.leadtime.prevsite")

# ---- LONG-TERM scenarios ----------------------------------------------------
build_scenario(optim.abund.leadtime.firstsite.trend.longterm,      "scenario.abund.leadtime.firstsite.longterm")
build_scenario(optim.abund.leadtime.previoussite.trend.longterm,   "scenario.abund.leadtime.previoussite.longterm")
build_scenario(optim.density.leadtime.firstsite.trend.longterm,    "scenario.density.leadtime.firstsite.longterm")
build_scenario(optim.density.leadtime.previoussite.trend.longterm, "scenario.density.leadtime.previoussite.longterm")
