#### ---- Long-term counterfactual (firstsite or previoussite) ---- ####
#
# Equivalent to optimloc_ecoregion_leadtime_trend.R but the counterfactual
# trajectory (ref.var) is taken from a long-term reference panel covering
# 1993 → last sampling year of the population, built by
# sp_optimloc_ecoreg_longterm_ref_build.R. Observed trajectory (obs.var)
# is the realised exposure at the yearly optimum sites during sampling
# years (sp.optim.abund.shift). obs and ref are stacked with their native
# year ranges intact and fit with the same interaction model.
#
# Args (6, longterm.parms.temp.firstsite.txt or .previoussite.txt):
#   var.type     Env | Abund
#   obs.env.var  e.g. shifted.cumtemp.above.mean   (column in sp.optim.abund.shift)
#   ref.env.var  same name as obs.env.var          (column in the long-term ref file)
#   autocor      no_autocor | ou | ar1
#   trans        NULL | Standardize
#   mod.family   gaussian | lognormal | nbinom2
#
# The choice of firstsite vs previoussite reference is controlled by which
# long-term reference file is loaded (see ref.mode arg below). Same script
# produces Optim_longterm_firstsite_leadtime_trend/ and
# Optim_longterm_previoussite_leadtime_trend/ outputs.
#
# Output: 6-row-per-population effect table (Ref.Int, Ref.Slope, Int.Dev,
# Slope.Dev, Int.Alone, Slope.Alone) — same shape as
# optim.abund.leadtime.firstsite.trend, so optimum_response_scenarios()
# can be applied directly.
#######################################################################

require(tidyr)
require(lubridate)
require(dplyr)
require(forcats)
require(glmmTMB)
require(datawizard)
require(foreach, quietly = TRUE)
require(doMC,    quietly = TRUE)
registerDoMC(cores = 20)
source("config.R")
# Set ref.mode and data.type at the top --------------------
ref.mode <- "firstsite" #c("firstsite", "previoussite")
data.type <- "abund"  # c("abund", "density")

# Parse args (5-column neutral parms, e.g. leadtime.parms.temp.longterm.txt)
# In the long-term version the ref panel uses the same exposure column
# names as the obs panel, so ref.env.var is set to obs.env.var (no
# contrast.var column needed).
args         <- commandArgs(trailingOnly = TRUE)
var.type     <- as.character(args[1])
obs.env.var  <- as.character(args[2])
autocor      <- as.character(args[3])
trans        <- as.character(args[4])
mod.family   <- gsub("\r", "", as.character(args[5]))
ref.env.var  <- obs.env.var

out.subdir <- paste0("Optim_leadtime_trend_", ref.mode, "_longterm")
if (!file.exists(file.path(dir_data, out.subdir)))
  dir.create(file.path(dir_data, out.subdir))

# Load data --------------------------------------------------------------
load(input_file("fish.traits.dat.RData"))

obs_obj_name <- paste0("sp.optim.", data.type, ".shift")
load(file.path(dir_data,
          paste0(obs_obj_name, ".RData")))
obs_obj <- get(obs_obj_name)

ref_obj_name <- paste0("sp.optim.", data.type, ".", ref.mode, ".longterm")
load(file.path(dir_data,
               paste0(ref_obj_name, ".RData")))
ref_panel <- get(ref_obj_name)

ftr <- fish.traits.dat |>
  select(SPECIES, FeedingType, sst_mean_TM) |>
  rename(temp.optim = sst_mean_TM)
attr(ftr$temp.optim, "names") <- NULL

obs_panel <- obs_obj |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON, ECOREGION = ORIG.ECOREGION) |>
  left_join(ftr, by = c("SPECIES.ORIG" = "SPECIES")) |>
  group_by(ECOREGION, SPECIES, SPECIES.ORIG) |>
  mutate(
    TEMP.OPTIM    = mean(temp.optim, na.rm = TRUE),
    THERMAL.GUILD = case_when(
      temp.optim <  23 ~ "Temperate",
      temp.optim >= 23 ~ "Tropical"),
    FEEDING.TYPE  = case_when(
      FeedingType %in% c("browsing on substrate") ~ "Microphages",
      FeedingType %in% c("selective plankton feeding", "filtering plankton",
                         "variable") ~ "Planktivores",
      FeedingType %in% c("grazing on aquatic plants") ~ "Grazers",
      FeedingType %in% c("hunting macrofauna (predator)",
                         "picking parasites off a host (cleaner)") ~ "Carnivores",
      FeedingType %in% c("other") ~ "Other",
      TRUE ~ as.character(FeedingType)),
    REALM = case_when(
      (LAT >=  40 & LAT <=  60) | (LAT >= -60 & LAT <= -40) ~ "Temperate",
      (LAT >  23.5 & LAT <  40) | (LAT > -40 & LAT < -23.5) ~ "Subtropical",
      (LAT <= 23.5 & LAT >= -23.5) ~ "Tropical"),
    THERMAL.GUILD = factor(THERMAL.GUILD),
    FEEDING.TYPE  = factor(FEEDING.TYPE),
    REALM         = factor(REALM),
    REALM         = fct_relevel(REALM, c("Tropical", "Subtropical", "Temperate"))) |>
  ungroup() |>
  relocate(REALM, .before = ECOREGION)

# Error/warning catcher (identical to original) --------------------------
myCatch <- function(expr, ...) {
  res  <- NULL
  diag <- NULL
  res <- tryCatch(
    expr  = withCallingHandlers(expr,
      warning = function(w) {
        parent <- parent.env(environment())
        parent$diag <- w
      }),
    error = function(e) NULL)
  diagnostic <- function(result, w) last.message <<- w$message
  warning <- diagnostic(res, diag)
  if (!is.null(res) & !is.null(warning)) res$warning <- warning
  return(res)
}

min.yrs <- 4
sp.id   <- unique(obs_panel$SPECIES)

shifted.trend.analysis.longterm <- foreach(i = seq_along(sp.id), .combine = rbind) %dopar% {

  cat("Doing Species ", i, " of ", length(sp.id), "\n", sep = "")

  tmp.obs <- obs_panel |> filter(SPECIES %in% sp.id[i])
  tmp.ref <- ref_panel |> filter(SPECIES %in% sp.id[i])

  n.yrs <- tmp.obs |> group_by(ECOREGION) |>
    summarise(nyrs = n(), .groups = "drop") |>
    filter(nyrs >= min.yrs)

  pop.ecoreg.df <- tmp.obs |> filter(ECOREGION %in% n.yrs$ECOREGION)

  if (nrow(pop.ecoreg.df) < min.yrs) return(NULL)

  model.out <- foreach(j = seq_len(nrow(n.yrs)), .combine = "rbind") %do% {

    eco <- unlist(n.yrs[j, "ECOREGION"])
    pop.obs <- pop.ecoreg.df |> filter(ECOREGION %in% eco)
    pop.ref <- tmp.ref |> filter(ECOREGION %in% eco)
    out <- NULL                       # default; remains NULL on early-skip

    if (var.type == "Env" &&
        obs.env.var %in% names(pop.obs) &&
        ref.env.var %in% names(pop.ref) &&
        nrow(pop.ref) >= min.yrs) {

      df.obs <- pop.obs |>
        transmute(REALM, ECOREGION, SPECIES, SPECIES.ORIG, FEEDING.TYPE, TEMP.OPTIM,
                  THERMAL.GUILD, YEAR,
                  ENV.VAR = .data[[obs.env.var]],
                  CONTR   = "obs.var")
      df.ref <- pop.ref |>
        transmute(ECOREGION, SPECIES, YEAR,
                  SPECIES.ORIG = df.obs$SPECIES.ORIG[1],
                  ENV.VAR = .data[[ref.env.var]],
                  CONTR   = "ref.var",
                  REALM         = df.obs$REALM[1],
                  FEEDING.TYPE  = df.obs$FEEDING.TYPE[1],
                  TEMP.OPTIM    = df.obs$TEMP.OPTIM[1],
                  THERMAL.GUILD = df.obs$THERMAL.GUILD[1])

      stacked <- bind_rows(df.ref, df.obs) |>
        drop_na(ENV.VAR) |>
        arrange(CONTR, YEAR) |>
        mutate(YEAR_c    = (YEAR - min(YEAR)) + 1,
               n.yrs     = n(),
               times.dec = as.numeric(paste(YEAR_c, 1, sep = ".")),
               cum.times = lead(cumsum(lag(times.dec, default = 0)),
                                default = sum(times.dec)),
               yrs       = glmmTMB::numFactor(cum.times),
               group     = 1,
               CONTR     = factor(CONTR, levels = c("ref.var", "obs.var")))

      if (trans == "Standardize") {
        stacked <- stacked |>
          mutate(YEAR    = as.vector(standardize(YEAR_c)),
                 ENV.VAR = as.vector(standardize(ENV.VAR)))
      } else if (trans == "NULL") {
        stacked <- stacked |> mutate(YEAR = YEAR_c)
      }

      if (autocor == "no_autocor")  my.fm1 <- formula(ENV.VAR ~ YEAR * CONTR + (1 | group))
      else if (autocor == "ou")     my.fm1 <- formula(ENV.VAR ~ YEAR * CONTR + ou(yrs + 0 | group))
      else if (autocor == "ar1")    my.fm1 <- formula(ENV.VAR ~ YEAR * CONTR + ar1(as.factor(YEAR) + 0 | group))

      m1 <- myCatch(glmmTMB(my.fm1, data = stacked,
                            ziformula = ~ 0, dispformula = ~ 1,
                            family = mod.family))
      if (any(names(m1) == "warning")) {
        mw1 <- m1$warning
        if (grepl("non-positive-definite Hessian matrix", mw1) |
            grepl("false convergence", mw1)) {
          m1 <- myCatch(glmmTMB(my.fm1, data = stacked,
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(optimizer = optim,
                                                        optArgs = list(method = "BFGS"))))
          if (any(names(m1) == "warning")) mw1 <- m1$warning
        }
        if (grepl("limit reached without convergence", mw1)) {
          m1 <- myCatch(glmmTMB(my.fm1, data = stacked,
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(
                                  optCtrl = list(iter.max = 1e3, eval.max = 1e3))))
        }
      }

      m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
                            data = stacked |> filter(CONTR == "obs.var"),
                            ziformula = ~ 0, dispformula = ~ 1,
                            family = mod.family))
      if (any(names(m2) == "warning")) {
        mw2 <- m2$warning
        if (grepl("non-positive-definite Hessian matrix", mw2) |
            grepl("false convergence", mw2)) {
          m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
                                data = stacked |> filter(CONTR == "obs.var"),
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(optimizer = optim,
                                                        optArgs = list(method = "BFGS"))))
          if (any(names(m2) == "warning")) mw2 <- m2$warning
        }
        if (grepl("limit reached without convergence", mw2)) {
          m2 <- myCatch(glmmTMB(update(my.fm1, ~ . - CONTR - YEAR:CONTR),
                                data = stacked |> filter(CONTR == "obs.var"),
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(
                                  optCtrl = list(iter.max = 1e3, eval.max = 1e3))))
        }
      }

      if (!is.null(m1) | !is.null(m2)) {
        if (!is.null(m1) & !is.null(m2)) {
          w1 <- if (any(names(m1) == "warning")) m1$warning else NA
          w2 <- if (any(names(m2) == "warning")) m2$warning else NA
          coefs.out1 <- data.frame(rbind(summary(m1)$coefficients$cond,
                                         summary(m2)$coefficients$cond),
                                   AIC     = c(rep(AIC(m1), 4), rep(AIC(m2), 2)),
                                   Warning = c(rep(w1, 4),     rep(w2, 2)))
        } else if (!is.null(m1) & is.null(m2)) {
          w1 <- if (any(names(m1) == "warning")) m1$warning else NA
          w2 <- NA
          coefs.out1 <- data.frame(rbind(summary(m1)$coefficients$cond,
                                         rep(NA, 4), rep(NA, 4)),
                                   AIC     = c(rep(AIC(m1), 4), rep(NA, 2)),
                                   Warning = c(rep(w1, 4),      rep(w2, 2)))
        } else {
          w1 <- NA
          w2 <- if (any(names(m2) == "warning")) m2$warning else NA
          coefs.out1 <- data.frame(rbind(rep(NA, 4), rep(NA, 4), rep(NA, 4), rep(NA, 4),
                                         summary(m2)$coefficients$cond),
                                   AIC     = c(rep(NA, 4),       rep(AIC(m2), 2)),
                                   Warning = c(rep(w1, 4),       rep(w2, 2)))
        }
        colnames(coefs.out1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
        rownames(coefs.out1) <- NULL

        out1 <- data.frame(
          REALM         = stacked[1, "REALM"],
          ECOREGION     = stacked[1, "ECOREGION"],
          SPECIES       = stacked[1, "SPECIES"],
          SPECIES.ORIG  = stacked[1, "SPECIES.ORIG"],
          FEEDING.TYPE  = stacked[1, "FEEDING.TYPE"],
          TEMP.OPTIM    = stacked[1, "TEMP.OPTIM"],
          THERMAL.GUILD = stacked[1, "THERMAL.GUILD"],
          FlagAnalysis  = "Trend_Environmental",
          NYEARS        = stacked[1, "n.yrs"],
          Effect        = c("Ref.Int", "Ref.Slope", "Int.Dev", "Slope.Dev",
                            "Int.Alone", "Slope.Alone"),
          coefs.out1)
      } else { out1 <- NULL }

      out <- out1
    }
    out
  }
  model.out
}

if (!is.null(shifted.trend.analysis.longterm)) {
  mod.parms <- paste(var.type, obs.env.var, ref.mode, autocor, trans, mod.family, sep = "_")
  shifted.trend.analysis.longterm <- shifted.trend.analysis.longterm |>
    mutate(mod.parms = mod.parms) |>
    relocate(mod.parms, .before = "Warning")

  outputName <- paste(mod.parms, ".RData", sep = "")
  outputPath <- file.path(dir_data, out.subdir,
                          outputName)
  assign(paste(mod.parms), value = shifted.trend.analysis.longterm,
         pos = 1, inherits = TRUE)
  save(list = paste(mod.parms), file = outputPath)
  cat("Saved:", outputPath, "\n")
}
