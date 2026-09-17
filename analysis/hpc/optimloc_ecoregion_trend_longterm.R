#### ---- Long-term reference-site trend (firstsite or previoussite) ---- ####
#
# Equivalent to optimloc_ecoregion_trend.R but the input data are exposures
# computed at the chosen REFERENCE optimum site of each population (first
# or second-to-last sampling year's optimum), extended back to 1993 and
# forward to the last sampling year. Assumption: the population would have
# been present at the reference site for the full length of the
# temperature time series.
#
# Args (parms.temp.longterm.txt, 6 columns):
#   var.type    Env | Abund | Lat | Lon
#   resp.env    shifted.cumtemp.above.{mean,q50,q70,q90,q95,q97.5} | ...
#               (column name in the long-term reference panel)
#   autocor     no_autocor | ou | ar1
#   trans       NULL | Standardize
#   mod.family  gaussian | lognormal | nbinom2
#   ref.mode    firstsite | previoussite
#
# Output saved to Optim_trend_longterm_<ref.mode>/<mod.parms>.RData with
# the same row structure as optim.abund.trend (REALM, ECOREGION, SPECIES,
# ..., Effect: Env.Int / Env.Trend, Estimate, SE, t.value, P.value, AIC,
# Warning, mod.parms).
#
# Expected input file (built by sp_optimloc_ecoreg_longterm_ref_build.R):
#   sp.optim.abund.<ref.mode>.longterm.RData
#######################################################################

require(tidyr)
require(lubridate)
require(dplyr)
require(forcats)
require(glmmTMB)
require(datawizard)

require(foreach, quietly = TRUE)
require(doMC,   quietly = TRUE)
registerDoMC(cores = 20)
source("config.R")
ref.mode <- "firstsite" #c("firstsite", "previoussite")
data.type <- "abund"  # c("abund", "density")

args        <- commandArgs(trailingOnly = TRUE)
var.type    <- as.character(args[1])
resp.env    <- as.character(args[2])
autocor     <- as.character(args[3])
trans       <- as.character(args[4])
mod.family  <- gsub("\r", "", as.character(args[5]))

out.subdir <- paste0("Optim_trend_", ref.mode, "_longterm")
out.dir    <- file.path(dir_data, out.subdir)
if (!file.exists(out.dir)) dir.create(out.dir)

# --- Load the long-term reference-site exposure panel ------------------
ref_obj_name <- paste0("sp.optim.", data.type, ".", ref.mode, ".longterm")
load(file.path(dir_data,
               paste0(ref_obj_name, ".RData")))
load(input_file("fish.traits.dat.RData"))

shift.df <- get(ref_obj_name)
if ("ORIG.ECOREGION" %in% names(shift.df))
  shift.df <- shift.df |> rename(ECOREGION = ORIG.ECOREGION)

ftr <- fish.traits.dat |>
  select(SPECIES, FeedingType, sst_mean_TM)

shifted.temp.df <- shift.df |>
  rename(LAT = SHIFTED.LAT, LON = SHIFTED.LON) |>
  left_join(ftr, by = "SPECIES") |>
  group_by(ECOREGION, SPECIES) |>
  mutate(
    ndays         = ifelse(leap_year(YEAR) == "TRUE", 366, 365),
    TEMP.OPTIM    = mean(sst_mean_TM, na.rm = TRUE),
    THERMAL.GUILD = case_when(
      TEMP.OPTIM <  23 ~ "Temperate",
      TEMP.OPTIM >= 23 ~ "Tropical"),
    FEEDING.TYPE = case_when(
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
  relocate(REALM, .before = ECOREGION)

# --- Error/warning catcher (identical to original scripts) -------------
myCatch <- function(expr, ...) {
  res  <- NULL
  diag <- NULL
  res <- tryCatch(
    expr = withCallingHandlers(expr,
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

pred.env <- resp.env
min.yrs  <- 4

sp.id <- unique(shifted.temp.df$SPECIES)

shifted.trend.analysis <- foreach(i = seq_along(sp.id), .combine = rbind) %dopar% {

  cat("Doing Species ", i, " of ", length(sp.id), "\n", sep = "")

  tmp.df <- shifted.temp.df |> filter(SPECIES %in% sp.id[i])

  n.yrs <- tmp.df |> group_by(ECOREGION) |>
    summarise(nyrs = n(), .groups = "drop") |>
    filter(nyrs >= min.yrs)

  pop.ecoreg.df <- tmp.df |> filter(ECOREGION %in% n.yrs$ECOREGION)

  if (nrow(pop.ecoreg.df) < min.yrs) return(NULL)

  model.out <- foreach(j = seq_len(nrow(n.yrs)), .combine = "rbind") %do% {

    pop.df <- pop.ecoreg.df |>
      filter(ECOREGION %in% unlist(n.yrs[j, "ECOREGION"]))

    if (var.type == "Env") {

      # Long-term first-site trend: no abund > 0 filter, the assumption is
      # that the population would have been present at the first site for
      # the full length of the temperature record.
      if (trans == "Standardize") {
        lm1.df <- pop.df |>
          mutate(
            RESP.ENV  = !!sym(resp.env),
            YEAR      = (YEAR - YEAR[1]) + 1,
            n.yrs     = n(),
            times.dec = as.numeric(paste(YEAR, 1, sep = ".")),
            cum.times = lead(cumsum(lag(times.dec, default = 0)),
                             default = sum(times.dec)),
            yrs       = glmmTMB::numFactor(cum.times),
            YEAR      = as.vector(standardize(YEAR)),
            group     = 1)
      } else if (trans == "NULL") {
        lm1.df <- pop.df |>
          mutate(
            RESP.ENV  = !!sym(resp.env),
            YEAR      = (YEAR - YEAR[1]) + 1,
            n.yrs     = n(),
            times.dec = as.numeric(paste(YEAR, 1, sep = ".")),
            cum.times = lead(cumsum(lag(times.dec, default = 0)),
                             default = sum(times.dec)),
            yrs       = glmmTMB::numFactor(cum.times),
            group     = 1)
      }

      if (autocor == "no_autocor")  my.fm1 <- formula(RESP.ENV ~ YEAR + (1 | group))
      else if (autocor == "ou")     my.fm1 <- formula(RESP.ENV ~ YEAR + ou(yrs + 0 | group))
      else if (autocor == "ar1")    my.fm1 <- formula(RESP.ENV ~ YEAR + ar1(as.factor(YEAR) + 0 | group))

      m1 <- myCatch(glmmTMB(my.fm1, data = lm1.df,
                            ziformula = ~ 0, dispformula = ~ 1,
                            family = mod.family))

      if (any(names(m1) == "warning")) {
        mw1 <- m1$warning
        if (grepl("non-positive-definite Hessian matrix", mw1) |
            grepl("false convergence", mw1)) {
          m1 <- myCatch(glmmTMB(my.fm1, data = lm1.df,
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(optimizer = optim,
                                                        optArgs = list(method = "BFGS"))))
          if (any(names(m1) == "warning")) mw1 <- m1$warning
        }
        if (grepl("limit reached without convergence", mw1)) {
          m1 <- myCatch(glmmTMB(my.fm1, data = lm1.df,
                                ziformula = ~ 0, dispformula = ~ 1,
                                family = mod.family,
                                control = glmmTMBControl(
                                  optCtrl = list(iter.max = 1e3, eval.max = 1e3))))
        }
      }

      if (!is.null(m1)) {
        w1 <- if (any(names(m1) == "warning")) m1$warning else NA
        coefs.out1.m1 <- data.frame(summary(m1)$coefficients$cond,
                                    AIC     = rep(AIC(m1), 2),
                                    Warning = rep(w1, 2))
        colnames(coefs.out1.m1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
        rownames(coefs.out1.m1) <- NULL
      } else {
        coefs.out1.m1 <- data.frame(t(rep(NA, 4)),
                                    AIC     = rep(NA, 2),
                                    Warning = rep(NA, 2))
        colnames(coefs.out1.m1) <- c("Estimate", "SE", "t.value", "P.value", "AIC", "Warning")
        rownames(coefs.out1.m1) <- NULL
      }

      out1 <- data.frame(
        REALM         = lm1.df[1, "REALM"],
        ECOREGION     = lm1.df[1, "ECOREGION"],
        SPECIES       = lm1.df[1, "SPECIES"],
        FEEDING.TYPE  = lm1.df[1, "FEEDING.TYPE"],
        TEMP.OPTIM    = lm1.df[1, "TEMP.OPTIM"],
        THERMAL.GUILD = lm1.df[1, "THERMAL.GUILD"],
        FlagAnalysis  = "Trend_Environmental",
        Nyears        = lm1.df[1, "n.yrs"],
        Effect        = c("Env.Int", "Env.Trend"),
        coefs.out1.m1)

      out <- out1
    }

    out
  }

  model.out
}

if (!is.null(shifted.trend.analysis)) {
  mod.parms <- paste(var.type, resp.env, ref.mode, autocor, trans, mod.family, sep = "_")
  shifted.trend.analysis$mod.parms <- mod.parms
  shifted.trend.analysis <- shifted.trend.analysis |>
    relocate(mod.parms, .before = "Warning")

  outputName <- paste(mod.parms, ".RData", sep = "")
  outputPath <- file.path(out.dir, outputName)
  assign(paste(mod.parms), value = shifted.trend.analysis, pos = 1, inherits = TRUE)
  save(list = paste(mod.parms), file = outputPath)
}
