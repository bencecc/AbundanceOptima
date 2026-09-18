# ==========================================================================
# unimodal_filter_helper.R
#
# Drops populations flagged by the modskurt unimodal-fit QC.
# Source this script then call unimodal_filter_helper() on whichever data
# frames you want filtered.
#
# USAGE
#   source(file.path(PROJ, "analysis/local/unimodal_filter_helper.R"))   # PROJ from config.R
#
#   # Abundance, all deviation classes, default 50% threshold
#   sp.optim.abund.shift <- unimodal_filter_helper(
#     sp.optim.abund.shift,
#     decision_csv = file.path(dir_data, "unimodal_fit_check_pop_decision.csv"))
#
#   # Bimodal-only (~5 pops), abundance
#   sp.optim.abund.shift <- unimodal_filter_helper(
#     sp.optim.abund.shift,
#     classes      = "bimodal",
#     decision_csv = file.path(dir_data, "unimodal_fit_check_pop_decision.csv"))
#
#   # Density data
#   sp.optim.density.shift <- unimodal_filter_helper(
#     sp.optim.density.shift,
#     decision_csv = file.path(dir_data, "unimodal_fit_check_pop_decision_density.csv"))
#
# ARGUMENTS
#   df            data frame with SPECIES + (ECOREGION or ORIG.ECOREGION)
#   classes       deviation classes counted toward exclusion. Subset of
#                 c("bimodal","edge_peak","flat"). Default = all three.
#   threshold     drop a population if sum(pct_<class>) >= threshold (default 50)
#   decision_csv  full path to the unimodal_fit_check_pop_decision csv
#                 (no default — supply abundance or density file explicitly)
#
# REQUIREMENT
#   sp.optim.abund.shift must be in memory; the QC csv stores ECOREGION_ID
#   (integer) and we map it back to the ECOREGION name using the same
#   `unique(test.dat$ECOREGION)` per-species order as modskurt_analysis.R.
# ==========================================================================

unimodal_filter_helper <- function(df,
                                   classes      = c("bimodal", "edge_peak", "flat"),
                                   threshold    = 50,
                                   decision_csv,
                                   sp_df        = NULL,
                                   map_df       = NULL) {
  # sp_df: the sampling data (sp.df). PREFERRED. The decision CSV's ECOREGION_ID
  #   is the per-species index into unique(sp.df$ECOREGION) in FIRST-APPEARANCE
  #   order (exactly how modskurt_analysis.R numbered the _ecoreg_<i>_ files), so
  #   sp.df is the only correct source for the id->name map. Bassian E/W split is
  #   handled via df$SPECIES.ORIG (suffixed species -> "Bassian", id 1).
  # map_df: DEPRECATED legacy fallback (used only if sp_df is NULL). Reconstructs
  #   the id from map_df's row order, which mis-maps multi-ecoregion species when
  #   map_df is sorted differently from sp.df (e.g. alphabetically by ECOREGION).
  # map_df: data frame (with SPECIES + ORIG.ECOREGION) used to map the QC's
  #   per-species ECOREGION_ID back to the ECOREGION name. Defaults to
  #   sp.optim.abund.shift in .GlobalEnv (back-compat). For DENSITY data pass
  #   map_df = sp.optim.density.shift, so the per-species ecoregion order matches
  #   the density fits rather than abundance.

  if (missing(decision_csv))
    stop("Provide decision_csv = '...unimodal_fit_check_pop_decision[...].csv'")
  if (!"SPECIES" %in% names(df)) {
    message("unimodal_filter_helper: SPECIES not in data frame - returning unchanged")
    return(df)
  }
  eco_col <- intersect(c("ECOREGION", "ORIG.ECOREGION"), names(df))[1]
  if (is.na(eco_col)) {
    message("unimodal_filter_helper: no ECOREGION column - returning unchanged")
    return(df)
  }
  # --- Build the REMOVE list from the QC csv -------------------------------
  pop_decision <- read.csv(decision_csv, stringsAsFactors = FALSE)
  pct_cols     <- paste0("pct_", classes)
  if (!all(pct_cols %in% names(pop_decision)))
    stop("classes must be a subset of c('bimodal','edge_peak','flat'); ",
         "missing column(s): ",
         paste(setdiff(pct_cols, names(pop_decision)), collapse = ", "))

  remove_list <- pop_decision |>
    dplyr::mutate(pct_excl = rowSums(dplyr::across(dplyr::all_of(pct_cols)))) |>
    dplyr::filter(pct_excl >= threshold) |>
    dplyr::select(SPECIES_dot = SPECIES, ECOREGION_ID)

  # --- Map (SPECIES, ECOREGION_ID) -> ECOREGION name -----------------------
  if (!is.null(sp_df)) {
    # CORRECT path: per-species ecoregion order = unique(sp.df$ECOREGION) in
    # first-appearance order (= modskurt_analysis.R's _ecoreg_<i>_ numbering).
    eco_lookup <- sp_df |>
      dplyr::distinct(SPECIES, ECOREGION) |>             # keeps first-appearance order
      dplyr::group_by(SPECIES) |>
      dplyr::mutate(ECOREGION_ID = dplyr::row_number()) |>
      dplyr::ungroup()
    # Bassian E/W split: suffixed species (SPECIES != SPECIES.ORIG in the data)
    # aren't in sp.df; each was fit on one coast -> ecoreg index 1 -> "Bassian".
    if ("SPECIES.ORIG" %in% names(df)) {
      sm <- df[df$SPECIES != df$SPECIES.ORIG, c("SPECIES", eco_col)]
      sm <- sm[!duplicated(sm), , drop = FALSE]
      if (nrow(sm)) {
        names(sm)[2] <- "ECOREGION"; sm$ECOREGION_ID <- 1L
        eco_lookup <- dplyr::bind_rows(eco_lookup, sm)
      }
    }
    eco_lookup$SPECIES_dot <- gsub(" ", ".", eco_lookup$SPECIES)
  } else {
    # DEPRECATED legacy path: reconstruct id from map_df row order (mis-maps
    # multi-ecoregion species when map_df is sorted unlike sp.df). Pass sp_df=.
    warning("unimodal_filter_helper: no sp_df= supplied; using legacy map_df ",
            "row-order id mapping, which can mis-assign multi-ecoregion species. ",
            "Pass sp_df = sp.df for the correct mapping.")
    if (is.null(map_df)) {
      if (!exists("sp.optim.abund.shift", envir = .GlobalEnv))
        stop("Provide sp_df=, map_df=, or load sp.optim.abund.shift.RData.")
      map_df <- get("sp.optim.abund.shift", envir = .GlobalEnv)
    }
    eco_lookup <- map_df |>
      dplyr::rename(ECOREGION = ORIG.ECOREGION) |>
      dplyr::group_by(SPECIES) |>
      dplyr::group_modify(~ tibble::tibble(
          ECOREGION    = unique(.x$ECOREGION),
          ECOREGION_ID = seq_along(unique(.x$ECOREGION)))) |>
      dplyr::ungroup() |>
      dplyr::mutate(SPECIES_dot = gsub(" ", ".", SPECIES))
  }

  remove_pops <- remove_list |>
    dplyr::inner_join(eco_lookup, by = c("SPECIES_dot", "ECOREGION_ID")) |>
    dplyr::select(SPECIES, ECOREGION)

  cat("[unimodal_filter_helper] classes = {",
      paste(classes, collapse = ", "), "} >= ",
      threshold, "%  =>  ", nrow(remove_pops),
      " populations to drop\n", sep = "")

  # --- Apply the filter ----------------------------------------------------
  keys <- remove_pops
  names(keys)[names(keys) == "ECOREGION"] <- eco_col
  n0   <- nrow(df)
  out  <- dplyr::anti_join(df, keys, by = c("SPECIES", eco_col))
  cat("  ", deparse(substitute(df)), ": dropped ", n0 - nrow(out),
      " of ", n0, " rows\n", sep = "")
  out
}
