# ============================================================================
# config.R — single source of paths for the AbundanceOptima analysis.
# Source this at the top of every script:   source("config.R")
# Every script reads/writes via the dir_* variables below — never hard-code
# an absolute path anywhere else.
# ============================================================================

# ---- EDIT THIS ONE LINE to point at your project root --------------------
# For released code, users create ~/abundance_optima and place the data there:
PROJ <- path.expand("~/abundance_optima")
# During development on this machine, point it at the staging copy instead:
# PROJ <- "G:/My Drive/workspace/BiodGlob/Fish_Thermal/AbundanceOptima"
# --------------------------------------------------------------------------

dir_data    <- file.path(PROJ, "data")      # inputs: provided .RData (Zenodo)
dir_packages<- file.path(PROJ, "packages")  # modskurt1 + CockR sources
dir_out     <- file.path(PROJ, "output")    # intermediate .RData generated along the pipeline
dir_results <- file.path(PROJ, "results")   # FINAL outputs (figures, tables)
dir_figures <- file.path(PROJ, "figures")   # figure scripts + their helpers (e.g. plot_sites.R)

# CONVENTION: a script writes its final outputs into its OWN named subfolder
# under dir_results, keeping the folder name it already uses in code, e.g.
#     outdir <- file.path(dir_results, "Fig1_abund_select"); dir.create(outdir, ...)
# dir_results is the common parent; the script owns the leaf folder name.
for (d in c(dir_out, dir_results))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Resolve an input .RData by name. Prefer a copy you generated in output/;
# otherwise use the provided (Zenodo) copy in data/. So the same script works
# whether you re-ran the pipeline or just downloaded the intermediates.
input_file <- function(name) {
  p_out <- file.path(dir_out, name)
  if (file.exists(p_out)) p_out else file.path(dir_data, name)
}

# ---- Large environmental grids (NOT redistributed; user downloads — see README)
# Point these at wherever the user stored the big files.
gebco_file <- file.path(dir_data, "GEBCO_2024.nc")   # GEBCO bathymetry
glorys_dir <- file.path(dir_data, "GLORYS")          # GLORYS thetao (subsurface temp)
meow_shapefile <- file.path(dir_data, "Marine_Ecoregions_Of_the_World__MEOW_.shp")  # MEOW ecoregions (+ .shx/.dbf/.prj)

# Fail early with a clear message if the root is missing.
if (!dir.exists(PROJ))
  stop("PROJ not found: ", PROJ,
       "\nCreate it (see README) or edit PROJ in config.R.")

# Convenience: R library path used for the analyses (edit or remove as needed).
# .libPaths(c("C:/R_libs/win-library/4.5", .libPaths()))
