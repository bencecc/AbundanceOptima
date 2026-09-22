# build_pooled_coasts.R
#
# SENSITIVITY to the ecoregion split (co-author query): builds ONE survey table in which species
# occurring in >= 2 adjacent ecoregions along a spatially continuous coast are treated as ONE
# population per coast, so the abundance optimum is fit over the whole latitudinal gradient of
# that coast instead of within each ecoregion. Three coasts qualify (contiguous coastline with a
# real latitudinal gradient); archipelagos (SE Asia / Coral Triangle), island arcs (Caribbean,
# Pacific islands), basins separated by a strait (W Mediterranean vs Atlantic shelf) and
# single-ecoregion coasts do not, and stay as in the main analysis.
#   East Australia : Torres Strait/N GBR > Central & S GBR > Tweed-Moreton > Manning-Hawkesbury
#                    > Cape Howe > Bassian (eastern group)
#   West/South Aus.: Exmouth to Broome > Ningaloo > Houtman > Leeuwin > South Australian Gulfs
#                    > Bassian (western + north-western groups)
#   California     : Northern California > Southern California Bight
# Bassian's site groups come from split_plan. All rows of a qualifying species from a coast are
# relabelled ECOREGION = "<coast> pooled"; ECOREGION.ORIG keeps the source ecoregion for the
# later per-species comparison with the main (per-ecoregion) results. A species can qualify on
# more than one coast (separate pooled populations).
#
# Input : sp.df.RData, split_plan.RData                                 [config.R paths]
# Output: data/sp.df.pooled.RData (object sp.df.pooled, same columns as sp.df + ECOREGION.ORIG)
#         analysis/hpc/spID_pooled.txt (species ids 1..n, in distinct(SPECIES) order)
# Then run modskurt_analysis(_singlenode).R with input alternative 2), and the downstream
# steps on the pooled outputs.

require(dplyr)
source("config.R")

load(input_file("sp.df.RData")); load(input_file("split_plan.RData"))
chains <- list(
  "East Australia pooled" = c("Torres Strait Northern Great Barrier Reef", "Central and Southern Great Barrier Reef",
                              "Tweed-Moreton", "Manning-Hawkesbury", "Cape Howe", "Bassian"),
  "West Australia pooled" = c("Exmouth to Broome", "Ningaloo", "Houtman", "Leeuwin", "South Australian Gulfs", "Bassian"),
  "California pooled"     = c("Northern California", "Southern California Bight"))
bassian_group <- list("East Australia pooled" = "E", "West Australia pooled" = c("W", "NW"), "California pooled" = character(0))

ckey <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
bs <- split_plan$sites |> filter(ECOREGION == "Bassian") |> mutate(k = ckey(LON, LAT))
stopifnot(all(c("E", "W", "NW") %in% bs$suffix))

pooled <- bind_rows(lapply(names(chains), function(nm) {
  d <- sp.df |> filter(ECOREGION %in% chains[[nm]])
  keep_bass <- bs |> filter(suffix %in% bassian_group[[nm]]) |> pull(k)
  d <- d |> filter(ECOREGION != "Bassian" | ckey(LON, LAT) %in% keep_bass)
  sp_multi <- d |> distinct(SPECIES, ECOREGION) |> count(SPECIES) |> filter(n >= 2) |> pull(SPECIES)
  out <- d |> filter(SPECIES %in% sp_multi) |> mutate(ECOREGION.ORIG = ECOREGION, ECOREGION = nm)
  cat(sprintf("%-22s ecoregions %d | species pooled %4d | rows %6d | sites %4d | lat %5.1f to %5.1f\n", nm,
              n_distinct(out$ECOREGION.ORIG), length(sp_multi), nrow(out), n_distinct(paste(out$LAT, out$LON)),
              min(out$LAT), max(out$LAT)))
  out
}))
sp.df.pooled <- pooled
cat("total pooled populations (species x coast):", nrow(distinct(sp.df.pooled, SPECIES, ECOREGION)),
    " distinct species:", n_distinct(sp.df.pooled$SPECIES), "\n")

save(sp.df.pooled, file = file.path(dir_data, "sp.df.pooled.RData"))
n_sp <- sp.df.pooled |> distinct(SPECIES) |> nrow()          # modskurt_analysis loops over distinct(SPECIES)
writeLines(as.character(seq_len(n_sp)), file.path(PROJ, "analysis/hpc/spID_pooled.txt"))
cat("written:", file.path(dir_data, "sp.df.pooled.RData"), "and analysis/hpc/spID_pooled.txt (", n_sp, "species )\n")
