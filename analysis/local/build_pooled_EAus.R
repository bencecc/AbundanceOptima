# build_pooled_EAus.R
#
# SENSITIVITY to the ecoregion split (co-author query): builds a survey table in which
# species occurring in >= 2 adjacent ecoregions along the continuous eastern Australian
# coast (Torres Strait -> Central/Southern GBR -> Tweed-Moreton -> Manning-Hawkesbury ->
# Cape Howe -> Bassian, eastern group only) are treated as ONE population each, so the
# abundance optimum is fit over the whole ~29-degree gradient instead of within each
# ecoregion. All rows of those species from the chain are relabelled
# ECOREGION = "East Australia pooled"; ECOREGION.ORIG keeps the source ecoregion for the
# later per-species comparison with the main (per-ecoregion) results.
#
# Input : sp.df.RData, split_plan.RData (Bassian site groups)      [config.R paths]
# Output: data/sp.df.pooled.EAus.RData (object sp.df.pooled.EAus, same columns as sp.df
#         + ECOREGION.ORIG) and analysis/hpc/spID_pooledEAus.txt (species ids 1..n).
# Then run modskurt_analysis.R with alternative 2) of its load block (see there), and the
# downstream steps on the pooled outputs.

require(dplyr)
source("config.R")

load(input_file("sp.df.RData")); load(input_file("split_plan.RData"))
chain <- c("Torres Strait Northern Great Barrier Reef", "Central and Southern Great Barrier Reef",
           "Tweed-Moreton", "Manning-Hawkesbury", "Cape Howe", "Bassian")
POOLED <- "East Australia pooled"

# Bassian: keep only its eastern site group (continuous with Cape Howe), from split_plan
bs <- split_plan$sites |> filter(ECOREGION == "Bassian")
stopifnot("E" %in% bs$suffix)
ckey <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
bass_E <- bs |> filter(suffix == "E") |> mutate(k = ckey(LON, LAT)) |> pull(k)

d <- sp.df |> filter(ECOREGION %in% chain) |>
  filter(ECOREGION != "Bassian" | ckey(LON, LAT) %in% bass_E)
sp_multi <- d |> distinct(SPECIES, ECOREGION) |> count(SPECIES) |> filter(n >= 2) |> pull(SPECIES)
sp.df.pooled.EAus <- d |> filter(SPECIES %in% sp_multi) |>
  mutate(ECOREGION.ORIG = ECOREGION, ECOREGION = POOLED)

cat("pooled species:", length(sp_multi), " rows:", nrow(sp.df.pooled.EAus),
    " sites:", n_distinct(paste(sp.df.pooled.EAus$LAT, sp.df.pooled.EAus$LON)),
    " lat:", paste(round(range(sp.df.pooled.EAus$LAT), 1), collapse = " to "), "\n")

save(sp.df.pooled.EAus, file = file.path(dir_data, "sp.df.pooled.EAus.RData"))
n_sp <- sp.df.pooled.EAus |> distinct(SPECIES) |> nrow()
writeLines(as.character(seq_len(n_sp)), file.path(PROJ, "analysis/hpc/spID_pooledEAus.txt"))
