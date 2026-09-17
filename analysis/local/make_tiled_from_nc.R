# make_tiled_from_nc.R — one-time: build the internally-TILED daily GLORYS tif
# from the NetCDF, so downstream AOI crops are ~10x faster and far lighter on RAM
# (lets you raise NCORES). terra carries the CF time axis into the tif, so
# time()/yr_layer keep working. Run ONCE (reused by all variants).
suppressMessages(require(terra))
source("config.R")

nc_in   <- glorys_daily_nc
# _TILED suffix marks it as the tiled build (distinct from the striped plain copy).
# Same pixel values; only the internal layout differs.
tif_out <- glorys_daily_tiled_tif

r <- rast(nc_in)                 # if it warns of multiple subdatasets: rast(nc_in, subds = "thetao")
cat("layers:", nlyr(r), " | time range:", format(range(time(r))), "\n")

terra::terraOptions(memmax = 8)  # cap RAM while writing; raise if the box is idle
terra::writeRaster(
  r, tif_out, overwrite = TRUE,
  gdal = c("TILED=YES", "BLOCKXSIZE=128", "BLOCKYSIZE=128",
           "COMPRESS=LZW", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"))

# verify time survived (yr_layer depends on it)
chk <- rast(tif_out)
cat("WROTE", basename(tif_out), "| layers:", nlyr(chk),
    "| time OK:", all(is.finite(as.numeric(time(chk)))), "\n")
