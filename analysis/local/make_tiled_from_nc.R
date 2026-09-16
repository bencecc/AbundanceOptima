# make_tiled_from_nc.R — one-time: build the internally-TILED daily GLORYS tif
# from the NetCDF, so AOI crops in slope_dev_nulldir_thebigone.r are ~10x faster
# and far lighter on RAM (lets you raise NCORES). terra carries the CF time axis
# into the tif, so time()/yr_layer keep working. Run ONCE (reused by all variants).
suppressMessages(require(terra))

env_dir <- "/home/lisandro/Lavori/MPA_timeseries/EnvData"
nc_in   <- file.path(env_dir, "glorys_daily_1993_2021_5_10_metres_mean.nc")
# _TILED suffix marks it as the tiled build (distinct from the striped OneDrive copy,
# which keeps the plain name). Same pixel values; only the internal layout differs.
tif_out <- file.path(env_dir, "raster_glorys_daily_1993_2021_5_10_metres_mean_TILED.tif")

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
