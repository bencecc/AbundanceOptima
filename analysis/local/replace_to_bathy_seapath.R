#' Relocates each out-of-range peak site onto reef habitat by a LEAST-COST SEA
#' PATH that avoids land, instead of a straight line. For each out-of-range point:
#'   - it is snapped to the nearest sea cell (so a land-origin point enters the
#'     sea graph at the coast);
#'   - least-cost sea-path distances (land = barrier, any depth passable, so the
#'     point can "escape" through deep water) are computed to every cell;
#'   - among 0-30 m reef cells reachable within `max_dist_km` BY SEA PATH, the
#'     destination preferentially shares the peak's temperature pixel
#'     (|delta_lat| <= lat_band_deg); if none, the minimum forced |delta_lat| is taken;
#'   - if no reef lies within the buffer by sea, the buffer is widened (x2),
#'     then latitude is allowed to shift, so a point always reaches a reachable
#'     reef rather than getting stuck.
#' The 0-30 m band is always the TARGET (destination); deep water is only
#' traversed to reach it. Movement never crosses land (sea-connected only).
#'
#' ESCAPE (resolution artifact): at the ~1 km working grid, reef-dense areas
#' fragment into isolated sea pockets, so a point's origin cell can have NO
#' sea path to ANY band cell. Such points are NOT left stranded at their origin
#' (that put points on land / in deep water). Instead the same destination rule
#' is applied with STRAIGHT-LINE distances (great-circle, land ignored) in place
#' of sea-path distances: candidates are the 0-30 m reef cells within
#' `max_dist_km` as the crow flies (buffer widened x2 if none), preferring the
#' same latitude pixel (|delta_lat| <= lat_band_deg), else the smallest
#' |delta_lat|, then the nearest. The point therefore lands on reef at (nearly)
#' its original latitude, but the route to it may cross land.
#' The sampled sites play only two indirect roles, in the sea-path and the
#' fallback case alike: (i) together with the peak sites they define the
#' rectangular bathymetry window (the "AOI", padded by `aoi_pad_deg`) from
#' which candidate reef cells are taken; (ii) when several reef cells qualify,
#' the choice among them weighs distance against how many sampled sites lie
#' near each cell (`density_weight`, default 0.5), so cells surrounded by
#' sampled sites are favoured. There is no interpolation between sampled
#' sites. The per-call message reports how many points used the fallback.
#'
#' Slow (one Dijkstra / `accCost` per relocated point) but the sea transition is
#' built once per call; intended for a many-core node.
#'
#' Returns a data structure of relocated coordinates + flags. Relocated points get
#'   `status = "to_bathy_seapath"` (whether reached by sea path or the
#'   straight-line fallback); `dist_km` is the sea-path distance (km) for
#'   sea-routed points, or the straight-line distance for fallback points.
#'   `status = "unplaced"` now occurs only if the AOI contains no band cell
#'   at all (essentially never).
replace_to_bathy_seapath <- function(
    modeled_df,
    sample_points,
    bathy_rast,
    target_depth      = 10,
    tol_m             = 1,
    target_range      = c(0, 30),
    sea_is_negative   = TRUE,
    depth_layer       = 1,
    aoi_pad_deg       = 0.75,
    agg_fact_sea      = 2,
    density_weight    = 0.5,
    k_neighbors       = 10,
    lat_band_deg      = 1/12,     # "same temperature pixel" = GLORYS lat-res
    max_dist_km       = 100,      # SEA-PATH search buffer (km)
    precomputed       = NULL
) {
  stopifnot(inherits(bathy_rast, "SpatRaster"))
  stopifnot(is.numeric(lat_band_deg), lat_band_deg > 0)
  stopifnot(is.numeric(max_dist_km), max_dist_km > 0)
  aoi_pad_deg <- max(aoi_pad_deg, max_dist_km / 111 * 1.1)   # crop must contain buffer

  find_coord_cols <- function(df, type = c("lon", "lat")) {
    type <- match.arg(type); nms <- names(df); low <- tolower(nms)
    patterns <- if (type == "lon") c("lon","longitude","x","long") else c("lat","latitude","y")
    for (p in patterns) { idx <- which(low == p); if (length(idx) == 1) return(nms[idx]) }
    NULL
  }

  # ---- inputs (as in the straight-line relocation) -------------------------
  if (inherits(modeled_df, "sf")) {
    P_xy <- sf::st_coordinates(sf::st_transform(modeled_df, 4326))
    attrs <- sf::st_drop_geometry(modeled_df)
  } else {
    modeled_df <- as.data.frame(modeled_df)
    lon_col <- find_coord_cols(modeled_df, "lon"); lat_col <- find_coord_cols(modeled_df, "lat")
    if (is.null(lon_col) || is.null(lat_col))
      stop("modeled_df must have LON/LAT columns. Found: ", paste(names(modeled_df), collapse = ", "))
    P_xy <- as.matrix(modeled_df[, c(lon_col, lat_col)]); colnames(P_xy) <- c("X","Y")
    attrs <- modeled_df; attrs[[lon_col]] <- NULL; attrs[[lat_col]] <- NULL
  }
  n_pts <- nrow(P_xy)

  if (inherits(sample_points, "sf")) {
    sample_xy <- sf::st_coordinates(sf::st_transform(sample_points, 4326))
  } else {
    sample_points <- as.data.frame(sample_points)
    lon_col_s <- find_coord_cols(sample_points, "lon"); lat_col_s <- find_coord_cols(sample_points, "lat")
    if (is.null(lon_col_s) || is.null(lat_col_s))
      stop("sample_points must have LON/LAT columns. Found: ", paste(names(sample_points), collapse = ", "))
    sample_xy <- as.matrix(sample_points[, c(lon_col_s, lat_col_s)])
  }
  sample_xy <- sample_xy[complete.cases(sample_xy), , drop = FALSE]
  use_density <- density_weight > 0 && nrow(sample_xy) >= 3

  # ---- band + ocean masks (band-build as in the straight-line relocation) --
  if (!is.null(precomputed)) {
    rZ <- precomputed$rZ; band_xy <- precomputed$band_xy
    band_density <- precomputed$band_density; bathy_line <- precomputed$bathy_line
    ocean_mask <- precomputed$ocean_mask
  } else {
    r_ll <- if (is.lonlat(bathy_rast)) bathy_rast else project(bathy_rast, "EPSG:4326", method = "bilinear")
    depth_layer <- max(1, min(depth_layer, nlyr(r_ll)))
    all_xy <- rbind(P_xy, sample_xy)
    ext_pad <- ext(min(all_xy[,1]) - aoi_pad_deg, max(all_xy[,1]) + aoi_pad_deg,
                   min(all_xy[,2]) - aoi_pad_deg, max(all_xy[,2]) + aoi_pad_deg)
    rZ0 <- crop(r_ll[[depth_layer]], ext_pad)
    rZ  <- if (agg_fact_sea > 1) aggregate(rZ0, fact = agg_fact_sea, fun = "mean", na.rm = TRUE) else rZ0
    r_abs <- if (sea_is_negative) -rZ else rZ
    r_sea <- if (sea_is_negative) (rZ < 0) else (rZ > 0)
    ocean_mask <- ifel(r_sea, 1, NA)                              # any sea (land-only barrier)
    band_mask <- ifel(!is.na(ocean_mask) &
                        (r_abs >= (target_depth - tol_m)) &
                        (r_abs <= (target_depth + tol_m)), 1, NA)
    bathy_line <- tryCatch(st_as_sf(as.contour(r_abs, levels = target_depth)), error = function(e) NULL)
    band_cells <- which(!is.na(values(band_mask)))
    if (length(band_cells) == 0) {
      P_dep <- terra::extract(rZ, P_xy, method = "simple")[, 1]
      P_abs <- if (sea_is_negative) -P_dep else P_dep
      P_sea <- if (sea_is_negative) P_dep < 0 else P_dep > 0
      ok_by_range <- is.finite(P_abs) & (P_sea %in% TRUE) & P_abs >= min(target_range) & P_abs <= max(target_range)
      out <- cbind(attrs, LON = P_xy[,1], LAT = P_xy[,2],
                   status = ifelse(ok_by_range, "ok", "unplaced"),
                   dist_km = 0, lat_ddeg = 0, lat_dkm = 0,
                   depth_m = as.numeric(P_dep), density_score = NA_real_)
      return(list(points_df = out, bathy_line = bathy_line, precomputed = NULL))
    }
    band_xy <- xyFromCell(band_mask, band_cells)
    if (use_density && requireNamespace("RANN", quietly = TRUE)) {
      k <- min(k_neighbors, nrow(sample_xy))
      nn <- RANN::nn2(sample_xy, band_xy, k = k)
      band_density <- rowMeans(1 / (nn$nn.dists + 0.01))
      band_density <- (band_density - min(band_density)) / (max(band_density) - min(band_density) + 1e-9)
    } else band_density <- rep(1, nrow(band_xy))
  }

  # ---- sea transition surface (built ONCE per call) -----------------------
  sea_r  <- raster::raster(ocean_mask)
  tr_sea <- tryCatch({
    t0 <- gdistance::transition(sea_r, transitionFunction = function(x) mean(x), directions = 8)
    gdistance::geoCorrection(t0, type = "c")
  }, error = function(e) NULL)
  sea_xy <- raster::rasterToPoints(sea_r)[, 1:2, drop = FALSE]
  snap_sea <- function(lon, lat) sea_xy[which.min((sea_xy[,1]-lon)^2 + (sea_xy[,2]-lat)^2), ]

  # ---- which points need relocation (identical) ---------------------------
  P_dep <- terra::extract(rZ, P_xy, method = "simple")[, 1]
  P_abs <- if (sea_is_negative) -P_dep else P_dep
  P_sea <- if (sea_is_negative) P_dep < 0 else P_dep > 0
  ok_by_range <- is.finite(P_abs) & (P_sea %in% TRUE) & P_abs >= min(target_range) & P_abs <= max(target_range)
  need <- which(!ok_by_range)

  out_xy <- P_xy
  status <- rep("ok", n_pts)
  moved_km <- rep(0, n_pts)
  density_at_dest <- rep(NA_real_, n_pts)
  band_lat <- band_xy[, 2]

  if (length(need) > 0 && nrow(band_xy) > 0 && !is.null(tr_sea) && nrow(sea_xy) >= 5) {
    # Parallel over need-points: ONE accCost (Dijkstra) per point. tr_sea / band_xy
    # / band_density live in this function's env and are forked to the workers
    # read-only (copy-on-write) — they are NOT duplicated per worker. Requires a
    # parallel backend registered by the caller (registerDoMC); runs serial (with
    # a foreach warning) if none. Each worker returns
    # c(i, lon, lat, dist_km, density, placed?) so we never write shared state.
    res <- foreach(j = seq_along(need), .combine = rbind) %dopar% {
      i  <- need[j]; plon <- P_xy[i, 1]; pa <- P_xy[i, 2]
      src <- snap_sea(plon, pa)                             # enter sea graph (coast for land pts)
      acc <- tryCatch(gdistance::accCost(tr_sea, src), error = function(e) NULL)
      dvec <- if (!is.null(acc)) raster::extract(acc, band_xy) / 1000 else rep(Inf, nrow(band_xy))
      via_sea <- any(is.finite(dvec))                       # any sea-reachable reef?

      # ESCAPE: at ~1 km the mask fragments reef-dense areas into isolated sea
      # pockets, so a point's origin cell has NO sea path to the 0-30 m band
      # (a resolution artifact, not a real barrier). Rather than leave it
      # stranded at its origin, apply the same destination rule below with
      # straight-line (land-ignoring) distances to the reef band - see header.
      # cos-lat great-circle km:
      if (!via_sea) {
        dvec <- sqrt(((band_xy[, 1] - plon) * cos(pa * pi / 180) * 111.32)^2 +
                     ((band_lat - pa) * 111.32)^2)
      }
      finite_d <- dvec[is.finite(dvec)]
      if (!length(finite_d)) return(c(i, NA, NA, NA, NA, 0)) # truly no reef in AOI
      dlat <- abs(band_lat - pa)

      # (1) reef cells within the buffer (sea-path km, or straight-line km on
      #     fallback); widen x2 ("escape") if none
      rad <- max_dist_km
      repeat {
        in_buf <- which(is.finite(dvec) & dvec <= rad)
        if (length(in_buf) || rad >= max(finite_d)) break
        rad <- rad * 2
      }
      if (!length(in_buf)) in_buf <- which(is.finite(dvec))

      # (2) prefer the same temperature pixel; else the minimum forced |Δlat|
      same_pix <- in_buf[dlat[in_buf] <= lat_band_deg]
      if (length(same_pix)) cand <- same_pix
      else { mind <- min(dlat[in_buf]); cand <- in_buf[dlat[in_buf] <= mind + 1e-9] }

      # (3) shortest path (density tie-break)
      if (use_density && length(cand) > 1) {
        ds <- dvec[cand]; dn <- band_density[cand]; denom <- max(ds)
        dscore <- if (denom > 0) 1 - ds / denom else rep(1, length(ds))
        best <- cand[which.max((1 - density_weight) * dscore + density_weight * dn)]
      } else best <- cand[which.min(dvec[cand])]

      # placed-code: 1 = via true sea path, 2 = via straight-line fallback
      c(i, band_xy[best, 1], band_xy[best, 2], dvec[best], band_density[best],
        if (via_sea) 1 else 2)
    }

    if (!is.null(res)) {
      if (is.null(dim(res))) res <- matrix(res, nrow = 1)
      placed <- res[res[, 6] %in% c(1, 2), , drop = FALSE]   # 1 = sea path, 2 = fallback
      unpl   <- res[res[, 6] == 0, , drop = FALSE]
      if (nrow(placed)) {
        ii <- placed[, 1]
        out_xy[ii, ] <- placed[, 2:3]
        moved_km[ii] <- placed[, 4]
        density_at_dest[ii] <- placed[, 5]
        status[ii] <- "to_bathy_seapath"                     # both land on the band, at sea
      }
      if (nrow(unpl)) status[unpl[, 1]] <- "unplaced"
      n_fb <- sum(res[, 6] == 2)
      if (n_fb > 0)
        message(sprintf("  seapath: %d/%d relocated points used the straight-line fallback (sea path unreachable at working resolution)",
                        n_fb, length(need)))
    }
  }

  # ---- output (as in the straight-line relocation) -------------------------
  lat_ddeg <- out_xy[, 2] - P_xy[, 2]
  lat_dkm  <- abs(lat_ddeg) * 111.32
  depth_fin <- terra::extract(rZ, out_xy, method = "simple")[, 1]
  out <- cbind(attrs, LON = out_xy[,1], LAT = out_xy[,2], status = status,
               dist_km = moved_km, lat_ddeg = lat_ddeg, lat_dkm = lat_dkm,
               depth_m = as.numeric(depth_fin), density_score = density_at_dest)
  rownames(out) <- NULL
  precomputed_out <- list(rZ = rZ, band_xy = band_xy, band_density = band_density,
                          bathy_line = bathy_line, ocean_mask = ocean_mask)
  list(points_df = out, bathy_line = bathy_line, precomputed = precomputed_out)
}
