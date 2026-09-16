# ============================================================================
# plot_sites() — map points (sites, or a year-to-year optimum trajectory) over
# an optional land / ocean basemap. A standalone plotting HELPER for the figures;
# it is NOT part of the modskurt1 package (that package's modskurt_plot() draws the
# fitted response curve, a different thing). In this release it is sourced by
# figures/Fig1_panels.R to draw element 3/3 of panel a (the optimum trajectory map).
#
# DEPENDENCIES
#   Packages — always required : sf, ggplot2
#   Packages — only for basemap = :
#       "ggom"        : ggOceanMaps (basemap + bathymetry); ggnewscale (separate
#                       bathy vs point fill scales); ggspatial (nicer point overlay)
#       "osm"/"esri"  : maptiles + ggspatial (downloads map tiles — needs internet)
#       "none"        : nothing beyond sf + ggplot2
#   NB: this file adds mapping deps (ggOceanMaps/maptiles/ggspatial) that the
#   analysis does not otherwise need — install only the branch you use.
#
# DATA / INPUTS (passed in by the caller — this helper loads no files itself)
#   plot.dat : data.frame/tibble of points; MUST have columns LON, LAT.
#              Optional: DATA.TYPE (groups points, sets colour + facet; defaults
#              to "points"), and the title.col column (default "ECOREGION").
#   plot.map : OPTIONAL sf polygons for land (e.g. rnaturalearth::ne_countries(...);
#              high-res land needs the rnaturalearthhires package). If NULL, no land
#              is drawn (basemap = "none"/tiles still work).
#   For basemap = "ggom", ggOceanMaps fetches its own bathymetry data on first use.
#
# RETURNS a ggplot object (facetted by DATA.TYPE when facet = TRUE).
# ============================================================================

plot_sites <- function(
		plot.dat,
		plot.map = NULL,               # sf polygons (e.g., rnaturalearth::ne_countries)
		title.col = "ECOREGION",
		mycols = c("original"="#1b9e77","relocated_sea"="#7570b3","relocated_land"="#d95f02",
				"relocated_sea_warn"="#e7298a","unresolved_sea"="#666666","points"="#2c7fb8"),
		basemap = c("none","osm","esri","ggom"),
		pad.km = 30,
		facet = TRUE,
		# point style
		point.size = 3, alpha = 0.8, shape = 21, col = "black",
		# ggOceanMaps options
		ggom.bathy = FALSE,          # TRUE to add ocean bathymetry
		ggom.bathy.style = NULL,     # e.g. "rcb", "contours", or a style list
		ggom.args = list()           # extra args to ggOceanMaps::basemap()
) {
	basemap <- match.arg(basemap)
	stopifnot(all(c("LON","LAT") %in% names(plot.dat)))
	if (!"DATA.TYPE" %in% names(plot.dat)) plot.dat$DATA.TYPE <- "points"

	# points → sf (EPSG:4326)
	pts <- sf::st_as_sf(plot.dat, coords = c("LON","LAT"), crs = 4326)
	
	# buffered bbox in km (stable across branches)
	buf <- sf::st_transform(pts, 3395) |>
			sf::st_buffer(dist = pad.km * 1000) |>
			sf::st_union() |>
			sf::st_transform(4326)
	bb <- sf::st_bbox(buf)
	x.range <- c(as.numeric(bb["xmin"]), as.numeric(bb["xmax"]))
	y.range <- c(as.numeric(bb["ymin"]), as.numeric(bb["ymax"]))
	
	# crop land if provided (for non-ggom branches)
	if (!is.null(plot.map)) {
		plot.map <- sf::st_make_valid(sf::st_transform(plot.map, 4326))
		plot.map <- suppressWarnings(suppressMessages(sf::st_crop(plot.map, bb)))
	}
	
	#main_title <- if (title.col %in% names(plot.dat)) paste0(unique(plot.dat[[title.col]])[1]) else ""
	main_title <- ""
	
	# -------------------- ggOceanMaps branch (no grid, no title) --------------------
	if (basemap == "ggom") {
		if (!requireNamespace("ggOceanMaps", quietly = TRUE)) {
			stop("basemap='ggom' requires the 'ggOceanMaps' package.")
		}
		has_spatial_point <- requireNamespace("ggspatial", quietly = TRUE)
		
		limits_vec <- c(x.range[1], x.range[2], y.range[1], y.range[2])
		
		# Sanitize user args: drop unnamed args and ones that would collide with 'limits'
		clean_args <- ggom.args
		if (length(clean_args)) {
			if (is.null(names(clean_args))) clean_args <- list()
			if (length(clean_args)) {
				bad <- which(is.na(names(clean_args)) | names(clean_args) == "")
				if (length(bad)) clean_args <- clean_args[-bad]
				# prevent duplicates with core args
				clean_args <- clean_args[!tolower(names(clean_args)) %in% c("limits","bathymetry","bathy","bathy.style")]
			}
		}
		
		# Minimal, safe basemap call (no 'grid' arg; grid removed via theme below)
		base <- try(
				do.call(ggOceanMaps::basemap, c(list(limits = limits_vec), clean_args)),
				silent = TRUE
		)
		if (inherits(base, "try-error")) {
			stop("ggOceanMaps::basemap failed: ", conditionMessage(attr(base, "condition")))
		}
		p <- base
		
		# Bathymetry (prefer geom_bathy; else LargeData-style; else older bathy=TRUE)
		if (isTRUE(ggom.bathy)) {
			has_geom_bathy <- "geom_bathy" %in% getNamespaceExports("ggOceanMaps")
			if (has_geom_bathy) {
				p <- p + if (is.null(ggom.bathy.style)) ggOceanMaps::geom_bathy()
						else ggOceanMaps::geom_bathy(style = ggom.bathy.style)
			} else {
				base2_args <- c(list(limits = limits_vec, bathymetry = TRUE), clean_args)
				if (!is.null(ggom.bathy.style)) base2_args$bathy.style <- ggom.bathy.style
				base2 <- try(do.call(ggOceanMaps::basemap, base2_args), silent = TRUE)
				if (!inherits(base2, "try-error")) {
					p <- base2
				} else {
					base3_args <- c(list(limits = limits_vec, bathy = TRUE), clean_args)
					if (!is.null(ggom.bathy.style)) base3_args$bathy.style <- ggom.bathy.style
					base3 <- try(do.call(ggOceanMaps::basemap, base3_args), silent = TRUE)
					if (!inherits(base3, "try-error")) {
						p <- base3
					} else {
						warning("Bathymetry layer not available for this ggOceanMaps build; proceeding without bathy.")
					}
				}
			}
		}
		
		# If bathymetry is drawn, reset fill scale for points to avoid continuous/discrete clash
		use_newscale <- FALSE
		if (isTRUE(ggom.bathy) && requireNamespace("ggnewscale", quietly = TRUE)) {
			p <- p + ggnewscale::new_scale_fill()
			use_newscale <- TRUE
		}
		
		# Overlay points (prefer fill mapping if we have a fresh fill scale; else map colour)
		if (has_spatial_point) {
			if (use_newscale) {
				p <- p +
						ggspatial::geom_spatial_point(
								data = plot.dat, ggplot2::aes(x = LON, y = LAT, fill = DATA.TYPE),
								shape = shape, color = col, size = point.size, alpha = alpha, crs = 4326
						) +
						ggplot2::scale_fill_manual(values = mycols)
			} else {
				p <- p +
						ggspatial::geom_spatial_point(
								data = plot.dat, ggplot2::aes(x = LON, y = LAT, colour = DATA.TYPE),
								shape = 19, size = point.size, alpha = alpha, crs = 4326
						) +
						ggplot2::scale_colour_manual(values = mycols)
			}
		} else {
			if (use_newscale) {
				p <- p +
						ggplot2::geom_sf(data = pts, ggplot2::aes(fill = DATA.TYPE),
								shape = shape, color = col, size = point.size, alpha = alpha) +
						ggplot2::scale_fill_manual(values = mycols)
			} else {
				p <- p +
						ggplot2::geom_sf(data = pts, ggplot2::aes(colour = DATA.TYPE),
								shape = 19, size = point.size, alpha = alpha) +
						ggplot2::scale_colour_manual(values = mycols)
			}
		}
		
		# Remove grid & title (ggplot panel gridlines; ggom doesn’t add a ggtitle by default here)
		p <- p + ggplot2::ggtitle(main_title) +
				ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
						panel.grid.minor = ggplot2::element_blank()) +
				ggplot2::guides(fill = "none", colour = "none") +
				ggplot2::theme(legend.position = "none")

		if (facet) p <- p + ggplot2::facet_wrap(~DATA.TYPE)
		return(p)
	}
	
	# -------------------- NONE (plain ggplot over sf land) --------------------
	if (basemap == "none") {
		p <- ggplot2::ggplot() +
				{ if (!is.null(plot.map) && nrow(plot.map) > 0)
						ggplot2::geom_sf(data = plot.map, fill = "wheat1", color = "gray30") } +
				ggplot2::geom_sf(data = pts, ggplot2::aes(fill = DATA.TYPE),
						shape = shape, color = col, size = point.size, alpha = alpha) +
				ggplot2::scale_fill_manual(values = mycols) +
				ggplot2::coord_sf(xlim = x.range, ylim = y.range, expand = FALSE) +
				ggplot2::ggtitle(main_title) +
				ggplot2::theme_bw() +
				ggplot2::guides(fill = "none", colour = "none") +
				ggplot2::theme(legend.position = "none") +
				ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
						panel.grid = ggplot2::element_blank())
		if (facet) p <- p + ggplot2::facet_wrap(~DATA.TYPE)
		return(p)
	}
	
	# -------------------- OSM / ESRI tiles (Web Mercator) --------------------
	if (!requireNamespace("maptiles", quietly = TRUE) ||
			!requireNamespace("ggspatial", quietly = TRUE)) {
		stop("basemap='", basemap, "' requires 'maptiles' and 'ggspatial'.")
	}
	provider <- if (basemap == "osm") "OpenStreetMap" else "Esri.WorldImagery"
	
	tiles <- try(maptiles::get_tiles(buf, provider = provider, crop = TRUE), silent = TRUE)
	if (inherits(tiles, "try-error")) {
		warning("Tile download failed (", basemap, "). Falling back to basemap='none'.")
		return(plot_sites(
						plot.dat, plot.map, title.col, mycols, basemap = "none",
						pad.km = pad.km, facet = facet,
						point.size = point.size, alpha = alpha, shape = shape, col = col,
						ggom.bathy = ggom.bathy, ggom.bathy.style = ggom.bathy.style, ggom.args = ggom.args
				))
	}
	
	pts_3857  <- sf::st_transform(pts, 3857)
	land_3857 <- if (!is.null(plot.map) && nrow(plot.map) > 0) sf::st_transform(plot.map, 3857) else NULL
	bb_3857   <- sf::st_bbox(sf::st_transform(buf, 3857))
	xlim_3857 <- c(as.numeric(bb_3857["xmin"]), as.numeric(bb_3857["xmax"]))
	ylim_3857 <- c(as.numeric(bb_3857["ymin"]), as.numeric(bb_3857["ymax"]))
	
	p <- ggplot2::ggplot() +
			ggspatial::annotation_spatial(tiles) +
			{ if (!is.null(land_3857) && nrow(land_3857) > 0)
					ggplot2::geom_sf(data = land_3857, fill = "wheat1", color = "gray30", alpha = 0.35) } +
			ggplot2::geom_sf(data = pts_3857, ggplot2::aes(fill = DATA.TYPE),
					shape = shape, color = col, size = point.size, alpha = alpha) +
			ggplot2::scale_fill_manual(values = mycols) +
			ggplot2::coord_sf(xlim = xlim_3857, ylim = ylim_3857, expand = FALSE, crs = 3857) +
			ggplot2::ggtitle(main_title) +
			ggplot2::theme_bw() +
			ggplot2::guides(fill = "none", colour = "none") +
			ggplot2::theme(legend.position = "none") +
			ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
					panel.grid = ggplot2::element_blank()) +
			
			if (facet) p <- p + ggplot2::facet_wrap(~DATA.TYPE)
	p
}

