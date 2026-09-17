# figures/ — main-text and Extended Data figure scripts

**Key convention: each main-figure script also generates the Extended Data (ED)
figures that share its analysis.** Rather than duplicate a script per panel, the
ED variants are produced by flipping a switch at the top of the same script and
re-running — e.g. `DATA_KIND <- "density"` regenerates the density-based ED
counterpart of an abundance main figure. Every script's header carries a
`PRODUCES:` block listing its main + ED outputs; the per-script list below is the
folder-level summary. (Journal terminology "Extended Data" is kept consistent with the
manuscript, which uses ED Table 1/2.)

**The rule (when an ED figure lives in a main-figure script):** an ED figure that
is a *variant or re-aggregation* of a main figure's analysis is embedded in that
figure's script (e.g. ED 4–7 in `Fig4.r`). An ED figure that poses a *distinct
question* would normally get its own script — the one current exception is
**ED Fig 8** (alternative explanations), which stays in `Fig4.r` because its panel a
reuses that script's counterfactual `Response` classification; its header flags it
as a separate question, not a Fig-4 variant.

## What each script produces

One block per script — Main figure, the Extended Data it also makes, and the
Switch(es) that select each output. (Same information as each script's own
`PRODUCES:` header.)

`Fig1_panels.R`

- Main: Fig 1 — single example species (Choerodon fasciatus): modskurt LAT/LON
  optimum + trajectory map (a), latitude vs year (b), observed vs no-shift
  temperature (c).
- Extended Data: none from this script. ED Fig 1 (analytical workflow + modskurt
  candidate-curve schematic) is a separate methods figure.
- Switch: none.

------------------------------------------------------------------------------

`Fig2.r`

- Main: Fig 2 — ecoregion warming map + directional roses + guild×range inset
  (abundance).
- Extended Data: ED Fig 2 — density sensitivity (the same map + roses + inset on
  density optima).
- Switch: `DATA_KIND` = `"abund"` / `"density"`.

------------------------------------------------------------------------------

`Fig3.r`

- Main: Fig 3 — P(poleward shift) ~ previous-year exposure, quadratic binary
  GLMM, faceted guild×range (abundance, 1-yr lag).
- Extended Data: ED Fig 3 — sensitivity, three panels: (a) abundance lag-2,
  (b) density lag-1, (c) density lag-2.
- Switch: `DATA.KIND` = `"abund"` / `"density"`; `LAG` = `1` / `2`.

------------------------------------------------------------------------------

`Fig4.r`

- Main: Fig 4a — counterfactual response by guild×range; Fig 4b — reachable/used
  best-path by guild×range.
- Extended Data:
  - ED Fig 4 & 6 — Fig 4a / 4b by ecoregion.
  - ED Fig 5 & 7 — counterfactual & best-path sensitivity (re-run Fig 4a/4b under
    the alternative `SCENARIO` settings: abund/previous-site, density/first-site,
    density/previous-site).
  - ED Fig 8 — "alternative explanations" (a distinct question, kept here only
    because panel a reuses this script's counterfactual `Response`): 8a STI-tolerance
    falsification (+ STI table), 8b depth-emigration proxy (self-contained).
- Switch: `SCENARIO` (kind abund/density × ref first/previous site × short/long
  term).

## Running a figure and its ED variant

1. Set the switch(es) at the top of the script (see above).
2. Run once per setting — outputs carry a kind/lag suffix so the main and ED
   versions do not overwrite each other.
3. Each run prints the population counts (Ns) quoted in that figure's legend and
   in the Results text, so the numbers in the manuscript can be checked directly
   against the console output.

## Conventions in these scripts

- **`ggsave()` calls are left commented** — uncomment and point at your own output
  folder (a single `dir_figs` once `config.R` path-generalisation is applied).
- **No `windows()` / `dev.off()`** interactive-device calls (Windows-only; removed
  for headless/Mac/Linux portability). Panels are written with `ggsave` / `png` /
  `cairo_pdf`.
- **Section separators use plain `----` / `====`** (keyboard-typeable).
- Palettes are always **CockR** (blue = poleward/cooling, red-orange =
  equatorward/warming); see the top-level `README.md` for setup and data.

## Helper scripts

- **`plot_sites.R`** — a standalone mapping helper (`plot_sites()`) used by
  `Fig1_panels.R` to draw the optimum-trajectory map (element 3/3 of panel a). It is
  *not* part of `modskurt1` (whose `modskurt_plot()` draws the fitted response curve —
  a different thing); it is bundled here and `source()`d beside the figure script. Its
  own header lists the packages (sf + ggplot2 always; ggOceanMaps/maptiles/ggspatial
  only for the ocean/tile basemaps) and the inputs it expects.

## Manual compositing

Some panels (e.g. the Fig 2 roses over the warming map) are exported as separate
layers and assembled in Adobe Illustrator — noted in the relevant script header.
