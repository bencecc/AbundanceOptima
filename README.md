# AbundanceOptima — code to reproduce *Geographic constraints and unused opportunities limit climate tracking by reef-fish abundance optima*

Reproduces the analyses, figures and tables in *[authors, year, journal]* on climate
tracking by reef-fish abundance optima. Optima are estimated with a Bayesian *modskurt*
framework; figures use *CockR* palettes.

---

## 1. What is and isn't provided
- **Provided (Zenodo [DOI]):** every intermediate `.RData` along the analysis path, so you
  can enter the pipeline at any step and reproduce everything downstream — including all
  figures and tables — **without** re-running the heavy or environmental-data steps.
- **Not provided (supply yourself, only to re-run from scratch):** the large **GLORYS**
  subsurface-temperature grids, **GEBCO** bathymetry, and the **MEOW** marine-ecoregions shapefile
  (`Marine_Ecoregions_Of_the_World__MEOW_.shp` + its `.shx`/`.dbf`/`.prj`; freely available from TNC).
  Point `config.R` (`glorys_dir`, `gebco_file`, `meow_shapefile`) at wherever you store them. GLORYS →
  steps 1–4; GEBCO → steps 4, 7; MEOW → step 2 (ecoregion temperature) and step 11 (Fig 2 map). Their
  outputs are in the Zenodo archive, so those stages are optional.
- **STI is provided, not reproduced:** `reef_fish_sti_glorys.RData` is the canonical input. Its
  generating scripts (stage 3) are included **for transparency only** — they query live
  OBIS/WoRMS/FishBase, so re-running will not reproduce the exact STI.

Survey input spans 1992–2021; the analysed span (with GLORYS exposure) is 1993–2021.

## 2. Setup
1. **Get the code (GitHub).** Clone this repository into your home dir as `~/abundance_optima`:
   ```sh
   git clone https://github.com/bencecc/AbundanceOptima.git ~/abundance_optima
   ```
   (or download the ZIP from the repo's green **Code** button and unzip it to that path). That
   gives you `config.R`, `analysis/`, `figures/`, `tables/` and `packages/` — the `data/`,
   `output/` and `results/` folders come empty (they are excluded from git).
2. **Get the data (Zenodo).** Download the Zenodo archive and unzip it into
   `~/abundance_optima/data/`. This is the large content that does not live in git.
3. (Optional, only to re-run stages 1/4) put GLORYS + GEBCO where `config.R` points.
4. Install the two bundled packages (neither is on CRAN — source tarballs are in `packages/`):
   ```r
   install.packages("packages/modskurt1_1.1.tar.gz", repos = NULL, type = "source")
   install.packages("packages/CockR_1.0.tar.gz",     repos = NULL, type = "source")
   ```
   - **modskurt1** (v1.1) fits the Bayesian *modskurt* abundance-response curves that locate each
     population's abundance optimum along a gradient (latitude, longitude, temperature). It compiles a
     Stan model, so it also needs **cmdstanr** with a working **CmdStan** install (mc-stan.org/cmdstanr).
   - **CockR** (v1.0) supplies the figure colour palettes — cocktail-inspired sequential, discrete and
     diverging scales for ggplot2. Needed for the figures, not the analysis.
5. Open `config.R`, confirm `PROJ = ~/abundance_optima`, and `source("config.R")` at the top
   of any script. All paths flow from `config.R` — nothing is hard-coded.

## 3. Two kinds of script (how they run)
- **Single session** — a plain R script (`Rscript foo.R`), optionally multicore (`foreach`/`doMC`);
  no `.sh`/`.txt`. Runs on a laptop **or** the 128-core single-node server ("thebigone").
- **HPC session** — heavy steps split into chunks by a parameter file (`*.txt`) submitted via a
  shell script (`*.sh`). HPC-produced outputs are reassembled with `cluster_summaries.R`. Templates live in `analysis/hpc/`;
  adapt the cluster-specific header (account, modules, queue) to your system.
  **Most users skip these entirely and load the provided outputs.**

## 4. Layout
```
abundance_optima/
  config.R            # single source of all paths (edit PROJ once)
  data/               # inputs (Zenodo) + GLORYS/GEBCO           [not in git]
  packages/           # modskurt1 + CockR
  analysis/local/     # single-session analysis scripts
  analysis/hpc/       # HPC-session scripts + their .sh / .txt params
  figures/            # Fig 1–4 (each also emits its Extended Data figs via a switch)
  tables/             # Extended Data table scripts
  output/             # intermediate .RData generated along the path  [not in git]
  results/            # FINAL figures & tables (each script's own subfolder) [not in git]
```
Save convention: a script writes into its **own named subfolder under `results/`**
(e.g. `results/Fig1_abund_select/`), created by the script; `config.R` gives the parent
(`dir_results`), the script owns the leaf folder name.

## 5. Quick start (skip-a-step)
- **Reproduce figures/tables only** (most users): fetch Zenodo data → run any script in
  `figures/` or `tables/` → outputs land in `results/…`.
- **Enter earlier**: start at whichever stage below you like; every stage's inputs are provided.
- **Full from scratch**: secure GLORYS + GEBCO, then run stages 1 → 11 in order.

## 6. Analysis path (in order — mirrors Extended Data Fig. 1a)

| # | Stage | Main script(s) | Runs on | Needs GLORYS/GEBCO | Input → output |
|---|---|---|---|---|---|
| 0 | Survey data input | *(assembled upstream)* | — | no | RLS + Reef Check → `sp.df.RData`; occurrences → `fish.occ.df`, `fish_names` |
| 1 | GLORYS tiling → ~1 km working grid | `make_tiled_from_nc.R` | single session | **GLORYS** | GLORYS `.nc` → tiled working grid |
| 2 | Site & ecoregion subsurface temperature | `all_sites_temperature.R`, `all_ecoregions_temperature.R` (+ their `.sh` + `.txt`) | HPC session | **GLORYS + MEOW** | → `temperature_all_sites`, `temperature_all_ecoregions` |
| 3 | Species Temperature Index (STI) — *transparency only* | `fish_env_distr_run.R`, `fish_env_distr.R`, `worms_validate.R` (+ `.sh` + `sti_indices.txt`) | HPC session | **GLORYS** | `fish.occ.df` + `fish_names` → `reef_fish_sti_glorys` — **provided; scripts query live OBIS/WoRMS/FishBase, will not reproduce it exactly** |
| 4 | Barrier identification & ecoregion split (Bassian, Hawaii) | `cluster_ecoregions.R` | single session | **GEBCO** | site coords + GEBCO → `split_plan.RData` (Bassian → `_W/_E/_NW`, Hawaii → `_SE/_NW`) |
| 5 | modskurt optima (lat & lon; abund & density) — **each split group fitted separately** | `modskurt_analysis.R` (+ `modskurt1.sh` + `spID_Bassian_Hawaii_Split.txt`; uses pkg `modskurt1`) | HPC session | no | `sp.df` + `split_plan` → `modskurt.optim.*` |
| 6 | Unimodal-fit QC filter | `unimodal_fit_check.R`, `unimodal_filter_helper.R` | single session | no | `modskurt.optim.*` → `…unimodal` |
| 7 | Relocation (SEA-PATH) | `replace_to_bathy_seapath.R` (via `cluster_summaries.R`) | single session | **GEBCO** | → `optim.*.relocated` |
| 8 | Shift + leadtime + long-term panels | `sp_optimloc_shift(.lead_time).R`, `…longterm_ref_build.R` | HPC session | no | → `sp.optim.*.shift(.leadtime/.longterm)` |
| 9 | Trends: directional + counterfactual scenarios | `optimloc_trend.R`, `optimloc_ecoregion_trend.R`, `optimum_response_scenarios.R` | HPC session | no | → `optim.*.trend`, `scenario.*` |
| 10 | Best-path search (+ sensitivity) | `bestpath_optimum_search.r` (+ `bestpath_greedy_search.r`, `slope_dev_nulldir_thebigone.r`) | single session | no | → `bestpath.opt.*` |
| 11 | Figures & tables | `figures/Fig1_panels.R`, `Fig2.r`, `Fig3.r`, `Fig4.r`; `tables/EDTables_1_and_2.R`, `Analysis/sensitivity_unrelocated.R` | single session | **MEOW** (Fig 2 map) | provided `.RData` → `results/…` |

*(HPC steps are 2, 3, 5, 8, 9 — each keeps its `.sh` and `.txt` together in `analysis/hpc/`; their
outputs are reassembled with `cluster_summaries.R`. HPC-produced `.RData` and the `.txt` parameter
files are **not** in `data/`.)*

*(Each stage lists the keeper script(s); the `_old`/`diag_`/`test_` variants in the source tree are
not released.)*

## 7. Software
R 4.5.2. Two packages are **bundled** in `packages/` (install per §2): **modskurt1** (v1.1) and
**CockR**; everything else is on CRAN. **Pin `glmmTMB 1.1.14` / `TMB 1.9.19`** — model selection is
AIC/P≈0.05-sensitive across versions. Full `sessionInfo()` is in `packages/`.

R packages used across the released scripts (install what your stage needs — a figures-only user
needs far fewer than a from-scratch re-run):

- **Data wrangling** — `tidyverse` (provides `dplyr`, `tidyr`, `ggplot2`, `stringr`, `forcats`,
  `purrr`, `tibble`), `dbplyr`.
- **Bayesian modskurt / Stan** — `modskurt1` (bundled), `cmdstanr` (+ a working **CmdStan**),
  `posterior`, `bayesplot`.
- **Mixed models & stats** — `glmmTMB`, `DHARMa`, `car`, `nnet`, `ggeffects`, `vegan`, `ape`,
  `cluster`.
- **Spatial & mapping** — `sf`, `terra`, `tidyterra`, `ggOceanMaps`, `rnaturalearth`,
  `rnaturalearthdata`, `rnaturalearthhires`; `ggspatial`, `ggnewscale`, `maptiles` (only for the
  optional tile / ocean basemaps in `figures/plot_sites.R`).
- **Plotting extras** — `CockR` (bundled), `patchwork`, `ggh4x`, `ggpmisc`, `scales`.
- **Extended Data tables (Word export)** — `flextable`, `officer`.
- **Parallel / HPC** — `foreach`, `doMC`, `doParallel`, `parallel`.
- **Online data APIs — stage 3 (STI) only, transparency** — `robis`, `rfishbase`, `worrms`
  (these query live services; STI is provided as `reef_fish_sti_glorys.RData`, so most users
  never need them).
