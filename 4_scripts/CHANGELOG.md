# Scripts Changelog

Consolidated history of features and fixes for the SDM pipeline scripts in this folder. Replaces the earlier per-feature README/SUMMARY files.

---

## Clamping (current + future, both clamped and unclamped)

The pipeline generates **both clamped and unclamped predictions** for current and future climate, plus diagnostic difference maps.

- **Clamped** (`*_clamped_*.tif`): environmental values constrained to training range. Use for publication.
- **Unclamped** (`*_unclamped_*.tif`): allows extrapolation beyond training range. Use for uncertainty assessment only.
- **Difference maps** (`clamp_effect_*.tif/png`): unclamped − clamped, diverging blue/white/red palette.
- **3-panel comparison plots** (`*_clamping_comparison_*.png`): clamped | unclamped | difference.

Console output reports % of area affected by clamping. Decision rules:
- <1% affected → trivial
- 1–10% → minor, use clamped
- 10–30% → moderate, use clamped + discuss
- \>30% → model may be inappropriate; consider alternative scenario or expanded training data

References: Elith et al. (2011), Merow et al. (2013), Phillips et al. (2006/2017).

---

## Variable filtering (BIO / SBIO / ENV)

Filtering strategy used by `variable_correlation_filter.R`, `quick_variable_filter.R`, `apply_variable_filter_to_sdm.R`, and integrated into `10_updated_sdmtune_loop.R` via `USE_VARIABLE_FILTER <- TRUE`:

1. **Keep ALL BIO variables** — climate is the primary driver, never removed at this stage.
2. **Remove SBIO/ENV** correlated with any BIO variable (Spearman |r| > 0.7).
3. Among SBIO/ENV, **SBIO has priority over ENV** when they correlate with each other.
4. Final variable selection (among BIO) is done by the SDM pipeline's `varSel` step.

Rationale: SBIO (soil bioclim) is largely a proxy for air temperature and creates redundancy with BIO. BIO variables that survive `varSel` are selected by permutation importance.

Outputs (in `./5_outputs/`): `selected_variable_names.txt`, `selected_variables.xlsx`, `variable_filtering_summary.png`. Filter is cached and reused on subsequent runs.

---

## GABI AntMaps integration

Local GABI AntMaps occurrence data is merged with online sources (GBIF, BISON, ALA, iNat, iDigBio).

- Loaded once at script start from `./2_processed_data/gabi_antmaps_data_clean.csv`
- Per species, matched by **name** AND **GBIF acceptedKey** (dual matching)
- Columns mapped to standard format: `valid_species_name → species`, `dec_long → decimalLongitude`, `dec_lat → decimalLatitude`, `source = "gabi_antmaps"`
- Records with missing coordinates are dropped

### GABI species name matching fix

Original exact-string match failed when loop names included author citations.

- Strip author citations: `gsub("\\s*\\([^)]+\\)\\s*$", "", species_name)`
- Trim whitespace and compare case-insensitively: `tolower(trimws(...))`
- Example: `"Camponotus cruentatus (Latreille, 1802)"` → `"Camponotus cruentatus"` → 853 records matched (was 0)

---

## Ensemble modeling (`ensamble.R` / `ensemble.R`)

Standalone example of multi-algorithm ensemble SDM with SDMtune.

- Trains **Maxnet**, **Random Forest**, and **Boosted Regression Trees** on the same SWD
- Combines predictions via **mean**, **median**, and **AUC-weighted mean**
- Outputs: ensemble + per-model rasters (current, future, change), evaluation table with AUC and weights, comparison plots
- Default test species: *Camponotus herculeanus* (GBIF key 1329074)
- Linear, no loops — intended for stepwise inspection

---

## Timing instrumentation

Per-species and overall pipeline timing.

- `overall_start_time <- Sys.time()` at pipeline start
- `species_start_time` / `species_end_time` per species, elapsed in minutes
- Written into the per-species `summary_*.xlsx`: `start_time`, `end_time`, `processing_time_mins`
- Console reports finish time and elapsed minutes per species

---

## SDMtune best-practice updates applied to the main loop

Applied to `updated_sdmtune_loop.R` / `10_updated_sdmtune_loop.R`:

1. **CV folds** raised from 2 to 4 (more robust evaluation)
2. **Spatial thinning** with `thinData()` — one occurrence per raster cell
3. **`addSamplesToBg()`** — mirrors MaxEnt's `addsamplestobackground=true`
4. **Initial variable importance** exported (Excel + bar chart via `plotVarImp()`)
5. **Correlation matrix** via `corVar()`, correlated pairs (>0.7) saved to Excel
6. **Selected variables** after correlation removal exported with their final importance

---

## Source files this changelog replaces

The following files were consolidated into this changelog and removed from the repo:

- `CLAMPING_IMPLEMENTATION_README.md`
- `ENSEMBLE_README.md`
- `GABI_INTEGRATION_SUMMARY.md`
- `GABI_MATCHING_FIX.md`
- `TIMING_FEATURE_README.md`
- `UPDATES_SUMMARY.md`
- `VARIABLE_FILTERING_CHANGELOG.md`
- `VARIABLE_FILTERING_README.md`
- `VARIABLE_FILTER_QUICK_START.md`
