# SDM Pipeline for Norwegian Ant Species
## Species Distribution Modeling using SDMtune & MaxEnt

---

# TABLE OF CONTENTS

1. [Project Overview](#1-project-overview)
2. [File Structure](#2-file-structure)
3. [Quick Start](#3-quick-start)
4. [Critical Technical Decisions](#4-critical-technical-decisions)
5. [Output Files Explained](#5-output-files-explained)
6. [Common Issues & Solutions](#6-common-issues--solutions)
7. [Key References](#7-key-references)
8. [Testing the Pipeline](#8-testing-the-pipeline)
9. [Extending the Pipeline](#9-extending-the-pipeline)
10. [Troubleshooting Checklist](#10-troubleshooting-checklist)
11. [Version History](#11-version-history)

---

# 1. PROJECT OVERVIEW

Production SDM pipeline for modeling ant species distributions in Norway using MaxEnt (SDMtune R package). Predicts current and future distributions under climate change scenarios.

## Design Principles

- **Occam's Razor**: Simple methods that work > complex methods that fail
- **MaxEnt Best Practices**: Merow et al. 2013, Elith et al. 2011, Radosavljevic & Anderson 2014
- **Minimal Complexity**: Only essential features, no over-engineering
- **Robustness**: Automatic fallbacks prevent pipeline failures

## Active Scripts

- **`10_updated_sdmtune_loop_merged.R`**: THE production pipeline — the only SDMtune script that counts. All other variants (`10_updated_sdmtune_loop_merged.R`, `updated_sdmtune_loop_2.R`, etc.) are obsolete and live in `4_scripts/Gamle skript/` for reference only.

---

# 2. FILE STRUCTURE

```
data/
├── 1_raw_data/
│   └── bioclim/
│       ├── current/        # Current climate (BIO1-BIO19)
│       ├── future/         # Future climate SSP585 2021-2040 (BIO1-BIO19)
│       ├── sbio/           # Soil bioclim layers
│       └── env/            # Additional environmental layers
├── 2_processed_data/
│   ├── complete_ant_data.xlsx
│   ├── gabi_antmaps_data_clean.csv
│   └── Arter_Norge.xlsx
├── 4_scripts/
│   ├── 10_updated_sdmtune_loop_merged.R   # PRODUCTION SCRIPT (only one that counts)
│   ├── Gamle skript/                       # Obsolete script versions
│   └── rename_bioclim_layers.R             # Fix internal raster names
├── 5_outputs/
│   └── screen_output/
│       └── prescreening_bio6_spocc.xlsx
├── species/
│   └── [speciesKey]_[species_name]/
│       └── SDM_maxnet/
│           ├── rasters/
│           │   ├── current_clamped_[key].tif
│           │   ├── future_clamped_[key].tif
│           │   ├── change_clamped_[key].tif
│           │   ├── uncertainty_[key].tif
│           │   └── mess_[key].tif
│           ├── tables/
│           │   ├── summary_[key].xlsx
│           │   ├── evaluation_[key].xlsx
│           │   ├── variable_importance_[key].xlsx
│           │   └── jackknife_[key].xlsx
│           └── [multiple visualization PNGs]
└── CLAUDE.md                         # This file
```

---

# 3. QUICK START

## 3.1 Prerequisites

```r
install.packages(c("spocc", "rgbif", "tidyverse", "CoordinateCleaner",
                   "SDMtune", "terra", "sf", "ggplot2", "rio", "rnaturalearth",
                   "raster", "dismo", "patchwork", "viridis", "tidyterra"))
```

## 3.2 Fix Raster Internal Names (One-time setup)

**Problem**: Bioclim .tif files may have renamed filenames but old internal metadata names, causing layer mismatches.

**Solution**:
```r
source("./4_scripts/rename_bioclim_layers.R")
```

Processes: current/, future/, sbio/, env/ folders
Uses temp file workflow to avoid "source and target cannot be the same" error

**Verify success**:
```r
library(terra)
r <- rast("./1_raw_data/bioclim/current/BIO1_Annual_Mean_Temperature.tif")
names(r)  # Should output: "BIO1_Annual_Mean_Temperature"
```

## 3.3 Run SDM Pipeline

```r
source("./4_scripts/10_updated_sdmtune_loop_merged.R")
```

**Configuration** (lines 13-18):
```r
MIN_RECORDS <- 50             # Minimum presence points required
DOWNLOAD_LIMIT <- 100000      # Max records per data source
CV_FOLDS <- 5                 # Cross-validation folds
RECENT_YEARS <- 1960          # Temporal filter (post-1960 records)
BUFFER_DISTANCE <- 3000000    # 3000 km buffer for background sampling
USE_VARIABLE_FILTER <- TRUE   # Pre-filter correlated SBIO/ENV layers
```

**Test subset** (comment out for full run):
```r
# test_species <- test_species[1:3, ]  # Test with 3 species first
```

---

# 4. CRITICAL TECHNICAL DECISIONS

## 4.1 Background Point Selection: Circle-Based with Fallback

### Primary Method: Biogeographically-Informed Circles

**Rationale**:
- Convex hull creates inappropriate background (temperate zones for tropical species)
- Caused instant overfitting for widespread species (*Monomorium floricola*)
- Circle method restricts background to biogeographically plausible regions

**Implementation** (`10_updated_sdmtune_loop_merged.R` lines 529-571):
1. Create 3000km buffer around each occurrence point (projected EPSG:3857)
2. Merge buffers with S2 disabled (prevents spherical geometry edge crossing errors)
3. Simplify geometry (1km tolerance, negligible at 3000km scale)
4. Transform to WGS84 while S2 still disabled
5. Sample 20,000 random points within merged buffer
6. Filter to land areas only

**S2 Spherical Geometry Handling**:
- S2 uses stricter validity rules for spherical geometry (enabled by default in sf ≥1.0)
- Overlapping buffers create edge intersections that violate S2 validity
- Solution: Disable S2 during union/transform, re-enable after geometry operations complete
- Reference: [r-spatial/sf Issue #1710](https://github.com/r-spatial/sf/issues/1710)

### Fallback Method: Random Global Background

**Trigger**: Any error during circle-based selection (geometry errors, S2 failures, etc.)

**Implementation**:
- Automatically switches to random global background from entire land surface
- Samples 20,000 points from world land polygons
- Less biogeographically informed but prevents pipeline failure
- Console messages clearly indicate which method was used

**Justification**: Robustness is critical for production pipelines processing hundreds of species unattended.

---

## 4.2 MaxEnt Hyperparameters: Conservative Approach

### Feature Classes: H, Q, HQ Only

**Production script** (`10_updated_sdmtune_loop_merged.R` lines 733-734):
```r
h <- list(reg = seq(2, 10, 2),    # 5 values: 2, 4, 6, 8, 10
          fc = c("h", "q", "hq"))  # Hinge, Quadratic, Hinge+Quadratic
```

**Rationale**:
- EPPO best practices: MaxEnt 3.4.0+ omits threshold features by default
- Hinge features provide smoother, ecologically interpretable response curves
- Quadratic features capture classic hump-shaped ecological responses
- Product features avoided (risk of overfitting with small sample sizes)
- Literature: Only 3.7% of SDM studies properly evaluate feature classes (Merow et al. 2013)

### Regularization: 2-10 by Steps of 2

**Rationale**:
- Radosavljevic & Anderson 2014: Higher regularization prevents overfitting in small samples
- Conservative range appropriate for 50-1000 occurrence records
- Genetic algorithm (pop=10, gen=10) searches 15 model configurations efficiently

### Fallback Grid (if insufficient models)

Expanded grid activates if optimization fails:
```r
reg = seq(0.5, 5, 0.5)  # 10 values
fc = c("l", "q", "h", "lq", "lh", "qh", "lqp", "lqh", "lph", "qph", "lqph")
```

Final fallback: gridSearch (exhaustive testing)

---

## 4.3 Variable Selection: Two-Stage Filter

### Stage 1: Permutation Importance (Top 15)

Reduces computational load by pre-filtering to most important variables.

### Stage 2: Correlation-Based Removal

**Primary threshold**: Spearman ρ > 0.7
**Fallback threshold**: ρ > 0.85 (if primary fails)
**Permutations**: 10 (primary), 5 (fallback)

**Rationale**:
- Removes redundant variables that inflate model complexity
- Spearman correlation robust to non-linear relationships
- 0.7 threshold standard in SDM literature (Dormann et al. 2013)
- Fallback ensures pipeline continues even with highly correlated predictor sets

**If both fail**: Species skipped (appropriate behavior - data quality issue)

---

## 4.4 Data Cleaning: CoordinateCleaner + Spatial Thinning

### CoordinateCleaner Pipeline

Sequential filters:
1. `cc_val`: Invalid coordinates (lat/lon out of range)
2. `cc_equ`: Equal lat/lon (often data entry errors)
3. `cc_cap`: Country capitals (museum/herbarium locations)
4. `cc_cen`: Country/province centroids (coarse georeferences)
5. `cc_sea`: Sea coordinates (land species only)
6. `cc_zero`: Null Island (0°N, 0°E - common error)
7. `cc_dupl`: Duplicate records
8. `cc_outl`: Geographic outliers (quantile method, 3× multiplier)

**Fallback**: If full cleaning fails, minimal cleaning (zeros + duplicates only)

### Spatial Thinning

**Method**: SDMtune::thinData (one point per raster cell)
**Rationale**: Removes spatial autocorrelation (violates independence assumption in MaxEnt)

---

## 4.5 Clamping: Both Clamped and Unclamped Predictions

### Clamped Predictions (Default for Publication)

Restricts environmental values to training range (prevents unrealistic extrapolation).

### Unclamped Predictions (Exploratory)

Allows extrapolation beyond training data (useful for novel climates under climate change).

### Outputs

- `current_clamped_[key].tif` / `current_unclamped_[key].tif`
- `future_clamped_[key].tif` / `future_unclamped_[key].tif`
- `clamp_effect_current_[key].tif` (difference map)
- `clamp_effect_future_[key].tif` (difference map)

**Justification**: Allows assessment of model extrapolation risk (critical for climate change projections).

---

## 4.6 MESS Analysis (Multivariate Environmental Similarity Surfaces)

**Method**: dismo::mess (Elith et al. 2010)
**Output**: `mess_[key].tif`

**Interpretation**:
- Negative values = extrapolation (novel climate, predictions less reliable)
- Positive values = interpolation (within training range, predictions reliable)

**Use case**: Identifies areas where future climate predictions are most uncertain.

---

# 5. OUTPUT FILES EXPLAINED

## 5.1 Rasters (`.tif` files in `rasters/`)

| File | Description | Values | Use |
|------|-------------|--------|-----|
| `current_clamped_[key].tif` | Current climate suitability (clamped) | 0-1 | Publication-ready |
| `current_unclamped_[key].tif` | Current climate suitability (unclamped) | 0-1 | Exploratory |
| `future_clamped_[key].tif` | Future climate suitability SSP585 2021-2040 (clamped) | 0-1 | Publication-ready |
| `future_unclamped_[key].tif` | Future climate suitability (unclamped) | 0-1 | Exploratory |
| `change_clamped_[key].tif` | Climate change impact (future - current, clamped) | -1 to +1 | Gain/loss analysis |
| `change_unclamped_[key].tif` | Climate change impact (unclamped) | -1 to +1 | Exploratory |
| `clamp_effect_current_[key].tif` | Clamping effect (unclamped - clamped) | -1 to +1 | Extrapolation assessment |
| `clamp_effect_future_[key].tif` | Clamping effect future | -1 to +1 | Extrapolation assessment |
| `uncertainty_[key].tif` | Prediction uncertainty (SD across CV folds) | ≥0 | Confidence assessment |
| `mess_[key].tif` | MESS analysis (extrapolation risk) | -100 to +100 | Novel climate detection |

**Projection**: WGS84 (EPSG:4326)
**Format**: GeoTIFF, LZW compression

---

## 5.2 Tables (`.xlsx` files in `tables/`)

### `summary_[key].xlsx`

Comprehensive model metadata:
- Species name, GBIF key, processing timestamp
- Data counts: raw, cleaned, final presence/background
- Data sources breakdown (GBIF, iNaturalist, ALA, GABI AntMaps)
- Variables: count, names
- Evaluation: test AUC/TSS, CV AUC/TSS
- Threshold: optimal (max TSS), sensitivity, specificity
- Confusion matrix: TP, FP, TN, FN
- Hyperparameters: feature class, regularization, CV folds
- Climate change: mean change, gain %, loss %

### `evaluation_[key].xlsx`

Model performance metrics:
- AUC (Area Under ROC Curve): 0-1, >0.7 acceptable, >0.8 good
- TSS (True Skill Statistic): -1 to +1, >0.4 acceptable, >0.6 good

### `initial_variable_importance_[key].xlsx`

Permutation importance before correlation filtering (all variables).

### `final_variable_importance_[key].xlsx`

Permutation importance after correlation filtering (selected variables only).

### `jackknife_[key].xlsx`

Variable contribution analysis:
- With only: AUC using only this variable
- Without: AUC excluding this variable
- Identifies dominant and redundant variables

### `optimization_results_[key].xlsx`

Hyperparameter tuning results:
- Tested configurations (fc, reg combinations)
- Train AUC, test AUC, difference AUC
- Sorted by test AUC (best first)

### `correlated_variables_[key].xlsx`

Pairs of correlated variables (Spearman ρ > 0.7) removed during selection.

### `mess_summary_[key].xlsx`

MESS analysis summary:
- Mean MESS value
- % extrapolation (MESS < 0)
- % interpolation (MESS ≥ 0)

---

## 5.3 Visualizations (`.png` files in `SDM_maxnet/`)

### Suitability Maps (Global)

- `current_clamped_[key].png`
- `current_unclamped_[key].png`
- `future_clamped_[key].png`
- `future_unclamped_[key].png`

### Suitability Maps (Europe Zoom)

- `current_europe_clamped_[key].png`
- `current_europe_unclamped_[key].png`
- `future_europe_clamped_[key].png`
- `future_europe_unclamped_[key].png`

### Analysis Maps

- `climate_change_impact_clamped_[key].png`: Gain/loss under climate change
- `uncertainty_[key].png`: Prediction uncertainty (SD across CV folds)
- `mess_analysis_[key].png`: Extrapolation risk map
- `clamp_effect_current_[key].png`: Clamping effect current
- `clamp_effect_future_[key].png`: Clamping effect future

### Clamping Comparison (3-Panel)

- `current_clamping_comparison_[key].png`: Clamped | Unclamped | Difference
- `future_clamping_comparison_[key].png`: Clamped | Unclamped | Difference

### Threshold Maps (Binary Presence/Absence)

- `th_current_binary_maxTSS_clamped.png`: Current distribution (max TSS threshold)
- `th_future_binary_maxTSS_clamped.png`: Future distribution (max TSS threshold)

### Data Inspection

- `data_inspection_[key].png`: Occurrence points + background points
- `distribution_[key].png`: Occurrence points only

### Variable Importance

- `initial_variable_importance_[key].png`: Before correlation filtering
- `final_variable_importance_[key].png`: After correlation filtering
- `jackknife_plot_[key].png`: Variable contributions

### Response Curves (in `response_curves/` subdirectory)

- `response_[varname].png`: Standard response curve (all data)
- `marginal_response_[varname].png`: Marginal response (presence only)

### Model Evaluation

- `roc_curve_[key].png`: ROC curve (AUC visualization)

---

# 6. COMMON ISSUES & SOLUTIONS

## 6.1 "Error: [rast] file does not exist"

**Cause**: Bioclim files missing or incorrect paths

**Solution**: Verify files exist
```r
list.files("./1_raw_data/bioclim/current", pattern = "\\.tif$")
list.files("./1_raw_data/bioclim/future", pattern = "\\.tif$")
```

---

## 6.2 "Error in prepareSWD ... different extents"

**Cause**: Internal raster names don't match between current/future stacks

**Solution**: Run `rename_bioclim_layers.R`

**Verify**:
```r
names(current_predictors)
names(future_predictors)
# Should be identical (order matters!)
```

---

## 6.3 "Insufficient data after cleaning: [N]"

**Cause**: Too few occurrence records after quality filters

**Options**:
1. Lower `MIN_RECORDS` threshold
2. Skip temporal filter (comment temporal filtering section)
3. Accept species will be skipped (appropriate for low-quality data)

---

## 6.4 "Variable selection failed completely"

**Cause**: High correlation among all variables OR too few background points

**Solution**: Automatic fallback to `cor_th = 0.85` (already implemented)

If still fails: Species skipped (appropriate - indicates data quality issue)

---

## 6.5 "Loop 0 is not valid: Edge 113 crosses edge 474"

**Cause**: S2 spherical geometry edge crossing error during circle union

**Solution**: Already implemented - automatic fallback to random global background

**Console indicates**:
```
Circle-based background selection FAILED: Loop 0 is not valid...
Falling back to random global background points...
Random global background selection successful: 18756 points
```

**No action needed** - pipeline continues automatically.

---

## 6.6 Background Points in Wrong Regions

**Cause**: Modified script to use convex hull instead of circles

**Solution**: Use circle-based method (`10_updated_sdmtune_loop_merged.R` lines 529-571)

**Do NOT use**: Convex hull (creates unrealistic background for widespread species)

---

# 7. KEY REFERENCES

## 7.1 MaxEnt Best Practices

**Merow, C., Smith, M.J., & Silander, J.A. (2013)**
*A practical guide to MaxEnt for modeling species' distributions: what it does, and why inputs and settings matter*
Ecography, 36(10), 1058-1069
Feature class selection, regularization importance

**Elith, J., Phillips, S.J., et al. (2011)**
*A statistical explanation of MaxEnt for ecologists*
Diversity and Distributions, 17(1), 43-57
MaxEnt algorithm explanation, output interpretation

**Radosavljevic, A. & Anderson, R.P. (2014)**
*Making better Maxent models of species distributions: complexity, overfitting and evaluation*
Journal of Biogeography, 41(4), 629-643
Regularization guidelines for small samples, model evaluation

**Phillips, S.J., Anderson, R.P., Dudík, M., Schapire, R.E., & Blair, M.E. (2017)**
*Opening the black box: an open-source release of Maxent*
Ecography, 40(7), 887-893
MaxEnt 3.4.0+ updates, threshold feature removal

---

## 7.2 Background Selection

**Barbet-Massin, M., Jiguet, F., Albert, C.H., & Thuiller, W. (2012)**
*Selecting pseudo-absences for species distribution models: how, where and how many?*
Methods in Ecology and Evolution, 3(2), 327-338
Background point number, sampling strategies

**VanDerWal, J., Shoo, L.P., Graham, C., & Williams, S.E. (2009)**
*Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know?*
Ecological Modelling, 220(4), 589-594
Geographic extent of background sampling

---

## 7.3 Model Evaluation

**Lobo, J.M., Jiménez-Valverde, A., & Real, R. (2008)**
*AUC: a misleading measure of the performance of predictive distribution models*
Global Ecology and Biogeography, 17(2), 145-151
AUC limitations, TSS as alternative

**Allouche, O., Tsoar, A., & Kadmon, R. (2006)**
*Assessing the accuracy of species distribution models: prevalence, kappa and the true skill statistic (TSS)*
Journal of Applied Ecology, 43(6), 1223-1232
TSS methodology

---

## 7.4 MESS Analysis

**Elith, J., Kearney, M., & Phillips, S. (2010)**
*The art of modelling range-shifting species*
Methods in Ecology and Evolution, 1(4), 330-342
MESS analysis methodology, extrapolation detection

---

## 7.5 Coordinate Cleaning

**Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., et al. (2019)**
*CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases*
Methods in Ecology and Evolution, 10(5), 744-751
CoordinateCleaner package methodology

---

## 7.6 Software

**Vignali, S., Barras, A.G., Arlettaz, R., & Braunisch, V. (2020)**
*SDMtune: An R package to tune and evaluate species distribution models*
Ecology and Evolution, 10(20), 11488-11506
SDMtune package

---

# 8. TESTING THE PIPELINE

## 8.1 Quick Test (3 Species)

Keep test subset line:
```r
test_species <- test_species[1:3, ]
```

**Expected runtime**: ~30-60 minutes (10-20 min per species)

**Verify**:
1. Console shows progress for 3 species
2. `species/` folder contains 3 species folders
3. Each has `SDM_maxnet/rasters/` with 10+ .tif files
4. Each has `SDM_maxnet/tables/` with 8+ .xlsx files
5. Each has `SDM_maxnet/` with 20+ .png visualizations
6. `summary_[key].xlsx` shows AUC > 0.7, TSS > 0.4 (acceptable models)

---

## 8.2 Full Run (All Species)

Comment out test subset:
```r
# test_species <- test_species[1:3, ]
```

**Expected runtime**: Varies by species count (~10-20 min per species)

**Monitoring**:
- Console messages show progress
- Each species completion shows: AUC, TSS, processing time
- Failures logged with reason (insufficient data, optimization failed, etc.)

**Expected failures**: 5-10% (normal - insufficient data, extreme distributions, etc.)

---

# 9. EXTENDING THE PIPELINE

## 9.1 Change Climate Scenario

Replace future bioclim files in `./1_raw_data/bioclim/future/` with:
- SSP245 (moderate emissions)
- SSP370 (high emissions)
- Different time periods (2041-2060, 2061-2080, 2081-2100)

Pipeline automatically uses whichever files are present.

---

## 9.2 Add Environmental Layers

1. Add .tif files to `./1_raw_data/bioclim/env/` or `sbio/`
2. Ensure filenames match internal raster names (use `rename_bioclim_layers.R`)
3. Set `USE_VARIABLE_FILTER = TRUE` to pre-filter correlated variables
4. Pipeline automatically includes new layers in modeling

---

## 9.3 Modify Species List

**Filter by geographic region**:
```r
# Only species occurring in specific countries
test_species <- test_species %>%
  filter(acceptedKey %in% c(1234567, 2345678, ...))
```

**Filter by taxonomy**:
```r
# Only specific genera
test_species <- test_species %>%
  filter(grepl("^Camponotus", name))
```

---

## 9.4 Adjust Background Sampling

**Change buffer distance** (line 18):
```r
BUFFER_DISTANCE <- 5000000  # 5000 km (wider sampling)
BUFFER_DISTANCE <- 1000000  # 1000 km (restricted sampling)
```

**Rationale**: Wider buffers for widespread species, narrower for restricted endemics.

---

# 10. TROUBLESHOOTING CHECKLIST

## 10.1 Before Running Pipeline

- [ ] All required packages installed (see section 3.1)
- [ ] Bioclim files exist in correct folders (current/, future/, sbio/, env/)
- [ ] Ran `rename_bioclim_layers.R` to fix internal names
- [ ] Verified raster names match: `names(terra::rast("path/to/BIO1.tif"))`
- [ ] Master data files exist: `complete_ant_data.xlsx`, `prescreening_bio6_spocc.xlsx`, `Arter_Norge.xlsx`
- [ ] GABI data exists: `gabi_antmaps_data_clean.csv` (optional but recommended)
- [ ] Working directory is `.../data/`: `getwd()`
- [ ] Configuration parameters set (MIN_RECORDS, BUFFER_DISTANCE, etc.)

---

## 10.2 If Species Fails

**Check console messages** - which step failed?

| Error Message | Meaning | Action |
|---------------|---------|--------|
| "Insufficient data after cleaning: [N]" | Too few records post-QC | Normal - skip species |
| "Variable selection failed completely" | All variables correlated | Normal - skip species |
| "SWD preparation failed" | Raster name mismatch | Run `rename_bioclim_layers.R` |
| "Circle-based background FAILED" | Geometry error | Automatic fallback (no action) |
| "Optimization failed: lewer than population" | Too few models | Automatic fallback (no action) |

**Most failures are expected** - data quality issues, extreme distributions, etc.

---

## 10.3 If Entire Script Crashes

1. **Check file paths** (Windows: backslashes → forward slashes in R)
2. **Verify working directory**: `getwd()` should end in `.../data/`
3. **Check cloud sync** (pause Proton Drive if files are locked)
4. **Run smaller batch**: `test_species <- test_species[1:5, ]`
5. **Check RAM** (large rasters require 8-16 GB)
6. **Update packages**: `update.packages(ask = FALSE)`

---

## 10.4 Performance Issues

**Slow species processing (>30 min per species)**:
- High occurrence density (>1000 records after thinning)
- Many environmental layers (>50 variables)
- Large raster extent/resolution

**Solutions**:
- Increase spatial thinning (modify `thinData` parameters)
- Pre-filter variables (`USE_VARIABLE_FILTER = TRUE`)
- Reduce raster resolution (resample environmental layers)

---

# 11. VERSION HISTORY

## v3.0 (2026-01-16) - Production Version

**Major changes**:
- Added robust background selection with automatic fallback (circle → random global)
- Fixed S2 spherical geometry edge crossing errors (Issue #1710)
- Simplified geometry with 1km tolerance before transform
- Automatic retry logic for hyperparameter optimization
- Comprehensive MESS analysis for extrapolation detection
- Clamped and unclamped predictions (both current and future)
- Full visualization suite (20+ plots per species)
- Processing time tracking per species

**Active script**:
- `10_updated_sdmtune_loop_merged.R` — the only SDMtune script that counts. All earlier variants are obsolete (kept under `4_scripts/Gamle skript/`).

**Status**: Production-ready, tested on 328 non-Norwegian ant species

---

## v2.0 (2025-12-06) - Simplified Version

**Changes**:
- Removed convex hull (replaced with circles)
- Single hyperparameter grid (LQH + reg 1-5)
- Removed widespread species detection
- Removed complex fallback chains
- Added comprehensive documentation

---

## v1.0 - Original Complex Version

**Features**:
- Multiple hyperparameter scenarios
- Convex hull background selection (problematic)
- Widespread species environmental filtering
- Many fallback mechanisms

**Issues**: Overfitting, complex failures, deprecated dependencies

---

**Last Updated**: 2026-01-16
**Active Script**: `10_updated_sdmtune_loop_merged.R`
**Status**: Production-ready, robust failsafe mechanisms
