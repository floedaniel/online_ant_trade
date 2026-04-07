# Variable Correlation Filtering

## Overview

This workflow tests correlations between environmental variables (BIO, SBIO, ENV) and automatically filters them based on this strategy:

**Filtering Strategy**:
1. **KEEP ALL BIO variables** (climate is primary driver - never removed)
2. **REMOVE SBIO/ENV variables** that are correlated with ANY BIO variable
3. **Among SBIO/ENV**: Apply priority SBIO > ENV when correlated with each other
4. **Further variable selection** happens in the SDM pipeline (varSel step)

**Rationale**: BIO variables (climate) are the primary drivers of species distributions. Correlated SBIO (soil climate) variables are redundant. The SDM pipeline will further select the most important BIO variables during modeling.

---

## Files

### 1. `variable_correlation_filter.R`
**Main analysis script** - analyzes correlations and creates filtered variable list

### 2. `apply_variable_filter_to_sdm.R`
**Helper functions** - load filtered variables for use in SDM pipeline

---

## Quick Start

### Step 1: Run Correlation Analysis

```r
source("./4_scripts/variable_correlation_filter.R")
```

**What it does:**
- Loads all environmental layers (BIO, SBIO, ENV)
- Samples 10,000 random points
- Calculates Spearman correlation matrix
- Identifies correlated pairs (|r| > 0.7)
- Removes lower-priority variables
- Creates visualizations and tables

**Outputs** (saved to `./5_outputs/`):
- `selected_variables.xlsx` - Variables to keep
- `removed_variables.xlsx` - Variables removed (with reasons)
- `correlation_heatmap_full.png` - All variables
- `correlation_heatmap_selected.png` - Filtered variables only
- `selected_variable_names.txt` - Simple list for scripting

### Step 2: Use Filtered Variables in SDM

#### Option A: Modify existing script

In `10_updated_sdmtune_loop.R`, replace lines 120-177 with:

```r
# Load filtered environmental data
source("./4_scripts/apply_variable_filter_to_sdm.R")
predictors <- load_filtered_predictors(use_filter = TRUE)
current_predictors <- predictors$current
future_predictors <- predictors$future
```

#### Option B: Manual loading

```r
# Read selected variable names
selected_vars <- readLines("./5_outputs/selected_variable_names.txt")

# Load all layers
current_predictors <- terra::rast(list.files("./1_raw_data/bioclim/current", full.names = TRUE))
# ... (load SBIO and ENV as needed)

# Filter to selected variables
current_predictors <- current_predictors[[selected_vars]]
```

---

## Configuration

Edit **line 13-14** in `variable_correlation_filter.R`:

```r
CORRELATION_THRESHOLD <- 0.7  # Lower = stricter filtering
SAMPLE_SIZE <- 10000          # More points = better estimate (slower)
```

### Recommended Thresholds:

- **0.7** (default) - Standard for ecological modeling
- **0.8** - More permissive, keeps more variables
- **0.9** - Very permissive, only removes extreme correlations
- **0.6** - Strict, aggressive filtering

---

## How Filtering Works

### Rule 1: BIO vs BIO - KEEP BOTH

```
BIO5 (Max Temp Warmest Month) ⟷ 0.85 ⟷ BIO10 (Mean Temp Warmest Quarter)
→ Keep BOTH (all BIO variables are kept regardless of correlation)
→ SDM pipeline will select most important via varSel
```

### Rule 2: BIO vs SBIO/ENV - ALWAYS KEEP BIO

```
BIO1 (Air Temperature) ⟷ 0.91 ⟷ SBIO1_0_5cm (Soil Temperature)
→ Keep BIO1, Remove SBIO1_0_5cm (SBIO redundant with BIO)
```

### Rule 3: SBIO vs ENV - KEEP SBIO

```
SBIO2 (Soil Diurnal Range) ⟷ 0.82 ⟷ clay_0_5cm (Clay content)
→ Keep SBIO2, Remove clay_0_5cm (SBIO has priority)
```

### Rule 4: Same Category (both SBIO or both ENV)

```
sand_0_5cm ⟷ -0.87 ⟷ clay_0_5cm
→ Keep clay_0_5cm (alphabetically first when same category)
```

---

## Understanding Outputs

### `selected_variables.xlsx`

| variable | category | priority |
|----------|----------|----------|
| BIO1_Annual_Mean_Temperature | BIO | 1 |
| BIO12_Annual_Precipitation | BIO | 1 |
| SBIO6_0_5cm_MinT_coldestMonth | SBIO | 2 |
| sand_0_5cm_mean_30s | ENV | 3 |

**Use these variables** in your SDM models.

### `removed_variables.xlsx`

| variable | category | priority | reason |
|----------|----------|----------|--------|
| SBIO1_0_5cm_Annual_Mean_Temperature | SBIO | 2 | Correlated with BIO1_Annual_Mean_Temperature (r=0.91); BIO1 has higher priority |
| clay_0_5cm_mean_30s | ENV | 3 | Correlated with sand_0_5cm_mean_30s (r=-0.87); sand has higher priority |

**These were removed** due to correlation with higher-priority variables.

### `high_correlation_pairs.xlsx`

Shows **all pairs** with |r| > threshold, sorted by correlation strength.

Use this to understand **why** variables were removed.

### Correlation Heatmaps

- **`correlation_heatmap_full.png`**: All variables before filtering
- **`correlation_heatmap_selected.png`**: Only selected variables

**Blue** = negative correlation
**White** = no correlation
**Red** = positive correlation

---

## Validation

After filtering, check:

```r
# Load selected correlation matrix
selected_cor <- rio::import("./5_outputs/selected_correlation_matrix.xlsx")

# Maximum absolute correlation
max_cor <- max(abs(selected_cor[upper.tri(selected_cor)]))
print(max_cor)
```

✅ Should be **< threshold** (e.g., < 0.7)

If not, some correlation chains remain. You can:
1. **Lower threshold** and re-run
2. **Manually inspect** `remaining_high_correlations.xlsx`
3. **Ignore** if chains are ecologically meaningful

---

## Integration with Existing Workflow

### Before filtering:
```
19 BIO variables  (climate)
22 SBIO variables (soil climate)
17 ENV variables  (soil properties, vegetation, etc.)
─────────────────
58 total variables
```

### After filtering (threshold 0.7):
```
19 BIO variables  (ALL KEPT - 100%)
~8 SBIO variables (reduced - many correlated with BIO)
~12 ENV variables (reduced - some correlated with BIO/SBIO)
─────────────────
~39 total variables (33% reduction)
```

**Benefits:**
- **All BIO variables preserved** (climate is primary driver)
- Removed redundant SBIO variables (soil climate proxy for air climate)
- Reduced multicollinearity between variable types
- Faster model training
- SDM pipeline's varSel will further optimize BIO variable selection
- Soil/environmental variables supplement where unique

---

## Troubleshooting

### Error: "No BIO files found"
→ Check paths in lines 15-17 of `variable_correlation_filter.R`

### Error: "Insufficient valid sample points"
→ Rasters have different extents/resolutions. Check raster alignment.

### Warning: "Some selected variables still have correlation > threshold"
→ Normal for correlation chains (A↔B↔C). Review `remaining_high_correlations.xlsx`.

### Variables missing after filtering
→ Expected! Check `removed_variables.xlsx` for reasons.

---

## Advanced: Custom Filtering Logic

### Current Logic (lines 181-268)

```r
# RULE 1: Keep ALL BIO variables (always)
# RULE 2: If BIO vs non-BIO, keep BIO
# RULE 3: Both BIO → keep both
# RULE 4: Neither BIO → SBIO > ENV priority
```

### Example Modification: Keep ALL BIO and ALL SBIO

If you want to keep ALL BIO and ALL SBIO, only remove ENV:

```r
# Replace RULE 2 (lines 210-228) with:
if (cat1 == "BIO" && cat2 == "ENV") {
  # Keep BIO, remove ENV
  selected_vars <- selected_vars[selected_vars != v2]
  # ... rest of removal logic
} else if (cat2 == "BIO" && cat1 == "ENV") {
  # Keep BIO, remove ENV
  selected_vars <- selected_vars[selected_vars != v1]
  # ... rest of removal logic
} else if (cat1 == "SBIO" && cat2 == "ENV") {
  # Keep SBIO, remove ENV
  selected_vars <- selected_vars[selected_vars != v2]
  # ... rest of removal logic
} else if (cat2 == "SBIO" && cat1 == "ENV") {
  # Keep SBIO, remove ENV
  selected_vars <- selected_vars[selected_vars != v1]
  # ... rest of removal logic
} else {
  # Keep both (BIO-BIO, BIO-SBIO, SBIO-SBIO)
  message("  Kept both: ", v1, " ↔ ", v2)
}
```

---

## Citation

If using this workflow in publications, cite:

- **Spearman correlation**: Spearman, C. (1904). *The proof and measurement of association between two things*
- **Multicollinearity in SDMs**: Dormann et al. (2013). *Collinearity: a review of methods to deal with it and a simulation study evaluating their performance*. Ecography, 36(1), 27-46.

---

## Next Steps

1. Run `variable_correlation_filter.R`
2. Review output tables and heatmaps
3. Modify threshold if needed and re-run
4. Integrate filtered variables into SDM pipeline
5. Compare model performance with/without filtering

**Tip**: Keep both filtered and unfiltered versions for comparison!
