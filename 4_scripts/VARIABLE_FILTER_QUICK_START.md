# Variable Filtering - Quick Start

## What It Does

**Keeps ALL BIO variables** (climate - primary drivers)
**Removes SBIO/ENV variables** correlated with BIO (redundant)

## Usage

### Integrated in Main Script (Automatic)

The variable filtering is now **built into** `10_updated_sdmtune_loop.R`:

```r
# Line 31 in 10_updated_sdmtune_loop.R
USE_VARIABLE_FILTER <- TRUE  # Set to FALSE to disable
```

When enabled:
1. On first run, automatically creates filter (saves to `./5_outputs/`)
2. On subsequent runs, loads cached filter
3. SDM runs with filtered variables

**No manual steps needed!**

---

### Standalone (Manual)

To run filtering separately:

```r
source("./4_scripts/quick_variable_filter.R")
```

**Outputs** (in `./5_outputs/`):
- `selected_variable_names.txt` - List of variables to use
- `selected_variables.xlsx` - Details table
- `variable_filtering_summary.png` - Before/after plot

**Time**: ~30 seconds for 58 variables

---

## Configuration

Edit line 13 in `quick_variable_filter.R`:

```r
CORRELATION_THRESHOLD <- 0.7  # Default (recommended)
```

**Higher** (0.8-0.9) = keeps more variables
**Lower** (0.5-0.6) = stricter filtering

---

## Typical Results

### Before: 58 variables
```
BIO1-BIO19   (19 climate)
SBIO1-SBIO11 (22 soil climate)
ENV          (17 soil/vegetation)
```

### After: ~39 variables
```
BIO1-BIO19   (19 ALL KEPT)
SBIO9, etc   (~8 non-correlated with BIO)
sand, pH     (~12 non-correlated with BIO)
```

### Removed: ~19 variables
```
SBIO1 → correlated with BIO1 (r=0.91)
SBIO5 → correlated with BIO5 (r=0.88)
...
```

---

## Why?

**SBIO** (soil temperature) is often a **proxy** for air temperature (BIO)
- Example: Soil mean temp ≈ Air mean temp
- Keeping both creates redundancy
- BIO is more important for species distributions

**Further selection** happens in SDM pipeline via `varSel`:
- Selects best BIO variables
- Removes correlated BIO variables
- Two-stage filtering = cleaner models

---

## Outputs Explained

### `variable_filtering_summary.png`

**Left panel**: Bar chart
- Gray = original count
- Blue = selected count
- Shows filtering effect

**Right panel**: Correlation heatmap
- Selected variables only
- Blue = negative correlation
- Red = positive correlation
- Should have low correlation overall

---

## Toggle On/Off

**To disable filtering:**

```r
# In 10_updated_sdmtune_loop.R, line 31:
USE_VARIABLE_FILTER <- FALSE
```

**To re-enable:**

```r
USE_VARIABLE_FILTER <- TRUE
```

**To regenerate filter:**

Delete old filter and re-run:
```bash
rm ./5_outputs/selected_variable_names.txt
```

---

## Troubleshooting

**Q: Filter file not found?**
A: Script will auto-generate on first run

**Q: Want to use all variables?**
A: Set `USE_VARIABLE_FILTER <- FALSE`

**Q: Filter seems wrong?**
A: Check `./5_outputs/variable_filtering_summary.png`

**Q: Script slow?**
A: Already optimized (5K samples, ~30 sec)

---

## Files

1. **`quick_variable_filter.R`** - Standalone filter script (simplified)
2. **`variable_correlation_filter.R`** - Full version with detailed outputs
3. **`10_updated_sdmtune_loop.R`** - Main SDM script (integrated)

Use **quick** for speed (default), **full** for detailed analysis.
