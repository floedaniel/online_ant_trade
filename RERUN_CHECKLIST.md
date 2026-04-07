# SDM Re-run Checklist

## What Was Fixed

✅ **Layer name alignment** - Fixed at data loading (lines 154-213)
- All current/future layer names now matched BEFORE modeling loop starts
- No redundant checks inside the loop

✅ **Overfitting prevention** - Adaptive hyperparameters (lines 635-655)
- Small samples (<200): Higher regularization, simpler features
- Medium samples (200-400): Moderate settings
- Large samples (>400): Full feature set

✅ **Variable selection stability** - glmnet fallback (lines 565-618)
- Added `addSamplesToBg(background, all = TRUE)` for stability
- Automatic fallback to cor_th=0.85 if 0.7 fails

✅ **Function parameter errors** - Removed invalid `permut` arguments
- Fixed `doJk()` call (line 780) - removed `permut = 5`
- Fixed `modelReport()` call (line 1007) - removed `permut = 5`
- These functions don't accept `permut` parameter

## Before Running

1. **Optional: Test layer alignment**
   ```r
   source("4_scripts/validate_layer_names.R")
   ```
   Look for: `✓✓✓ SUCCESS: Layer names match perfectly!`

## Run the Script

```r
source("4_scripts/updated_sdmtune_loop.R")
```

## What to Monitor

### At Start (Data Loading)
Look for these messages:
```
=== Standardizing layer names ===
Fixed hyphens: Current (X), Future (Y)
Standardized N future bioclim layers to match current naming
✓ Layer names successfully aligned
  Current predictors: 58 layers
  Future predictors: 58 layers
=== Layer standardization complete ===
```

### During Species Processing
For each species, check for:

**Sample size messages**:
- `Using conservative hyperparameters for small sample size (n=XXX)`
- `Using moderate hyperparameters for medium sample size (n=XXX)`
- `Using full hyperparameters for large sample size (n=XXX)`

**Variable selection**:
- If it says `Attempting variable selection with relaxed correlation threshold...` - this is OK, it's the fallback working

**NO MORE of these errors**:
- ❌ `Optimization algorithm interrupted at generation 0 because it overfits`
- ❌ `glmnet failed to complete regularization path`
- ❌ `Error: [subset] invalid name(s)`

## Expected Results

All 4 species should now complete successfully:

1. ✅ **Camponotus cruentatus** - Already worked
2. ✅ **Camponotus fallax** - Should now complete (overfitting fixed)
3. ✅ **Camponotus lateralis** - Should now complete (glmnet fixed)
4. ✅ **Camponotus pilicornis** - Should now complete (layer names fixed)

## Output Files per Species

Each species folder (`species/KEY_NAME/SDM_maxnet/`) should contain:
- `current_KEY.tif` - Current habitat suitability
- `future_KEY.tif` - Future habitat suitability
- `change_KEY.tif` - Climate change impact
- `uncertainty_KEY.tif` - Model uncertainty
- `selected_variables_KEY.xlsx` - Variables used
- `summary_KEY.xlsx` - Model statistics
- Various PNG plots

## If Issues Persist

### Check the error log:
```r
# List species that failed
list.files("./species", pattern = "error_log.txt", recursive = TRUE, full.names = TRUE)
```

### For debugging a specific species:
1. Check `species/KEY_NAME/SDM_maxnet/error_log.txt`
2. Look at the console output for that species
3. Check if the issue is different from the original three errors

## Notes

- The script now runs with test_species[1:4] - expand to full dataset when all 4 work
- Layer alignment happens ONCE at the start (more efficient)
- Adaptive methods prevent overfitting for species with few records
- Fallback mechanisms provide graceful degradation instead of failures
