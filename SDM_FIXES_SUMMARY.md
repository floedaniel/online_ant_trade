# SDM Script Fixes - Summary

Date: 2025-12-06
Script: `updated_sdmtune_loop.R`

## Issues Identified and Fixed

### Issue 1: Species 2 (Camponotus fallax) - Hyperparameter Optimization Overfitting

**Error Message:**
```
Error during hyperparameter optimization: Optimization algorithm interrupted at generation 0 because it overfits validation dataset.
```

**Root Cause:**
- Species had only 372 records after spatial thinning
- Hyperparameter grid was too aggressive (low regularization, complex features)
- Model overfitted immediately on first generation

**Fix Applied:**
Added adaptive hyperparameter selection based on sample size (lines 635-655):
- **Small samples (<200)**: reg = 1-10, simpler features (lq, lp, lh, qp, qh, lqp, lqh)
- **Medium samples (200-400)**: reg = 0.5-10 (step=1), moderate features
- **Large samples (>400)**: reg = 0.5-10 (step=0.5), full feature set

This prevents overfitting in species with limited data by forcing higher regularization and simpler feature combinations.

---

### Issue 2: Species 3 (Camponotus lateralis) - Variable Selection Failure

**Error Message:**
```
Error during variable selection: Error: glmnet failed to complete regularization path. Model may be infeasible. Try re-running with addsamplestobackground=T.
```

**Root Cause:**
- glmnet algorithm could not find a valid regularization path
- Insufficient background samples relative to predictor variables
- High multicollinearity in predictors

**Fixes Applied:**

1. **Added samples to background** (line 566):
   ```r
   background <- addSamplesToBg(background, all = TRUE)
   ```
   This increases the background sample size by adding all presence points, improving glmnet stability.

2. **Fallback with relaxed correlation threshold** (lines 597-618):
   - If varSel fails with cor_th=0.7, retry with cor_th=0.85
   - Reduces permutations from 10 to 5 in fallback
   - Provides graceful degradation instead of complete failure

---

### Issue 3: Species 4 (Camponotus pilicornis) - Future Prediction Name Mismatch

**Error Message:**
```
Error: [subset] invalid name(s)
```

**Root Cause:**
Layer name mismatch between current and future climate rasters:
- **Current**: `wc2.1_2.5m_bio_1.tif` ... `wc2.1_2.5m_bio_19.tif`
- **Future**: `wc2.1_2.5m_bioc_MPI_ESM1_2_HR_ssp585_2021_2040_processed_layer1.tif` ... `layer19.tif`

The script attempted to standardize names but didn't verify success. When subsetting future predictors at line 752, it failed because selected variable names didn't exist in the future stack.

**Fixes Applied:**

1. **Layer name standardization at data loading** (lines 154-213):
   ```r
   # Fix hyphens
   current_names_fixed <- gsub("-", "_", current_names)
   future_names_fixed <- gsub("-", "_", future_names)

   # Standardize future bioclim naming (layer1 -> wc2.1_2.5m_bio_1)
   for (i in seq_along(future_names_std)) {
     if (grepl("layer\\d+", future_names_std[i])) {
       layer_num <- as.numeric(sub(".*layer(\\d+).*", "\\1", future_names_std[i]))
       future_names_std[i] <- paste0("wc2.1_2.5m_bio_", layer_num)
     }
   }

   # Keep only common layers
   common_layers <- intersect(names(current_predictors), names(future_predictors))
   current_predictors <- current_predictors[[common_layers]]
   future_predictors <- future_predictors[[common_layers]]

   # Final verification
   if (!identical(names(current_predictors), names(future_predictors))) {
     stop("CRITICAL ERROR: Layer names still don't match!")
   }
   ```
   **Key improvement**: All layer name fixes happen ONCE at the beginning, before the modeling loop starts. This ensures perfect alignment for all species and avoids redundant checks inside the loop.

2. **Simplified future subsetting** (lines 807-812):
   ```r
   # Subset future predictors to match final_predictors
   # Note: Layer names are already aligned at the beginning of the script
   future_predictors_subset <- future_predictors[[names(final_predictors)]]

   future_pred <- predict(cv_model, data = future_predictors_subset, ...)
   ```
   Simple direct subsetting with no error checking needed, because layer alignment is guaranteed by the initialization code.

---

## Environmental Data Structure

### Current Predictors:
- 19 WorldClim bioclim variables: `wc2.1_2.5m_bio_1` through `wc2.1_2.5m_bio_19`
- 39 additional environmental layers (soil, land cover, PET, etc.)
- **Total**: 58 layers

### Future Predictors:
- 19 WorldClim bioclim variables (renamed from `layer1`-`layer19` to match current)
- 39 additional environmental layers (shared with current)
- **Total**: 58 layers

### Variable Selection Examples:

**C. cruentatus** (successful):
- Selected: `pet_pm_sr_yr`, `SBIO4_5_15cm_Temperature_Seasonality`, `clay_0_5cm_mean`
- All 3 variables exist in both current and future stacks

**C. pilicornis** (failed before fix):
- Selected: `pet_pm_sr_yr`, `SBIO4_5_15cm_Temperature_Seasonality`, `wc2.1_2.5m_bio_4`, `Shrubs`, `wc2.1_2.5m_bio_15`
- All 5 variables should exist in both stacks after name standardization

---

## Testing Recommendations

1. **Re-run the failed species** to verify fixes:
   - C. fallax (overfitting issue)
   - C. lateralis (glmnet issue)
   - C. pilicornis (name mismatch - should now complete)

2. **Monitor output messages** for:
   - "Keeping X common layers between current and future"
   - Sample size-based hyperparameter selection messages
   - Variable selection fallback attempts

3. **Check for warnings**:
   - Layers only in current or future (should be minimal or zero)
   - Missing variables in future projections (should be zero after layer alignment fix)

4. **Verify outputs**:
   - Future projection TIFs created successfully
   - Climate change difference maps computed
   - No error logs written

---

## Additional Improvements Made

1. **Better error messages**: More informative messages about what went wrong and why
2. **Defensive programming**: Multiple checks before operations that might fail
3. **Graceful degradation**: Fallback options instead of immediate failures
4. **Adaptive strategies**: Different approaches for different data characteristics

---

## Files Modified

- `updated_sdmtune_loop.R`
  - Lines 154-213: Layer name standardization (at data loading)
  - Lines 565-618: Variable selection with glmnet fallback
  - Lines 635-655: Adaptive hyperparameter selection
  - Lines 807-812: Simplified future subsetting (no redundant checks)

## Files Created

- `diagnose_sdm_issues.R` - Diagnostic script to inspect layer names and identify mismatches
- `SDM_FIXES_SUMMARY.md` - This summary document

---

## References

- **Overfitting prevention**: Warren & Seifert (2011) - Regularization in Maxent
- **Background selection**: Barve et al. (2011) - The crucial role of the accessible area
- **Variable selection**: Dormann et al. (2013) - Collinearity: a review of methods
- **Sample size**: van Proosdij et al. (2016) - Minimum required number of specimen records
