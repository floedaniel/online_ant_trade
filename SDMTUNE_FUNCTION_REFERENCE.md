# SDMtune Function Parameter Reference

Quick reference for commonly used SDMtune functions to avoid parameter errors.

## Variable Importance Functions

### `varImp()` - Variable Importance
**HAS** `permut` parameter ✓
```r
varImp(model, permut = 10)
```
- `model`: SDMmodel object
- `permut`: Number of permutations (default = 10)

### `doJk()` - Jackknife Test
**NO** `permut` parameter ✗
```r
doJk(model, metric = "auc", test = test_data, with_only = TRUE)
```
- `model`: SDMmodel or SDMmodelCV object
- `metric`: "auc", "tss", or "aicc"
- `test`: SWD object for testing (optional)
- `with_only`: Run test for each variable in isolation (default = TRUE)
- `env`: Raster for environment (used with aicc)
- `return_models`: Return all models (default = FALSE)
- `progress`: Show progress bar (default = TRUE)

**NOTE**: Variable importance within Jackknife is computed internally, no permutation parameter needed.

---

## Variable Selection Functions

### `varSel()` - Variable Selection
**HAS** `permut` parameter ✓ (passed to internal `varImp` calls)
```r
varSel(
  model = model,
  metric = "auc",
  test = test_data,
  bg4cor = background,
  method = "spearman",
  cor_th = 0.7,
  env = predictors,
  use_pc = FALSE,
  progress = TRUE,
  permut = 10
)
```
- `permut`: Passed to `varImp()` for importance calculation

---

## Model Reporting Functions

### `modelReport()` - Generate Model Report
**NO** `permut` parameter ✗
```r
modelReport(
  model,
  type = "cloglog",
  folder = "output_folder",
  test = test_data,
  response_curves = TRUE,
  only_presence = TRUE,
  jk = TRUE,
  env = predictors
)
```
- `model`: SDMmodel object
- `type`: Output type ("cloglog", "logistic", "raw")
- `folder`: Output folder path
- `test`: SWD object for testing (optional)
- `response_curves`: Generate response curves (default = FALSE)
- `only_presence`: Use only presence for response curves (default = FALSE)
- `jk`: Run Jackknife test (default = FALSE)
- `env`: Raster for environment (required if jk = TRUE)

**NOTE**: Jackknife is run internally if `jk = TRUE`, with default parameters.

---

## Summary

| Function | Has `permut`? | Notes |
|----------|---------------|-------|
| `varImp()` | ✓ YES | Number of permutations for importance calculation |
| `varSel()` | ✓ YES | Passed to internal `varImp()` calls |
| `doJk()` | ✗ NO | Jackknife has no permutation parameter |
| `modelReport()` | ✗ NO | Calls Jackknife internally with defaults |

---

## Fixed Issues in Script

### Line 780: `doJk()` call
**Before** (incorrect):
```r
doJk(final_model, metric = "auc", test = test_final, with_only = TRUE, permut = 5)
```

**After** (correct):
```r
doJk(final_model, metric = "auc", test = test_final, with_only = TRUE)
```

### Line 1007: `modelReport()` call
**Before** (incorrect):
```r
modelReport(final_model, type = "cloglog", folder = folder_path,
            test = test_final, response_curves = TRUE, only_presence = TRUE,
            jk = TRUE, env = final_predictors, permut = 5)
```

**After** (correct):
```r
modelReport(final_model, type = "cloglog", folder = folder_path,
            test = test_final, response_curves = TRUE, only_presence = TRUE,
            jk = TRUE, env = final_predictors)
```

---

## Additional Notes

- The `permut` parameter in `varImp()` and `varSel()` controls how many random permutations are used to assess variable importance
- Higher `permut` values give more stable estimates but take longer
- For `doJk()`, the importance is calculated by actually removing variables, not through permutation
- `modelReport()` will run Jackknife internally if `jk = TRUE`, but you cannot control the Jackknife parameters
