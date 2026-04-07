# Variable Filtering - Changelog

## Updated Strategy (Current Version)

**Date**: 2025-12-12

### New Filtering Rules

1. ✅ **Keep ALL BIO variables** - Never removed, regardless of correlation
2. ✅ **Remove SBIO** correlated with ANY BIO variable
3. ✅ **Remove ENV** correlated with ANY BIO variable
4. ✅ **Among SBIO/ENV**: SBIO > ENV priority when correlated
5. ✅ **Further selection**: SDM pipeline varSel will optimize BIO variables

### Rationale

**Why keep all BIO variables?**
- BIO variables (WorldClim bioclimatic variables) are the **primary drivers** of species distributions
- Climate is more important than soil properties for most species
- The SDM pipeline already has a **robust variable selection step** (`varSel`) that will:
  - Remove correlated BIO variables
  - Select the most important BIO variables via permutation importance
  - Use correlation threshold of 0.7 (Spearman)

**Why remove correlated SBIO variables?**
- SBIO variables (soil bioclim) are essentially **soil temperature** metrics
- Many SBIO variables are highly correlated with air temperature (BIO) variables
  - Example: SBIO1 (soil annual mean temp) ↔ BIO1 (air annual mean temp) ~ r=0.9
- Soil temperature is often a **proxy** for air temperature, not an independent predictor
- Keeping both creates redundancy and multicollinearity

**Why this is better:**
- Cleaner input to SDM pipeline
- Fewer redundant variables to process
- SDM varSel focuses on selecting among BIO variables (climate)
- SBIO/ENV variables that are NOT correlated with BIO are kept (unique information)

---

## Example Output

### Before Filtering (58 variables):
```
BIO1, BIO2, ..., BIO19                    (19 BIO)
SBIO1_0_5cm, SBIO1_5_15cm, ..., SBIO11   (22 SBIO)
sand, clay, nitrogen, vegetation types     (17 ENV)
```

### After Filtering (~39 variables):

**Kept:**
```
BIO1, BIO2, ..., BIO19                    (ALL 19 BIO kept)
SBIO9_0_5cm, SBIO9_5_15cm                 (Only SBIO NOT correlated with BIO)
sand, nitrogen, vegetation types           (Only ENV NOT correlated with BIO)
```

**Removed:**
```
SBIO1_0_5cm    → correlated with BIO1 (r=0.91)
SBIO10_0_5cm   → correlated with BIO10 (r=0.88)
clay           → correlated with sand (r=-0.87)
```

---

## Expected Correlation Patterns

Based on ecological theory, we expect:

### High Correlation (SBIO removed):
- **SBIO1** (soil mean temp) ↔ **BIO1** (air mean temp)
- **SBIO5** (soil max temp) ↔ **BIO5** (air max temp)
- **SBIO10** (soil warmest quarter) ↔ **BIO10** (air warmest quarter)

### Low Correlation (SBIO kept):
- **SBIO moisture** metrics ↔ BIO temperature (different dimensions)
- **Soil pH** ↔ BIO variables (independent property)
- **Soil nutrients** ↔ BIO variables (independent property)

### Result:
Most temperature-related SBIO variables removed, moisture/chemical SBIO kept.

---

## Validation

After running the filtering script, check:

```r
# Load results
selected <- rio::import("./5_outputs/selected_variables.xlsx")
removed <- rio::import("./5_outputs/removed_variables.xlsx")

# Check: ALL BIO should be selected
table(selected$category)
# BIO: 19 ✓

# Check: No BIO should be removed
table(removed$category)
# BIO: 0 ✓

# View removed SBIO and their reasons
removed %>% filter(category == "SBIO") %>% select(variable, reason)
```

---

## Integration with SDM Pipeline

### Current Workflow (unchanged):

1. **Pre-filtering** (this script) → Removes redundant SBIO/ENV
2. **SDM pipeline loads filtered variables**
3. **Variable importance** → Ranks all variables (permutation)
4. **varSel (correlation-based)** → Removes correlated variables (cor > 0.7)
   - Now focuses on BIO variables (most important)
   - SBIO/ENV already pre-filtered
5. **Hyperparameter optimization** → Final model with selected variables

### Benefits of New Approach:

✅ Cleaner variable set from the start
✅ varSel runs faster (fewer variables to compare)
✅ Focus on climate variables (primary drivers)
✅ Soil/env variables supplement where unique
✅ More ecologically meaningful variable selection

---

## Technical Changes

### File: `variable_correlation_filter.R`

**Lines 181-268**: Filtering logic replaced

**Old logic:**
```r
if (p1 < p2) remove v2
else if (p2 < p1) remove v1
else remove alphabetically
```

**New logic:**
```r
if (cat1 == "BIO" && cat2 != "BIO") remove v2, keep BIO
else if (cat2 == "BIO" && cat1 != "BIO") remove v1, keep BIO
else if (both BIO) KEEP BOTH
else apply priority SBIO > ENV
```

### File: `VARIABLE_FILTERING_README.md`

**Updated sections:**
- Overview → New filtering strategy
- How Filtering Works → 4 rules instead of 3
- Examples → Updated to reflect new logic
- Integration → Emphasize all BIO kept

---

## Testing

Run the script and verify:

```bash
Rscript ./4_scripts/variable_correlation_filter.R
```

**Expected console output:**
```
Protecting 19 BIO variables (never removed)
  Removed: SBIO1_0_5cm (correlated with BIO variable BIO1, r=0.91)
  Removed: SBIO10_5_15cm (correlated with BIO variable BIO10, r=0.88)
  Kept both: BIO5 ↔ BIO10 (r=0.82) - both are BIO variables
...

FILTERING COMPLETE
Original variables: 58
Selected variables: 39
Removed variables: 19
  BIO: 0 (should be 0) ✓
  SBIO: 14
  ENV: 5
```

---

## References

**Why climate > soil for SDMs:**
- Austin, M. (2002). *Spatial prediction of species distribution: an interface between ecological theory and statistical modelling*. Ecological Modelling, 157(2-3), 101-118.
- Elith, J., & Leathwick, J.R. (2009). *Species distribution models: ecological explanation and prediction across space and time*. Annual Review of Ecology, Evolution, and Systematics, 40, 677-697.

**Variable selection in SDMs:**
- Dormann et al. (2013). *Collinearity: a review of methods to deal with it and a simulation study evaluating their performance*. Ecography, 36(1), 27-46.
- Merow et al. (2013). *A practical guide to MaxEnt for modeling species' distributions*. Ecography, 36(10), 1058-1069.

---

## Future Improvements (Optional)

Potential enhancements (not implemented):

1. **Weighted ensemble**: Weight BIO variables higher in correlation checks
2. **Ecological groups**: Group BIO by temperature/precipitation, keep best from each
3. **PCA alternative**: Use PCA on BIO variables instead of correlation filtering
4. **Hierarchical filtering**: Multi-step filtering (BIO→SBIO→ENV sequentially)

Current approach is simpler and works well with existing SDM pipeline.
