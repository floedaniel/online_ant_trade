# SDMtune Script Updates Summary

## Overview
This document summarizes all improvements made to `updated_sdmtune_loop.R` based on SDMtune best practices.

## Updates Implemented

### 1. ✅ Increased Cross-Validation Folds
**Line 22**
- Changed `CV_FOLDS` from 2 to 4
- **Benefit**: More robust model evaluation with better statistical reliability

### 2. ✅ Spatial Thinning of Occurrence Data
**Lines 279-310**
- Added `thinData()` function to remove spatial duplicates
- Removes multiple occurrences in the same raster cell
- **Benefit**: Reduces spatial autocorrelation and sampling bias

### 3. ✅ Add Presence to Background
**Lines 377-379**
- Implemented `addSamplesToBg()` function
- Mirrors Maxent's `addsamplestobackground=true` behavior
- **Benefit**: Incorporates presence data diversity into background sampling

### 4. ✅ Variable Importance Visualization & Export
**Lines 403-414**
- Save initial variable importance as Excel file
- Generate bar chart visualization using `plotVarImp()`
- **Benefit**: Clear documentation of which variables drive model predictions

### 5. ✅ Correlation Matrix Export
**Lines 426-437**
- Use `corVar()` to identify correlated variables (>0.7)
- Save correlation pairs to Excel file
- **Benefit**: Document which variables were removed due to multicollinearity

### 6. ✅ Selected Variables Documentation
**Lines 472-491**
- Save list of selected variables after correlation removal
- Export final variable importance after selection
- Generate final variable importance plot
- **Benefit**: Complete audit trail of variable selection process

### 7. ✅ Optimization Results Export
**Lines 532-540**
- Save hyperparameter optimization results
- Documents all tested configurations with AUC scores
- **Benefit**: Track model tuning process and alternative configurations

### 8. ✅ Confusion Matrix Statistics
**Lines 554-564**
- Calculate confusion matrix at optimal threshold
- Compute sensitivity and specificity
- Display true/false positives/negatives
- **Benefit**: Detailed performance metrics beyond AUC/TSS

### 9. ✅ ROC Curve Visualization
**Lines 594-601**
- Generate ROC curve plot using `plotROC()`
- Display with cross-validated AUC value
- **Benefit**: Visual representation of model discrimination ability

### 10. ✅ Jackknife Test for Variable Importance
**Lines 603-624**
- Perform leave-one-out jackknife test with `doJk()`
- Export results to Excel
- Generate jackknife visualization plot
- **Benefit**: Alternative measure of variable importance via model performance changes

### 11. ✅ Response Curves (Previously Added)
**Lines 747-813**
- Generate standard and marginal response curves for each variable
- Uses `plotResponse()` function from SDMtune
- **Benefit**: Understand how each environmental variable affects species suitability

### 12. ✅ Enhanced Summary Statistics
**Lines 825-853**
- Added confusion matrix statistics to summary
- Added number of variables used
- Added sensitivity and specificity metrics
- **Benefit**: Comprehensive model performance documentation in single file

## Output Files Generated

### Per Species, the script now generates:

**Data Documentation:**
- `initial_variable_importance_[key].xlsx` - Variable importance before selection
- `initial_variable_importance_[key].png` - Bar chart visualization
- `correlated_variables_[key].xlsx` - Pairs of correlated variables
- `selected_variables_[key].xlsx` - Final variables after selection
- `final_variable_importance_[key].xlsx` - Importance of selected variables
- `final_variable_importance_[key].png` - Bar chart of final variables

**Model Tuning:**
- `optimization_results_[key].xlsx` - All hyperparameter configurations tested

**Model Evaluation:**
- `roc_curve_[key].png` - ROC curve visualization
- `jackknife_results_[key].xlsx` - Jackknife test results
- `jackknife_plot_[key].png` - Jackknife visualization

**Response Analysis:**
- `response_curves/response_[variable].png` - Standard response curves (blue)
- `response_curves/marginal_response_[variable].png` - Marginal response curves (red)

**Enhanced Summary:**
- `summary_[key].xlsx` - Now includes confusion matrix stats, sensitivity, specificity

**Existing Files (retained):**
- All spatial prediction maps (current, future, change, uncertainty)
- Binary threshold maps
- Data inspection plots
- Model report HTML

## Performance Improvements

1. **Better Statistical Rigor**: 4-fold CV instead of 2-fold
2. **Reduced Bias**: Spatial thinning removes duplicate locations
3. **Better Documentation**: Complete audit trail of all modeling decisions
4. **Enhanced Evaluation**: Multiple complementary metrics (AUC, TSS, sensitivity, specificity)
5. **Transparency**: All intermediate results saved for review

## Best Practices Implemented

✅ Data preparation with spatial thinning
✅ Correlation-based variable selection
✅ Importance-based variable filtering
✅ Genetic algorithm hyperparameter optimization
✅ Cross-validation for robust evaluation
✅ Multiple evaluation metrics
✅ Response curve analysis
✅ Jackknife variable importance testing
✅ ROC curve visualization
✅ Comprehensive result documentation

## Configuration Recommendations

Current settings are optimized for:
- Minimum 100 occurrence records
- 4-fold cross-validation
- Correlation threshold: 0.7 (Spearman)
- Top 15 variables by importance (before correlation removal)
- 2000 km buffer for background selection
- Post-1950 temporal filtering

## Next Steps

1. Run the updated script on your species list
2. Review the enhanced outputs for each species
3. Use variable importance plots to interpret ecological relationships
4. Compare jackknife vs permutation importance results
5. Examine response curves to understand species-environment relationships

---
**Script Version**: Enhanced with all SDMtune best practices
**Date**: 2025-12-04
**Documentation**: Complete
