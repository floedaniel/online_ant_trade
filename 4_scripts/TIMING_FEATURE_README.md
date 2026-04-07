# Timing Feature - README

## Overview

The SDM pipeline now tracks processing time at two levels:
1. **Per-species timing** - How long each species takes
2. **Overall timing** - Total pipeline runtime

## What Was Added

### Start of Pipeline (Line 236-240)
```r
overall_start_time <- Sys.time()
message("STARTING SDM PIPELINE")
message("Start time: ", format(overall_start_time, "%Y-%m-%d %H:%M:%S"))
```

### Start of Each Species (Line 248)
```r
species_start_time <- Sys.time()
message("Started: ", format(species_start_time, "%H:%M:%S"))
```

### End of Each Species (Line 1382-1383)
```r
species_end_time <- Sys.time()
species_elapsed <- as.numeric(difftime(species_end_time, species_start_time, units = "mins"))
```

### Added to Summary File (Lines 1390-1392)
```r
start_time = format(species_start_time, "%Y-%m-%d %H:%M:%S"),
end_time = format(species_end_time, "%Y-%m-%d %H:%M:%S"),
processing_time_mins = round(species_elapsed, 2),
```

### Console Output Per Species (Lines 1430-1431)
```r
message("Finished: ", format(species_end_time, "%H:%M:%S"))
message("Processing time: ", round(species_elapsed, 2), " minutes")
```

### End of Pipeline (Lines 1437-1453)
```r
overall_end_time <- Sys.time()
overall_elapsed <- as.numeric(difftime(overall_end_time, overall_start_time, units = "mins"))

message("TIMING SUMMARY:")
message("  Start time: ", format(overall_start_time, "%Y-%m-%d %H:%M:%S"))
message("  End time:   ", format(overall_end_time, "%Y-%m-%d %H:%M:%S"))
message("  Total time: ", round(overall_elapsed, 2), " minutes")
message("  Average per species: ", round(overall_elapsed / nrow(test_species), 2), " minutes")
```

---

## Console Output Examples

### Start of Pipeline
```
================================================================================
STARTING SDM PIPELINE
Start time: 2025-12-12 14:30:45
================================================================================
```

### During Processing
```
================================================================================
Processing species 1 of 10
Species: Formica fusca | Key: 1234567
Started: 14:30:50
...
============================================================
COMPLETED: Formica fusca
============================================================
Results: ./species/1234567_Formica_fusca/SDM_maxnet
Test - AUC: 0.876 | TSS: 0.723
CV   - AUC: 0.864 | TSS: 0.708
Finished: 14:45:23
Processing time: 14.55 minutes
```

### End of Pipeline
```
================================================================================
ALL SPECIES PROCESSING COMPLETED!
================================================================================
Total species processed: 10
Results directory: ./species

TIMING SUMMARY:
  Start time: 2025-12-12 14:30:45
  End time:   2025-12-12 17:15:32
  Total time: 164.78 minutes (2.75 hours)
  Average per species: 16.48 minutes
================================================================================
```

---

## Summary File Fields

The `summary_[speciesKey].xlsx` file now includes:

| Field | Description | Example |
|-------|-------------|---------|
| `start_time` | When species processing started | "2025-12-12 14:30:50" |
| `end_time` | When species processing finished | "2025-12-12 14:45:23" |
| `processing_time_mins` | Total minutes for this species | 14.55 |

These fields appear **after** `date` and **before** `raw_records` in the summary table.

---

## Uses

### 1. Performance Monitoring
- Track which species take longest
- Identify bottlenecks in processing
- Estimate time for future runs

### 2. Resource Planning
```r
# Load all summaries
summaries <- list.files("./species", pattern = "summary_.*\\.xlsx",
                        recursive = TRUE, full.names = TRUE)
timing_data <- map_dfr(summaries, rio::import)

# Analyze timing
mean(timing_data$processing_time_mins)    # Average time
max(timing_data$processing_time_mins)     # Longest species
sum(timing_data$processing_time_mins)     # Total time

# Predict time for N species
avg_time <- mean(timing_data$processing_time_mins)
n_species <- 100
estimated_hours <- (avg_time * n_species) / 60
```

### 3. Identify Slow Species
```r
# Which species took longest?
timing_data %>%
  arrange(desc(processing_time_mins)) %>%
  select(species, processing_time_mins, final_records, num_variables) %>%
  head(10)
```

### 4. Correlate with Data
```r
# Does more data = longer processing?
library(ggplot2)

ggplot(timing_data, aes(x = final_records, y = processing_time_mins)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Processing Time vs. Data Size",
       x = "Number of Presence Records",
       y = "Processing Time (minutes)")
```

---

## Expected Timing

Typical processing times per species (with ensemble models):

| Component | Time (minutes) | Notes |
|-----------|----------------|-------|
| Data download | 0.5 - 2 | Depends on network |
| Data cleaning | 0.2 - 0.5 | |
| Variable selection | 1 - 3 | Depends on # variables |
| Hyperparameter optimization | 3 - 8 | Most time-intensive |
| Cross-validation | 2 - 5 | Depends on CV_FOLDS |
| Predictions & plots | 1 - 2 | |
| Ensemble models (RF+GLM+GAM) | 3 - 6 | If enabled |
| **Total** | **10 - 25 minutes** | Average ~15 min |

**Factors affecting time:**
- Number of presence records (more = longer)
- Number of variables (more = longer)
- CV_FOLDS setting (higher = longer)
- Ensemble enabled/disabled (+5-10 min if enabled)
- Computer performance

---

## Optimization Tips

If processing is too slow:

### 1. Reduce CV Folds (Line 24)
```r
CV_FOLDS <- 3  # Instead of 5 or 10
```
**Impact**: Faster but less robust validation

### 2. Disable Ensemble Models (Line 31)
```r
USE_VARIABLE_FILTER <- TRUE
```
Add after line 31:
```r
USE_ENSEMBLE <- FALSE  # Skip ensemble models
```
**Impact**: -5 to -10 minutes per species

### 3. Reduce Hyperparameter Grid (Lines 693-696)
```r
h <- list(reg = c(1, 3, 5),  # Fewer values
          fc = c("lq", "lqh", "lqpt"))  # Fewer feature classes
```
**Impact**: Faster but may miss optimal settings

### 4. Use Variable Filter (Already Enabled)
```r
USE_VARIABLE_FILTER <- TRUE  # Already line 31
```
**Impact**: -2 to -5 minutes per species

### 5. Parallel Processing (Advanced)
Not implemented but possible with `future` package for multiple species simultaneously.

---

## Troubleshooting

### Issue: Processing time seems stuck

**Check:**
1. Variable optimization step - can take 5-10 minutes
2. Ensemble RF training - can take 3-5 minutes
3. Network issues during data download

**Solution:**
- Console shows current step
- Wait for current step to complete
- Check Task Manager for CPU usage

### Issue: One species takes 2+ hours

**Likely causes:**
- Very large dataset (>5000 records)
- Many variables (>50)
- Optimization not converging

**Solution:**
```r
# Check the summary file:
summary <- rio::import("./species/[key]_[name]/SDM_maxnet/tables/summary_[key].xlsx")
summary$final_records  # How many records?
summary$num_variables  # How many variables?

# If excessive, consider:
# - Reducing records via spatial thinning
# - Reducing variables via filtering
```

---

## Future Enhancements

Possible additions:

1. **Per-step timing** - Track each major step separately
2. **Logging to file** - Save timing data to CSV
3. **Progress bar** - Visual progress indicator
4. **Email notification** - Alert when completed
5. **Parallel processing** - Multiple species at once

Current implementation is simple and sufficient for most use cases.
