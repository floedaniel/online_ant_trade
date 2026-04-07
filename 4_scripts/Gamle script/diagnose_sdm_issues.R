# Diagnostic script for SDM issues
library(terra)
library(rio)

# Load environmental layers
current_dir <- "./1_raw_data/bioclim/current"
future_dir <- "./1_raw_data/bioclim/future"
env_dir <- "./1_raw_data/bioclim/env"

current_files <- list.files(current_dir, pattern = "\\.tif$", full.names = TRUE)
future_files <- list.files(future_dir, pattern = "\\.tif$", full.names = TRUE)
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

current_predictors <- terra::rast(current_files)
future_predictors <- terra::rast(future_files)
env_stack <- terra::rast(env_files)

# Combine layers
current_predictors <- c(current_predictors, env_stack)
future_predictors <- c(future_predictors, env_stack)

# Fix layer names (from script)
current_names <- names(current_predictors)
current_names_fixed <- gsub("-", "_", current_names)
names(current_predictors) <- current_names_fixed

future_names <- names(future_predictors)
future_names_fixed <- gsub("-", "_", future_names)
names(future_predictors) <- future_names_fixed

# Standardize future bioclim layer names
future_names_std <- names(future_predictors)
for (i in seq_along(future_names_std)) {
  if (grepl("layer\\d+", future_names_std[i])) {
    layer_num <- as.numeric(sub(".*layer(\\d+).*", "\\1", future_names_std[i]))
    future_names_std[i] <- paste0("wc2.1_2.5m_bio_", layer_num)
  }
}
names(future_predictors) <- future_names_std

cat("\n=== CURRENT PREDICTORS ===\n")
cat("Number of layers:", terra::nlyr(current_predictors), "\n")
cat("\nLayer names:\n")
print(names(current_predictors))

cat("\n=== FUTURE PREDICTORS ===\n")
cat("Number of layers:", terra::nlyr(future_predictors), "\n")
cat("\nLayer names:\n")
print(names(future_predictors))

cat("\n=== NAME COMPARISON ===\n")
cat("Layers in CURRENT but not in FUTURE:\n")
print(setdiff(names(current_predictors), names(future_predictors)))

cat("\nLayers in FUTURE but not in CURRENT:\n")
print(setdiff(names(future_predictors), names(current_predictors)))

cat("\n=== SPECIES 4 (C. pilicornis) SELECTED VARIABLES ===\n")
if (file.exists("./species/1312770_Camponotus_pilicornis_(Roger,_1859)/SDM_maxnet/selected_variables_1312770.xlsx")) {
  selected_vars <- rio::import("./species/1312770_Camponotus_pilicornis_(Roger,_1859)/SDM_maxnet/selected_variables_1312770.xlsx")
  cat("\nSelected variables:\n")
  print(selected_vars)

  cat("\n=== CHECKING IF SELECTED VARIABLES EXIST IN FUTURE ===\n")
  missing_in_future <- setdiff(selected_vars$variable, names(future_predictors))
  if (length(missing_in_future) > 0) {
    cat("ERROR: The following selected variables are MISSING in future predictors:\n")
    print(missing_in_future)
  } else {
    cat("OK: All selected variables exist in future predictors\n")
  }
} else {
  cat("File not found - species may have failed before variable selection\n")
}

cat("\n=== SPECIES 1 (C. cruentatus) SELECTED VARIABLES ===\n")
if (file.exists("./species/1312486_Camponotus_cruentatus_(Latreille,_1802)/SDM_maxnet/selected_variables_1312486.xlsx")) {
  selected_vars <- rio::import("./species/1312486_Camponotus_cruentatus_(Latreille,_1802)/SDM_maxnet/selected_variables_1312486.xlsx")
  cat("\nSelected variables:\n")
  print(selected_vars)
}

cat("\n=== SPECIES 2 (C. fallax) SELECTED VARIABLES ===\n")
if (file.exists("./species/1312649_Camponotus_fallax_(Nylander,_1856)/SDM_maxnet/selected_variables_1312649.xlsx")) {
  selected_vars <- rio::import("./species/1312649_Camponotus_fallax_(Nylander,_1856)/SDM_maxnet/selected_variables_1312649.xlsx")
  cat("\nSelected variables:\n")
  print(selected_vars)
}

cat("\n=== RECOMMENDATIONS ===\n")
cat("1. Ensure all layer names match between current and future rasters\n")
cat("2. Check for duplicate layer names\n")
cat("3. Verify that variable selection doesn't pick layers unique to current\n")
cat("4. For overfitting issues, try:\n")
cat("   - Increase regularization minimum (reg = seq(1, 10, 0.5))\n")
cat("   - Reduce feature complexity (fewer fc combinations)\n")
cat("   - Add more background points\n")
cat("5. For glmnet issues, try:\n")
cat("   - addsamplestobackground=TRUE in prepareSWD\n")
cat("   - Reduce correlation threshold (cor_th = 0.8 or 0.9)\n")
cat("   - Skip varSel if data is too sparse\n")
