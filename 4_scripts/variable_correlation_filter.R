# ========================================
# Variable Correlation Analysis & Filtering
# Prioritizes BIO > SBIO > ENV variables
# ========================================

library(terra)
library(tidyverse)
library(corrplot)
library(rio)
library(viridis)

# ---------- Configuration ----------
CORRELATION_THRESHOLD <- 0.7  # Spearman correlation threshold
SAMPLE_SIZE <- 10000          # Number of points to sample for correlation analysis

# ---------- Load Environmental Data ----------
message("Loading environmental layers...")

# Load current bioclim (BIO variables)
current_predictors_dir <- "./1_raw_data/bioclim/current"
current_predictors_files <- list.files(current_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(current_predictors_files), " BIO variables")
if (length(current_predictors_files) == 0) {
  stop("ERROR: No BIO files found in ", current_predictors_dir)
}
bio_stack <- terra::rast(current_predictors_files)
message("Loaded BIO variables: ", terra::nlyr(bio_stack), " layers")

# Load SBIO layers (soil bioclim)
sbio_dir <- "./1_raw_data/bioclim/sbio"
sbio_files <- list.files(sbio_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(sbio_files), " SBIO files")
if (length(sbio_files) == 0) {
  warning("WARNING: No SBIO files found in ", sbio_dir)
  sbio_stack <- NULL
} else {
  sbio_stack <- terra::rast(sbio_files)
  message("Loaded SBIO variables: ", terra::nlyr(sbio_stack), " layers")
}

# Load environmental layers
env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
message("Found ", length(env_files), " ENV files")
if (length(env_files) == 0) {
  warning("WARNING: No ENV files found in ", env_dir)
  env_stack <- NULL
} else {
  env_stack <- terra::rast(env_files)
  message("Loaded ENV variables: ", terra::nlyr(env_stack), " layers")
}

# Combine all layers
if (!is.null(sbio_stack) && !is.null(env_stack)) {
  all_predictors <- c(bio_stack, sbio_stack, env_stack)
} else if (!is.null(sbio_stack)) {
  all_predictors <- c(bio_stack, sbio_stack)
} else if (!is.null(env_stack)) {
  all_predictors <- c(bio_stack, env_stack)
} else {
  all_predictors <- bio_stack
}

message("Total layers loaded: ", terra::nlyr(all_predictors))
message("Variable names: ", paste(names(all_predictors), collapse = ", "))

# ---------- Create Variable Categories ----------
message("\nCategorizing variables by priority...")

var_names <- names(all_predictors)

# Categorize each variable
var_category <- sapply(var_names, function(vname) {
  if (grepl("^BIO[0-9]+", vname)) {
    return("BIO")
  } else if (grepl("^SBIO[0-9]+", vname)) {
    return("SBIO")
  } else {
    return("ENV")
  }
})

# Priority levels: BIO = 1 (highest), SBIO = 2, ENV = 3 (lowest)
var_priority <- ifelse(var_category == "BIO", 1,
                       ifelse(var_category == "SBIO", 2, 3))

var_info <- data.frame(
  variable = var_names,
  category = var_category,
  priority = var_priority,
  stringsAsFactors = FALSE
)

message("BIO variables: ", sum(var_category == "BIO"))
message("SBIO variables: ", sum(var_category == "SBIO"))
message("ENV variables: ", sum(var_category == "ENV"))

# Export variable categorization
rio::export(var_info, "./5_outputs/variable_categorization.xlsx")

# ---------- Sample Points for Correlation Analysis ----------
message("\nSampling ", SAMPLE_SIZE, " random points for correlation analysis...")

# Create random sample points
set.seed(123)
sample_points <- spatSample(all_predictors, size = SAMPLE_SIZE,
                            method = "random", na.rm = TRUE, as.df = TRUE)

# Remove any rows with NA values
sample_points <- na.omit(sample_points)
message("Valid sample points after removing NAs: ", nrow(sample_points))

if (nrow(sample_points) < 100) {
  stop("ERROR: Insufficient valid sample points (", nrow(sample_points), "). Check raster coverage.")
}

# ---------- Calculate Correlation Matrix ----------
message("\nCalculating Spearman correlation matrix...")

# Calculate Spearman correlation
cor_matrix <- cor(sample_points, method = "spearman", use = "pairwise.complete.obs")

# Save full correlation matrix
rio::export(as.data.frame(cor_matrix), "./5_outputs/full_correlation_matrix.xlsx")
message("Full correlation matrix saved to ./5_outputs/full_correlation_matrix.xlsx")

# ---------- Visualize Correlation Matrix ----------
message("\nCreating correlation heatmap...")

png("./5_outputs/correlation_heatmap_full.png", width = 2400, height = 2400, res = 150)
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.cex = 0.7,
         tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = NULL,  # Don't add coefficients (too many variables)
         cl.lim = c(-1, 1),
         title = "Full Correlation Matrix (Spearman)",
         mar = c(0, 0, 2, 0))
dev.off()

message("Full correlation heatmap saved to ./5_outputs/correlation_heatmap_full.png")

# ---------- Identify Highly Correlated Variable Pairs ----------
message("\nIdentifying highly correlated variable pairs (|r| > ", CORRELATION_THRESHOLD, ")...")

# Get upper triangle of correlation matrix
cor_upper <- cor_matrix
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA

# Convert to long format
cor_long <- as.data.frame(as.table(cor_upper))
names(cor_long) <- c("var1", "var2", "correlation")

# Filter for high correlations
high_cor <- cor_long %>%
  filter(!is.na(correlation)) %>%
  filter(abs(correlation) > CORRELATION_THRESHOLD) %>%
  arrange(desc(abs(correlation)))

message("Found ", nrow(high_cor), " variable pairs with |correlation| > ", CORRELATION_THRESHOLD)

if (nrow(high_cor) > 0) {
  # Add category and priority info
  high_cor <- high_cor %>%
    left_join(var_info %>% rename(var1 = variable, cat1 = category, pri1 = priority), by = "var1") %>%
    left_join(var_info %>% rename(var2 = variable, cat2 = category, pri2 = priority), by = "var2")

  # Export high correlation pairs
  rio::export(high_cor, "./5_outputs/high_correlation_pairs.xlsx")
  message("High correlation pairs saved to ./5_outputs/high_correlation_pairs.xlsx")

  # Show top 20 correlations
  message("\nTop 20 correlated pairs:")
  print(head(high_cor %>% select(var1, var2, correlation, cat1, cat2), 20))
}

# ---------- Filter Variables Based on Priority ----------
message("\nFiltering variables based on priority:")
message("  - KEEP: All BIO variables (no filtering)")
message("  - REMOVE: SBIO/ENV variables correlated with ANY BIO variable")
message("  - REMOVE: SBIO/ENV variables correlated with each other")

# Start with all variables as candidates
selected_vars <- var_names
removed_vars <- character(0)
removal_reasons <- character(0)

# RULE 1: Keep ALL BIO variables (always)
bio_vars <- var_names[var_category == "BIO"]
message("\nProtecting ", length(bio_vars), " BIO variables (never removed)")

# Process high correlation pairs
if (nrow(high_cor) > 0) {
  for (i in seq_len(nrow(high_cor))) {
    v1 <- as.character(high_cor$var1[i])
    v2 <- as.character(high_cor$var2[i])
    cat1 <- as.character(high_cor$cat1[i])
    cat2 <- as.character(high_cor$cat2[i])
    cor_val <- round(high_cor$correlation[i], 3)

    # Skip if both variables already removed
    if (!(v1 %in% selected_vars) || !(v2 %in% selected_vars)) {
      next
    }

    # RULE 2: If one variable is BIO, always keep BIO and remove the other
    if (cat1 == "BIO" && cat2 != "BIO") {
      # Keep BIO (v1), remove non-BIO (v2)
      selected_vars <- selected_vars[selected_vars != v2]
      removed_vars <- c(removed_vars, v2)
      removal_reasons <- c(removal_reasons,
                          paste0("Correlated with BIO variable ", v1, " (r=", cor_val, ")"))
      message("  Removed: ", v2, " (correlated with BIO variable ", v1, ", r=", cor_val, ")")
    } else if (cat2 == "BIO" && cat1 != "BIO") {
      # Keep BIO (v2), remove non-BIO (v1)
      selected_vars <- selected_vars[selected_vars != v1]
      removed_vars <- c(removed_vars, v1)
      removal_reasons <- c(removal_reasons,
                          paste0("Correlated with BIO variable ", v2, " (r=", cor_val, ")"))
      message("  Removed: ", v1, " (correlated with BIO variable ", v2, ", r=", cor_val, ")")
    } else if (cat1 == "BIO" && cat2 == "BIO") {
      # RULE 3: Both are BIO - KEEP BOTH (no removal)
      message("  Kept both: ", v1, " ↔ ", v2, " (r=", cor_val, ") - both are BIO variables")
    } else {
      # RULE 4: Neither is BIO - apply priority SBIO > ENV
      p1 <- high_cor$pri1[i]
      p2 <- high_cor$pri2[i]

      if (p1 < p2) {
        # SBIO has priority over ENV
        selected_vars <- selected_vars[selected_vars != v2]
        removed_vars <- c(removed_vars, v2)
        removal_reasons <- c(removal_reasons,
                            paste0("Correlated with ", cat1, " variable ", v1, " (r=", cor_val, "); ",
                                   cat1, " has higher priority than ", cat2))
        message("  Removed: ", v2, " (", cat2, " correlated with ", cat1, " variable ", v1, ", r=", cor_val, ")")
      } else if (p2 < p1) {
        selected_vars <- selected_vars[selected_vars != v1]
        removed_vars <- c(removed_vars, v1)
        removal_reasons <- c(removal_reasons,
                            paste0("Correlated with ", cat2, " variable ", v2, " (r=", cor_val, "); ",
                                   cat2, " has higher priority than ", cat1))
        message("  Removed: ", v1, " (", cat1, " correlated with ", cat2, " variable ", v2, ", r=", cor_val, ")")
      } else {
        # Same priority (both SBIO or both ENV) - keep first alphabetically
        if (v1 < v2) {
          selected_vars <- selected_vars[selected_vars != v2]
          removed_vars <- c(removed_vars, v2)
          removal_reasons <- c(removal_reasons,
                              paste0("Correlated with ", v1, " (r=", cor_val, "); ",
                                     "same category (", cat1, "), ", v1, " kept alphabetically"))
          message("  Removed: ", v2, " (same category as ", v1, ", r=", cor_val, ")")
        } else {
          selected_vars <- selected_vars[selected_vars != v1]
          removed_vars <- c(removed_vars, v1)
          removal_reasons <- c(removal_reasons,
                              paste0("Correlated with ", v2, " (r=", cor_val, "); ",
                                     "same category (", cat1, "), ", v2, " kept alphabetically"))
          message("  Removed: ", v1, " (same category as ", v2, ", r=", cor_val, ")")
        }
      }
    }
  }
}

message("\n", paste(rep("=", 80), collapse = ""))
message("FILTERING COMPLETE")
message(paste(rep("=", 80), collapse = ""))
message("Original variables: ", length(var_names))
message("Selected variables: ", length(selected_vars))
message("Removed variables: ", length(removed_vars))

# ---------- Summary of Selected Variables ----------
selected_info <- var_info %>%
  filter(variable %in% selected_vars) %>%
  arrange(priority, variable)

message("\nSelected variables by category:")
message("  BIO: ", sum(selected_info$category == "BIO"))
message("  SBIO: ", sum(selected_info$category == "SBIO"))
message("  ENV: ", sum(selected_info$category == "ENV"))

# Export selected variables
rio::export(selected_info, "./5_outputs/selected_variables.xlsx")
message("\nSelected variables saved to ./5_outputs/selected_variables.xlsx")

# ---------- Summary of Removed Variables ----------
if (length(removed_vars) > 0) {
  removed_info <- data.frame(
    variable = removed_vars,
    category = var_category[match(removed_vars, var_names)],
    priority = var_priority[match(removed_vars, var_names)],
    reason = removal_reasons,
    stringsAsFactors = FALSE
  ) %>%
    arrange(priority, variable)

  message("\nRemoved variables by category:")
  message("  BIO: ", sum(removed_info$category == "BIO"))
  message("  SBIO: ", sum(removed_info$category == "SBIO"))
  message("  ENV: ", sum(removed_info$category == "ENV"))

  # Export removed variables
  rio::export(removed_info, "./5_outputs/removed_variables.xlsx")
  message("\nRemoved variables saved to ./5_outputs/removed_variables.xlsx")

  # Show removed variables
  message("\nRemoved variables:")
  print(removed_info)
}

# ---------- Correlation Matrix of Selected Variables ----------
message("\nCalculating correlation matrix for selected variables...")

selected_cor_matrix <- cor_matrix[selected_vars, selected_vars]

# Export selected correlation matrix
rio::export(as.data.frame(selected_cor_matrix), "./5_outputs/selected_correlation_matrix.xlsx")

# Visualize selected correlation matrix
png("./5_outputs/correlation_heatmap_selected.png",
    width = max(1200, length(selected_vars) * 40),
    height = max(1200, length(selected_vars) * 40),
    res = 150)
corrplot(selected_cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black",
         number.cex = 0.6,
         cl.lim = c(-1, 1),
         title = "Selected Variables Correlation Matrix (Spearman)",
         mar = c(0, 0, 2, 0))
dev.off()

message("Selected correlation heatmap saved to ./5_outputs/correlation_heatmap_selected.png")

# ---------- Check Maximum Correlation in Selected Set ----------
selected_cor_upper <- selected_cor_matrix
selected_cor_upper[lower.tri(selected_cor_upper, diag = TRUE)] <- NA

max_selected_cor <- max(abs(selected_cor_upper), na.rm = TRUE)
message("\nMaximum absolute correlation in selected variables: ", round(max_selected_cor, 3))

if (max_selected_cor > CORRELATION_THRESHOLD) {
  warning("WARNING: Some selected variables still have correlation > ", CORRELATION_THRESHOLD)
  warning("This can happen when multiple variables form correlation chains")

  # Find which pairs
  selected_cor_long <- as.data.frame(as.table(selected_cor_upper))
  names(selected_cor_long) <- c("var1", "var2", "correlation")

  remaining_high_cor <- selected_cor_long %>%
    filter(!is.na(correlation)) %>%
    filter(abs(correlation) > CORRELATION_THRESHOLD) %>%
    arrange(desc(abs(correlation)))

  if (nrow(remaining_high_cor) > 0) {
    message("\nRemaining high correlations:")
    print(remaining_high_cor)
    rio::export(remaining_high_cor, "./5_outputs/remaining_high_correlations.xlsx")
  }
}

# ---------- Create Filtered Raster Stack ----------
message("\nCreating filtered raster stack with selected variables...")

filtered_stack <- all_predictors[[selected_vars]]

# Export list of selected variable names for use in other scripts
writeLines(selected_vars, "./5_outputs/selected_variable_names.txt")
message("Selected variable names saved to ./5_outputs/selected_variable_names.txt")

# ---------- Summary Report ----------
message("\n", paste(rep("=", 80), collapse = ""))
message("VARIABLE CORRELATION ANALYSIS - SUMMARY REPORT")
message(paste(rep("=", 80), collapse = ""))
message("Filtering Strategy: KEEP ALL BIO, REMOVE correlated SBIO/ENV")
message("Correlation threshold: ", CORRELATION_THRESHOLD)
message("Sample size: ", nrow(sample_points), " points")
message("")
message("ORIGINAL VARIABLES:")
message("  BIO: ", sum(var_category == "BIO"))
message("  SBIO: ", sum(var_category == "SBIO"))
message("  ENV: ", sum(var_category == "ENV"))
message("  Total: ", length(var_names))
message("")
message("SELECTED VARIABLES (after filtering):")
message("  BIO: ", sum(selected_info$category == "BIO"), " (ALL BIO KEPT)")
message("  SBIO: ", sum(selected_info$category == "SBIO"))
message("  ENV: ", sum(selected_info$category == "ENV"))
message("  Total: ", length(selected_vars))
message("")
if (length(removed_vars) > 0) {
  message("REMOVED VARIABLES:")
  message("  BIO: ", sum(removed_info$category == "BIO"), " (should be 0)")
  message("  SBIO: ", sum(removed_info$category == "SBIO"))
  message("  ENV: ", sum(removed_info$category == "ENV"))
  message("  Total: ", length(removed_vars))
} else {
  message("REMOVED VARIABLES: None")
}
message("")
message("FILTERING RULES APPLIED:")
message("  1. ALL BIO variables kept (no BIO removed)")
message("  2. SBIO/ENV removed if correlated with ANY BIO variable")
message("  3. Among SBIO/ENV: SBIO > ENV priority")
message("  4. Further variable selection happens in SDM pipeline")
message("")
message("OUTPUT FILES:")
message("  ./5_outputs/variable_categorization.xlsx")
message("  ./5_outputs/full_correlation_matrix.xlsx")
message("  ./5_outputs/correlation_heatmap_full.png")
message("  ./5_outputs/high_correlation_pairs.xlsx")
message("  ./5_outputs/selected_variables.xlsx")
message("  ./5_outputs/removed_variables.xlsx")
message("  ./5_outputs/selected_correlation_matrix.xlsx")
message("  ./5_outputs/correlation_heatmap_selected.png")
message("  ./5_outputs/selected_variable_names.txt")
message(paste(rep("=", 80), collapse = ""))
message("\nAnalysis complete!")

# Cleanup
gc()
