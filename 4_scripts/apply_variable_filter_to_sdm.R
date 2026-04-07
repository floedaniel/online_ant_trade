# ========================================
# Apply Variable Filter to SDM Pipeline
# Helper script to use filtered variables
# ========================================

# This script provides a function to load and apply the filtered variable set
# from the variable_correlation_filter.R output

#' Load filtered environmental variables
#'
#' @param filter_file Path to selected_variable_names.txt (default: "./5_outputs/selected_variable_names.txt")
#' @param use_filter Logical, whether to use filtering (TRUE) or all variables (FALSE)
#' @return A list with filtered current and future predictor stacks
#' @export
load_filtered_predictors <- function(filter_file = "./5_outputs/selected_variable_names.txt",
                                     use_filter = TRUE) {

  message("Loading environmental data...")

  # Load current bioclim
  current_predictors_dir <- "./1_raw_data/bioclim/current"
  current_predictors_files <- list.files(current_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
  message("Found ", length(current_predictors_files), " current bioclim files")
  if (length(current_predictors_files) == 0) {
    stop("ERROR: No current bioclim files found in ", current_predictors_dir)
  }
  current_predictors <- terra::rast(current_predictors_files)
  message("Loaded current bioclim: ", terra::nlyr(current_predictors), " layers")

  # Load future bioclim
  future_predictors_dir <- "./1_raw_data/bioclim/future"
  future_predictors_files <- list.files(future_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
  message("Found ", length(future_predictors_files), " future bioclim files")
  if (length(future_predictors_files) == 0) {
    stop("ERROR: No future bioclim files found in ", future_predictors_dir)
  }
  future_predictors <- terra::rast(future_predictors_files)
  message("Loaded future bioclim: ", terra::nlyr(future_predictors), " layers")

  # Load SBIO layers
  sbio_dir <- "./1_raw_data/bioclim/sbio"
  sbio_files <- list.files(sbio_dir, pattern = "\\.tif$", full.names = TRUE)
  message("Found ", length(sbio_files), " SBIO files")
  if (length(sbio_files) == 0) {
    warning("WARNING: No SBIO files found in ", sbio_dir)
    sbio_stack <- NULL
  } else {
    sbio_stack <- terra::rast(sbio_files)
    message("Loaded SBIO layers: ", terra::nlyr(sbio_stack), " layers")
  }

  # Load environmental layers
  env_dir <- "./1_raw_data/bioclim/env"
  env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
  message("Found ", length(env_files), " environmental files")
  if (length(env_files) == 0) {
    warning("WARNING: No environmental files found in ", env_dir)
    env_stack <- NULL
  } else {
    env_stack <- terra::rast(env_files)
    message("Loaded environmental layers: ", terra::nlyr(env_stack), " layers")
  }

  # Combine all layers (for current)
  if (!is.null(sbio_stack) && !is.null(env_stack)) {
    current_predictors <- c(current_predictors, sbio_stack, env_stack)
    # Note: SBIO and ENV are typically static, so they're the same for current/future
    future_predictors <- c(future_predictors, sbio_stack, env_stack)
  } else if (!is.null(sbio_stack)) {
    current_predictors <- c(current_predictors, sbio_stack)
    future_predictors <- c(future_predictors, sbio_stack)
  } else if (!is.null(env_stack)) {
    current_predictors <- c(current_predictors, env_stack)
    future_predictors <- c(future_predictors, env_stack)
  }

  message("Total layers before filtering - Current: ", terra::nlyr(current_predictors),
          ", Future: ", terra::nlyr(future_predictors))

  # Apply filtering if requested
  if (use_filter) {
    if (!file.exists(filter_file)) {
      stop("ERROR: Filter file not found: ", filter_file,
           "\nPlease run variable_correlation_filter.R first!")
    }

    message("\nApplying variable filter from: ", filter_file)
    selected_vars <- readLines(filter_file)
    message("Selected variables: ", length(selected_vars))

    # Check which selected variables exist in the current stack
    available_vars <- names(current_predictors)
    missing_vars <- setdiff(selected_vars, available_vars)

    if (length(missing_vars) > 0) {
      warning("WARNING: ", length(missing_vars), " selected variables not found in current predictors:")
      warning("  ", paste(missing_vars, collapse = ", "))
      selected_vars <- intersect(selected_vars, available_vars)
      message("Using ", length(selected_vars), " available selected variables")
    }

    # Filter to selected variables only
    current_predictors <- current_predictors[[selected_vars]]
    future_predictors <- future_predictors[[selected_vars]]

    message("Filtered to ", length(selected_vars), " variables")
    message("Variables: ", paste(selected_vars, collapse = ", "))
  } else {
    message("\nNo filtering applied - using all variables")
  }

  return(list(
    current = current_predictors,
    future = future_predictors
  ))
}

# Example usage:
# predictors <- load_filtered_predictors(use_filter = TRUE)
# current_predictors <- predictors$current
# future_predictors <- predictors$future
