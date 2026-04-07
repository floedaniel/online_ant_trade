# ========================================
# Quick Variable Filter - Simplified Version
# KEEP ALL BIO, REMOVE correlated SBIO/ENV
# ========================================

library(terra)
library(tidyverse)
library(rio)

message("Quick Variable Filtering...")

# ---------- Configuration ----------
CORRELATION_THRESHOLD <- 0.7
SAMPLE_SIZE <- 5000  # Reduced for speed

# ---------- Load Environmental Data ----------
bio_stack <- terra::rast(list.files("./1_raw_data/bioclim/current", pattern = "\\.tif$", full.names = TRUE))

sbio_files <- list.files("./1_raw_data/bioclim/sbio", pattern = "\\.tif$", full.names = TRUE)
sbio_stack <- if (length(sbio_files) > 0) terra::rast(sbio_files) else NULL

env_files <- list.files("./1_raw_data/bioclim/env", pattern = "\\.tif$", full.names = TRUE)
env_stack <- if (length(env_files) > 0) terra::rast(env_files) else NULL

# Combine layers
all_predictors <- bio_stack
if (!is.null(sbio_stack)) all_predictors <- c(all_predictors, sbio_stack)
if (!is.null(env_stack)) all_predictors <- c(all_predictors, env_stack)

message("  Loaded ", terra::nlyr(all_predictors), " variables")

# ---------- Categorize Variables ----------
var_names <- names(all_predictors)
var_category <- ifelse(grepl("^BIO[0-9]+", var_names), "BIO",
                       ifelse(grepl("^SBIO[0-9]+", var_names), "SBIO", "ENV"))

bio_vars <- var_names[var_category == "BIO"]
message("  ", length(bio_vars), " BIO variables (all kept)")

# ---------- Sample & Calculate Correlation ----------
set.seed(123)
sample_points <- spatSample(all_predictors, size = SAMPLE_SIZE, method = "random", na.rm = TRUE, as.df = TRUE)
sample_points <- na.omit(sample_points)

cor_matrix <- cor(sample_points, method = "spearman", use = "pairwise.complete.obs")

# ---------- Filter Variables ----------
selected_vars <- var_names
removed_vars <- character(0)

# Get upper triangle
cor_upper <- cor_matrix
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA

# Find high correlations
high_cor <- as.data.frame(as.table(cor_upper)) %>%
  filter(!is.na(Freq), abs(Freq) > CORRELATION_THRESHOLD) %>%
  rename(var1 = Var1, var2 = Var2, cor = Freq)

# Add categories
high_cor$cat1 <- var_category[match(high_cor$var1, var_names)]
high_cor$cat2 <- var_category[match(high_cor$var2, var_names)]

# Apply filtering rules
for (i in seq_len(nrow(high_cor))) {
  v1 <- as.character(high_cor$var1[i])
  v2 <- as.character(high_cor$var2[i])
  cat1 <- high_cor$cat1[i]
  cat2 <- high_cor$cat2[i]

  if (!(v1 %in% selected_vars) || !(v2 %in% selected_vars)) next

  # Rule: Keep BIO, remove non-BIO
  if (cat1 == "BIO" && cat2 != "BIO") {
    selected_vars <- selected_vars[selected_vars != v2]
    removed_vars <- c(removed_vars, v2)
  } else if (cat2 == "BIO" && cat1 != "BIO") {
    selected_vars <- selected_vars[selected_vars != v1]
    removed_vars <- c(removed_vars, v1)
  } else if (cat1 != "BIO" && cat2 != "BIO") {
    # Neither BIO: SBIO > ENV
    if (cat1 == "SBIO" && cat2 == "ENV") {
      selected_vars <- selected_vars[selected_vars != v2]
      removed_vars <- c(removed_vars, v2)
    } else if (cat2 == "SBIO" && cat1 == "ENV") {
      selected_vars <- selected_vars[selected_vars != v1]
      removed_vars <- c(removed_vars, v1)
    } else if (v1 < v2) {
      selected_vars <- selected_vars[selected_vars != v2]
      removed_vars <- c(removed_vars, v2)
    } else {
      selected_vars <- selected_vars[selected_vars != v1]
      removed_vars <- c(removed_vars, v1)
    }
  }
}

removed_vars <- unique(removed_vars)

# ---------- Summary ----------
selected_info <- data.frame(
  variable = selected_vars,
  category = var_category[match(selected_vars, var_names)]
)

message("  Selected: ", length(selected_vars), " (BIO=", sum(selected_info$category == "BIO"),
        ", SBIO=", sum(selected_info$category == "SBIO"),
        ", ENV=", sum(selected_info$category == "ENV"), ")")
message("  Removed: ", length(removed_vars))

# ---------- Save ----------
dir.create("./5_outputs", showWarnings = FALSE, recursive = TRUE)
writeLines(selected_vars, "./5_outputs/selected_variable_names.txt")
rio::export(selected_info, "./5_outputs/selected_variables.xlsx")

# Quick summary plot
png("./5_outputs/variable_filtering_summary.png", width = 800, height = 600)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Before/After bar plot
counts_before <- table(var_category)
counts_after <- table(selected_info$category)
barplot(rbind(counts_before, counts_after),
        beside = TRUE,
        col = c("gray70", "steelblue"),
        legend = c("Original", "Selected"),
        main = "Variable Filtering Results",
        ylab = "Count",
        xlab = "Category",
        ylim = c(0, max(counts_before) * 1.2))

# Correlation heatmap (selected only)
sel_cor <- cor_matrix[selected_vars, selected_vars]
max_vars_plot <- min(30, length(selected_vars))  # Limit to 30 for readability
if (length(selected_vars) > max_vars_plot) {
  # Sample evenly across categories
  sample_idx <- c(1:length(bio_vars),
                  sample(which(selected_info$category == "SBIO"), min(5, sum(selected_info$category == "SBIO"))),
                  sample(which(selected_info$category == "ENV"), min(5, sum(selected_info$category == "ENV"))))
  sample_idx <- unique(sample_idx[1:min(max_vars_plot, length(sample_idx))])
  sel_cor <- sel_cor[sample_idx, sample_idx]
}

image(1:nrow(sel_cor), 1:ncol(sel_cor), sel_cor,
      col = colorRampPalette(c("blue", "white", "red"))(50),
      xlab = "", ylab = "", main = "Selected Variables Correlation",
      axes = FALSE)
axis(1, at = 1:nrow(sel_cor), labels = rownames(sel_cor), las = 2, cex.axis = 0.6)
axis(2, at = 1:ncol(sel_cor), labels = colnames(sel_cor), las = 2, cex.axis = 0.6)
abline(h = (0:nrow(sel_cor)) + 0.5, col = "gray90", lty = 1)
abline(v = (0:ncol(sel_cor)) + 0.5, col = "gray90", lty = 1)

dev.off()

message("  Outputs saved to ./5_outputs/")
message("  - selected_variable_names.txt")
message("  - selected_variables.xlsx")
message("  - variable_filtering_summary.png")

gc()
