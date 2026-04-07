library(SDMtune)
library(terra)
library(mgcv)  # For GAM

# ===== STEP 1: CREATE FOLDS =====
message("Creating folds...")
folds <- randomFolds(train_final, k = CV_FOLDS, only_presence = TRUE, seed = 124)

# ===== STEP 2: TRAIN RF WITH CV =====
message("Training RF...")

n_vars <- ncol(train_final@data) - 1
default_mtry <- floor(sqrt(n_vars))

h_rf <- list(
  ntree = c(300, 500, 1000),
  nodesize = c(1, 3, 5, 10)
)

rf_model <- train(
  method = "RF",
  data = train_final,
  folds = folds,
  mtry = default_mtry,
  ntree = 500,
  nodesize = 1,
  progress = TRUE
)

rf_optimized <- optimizeModel(
  model = rf_model,
  hypers = h_rf,
  metric = "auc",
  test = test_final,
  pop = 10,
  seed = 124,
  progress = TRUE
)

rf_best <- rf_optimized@models[[1]]
pred_rf <- predict(rf_best, data = final_predictors, type = "cloglog")

# ===== STEP 3: PREPARE DATA FOR GLM/GAM =====
message("Preparing data for GLM and GAM...")

# Extract training data as data frame
train_df <- data.frame(
  pa = train_final@pa,  # Presence/absence (1/0)
  train_final@data
)

# Extract prediction data
pred_df <- as.data.frame(final_predictors, xy = TRUE, na.rm = FALSE)

# Get variable names (exclude coordinates)
var_names <- setdiff(names(pred_df), c("x", "y"))

# ===== STEP 4: TRAIN GLM =====
message("Training GLM...")

# Create formula
glm_formula <- as.formula(paste("pa ~", paste(var_names, collapse = " + ")))

# Fit GLM with binomial family
glm_model <- glm(glm_formula, 
                 data = train_df, 
                 family = binomial(link = "logit"))

# Predict
pred_glm_df <- predict(glm_model, 
                       newdata = pred_df, 
                       type = "response",
                       na.action = na.pass)

# Convert to raster
pred_df$glm <- pred_glm_df
pred_glm <- rast(pred_df[, c("x", "y", "glm")], type = "xyz", crs = crs(final_predictors))
names(pred_glm) <- "GLM"

# ===== STEP 5: TRAIN GAM =====
message("Training GAM...")

# Create GAM formula with smooth terms
# Use s() for smoothing splines on continuous variables
gam_formula <- as.formula(paste("pa ~", 
                                paste0("s(", var_names, ")", collapse = " + ")))

# Fit GAM
gam_model <- gam(gam_formula,
                 data = train_df,
                 family = binomial(link = "logit"),
                 method = "REML")

# Predict
pred_gam_df <- predict(gam_model,
                       newdata = pred_df,
                       type = "response",
                       na.action = na.pass)

# Convert to raster
pred_df$gam <- pred_gam_df
pred_gam <- rast(pred_df[, c("x", "y", "gam")], type = "xyz", crs = crs(final_predictors))
names(pred_gam) <- "GAM"

# ===== STEP 6: CREATE ENSEMBLE =====
message("Creating ensemble...")

# Stack all predictions
ensemble_stack <- c(pred_rf, pred_glm, pred_gam)
names(ensemble_stack) <- c("RF", "GLM", "GAM")

# Mean ensemble
ensemble_mean <- mean(ensemble_stack, na.rm = TRUE)
names(ensemble_mean) <- "ensemble_mean"

# ===== STEP 7: PLOT =====

# Alternative: ggplot
library(ggplot2)
library(tidyterra)

ggplot() +
  tidyterra::geom_spatraster(data = ensemble_mean, maxcell = 5e+07) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", 
                       limits = c(0, 1), na.value = "white") +
  labs(title = paste("RF+GLM+GAM Ensemble -", species_name),
       x = "Longitude", y = "Latitude") +
  map_theme()

ggsave(file.path(sdm_folder, paste0("Ensemble_", species_key, ".png")),
       current_plot, width = 14, height = 10, dpi = 300, bg = "white")

# ===== STEP 8: SAVE =====
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_RF_", species_key, ".tif")), overwrite = TRUE)
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_GLM_", species_key, ".tif")), overwrite = TRUE)
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_GAM_", species_key, ".tif")), overwrite = TRUE)
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_ensemble_", species_key, ".tif")), overwrite = TRUE)
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_ensemble_sd_", species_key, ".tif")), overwrite = TRUE)
terra::writeRaster(climate_change_diff, file.path(raster_folder, paste0("prediction_ensemble_cv_", species_key, ".tif")), overwrite = TRUE)
