# ---------- Load Required Libraries ----------
library(spocc)
library(tidyverse)
library(CoordinateCleaner)
library(SDMtune)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(viridis)
library(zeallot)
library(rio)
library(rnaturalearth)
library(jsonlite)

# ---------- Configuration ----------
species_name <- " Monomorium floricola"  

base_dir <- "./Species"
species_folder <- file.path(base_dir, eppocode)
sdm_folder <- file.path(species_folder, "SDM_maxnet")

if (!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)
if (!dir.exists(sdm_folder)) dir.create(sdm_folder, recursive = TRUE)

MIN_RECORDS <- 100           
CV_FOLDS <- 2               
RECENT_YEARS <- 1950        
BUFFER_DISTANCE <- 2000000  # 2000 km

# ---------- Load Environmental Data ----------
message("Loading environmental data...")

current_dir <- "./1_raw_data/bioclim/current"
current_files <- list.files(current_dir, pattern = "\\.tif$", full.names = TRUE)
current <- terra::rast(current_files)

future_dir <- "./1_raw_data/bioclim/future"
future_files <- list.files(future_dir, pattern = "\\.tif$", full.names = TRUE)
future <- terra::rast(future_files)

env_dir <- "./1_raw_data/bioclim/env"
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)

if (length(env_files) > 0) {
  env_stack <- terra::rast(env_files)
  current_predictors <- c(current, env_stack)
  future_predictors <- c(future, env_stack)
}

raster_extent <- terra::ext(-180, 180, -60, 90)
current_predictors <- terra::crop(current_predictors, raster_extent)
future_predictors <- terra::crop(future_predictors, raster_extent)

message("Loaded ", terra::nlyr(current_predictors), " environmental layers")

# ---------- Download Occurrence Data ----------
message("Downloading occurrence data for: ", species_name)

occ_data <- spocc::occ(
  query = species_name, 
  from = c("gbif", "bison", "ala", "inat", "idigbio"), 
  limit = 1000
)

occ_df <- occ2df(occ_data)
occ_df <- dplyr::select(occ_df, name, longitude, latitude, date, prov)
names(occ_df) <- c("species", "decimalLongitude", "decimalLatitude", "date", "source")
occ_df$species <- species_name

occ_df$year <- NA

if (!all(is.na(occ_df$date))) {
  occ_df$year <- as.numeric(substr(as.character(occ_df$date), 1, 4))
}

message("Raw occurrence records: ", nrow(occ_df))

if (nrow(occ_df) < 100) {
  stop("Insufficient occurrence data (< 100 records).")
}

# ---------- Temporal Filtering (Fithian et al. 2015) ----------
message("Applying temporal filtering...")

if (sum(!is.na(occ_df$year)) > 0.5 * nrow(occ_df)) {
  occ_df_temporal <- occ_df %>% filter(is.na(year) | year >= RECENT_YEARS)
  
  recent_with_year <- occ_df_temporal %>% filter(!is.na(year))
  if (nrow(recent_with_year) > 20) {
    temporal_bias_lat <- lm(decimalLatitude ~ year, data = recent_with_year)
    temporal_bias_lon <- lm(decimalLongitude ~ year, data = recent_with_year)
    
    if (summary(temporal_bias_lat)$coefficients[2,4] < 0.05 | 
        summary(temporal_bias_lon)$coefficients[2,4] < 0.05) {
      warning("Temporal bias detected - interpret with caution")
    }
  }
} else {
  occ_df_temporal <- occ_df
}

message("Records after temporal filtering: ", nrow(occ_df_temporal))

# ---------- Clean Occurrence Data ----------
message("Cleaning occurrence data...")

occ_df_clean <- occ_df_temporal %>%
  cc_val() %>%
  cc_equ() %>%
  cc_cap() %>%
  cc_cen() %>%
  cc_sea() %>%
  cc_zero() %>%
  cc_dupl() %>%
  cc_outl(method = "quantile", mltpl = 3, verbose = TRUE, value = "clean")

message("Records after cleaning: ", nrow(occ_df_clean))

if (nrow(occ_df_clean) < MIN_RECORDS) {
  stop("Insufficient data (< ", MIN_RECORDS, " records).")
}

#rio::export(occ_df, file.path(sdm_folder, paste0("raw_occurrence_", eppocode, ".csv")))
#rio::export(occ_df_clean, file.path(sdm_folder, paste0("cleaned_occurrence_", eppocode, ".csv")))

# ---------- Improved Background Selection (Barve et al. 2011) ----------
message("Generating background points...")

p_coords <- occ_df_clean[, c("decimalLongitude", "decimalLatitude")]

presence_sf <- sf::st_as_sf(p_coords, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Convex hull + buffer approach
presence_projected <- sf::st_transform(presence_sf, crs = 3857)
convex_hull <- sf::st_convex_hull(sf::st_union(presence_projected))
study_area <- sf::st_buffer(convex_hull, dist = BUFFER_DISTANCE)
study_area_geo <- sf::st_transform(study_area, crs = 4326)

bg_points <- sf::st_sample(study_area_geo, size = 10000, type = "random")
bg_coords_sf <- sf::st_sf(geometry = bg_points)

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
world <- sf::st_make_valid(world)
bg_coords_land <- sf::st_filter(bg_coords_sf, world)

bg_coords <- sf::st_coordinates(bg_coords_land)
bg_coords <- as.data.frame(bg_coords[, c("X", "Y")])
names(bg_coords) <- c("x", "y")

valid_presence <- terra::extract(current_predictors, as.matrix(p_coords), cells = TRUE)
p_coords <- p_coords[!is.na(valid_presence[, 1]), ]

valid_bg <- terra::extract(current_predictors, as.matrix(bg_coords), cells = TRUE)
bg_coords <- bg_coords[!is.na(valid_bg[, 1]), ]

message("Presence: ", nrow(p_coords), " | Background: ", nrow(bg_coords))

#rio::export(p_coords, file.path(sdm_folder, paste0("final_occurrence_", eppocode, ".csv")))

# ---------- Prepare SWD and Split Data ----------
message("Preparing data...")

data_swd <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors)
c(train_data, test_data) %<-% trainValTest(data_swd, test = 0.2, only_presence = TRUE, seed = 25)

# ---------- Variable Selection ----------
message("Variable selection...")

# Step 1: Variable importance
initial_model <- train(method = "Maxnet", data = train_data, fc = "lqpth", reg = 1)
vi <- varImp(initial_model, permut = 10)

top_vars <- vi %>% 
  arrange(desc(Permutation_importance)) %>% 
  slice_head(n = min(15, nrow(vi))) %>%
  pull(Variable)

current_predictors_subset <- current_predictors[[top_vars]]

# Step 2: Remove correlated variables using varSel
background <- prepareSWD(species = species_name, a = bg_coords, env = current_predictors_subset)
swd_subset <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = current_predictors_subset)
c(train_subset, test_subset) %<-% trainValTest(swd_subset, test = 0.2, only_presence = TRUE, seed = 25)

selected_variables_model <- varSel(
  model = train(method = "Maxnet", data = train_subset),
  metric = "auc",
  test = test_subset,
  bg4cor = background,
  method = "spearman",
  cor_th = 0.7,
  env = current_predictors_subset,
  use_pc = FALSE,
  progress = TRUE,
  permut = 10
)

selected_vars <- names(selected_variables_model@data@data)
final_predictors <- current_predictors_subset[[selected_vars]]

message("Final variables (", length(selected_vars), "): ", paste(selected_vars, collapse = ", "))

# Prepare final SWD
final_swd <- prepareSWD(species = species_name, p = p_coords, a = bg_coords, env = final_predictors)
c(train_final, test_final) %<-% trainValTest(final_swd, test = 0.2, only_presence = TRUE, seed = 25)

# ---------- Hyperparameter Optimization ----------
message("Optimizing hyperparameters...")

h <- list(
  reg = seq(0.5, 10, 0.5),
  fc = c("lq", "lp", "lh", "qp", "qh", "ph", "lqp", "lqh", "lph", "qph", "lqph")
)

optimized_model <- optimizeModel(
  model = train(method = "Maxnet", data = train_final),
  hypers = h,
  metric = "auc",
  test = test_final,
  keep_best = 0.5,
  seed = 124,
  interactive = TRUE,
  progress = TRUE
)

final_model <- optimized_model@models[[which.max(optimized_model@results$test_AUC)]]

message("Best FC: ", final_model@model@fc, " | Reg: ", final_model@model@reg)

# ---------- Evaluation Metrics (Lobo et al. 2008) ----------
message("Calculating evaluation metrics...")

auc_val <- auc(final_model, test = test_final)
tss_val <- tss(final_model, test = test_final)

ths <- thresholds(final_model, type = "cloglog")
optimal_threshold <- ths[3, 2]

message("AUC: ", round(auc_val, 3), " | TSS: ", round(tss_val, 3), 
        " | Threshold: ", round(optimal_threshold, 3))

# ---------- Cross-Validation ----------
message("Performing ", CV_FOLDS, "-fold cross-validation...")

folds <- randomFolds(final_swd, k = CV_FOLDS, only_presence = TRUE, seed = 321)

cv_model <- train("Maxnet",
                  data = final_swd,
                  fc = final_model@model@fc,
                  reg = final_model@model@reg,
                  folds = folds)

cv_auc_test <- auc(cv_model, test = TRUE)
cv_tss_test <- tss(cv_model, test = TRUE)

message("CV Testing - AUC: ", round(cv_auc_test, 3), " | TSS: ", round(cv_tss_test, 3))

# Predictions aggregated across folds
pred <- predict(cv_model, data = final_predictors, type = "cloglog", 
                fun = c("mean", "max", "sd"), progress = TRUE)

pred_mean <- pred$mean
pred_max <- pred$max
pred_sd <- pred$sd  # Standard deviation across CV folds

# ---------- Future Projections ----------
message("Future projections...")

future_predictors_subset <- future_predictors[[names(final_predictors)]]

future_pred <- predict(cv_model, data = future_predictors_subset, type = "cloglog", 
                       fun = c("mean", "max"), progress = TRUE)

future_pred_mean <- future_pred$mean

# Ensure same geometry
if (!terra::compareGeom(pred_mean, future_pred_mean, stopOnError = FALSE)) {
  future_pred_mean <- terra::resample(future_pred_mean, pred_mean, method = "bilinear")
}

climate_change_diff <- future_pred_mean - pred_mean

diff_values <- terra::values(climate_change_diff, na.rm = TRUE)
mean_change <- mean(diff_values, na.rm = TRUE)
pos_change_area <- sum(diff_values > 0, na.rm = TRUE) / length(diff_values) * 100
neg_change_area <- sum(diff_values < 0, na.rm = TRUE) / length(diff_values) * 100

message("Climate change - Mean: ", round(mean_change, 4), " | Gain: ", round(pos_change_area, 1), 
        "% | Loss: ", round(neg_change_area, 1), "%")

# ---------- Save Rasters ----------
#terra::writeRaster(pred_mean, file.path(sdm_folder, paste0("current_", eppocode, ".tif")), overwrite = TRUE)
#terra::writeRaster(pred_sd, file.path(sdm_folder, paste0("uncertainty_", eppocode, ".tif")), overwrite = TRUE)
#terra::writeRaster(future_pred_mean, file.path(sdm_folder, paste0("future_", eppocode, ".tif")), overwrite = TRUE)
#terra::writeRaster(climate_change_diff, file.path(sdm_folder, paste0("change_", eppocode, ".tif")), overwrite = TRUE)

# ---------- Visualizations ----------
message("Creating visualizations...")

map_theme <- function() {
  theme_minimal() + theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_line(color = "white", size = 0.3),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  )
}

# Current prediction
current_plot <- ggplot() +
  tidyterra::geom_spatraster(data = pred_mean) +
  geom_point(data = p_coords, aes(x = decimalLongitude, y = decimalLatitude), 
             color = "red", size = 0.5, alpha = 0.5) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1)) +
  labs(title = paste("Current Habitat Suitability -", species_name),
       subtitle = paste("Test AUC:", round(auc_val, 3), "| CV AUC:", round(cv_auc_test, 3), 
                        "| TSS:", round(tss_val, 3))) +
  map_theme()

current_plot

ggsave(file.path(sdm_folder, paste0("current_", eppocode, ".png")), 
       current_plot, width = 14, height = 10, dpi = 300, bg = "white")

# Uncertainty (from CV folds)
uncertainty_plot <- ggplot() +
  tidyterra::geom_spatraster(data = pred_sd) +
  scale_fill_viridis_c(name = "SD", option = "magma") +
  labs(title = paste("Prediction Uncertainty -", species_name),
       subtitle = paste("Standard deviation across", CV_FOLDS, "CV folds")) +
  map_theme()

uncertainty_plot

ggsave(file.path(sdm_folder, paste0("uncertainty_", eppocode, ".png")), 
       uncertainty_plot, width = 14, height = 10, dpi = 300, bg = "white")

# Future prediction
future_plot <- ggplot() +
  tidyterra::geom_spatraster(data = future_pred_mean) +
  scale_fill_viridis_c(name = "Suitability", option = "plasma", limits = c(0, 1)) +
  labs(title = paste("Future Habitat Suitability -", species_name),
       subtitle = "Future climate scenario") +
  map_theme()

future_plot

ggsave(file.path(sdm_folder, paste0("future_", eppocode, ".png")), 
       future_plot, width = 14, height = 10, dpi = 300, bg = "white")

# Climate change impact
change_plot <- ggplot() +
  tidyterra::geom_spatraster(data = climate_change_diff) +
  scale_fill_gradient2(name = "Change", low = "red", mid = "white", high = "blue", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = paste("Climate Change Impact -", species_name),
       subtitle = paste("Gain:", round(pos_change_area, 1), "% | Loss:", 
                        round(neg_change_area, 1), "%")) +
  map_theme()

change_plot

ggsave(file.path(sdm_folder, paste0("change_", eppocode, ".png")), 
       change_plot, width = 14, height = 10, dpi = 300, bg = "white")


#  Treshold ---------------------------------------------------------------
# ---------- Threshold Maps ----------
message("Creating threshold maps...")

# Get thresholds
ths <- thresholds(final_model, type = "cloglog")

# Display threshold information
message("Thresholds:")
message("  Minimum training presence: ", round(ths[1, 2], 4))
message("  Equal sens/spec: ", round(ths[2, 2], 4))
message("  Max TSS: ", round(ths[3, 2], 4))

# Create binary maps using plotPA
# Max TSS threshold (recommended)
plotPA(pred_mean, 
       th = ths[3, 2], 
       filename = "binary_maxTSS.tif", 
       format = "GTiff")





# ---------- Model Report ----------
tryCatch({
  modelReport(final_model, type = "cloglog", 
              folder = file.path(sdm_folder, paste0("report_", eppocode)),
              test = test_final, response_curves = TRUE, only_presence = TRUE, 
              jk = TRUE, env = final_predictors, permut = 5)
}, error = function(e) warning("Model report failed: ", e$message))

# ---------- Summary ----------
model_summary <- list(
  species = species_name,
  code = eppocode,
  date = as.character(Sys.Date()),
  data = list(
    raw = nrow(occ_df),
    cleaned = nrow(occ_df_clean),
    final = nrow(p_coords),
    background = nrow(bg_coords)
  ),
  variables = names(final_predictors),
  performance = list(
    test_auc = round(auc_val, 3),
    test_tss = round(tss_val, 3),
    cv_auc = round(cv_auc_test, 3),
    cv_tss = round(cv_tss_test, 3),
    threshold = round(optimal_threshold, 3)
  ),
  parameters = list(
    fc = final_model@model@fc,
    reg = final_model@model@reg,
    cv_folds = CV_FOLDS
  ),
  climate_change = list(
    mean = round(mean_change, 4),
    gain_pct = round(pos_change_area, 1),
    loss_pct = round(neg_change_area, 1)
  )
)

jsonlite::write_json(model_summary, file.path(sdm_folder, paste0("summary_", eppocode, ".json")), 
                     pretty = TRUE, auto_unbox = TRUE)

message("\n", paste(rep("=", 60), collapse = ""))
message("COMPLETED!")
message(paste(rep("=", 60), collapse = ""))
message("Results: ", sdm_folder)
message("Test - AUC: ", round(auc_val, 3), " | TSS: ", round(tss_val, 3))
message("CV   - AUC: ", round(cv_auc_test, 3), " | TSS: ", round(cv_tss_test, 3))

gc()