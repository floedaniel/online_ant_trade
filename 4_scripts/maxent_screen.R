library(tidyverse)
library(terra)
library(geodata)
library(rgbif)
library(spocc)
library(sf)
library(rnaturalearth)
library(rio)

# =============================================================================
# SIMPLE COLD-TOLERANCE PRESCREENING FOR NORWAY
# Using bio1 (Mean Annual Temp) + bio6 (Min Temp of Coldest Month)
# Data from multiple sources via spocc
# =============================================================================

# Parameters
MIN_OCCURRENCES <- 100
SPOCC_LIMIT <- 100000

# INCLUDE Species has ≥5 occurrences where bio6 ≤ Norway's mildest winter. Strong evidence of cold tolerance - multiple records from cold climates.
# INCLUDE_MARGINALSpecies has at least 1 occurrence cold enough, but <5 total. Could be outlier/error, so less confident.
MIN_COLD_OCCURRENCES <- 5

dir.create("./5_outputs/screen_output", showWarnings = FALSE)
dir.create("./5_outputs/screen_output/species_histograms", showWarnings = FALSE)

# =============================================================================
# 0. LOAD MASTER SPECIES LIST AND FILTER OUT NORWEGIAN NATIVES
# =============================================================================

master_file <- "./2_processed_data/complete_ant_data.xlsx"
master_data <- rio::import(master_file)

# Norwegian native ant species (to exclude)
norwegian_ant_species <- c(
  "Hypoponera punctatissima", "Dolichoderus quadripunctatus", "Camponotus herculeanus",
  "Camponotus ligniperda", "Camponotus vagus", "Polyergus rufescens",
  "Lasius bicornis", "Lasius brunneus", "Lasius carniolicus", "Lasius citrinus",
  "Lasius flavus", "Lasius fuliginosus", "Lasius meridionalis", "Lasius mixtus",
  "Lasius niger", "Lasius platythorax", "Lasius psammophilus", "Lasius sabularum",
  "Lasius umbratus", "Formica cinerea", "Formica fusca", "Formica gagatoides",
  "Formica lemani", "Formica rufibarbis", "Formica transkaucasica", "Formica aquilonia",
  "Formica lugubris", "Formica polyctena", "Formica pratensis", "Formica rufa",
  "Formica truncorum", "Formica uralensis", "Formica sanguinea", "Formica exsecta",
  "Formica foreli", "Formica forsslundi", "Formica pressilabris", "Formica suecica",
  "Formica cunicularia", "Myrmica lobicornis", "Myrmica lonae", "Myrmica microrubra",
  "Myrmica rubra", "Myrmica ruginodis", "Myrmica rugulosa", "Myrmica sabuleti",
  "Myrmica scabrinodis", "Myrmica schencki", "Myrmica specioides", "Myrmica sulcinodis",
  "Symbiomyrma karavajevi", "Leptothorax acervorum", "Leptothorax goesswaldi",
  "Leptothorax kutteri", "Leptothorax muscorum", "Leptothorax tuberum",
  "Formicoxenus nitidulus", "Harpagoxenus sublaevis", "Stenamma debile",
  "Temnothorax interruptus", "Temnothorax nylanderi", "Tetramorium caespitum",
  "Anergates atratulus", "Myrmecina graminicola"
)

message("Fetching GBIF keys for Norwegian native species...")

norwegian_keys <- map_dfr(norwegian_ant_species, function(sp) {
  match <- tryCatch({
    name_backbone(name = sp, rank = "SPECIES")
  }, error = function(e) NULL)
  
  if (!is.null(match) && !is.null(match$usageKey)) {
    tibble(
      input_name = sp,
      gbif_name = match$canonicalName %||% NA_character_,
      usageKey = match$usageKey,
      acceptedKey = match$acceptedUsageKey %||% match$usageKey,
      status = match$status
    )
  } else {
    tibble(
      input_name = sp,
      gbif_name = NA_character_,
      usageKey = NA_integer_,
      acceptedKey = NA_integer_,
      status = "NOT_FOUND"
    )
  }
})

message(sprintf("Found GBIF keys for %d / %d Norwegian species", 
                sum(!is.na(norwegian_keys$usageKey)), length(norwegian_ant_species)))

norwegian_accepted_keys <- norwegian_keys %>%
  filter(!is.na(acceptedKey)) %>%
  pull(acceptedKey) %>%
  unique()

# Filter master list
test_species <- master_data %>%
  filter(status == "ACCEPTED") %>%
  filter(!acceptedKey %in% norwegian_accepted_keys) %>%
  pull(name)

# For testing - remove this line for full run
# test_species <- test_species[1:250]

message(sprintf("Species for analysis after excluding Norwegian natives: %d", length(test_species)))

export(norwegian_keys, "./5_outputs//norwegian_species_gbif_keys.xlsx")

# =============================================================================
# 1. LOAD CLIMATE DATA (bio1 and bio6)
# =============================================================================

message("\nLoading WorldClim bio1 and bio6...")

bioclim <- worldclim_global(var = "bio", res = 2.5, path = "data/worldclim")
bio1 <- bioclim[[1]]
bio6 <- bioclim[[6]]
names(bio1) <- "bio1"
names(bio6) <- "bio6"

# =============================================================================
# 2. EXTRACT NORWAY'S CLIMATE DISTRIBUTION
# =============================================================================

message("Extracting Norway climate reference...")

norway_sf <- ne_countries(scale = "medium", country = "Norway", returnclass = "sf") %>%
  st_crop(xmin = -10, xmax = 35, ymin = 56, ymax = 72)

norway_poly <- vect(norway_sf)

# Extract both bio1 and bio6
norway_bio1 <- terra::extract(bio1, norway_poly, na.rm = TRUE) %>%
  as_tibble() %>%
  pull(bio1) %>%
  na.omit()

norway_bio6 <- terra::extract(bio6, norway_poly, na.rm = TRUE) %>%
  as_tibble() %>%
  pull(bio6) %>%
  na.omit()

# Norway climate summaries
norway_bio1_summary <- tibble(
  variable = "bio1",
  metric = c("min", "q05", "q25", "median", "q75", "q95", "max"),
  value = c(min(norway_bio1), quantile(norway_bio1, 0.05), quantile(norway_bio1, 0.25),
            median(norway_bio1), quantile(norway_bio1, 0.75), quantile(norway_bio1, 0.95),
            max(norway_bio1))
)

norway_bio6_summary <- tibble(
  variable = "bio6",
  metric = c("min", "q05", "q25", "median", "q75", "q95", "max"),
  value = c(min(norway_bio6), quantile(norway_bio6, 0.05), quantile(norway_bio6, 0.25),
            median(norway_bio6), quantile(norway_bio6, 0.75), quantile(norway_bio6, 0.95),
            max(norway_bio6))
)

norway_climate_summary <- bind_rows(norway_bio1_summary, norway_bio6_summary)

message("\nNorway climate summary (°C × 10):")
print(norway_climate_summary)

# Key thresholds
NORWAY_BIO1_MIN <- min(norway_bio1)
NORWAY_BIO1_MAX <- max(norway_bio1)
NORWAY_BIO6_MIN <- min(norway_bio6)
NORWAY_BIO6_MAX <- max(norway_bio6)
NORWAY_BIO6_Q95 <- quantile(norway_bio6, 0.95)
NORWAY_LAT_MIN <- 58

message(sprintf("\nNorway bio1 range: %.1f to %.1f", NORWAY_BIO1_MIN, NORWAY_BIO1_MAX))
message(sprintf("Norway bio6 range: %.1f to %.1f", NORWAY_BIO6_MIN, NORWAY_BIO6_MAX))

# =============================================================================
# 3. FETCH OCCURRENCES FUNCTION (using spocc) - FIXED
# =============================================================================

fetch_occurrences_spocc <- function(species_name, limit) {
  
  message(sprintf("  Fetching from multiple databases: %s", species_name))
  
  occ_data <- tryCatch({
    spocc::occ(
      query = species_name,
      from = c("gbif", "bison", "ala", "inat", "idigbio"),
      limit = limit,
      has_coords = TRUE
    )
  }, error = function(e) {
    message(sprintf("    ERROR: %s", e$message))
    return(NULL)
  })
  
  if (is.null(occ_data)) {
    return(NULL)
  }
  
  # Convert to dataframe
  occ_df <- tryCatch({
    spocc::occ2df(occ_data)
  }, error = function(e) {
    message(sprintf("    ERROR converting to df: %s", e$message))
    return(NULL)
  })
  
  if (is.null(occ_df) || nrow(occ_df) == 0) {
    return(NULL)
  }
  
  # Check which columns exist and select only those
  required_cols <- c("name", "longitude", "latitude", "prov")
  optional_cols <- c("date")
  
  available_cols <- intersect(required_cols, names(occ_df))
  available_optional <- intersect(optional_cols, names(occ_df))
  
  # Make sure we have the minimum required columns
  if (!all(c("name", "longitude", "latitude", "prov") %in% names(occ_df))) {
    message("    WARNING: Missing required columns, skipping")
    return(NULL)
  }
  
  # Select available columns
  occ_df <- occ_df %>%
    dplyr::select(all_of(c(available_cols, available_optional)))
  
  # Add date column if missing
  if (!"date" %in% names(occ_df)) {
    occ_df$date <- NA
  }
  
  # Standardize column names
  occ_df <- occ_df %>%
    dplyr::rename(
      species = name,
      decimalLongitude = longitude,
      decimalLatitude = latitude,
      source = prov
    ) %>%
    mutate(
      decimalLongitude = as.numeric(decimalLongitude),
      decimalLatitude = as.numeric(decimalLatitude)
    )
  
  # Report sources
  source_counts <- occ_df %>% count(source)
  message(sprintf("    Retrieved %d records from: %s", 
                  nrow(occ_df),
                  paste(sprintf("%s(%d)", source_counts$source, source_counts$n), collapse = ", ")))
  
  return(occ_df)
}
# =============================================================================
# 4. PRESCREENING FUNCTION
# =============================================================================

prescreen_species <- function(species_name, bio1_raster, bio6_raster, 
                              norway_bio6_max, norway_lat_min,
                              min_occ, limit, min_cold_occ) {
  
  result <- tibble(
    species = species_name,
    status = "pending",
    n_total = NA_integer_,
    n_with_climate = NA_integer_,
    n_cold = NA_integer_,
    bio1_min = NA_real_,
    bio1_median = NA_real_,
    bio1_max = NA_real_,
    bio6_min = NA_real_,
    bio6_q05 = NA_real_,
    bio6_median = NA_real_,
    max_latitude = NA_real_,
    sources = NA_character_,
    decision = NA_character_,
    reason = NA_character_
  )
  
  # Fetch occurrences using spocc
  occ_df <- fetch_occurrences_spocc(species_name, limit = limit)
  
  if (is.null(occ_df) || nrow(occ_df) == 0) {
    result$status <- "no_data"
    result$n_total <- 0
    result$decision <- "SKIP"
    result$reason <- "No occurrences found in any database"
    return(list(result = result, bio1_vals = NULL, bio6_vals = NULL, coords = NULL))
  }
  
  # Clean coordinates
  coords <- occ_df %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    filter(abs(decimalLatitude) <= 90, abs(decimalLongitude) <= 180) %>%
    filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  result$n_total <- nrow(coords)
  result$sources <- paste(unique(coords$source), collapse = ", ")
  
  if (result$n_total < min_occ) {
    result$status <- "insufficient_data"
    result$decision <- "REVIEW"
    result$reason <- sprintf("Only %d occurrences (need %d)", result$n_total, min_occ)
    return(list(result = result, bio1_vals = NULL, bio6_vals = NULL, coords = NULL))
  }
  
  # Extract bio1 and bio6
  pts <- vect(coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  bio1_vals <- terra::extract(bio1_raster, pts, ID = FALSE)[[1]]
  bio6_vals <- terra::extract(bio6_raster, pts, ID = FALSE)[[1]]
  
  # Keep only points with both values
  valid_idx <- !is.na(bio1_vals) & !is.na(bio6_vals)
  bio1_vals <- bio1_vals[valid_idx]
  bio6_vals <- bio6_vals[valid_idx]
  coords_valid <- coords[valid_idx, ]
  
  result$n_with_climate <- length(bio6_vals)
  
  if (result$n_with_climate < 10) {
    result$status <- "insufficient_climate"
    result$decision <- "REVIEW"
    result$reason <- "Too few points with valid climate data"
    return(list(result = result, bio1_vals = NULL, bio6_vals = NULL, coords = NULL))
  }
  
  # Compute statistics
  result$bio1_min <- min(bio1_vals)
  result$bio1_median <- median(bio1_vals)
  result$bio1_max <- max(bio1_vals)
  result$bio6_min <- min(bio6_vals)
  result$bio6_q05 <- quantile(bio6_vals, 0.05)
  result$bio6_median <- median(bio6_vals)
  result$max_latitude <- max(coords_valid$decimalLatitude)
  result$n_cold <- sum(bio6_vals <= norway_bio6_max)
  result$status <- "processed"
  result$bio6_max <- max(bio6_vals)
  
  # Decision logic
  if (result$bio6_min > 150) {
    result$decision <- "EXCLUDE"
    result$reason <- sprintf("Tropical: coldest record %.1f°C, max lat %.1f°", 
                             result$bio6_min/10, result$max_latitude)
  } else if (result$bio6_min > 100 & result$max_latitude < 40) {
    result$decision <- "EXCLUDE"
    result$reason <- sprintf("Warm climate only: bio6_min=%.1f°C, max_lat=%.1f°",
                             result$bio6_min/10, result$max_latitude)
  } else if (result$bio6_min > 50 & result$max_latitude < 35) {
    result$decision <- "EXCLUDE"
    result$reason <- sprintf("Subtropical: bio6_min=%.1f°C, max_lat=%.1f°",
                             result$bio6_min/10, result$max_latitude)
  } else if (result$n_cold >= min_cold_occ) {
    result$decision <- "INCLUDE"
    result$reason <- sprintf("%d occurrences in Norway-compatible climate", result$n_cold)
  } else if (result$bio6_min <= norway_bio6_max) {
    result$decision <- "INCLUDE_MARGINAL"
    result$reason <- sprintf("Some cold tolerance: bio6_min=%.1f°C", result$bio6_min/10)
  } else if (result$max_latitude >= 50) {
    result$decision <- "REVIEW"
    result$reason <- sprintf("Northern range but warmer climate: lat=%.1f°", result$max_latitude)
  } else {
    result$decision <- "EXCLUDE"
    result$reason <- sprintf("No cold tolerance: bio6_min=%.1f°C, max_lat=%.1f°",
                             result$bio6_min/10, result$max_latitude)
  }
  
  return(list(result = result, bio1_vals = bio1_vals, bio6_vals = bio6_vals, coords = coords_valid))
}

# =============================================================================
# 5. RUN PRESCREENING (PARALLEL - FIXED)
# =============================================================================

library(furrr)

# Set up parallel workers
n_cores <- parallelly::availableCores() - 1
plan(multisession, workers = n_cores)

message("\n", paste(rep("=", 60), collapse = ""))
message(sprintf("Running prescreening on %d species using %d cores...", 
                length(test_species), n_cores))
message("Using spocc to query: GBIF, BISON, ALA, iNaturalist, iDigBio")
message(paste(rep("=", 60), collapse = ""))

# Parallel processing - load rasters inside each worker
all_results <- future_map(test_species, function(sp) {
  
  # Load rasters fresh inside each worker (required for terra + parallel)
  bioclim_worker <- geodata::worldclim_global(var = "bio", res = 2.5, path = "data/worldclim")
  bio1_worker <- bioclim_worker[[1]]
  bio6_worker <- bioclim_worker[[6]]
  
  # Run prescreening
  res <- prescreen_species(sp, bio1_worker, bio6_worker, NORWAY_BIO6_MAX, NORWAY_LAT_MIN,
                           min_occ = MIN_OCCURRENCES, 
                           limit = SPOCC_LIMIT,
                           min_cold_occ = MIN_COLD_OCCURRENCES)
  
  return(res)
  
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Extract results
results <- map_dfr(all_results, ~ .x$result)

species_bio1_data <- list()
species_bio6_data <- list()
species_coords_data <- list()

for (i in seq_along(test_species)) {
  sp <- test_species[i]
  if (!is.null(all_results[[i]]$bio6_vals)) {
    species_bio1_data[[sp]] <- all_results[[i]]$bio1_vals
    species_bio6_data[[sp]] <- all_results[[i]]$bio6_vals
    species_coords_data[[sp]] <- all_results[[i]]$coords
  }
}

plan(sequential)

# =============================================================================
# 6. SUMMARY
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PRESCREENING SUMMARY")
message(paste(rep("=", 60), collapse = ""))

summary_table <- results %>% count(decision)
print(summary_table)

# Source summary
source_summary <- results %>%
  filter(!is.na(sources)) %>%
  separate_rows(sources, sep = ", ") %>%
  count(sources, name = "n_species") %>%
  arrange(desc(n_species))

message("\nData sources used:")
print(source_summary)

species_to_model <- results %>%
  filter(decision %in% c("INCLUDE", "INCLUDE_MARGINAL", "REVIEW")) %>%
  pull(species)

message(sprintf("\nSpecies excluded: %d", sum(results$decision == "EXCLUDE", na.rm = TRUE)))
message(sprintf("Species for MaxEnt modeling: %d", length(species_to_model)))
message(sprintf("Reduction: %.1f%%", 100 * (1 - length(species_to_model)/nrow(results))))

# =============================================================================
# 7. HISTOGRAM PLOTS
# =============================================================================

message("\nGenerating individual species histograms...")

decision_colors <- c(
  "EXCLUDE" = "#d73027",
  "REVIEW" = "#fee090", 
  "INCLUDE_MARGINAL" = "#a6d96a",
  "INCLUDE" = "#1a9850",
  "SKIP" = "gray70"
)

norway_bio6_df <- tibble(bio6 = norway_bio6)

for (sp in names(species_bio6_data)) {
  
  sp_bio6 <- species_bio6_data[[sp]]
  sp_result <- results %>% filter(species == sp)
  sp_decision <- sp_result$decision
  sp_color <- decision_colors[sp_decision]
  
  p <- ggplot() +
    geom_histogram(data = norway_bio6_df,
                   aes(x = bio6, y = after_stat(density)),
                   fill = "steelblue", alpha = 0.4, bins = 50) +
    geom_histogram(data = tibble(bio6 = sp_bio6),
                   aes(x = bio6, y = after_stat(density)),
                   fill = sp_color, alpha = 0.6, bins = 50) +
    geom_vline(xintercept = NORWAY_BIO6_MAX, linetype = "dashed", 
               color = "steelblue", linewidth = 1) +
    annotate("text", x = NORWAY_BIO6_MAX, y = Inf, label = "Norway max", 
             hjust = -0.1, vjust = 2, color = "steelblue", size = 3) +
    geom_vline(xintercept = min(sp_bio6), linetype = "dashed",
               color = sp_color, linewidth = 1) +
    labs(
      x = "Bio6: Min Temperature of Coldest Month (°C × 10)",
      y = "Density",
      title = sp,
      subtitle = sprintf("%s | n=%d | bio6_min=%.0f",
                         sp_decision, length(sp_bio6), min(sp_bio6))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  
  filename_safe <- str_replace_all(sp, " ", "_")
  ggsave(sprintf("./5_outputs/screen_output/species_histograms/%s.png", filename_safe),
         p, width = 10, height = 6, dpi = 150)
}

message(sprintf("Saved %d histograms", length(species_bio6_data)))

# =============================================================================
# 8. OVERVIEW BOXPLOT
# =============================================================================

message("\nGenerating overview boxplot...")

all_bio6_data <- map_dfr(names(species_bio6_data), function(sp) {
  tibble(
    species = sp,
    bio6 = species_bio6_data[[sp]],
    decision = results %>% filter(species == sp) %>% pull(decision)
  )
})

all_bio6_data <- bind_rows(
  all_bio6_data,
  tibble(species = "NORWAY (reference)", bio6 = norway_bio6, decision = "NORWAY")
)

# Order species by bio6_min (coldest first)
species_order <- results %>%
  filter(status == "processed") %>%
  arrange(bio6_min) %>%
  pull(species)
species_order <- c("NORWAY (reference)", species_order)

all_bio6_data <- all_bio6_data %>%
  mutate(species = factor(species, levels = species_order))

plot_colors <- c(
  "NORWAY" = "steelblue",
  "EXCLUDE" = "#d73027",
  "REVIEW" = "#fee090",
  "INCLUDE_MARGINAL" = "#a6d96a", 
  "INCLUDE" = "#1a9850"
)

p_overview <- ggplot(all_bio6_data, aes(x = bio6, y = species, fill = decision)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, alpha = 0.7,  outlier.colour = "red") +
  geom_vline(xintercept = NORWAY_BIO6_MAX, linetype = "dashed", 
             color = "steelblue", linewidth = 1) +
  annotate("text", x = NORWAY_BIO6_MAX, y = 0.5, 
           label = sprintf("Norway max (%.0f)", NORWAY_BIO6_MAX),
           hjust = -0.05, color = "steelblue", size = 3, fontface = "bold") +
  scale_fill_manual(values = plot_colors) +
  labs(
    x = "Bio6: Min Temperature of Coldest Month",
    y = NULL,
    title = "Climate Prescreening: Species vs Norway",
    subtitle = "Species ordered by coldest occurrence | Dashed line = Norway's mildest winter",
    fill = "Decision"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 8),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

print(p_overview)

ggsave("./5_outputs/screen_output/prescreening_overview.png", p_overview, 
       width = 12, height = max(8, length(species_order) * 0.35), dpi = 150, limitsize = FALSE)

message("Saved: ./5_outputs/screen_output/prescreening_overview.png")

# =============================================================================
# 9. SCATTER PLOT: BIO1 vs BIO6 (Climate Space) WITH ERROR BARS
# =============================================================================


message("Generating climate space scatter plot...")

scatter_data <- results %>%
  filter(status == "processed") %>%
  mutate(decision = factor(decision, levels = c("EXCLUDE", "REVIEW", "INCLUDE_MARGINAL", "INCLUDE")))

# Norway climate box boundaries (convert to °C)
norway_box <- tibble(
  xmin = NORWAY_BIO1_MIN / 10,
  xmax = NORWAY_BIO1_MAX / 10,
  ymin = NORWAY_BIO6_MIN / 10,
  ymax = NORWAY_BIO6_MAX / 10
)

p_scatter <- ggplot(scatter_data, aes(x = bio1_min/10, y = bio6_min/10, color = decision)) +
  # Norway climate envelope (rectangle)
  annotate("rect", 
           xmin = norway_box$xmin, xmax = norway_box$xmax,
           ymin = norway_box$ymin, ymax = norway_box$ymax,
           fill = "steelblue", alpha = 0.15, color = "steelblue", 
           linewidth = 1.2, linetype = "solid") +
  # Horizontal error bars (bio1 range: min to max)
  geom_errorbar(aes(xmin = bio1_min/10, xmax = bio1_max/10), height = 0.5, alpha = 0.4, linewidth = 0.5) +
  # Vertical error bars (bio6 range: min to max)
  geom_errorbar(aes(ymin = bio6_min/10, ymax = bio6_max/10), width = 0.5, alpha = 0.4, linewidth = 0.5) +
  # Species points (at minimum values)
  geom_point(aes(size = n_with_climate), alpha = 0.8) +
  # Species labels
  geom_text(aes(label = species, size = 2.5, vjust = -1.2, check_overlap = TRUE)) +
  # Colors
  scale_color_manual(values = c(
    "EXCLUDE" = "#d73027", 
    "REVIEW" = "#fee090",
    "INCLUDE_MARGINAL" = "#a6d96a", 
    "INCLUDE" = "#1a9850"
  )) +
  scale_size_continuous(range = c(3, 10), name = "N occurrences") +
  labs(
    x = "Mean Annual Temperature (bio1, C)",
    y = "Coldest Month Temperature (bio6, C)",
    title = "Climate Space Prescreening: Species vs Norway",
    subtitle = "Error bars show full climate range per species | Blue box = Norway's climate range",
    color = "Decision"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off")

print(p_scatter)

ggsave("./5_outputs/screen_output/prescreening_climate_scatter.png", p_scatter, width = 12, height = 10, dpi = 150)

message("Saved: ./5_outputs/screen_output/prescreening_climate_scatter.png")
# =============================================================================
# 10. SPECIES OCCURRENCE MAPS
# =============================================================================

message("\nGenerating species occurrence maps...")

dir.create("./5_outputs/screen_output/species_maps", showWarnings = FALSE)

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Loop over species with occurrence data
for (sp in names(species_coords_data)) {
  
  message(sprintf("  Mapping: %s", sp))
  
  sp_coords <- species_coords_data[[sp]]
  sp_result <- results %>% filter(species == sp)
  sp_decision <- sp_result$decision
  sp_n <- nrow(sp_coords)
  sp_color <- decision_colors[sp_decision]
  
  # Create map
  p_map <- ggplot() +
    # World background
    geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
    # Occurrence points
    geom_point(data = sp_coords,
               aes(x = decimalLongitude, y = decimalLatitude),
               color = sp_color, alpha = 0.6, size = 1.5) +
    # Labels
    labs(
      title = sp,
      subtitle = sprintf("%s | n = %d | bio6_min = %.1f°C | Sources: %s",
                         sp_decision, sp_n, sp_result$bio6_min/10, sp_result$sources)
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Save
  filename_safe <- str_replace_all(sp, " ", "_")
  ggsave(sprintf("./5_outputs/screen_output/species_maps/%s.png", filename_safe),
         p_map, width = 12, height = 7, dpi = 150)
}

message(sprintf("Saved %d occurrence maps to ./5_outputs/screen_output/species_maps/", length(species_coords_data)))

# =============================================================================
# 11. EXPORT RESULTS
# =============================================================================

message("\nExporting results...")

# Add Norway reference values to results for comparison
results_with_norway <- results %>%
  mutate(
    norway_bio1_min = NORWAY_BIO1_MIN,
    norway_bio1_max = NORWAY_BIO1_MAX,
    norway_bio6_min = NORWAY_BIO6_MIN,
    norway_bio6_max = NORWAY_BIO6_MAX
  ) %>%
  arrange(decision, bio6_min)

export(
  list(
    "Prescreening_Results" = results_with_norway,
    "Species_To_Model" = tibble(species = species_to_model),
    "Excluded_Species" = results %>% filter(decision == "EXCLUDE") %>% 
      select(species, bio1_min, bio6_min, sources, reason),
    "Norwegian_Natives" = norwegian_keys,
    "Norway_Climate" = norway_climate_summary,
    "Source_Summary" = source_summary,
    "Summary" = summary_table
  ),
  "./5_outputs/screen_output/prescreening_bio6_spocc.xlsx"
)