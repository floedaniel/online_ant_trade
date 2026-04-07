# =============================================================================
# HOLDRIDGE TRIANGLE PLOTS FOR ANT SPECIES
# =============================================================================

library(Ternary)
library(tidyverse)
library(terra)
library(geodata)
library(spocc)

data(holdridgeLifeZonesUp, package = "Ternary")

dir.create("./holdridge_output/holdridge_plots", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. PARAMETERS
# =============================================================================

SPECIES_LIST <- c(
  "Formica rufa",
  "Lasius niger", 
  "Wasmannia auropunctata"
)

OCC_LIMIT <- 500

# =============================================================================
# 2. LOAD CLIMATE DATA
# =============================================================================

message("Loading WorldClim data...")

bioclim <- worldclim_global(var = "bio", res = 2.5, path = "data/worldclim")
bio1 <- bioclim[[1]]   # Mean annual temperature
bio12 <- bioclim[[12]] # Annual precipitation

# =============================================================================
# 3. FETCH OCCURRENCES
# =============================================================================

fetch_occurrences <- function(species_name, limit = 500) {
  
  message(sprintf("  Fetching: %s", species_name))
  
  occ_data <- tryCatch({
    spocc::occ(query = species_name, from = "gbif", limit = limit, has_coords = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(occ_data)) return(NULL)
  
  occ_df <- tryCatch({
    spocc::occ2df(occ_data)
  }, error = function(e) NULL)
  
  if (is.null(occ_df) || nrow(occ_df) == 0) return(NULL)
  
  occ_df <- occ_df %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    filter(abs(latitude) <= 90, abs(longitude) <= 180) %>%
    filter(!(latitude == 0 & longitude == 0)) %>%
    mutate(
      species = species_name,
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude)
    ) %>%
    select(species, longitude, latitude) %>%
    distinct()
  
  message(sprintf("    Retrieved %d records", nrow(occ_df)))
  return(occ_df)
}

message("\nFetching occurrence data...")
all_occurrences <- map_dfr(SPECIES_LIST, ~fetch_occurrences(.x, OCC_LIMIT))
message(sprintf("Total occurrences: %d", nrow(all_occurrences)))

# =============================================================================
# 4. EXTRACT CLIMATE VALUES
# =============================================================================

message("\nExtracting climate values...")

coords <- as.matrix(all_occurrences[, c("longitude", "latitude")])
all_occurrences$mat <- terra::extract(bio1, coords)[[1]]
all_occurrences$precipitation <- terra::extract(bio12, coords)[[1]]

# Calculate Holdridge variables
all_occurrences <- all_occurrences %>%
  filter(!is.na(mat), !is.na(precipitation), precipitation > 0) %>%
  mutate(
    bioT = pmax(0, pmin(30, mat)),
    pet_ratio = (bioT * 58.93) / precipitation
  ) %>%
  filter(pet_ratio > 0)

message(sprintf("Valid points with climate: %d", nrow(all_occurrences)))

# =============================================================================
# 5. PLOT FUNCTION
# =============================================================================

plot_holdridge <- function(occ_data, species_name, output_dir = "./holdridge_output/holdridge_plots") {
  
  sp <- occ_data %>% filter(species == species_name)
  
  if (nrow(sp) < 5) {
    message(sprintf("  Skipping %s (n=%d)", species_name, nrow(sp)))
    return(NULL)
  }
  
  filename <- str_replace_all(species_name, "[^A-Za-z0-9]", "_")
  filepath <- file.path(output_dir, paste0(filename, "_holdridge.png"))
  
  png(filepath, width = 1200, height = 1000, res = 120)
  par(mar = c(0, 0, 2, 0))
  
  HoldridgePlot(hex.labels = holdridgeLifeZonesUp, hex.cex = 0.55)
  HoldridgeBelts()
  
  HoldridgePoints(
    pet = sp$pet_ratio,
    prec = sp$precipitation,
    col = adjustcolor("gray50", 0.2),
    pch = 16,
    cex = 1
  )
  
  title(main = bquote(italic(.(species_name)) ~ "(n =" ~ .(nrow(sp)) ~ ")"),
        line = 0.5, cex.main = 1.3)
  
  dev.off()
  
  message(sprintf("  Saved: %s (n=%d)", species_name, nrow(sp)))
  return(tibble(species = species_name, n = nrow(sp), filepath = filepath))
}

# =============================================================================
# 6. GENERATE ALL PLOTS
# =============================================================================

message("\nGenerating Holdridge plots...")

species_list <- unique(all_occurrences$species)

results <- map_dfr(species_list, ~plot_holdridge(all_occurrences, .x))

message(sprintf("\nDone! Generated %d plots", nrow(results)))

print(results)
