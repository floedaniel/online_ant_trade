# =============================================================================
# ANT SPECIES CLASSIFICATION BY GLOBAL ECOZONE
# =============================================================================

library(tidyverse)
library(sf)
library(spocc)
library(rgbif)
library(rio)

# =============================================================================
# 1. PARAMETERS
# =============================================================================

# Path to ecozone shapefile (adjust for your system)
ECOZONE_PATH <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/Daniel jobb/Climat layers files/Present/eco_zone_shp/eco_zone.shp"

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
SPECIES_LIST <- master_data %>%
  filter(status == "ACCEPTED") %>%
  filter(!acceptedKey %in% norwegian_accepted_keys) %>%
  pull(name)



# Occurrence fetch limit per species
OCC_LIMIT <- 100000

# Output directory
dir.create("./ecozone_output", showWarnings = FALSE)

# =============================================================================
# 2. LOAD ECOZONE SHAPEFILE
# =============================================================================

message("Loading ecozone shapefile...")

ecozones <- st_read(ECOZONE_PATH)

message(sprintf("Loaded %d ecozone polygons", nrow(ecozones)))

# Use GEZ_TERM as the zone identifier
zone_col <- "GEZ_TERM"

message(sprintf("\nUnique ecozones:"))
print(unique(ecozones[[zone_col]]))

# Ensure valid CRS (WGS84)
if (is.na(st_crs(ecozones))) {
  ecozones <- st_set_crs(ecozones, 4326)
} else if (st_crs(ecozones)$epsg != 4326) {
  ecozones <- st_transform(ecozones, 4326)
}

# Fix any invalid geometries
ecozones <- st_make_valid(ecozones)

message("Ecozone data ready")

# =============================================================================
# 3. FETCH OCCURRENCE FUNCTION
# =============================================================================

fetch_occurrences <- function(species_name, limit = 5000) {
  
  message(sprintf("  Fetching: %s", species_name))
  
  occ_data <- tryCatch({
    spocc::occ(
      query = species_name,
      from = c("gbif"),
      limit = limit,
      has_coords = TRUE
    )
  }, error = function(e) {
    message(sprintf("    ERROR: %s", e$message))
    return(NULL)
  })
  
  if (is.null(occ_data)) return(NULL)
  
  occ_df <- tryCatch({
    spocc::occ2df(occ_data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(occ_df) || nrow(occ_df) == 0) return(NULL)
  
  # Clean and standardize
  occ_df <- occ_df %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    filter(abs(latitude) <= 90, abs(longitude) <= 180) %>%
    filter(!(latitude == 0 & longitude == 0)) %>%
    mutate(
      species = species_name,
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude)
    ) %>%
    select(species, longitude, latitude)
  
  message(sprintf("    Retrieved %d records", nrow(occ_df)))
  
  return(occ_df)
}

# =============================================================================
# 4. FETCH ALL SPECIES OCCURRENCES
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("Fetching occurrence data for ", length(SPECIES_LIST), " species...")
message(paste(rep("=", 60), collapse = ""))

all_occurrences <- map_dfr(SPECIES_LIST, function(sp) {
  fetch_occurrences(sp, limit = OCC_LIMIT)
})

message(sprintf("\nTotal occurrences: %d", nrow(all_occurrences)))
message(sprintf("Species with data: %d", n_distinct(all_occurrences$species)))

# =============================================================================
# 5. SPATIAL JOIN: ASSIGN ECOZONES TO OCCURRENCES
# =============================================================================

message("\nAssigning ecozones to occurrences...")

# Convert occurrences to sf object
occ_sf <- st_as_sf(all_occurrences, 
                   coords = c("longitude", "latitude"), 
                   crs = 4326)

# Spatial join - find which ecozone each point falls in
occ_with_zone <- st_join(occ_sf, ecozones[, zone_col], left = TRUE)

# Convert back to dataframe with coordinates
occ_results <- occ_with_zone %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  rename(ecozone = all_of(zone_col))

# Handle points outside any ecozone (ocean, etc.)
occ_results <- occ_results %>%
  mutate(ecozone = ifelse(is.na(ecozone), "Outside/Ocean", ecozone))

message(sprintf("Assigned ecozones to %d occurrences", nrow(occ_results)))

# =============================================================================
# 6. SUMMARIZE BY SPECIES AND ECOZONE
# =============================================================================

message("\nGenerating summaries...")

# Summary: count by species and ecozone
species_ecozone_summary <- occ_results %>%
  group_by(species, ecozone) %>%
  summarise(n_occurrences = n(), .groups = "drop") %>%
  arrange(species, desc(n_occurrences))

# Wide format: species × ecozone matrix
species_ecozone_wide <- species_ecozone_summary %>%
  pivot_wider(
    names_from = ecozone,
    values_from = n_occurrences,
    values_fill = 0
  )

# Calculate dominant ecozone per species
dominant_ecozone <- species_ecozone_summary %>%
  group_by(species) %>%
  slice_max(n_occurrences, n = 1) %>%
  ungroup() %>%
  rename(dominant_ecozone = ecozone, n_in_dominant = n_occurrences)

# Species summary with total counts
species_summary <- occ_results %>%
  group_by(species) %>%
  summarise(
    n_total = n(),
    n_ecozones = n_distinct(ecozone),
    .groups = "drop"
  ) %>%
  left_join(dominant_ecozone, by = "species") %>%
  mutate(pct_in_dominant = 100 * n_in_dominant / n_total) %>%
  arrange(desc(n_total))

# Ecozone summary
ecozone_summary <- occ_results %>%
  group_by(ecozone) %>%
  summarise(
    n_occurrences = n(),
    n_species = n_distinct(species),
    species_list = paste(unique(species), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_occurrences))

# =============================================================================
# 7. PRINT RESULTS
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("RESULTS: SPECIES BY ECOZONE")
message(paste(rep("=", 60), collapse = ""))

message("\n--- Species Summary ---")
print(species_summary, n = 20)

message("\n--- Ecozone Summary ---")
print(ecozone_summary, n = 20)

message("\n--- Detailed Species × Ecozone Counts ---")
print(species_ecozone_summary, n = 50)

# =============================================================================
# 8. VISUALIZATION
# =============================================================================

message("\nGenerating plots...")

library(rnaturalearth)

# Get world basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# Color palette for ecozones (consistent across all plots)
all_zones <- unique(ecozones[[zone_col]])
n_zones <- length(all_zones)
zone_colors <- setNames(scales::hue_pal()(n_zones), all_zones)

# Create individual maps per species with only relevant biomes
dir.create("./ecozone_output/species_maps", showWarnings = FALSE)

for (sp in unique(occ_results$species)) {
  
  message(sprintf("  Mapping: %s", sp))
  
  # Get occurrences for this species
  sp_occ <- occ_results %>% filter(species == sp)
  
  # Get unique ecozones for this species (excluding Outside/Ocean)
  sp_zones <- sp_occ %>%
    filter(ecozone != "Outside/Ocean") %>%
    pull(ecozone) %>%
    unique()
  
  # Subset ecozone polygons to only those where species occurs
  sp_ecozones <- ecozones %>%
    filter(.data[[zone_col]] %in% sp_zones)
  
  # Count per zone for legend
  zone_counts <- sp_occ %>%
    filter(ecozone != "Outside/Ocean") %>%
    count(ecozone) %>%
    mutate(label = sprintf("%s (n=%d)", ecozone, n))
  
  # Create color mapping for this species' zones
  sp_zone_colors <- zone_colors[sp_zones]
  
  # Build the map
  p_map <- ggplot() +
    # World background
    geom_sf(data = world, fill = "gray95", color = "gray80", linewidth = 0.1) +
    # Ecozone polygons (only species-relevant ones)
    geom_sf(data = sp_ecozones, aes(fill = .data[[zone_col]]),  alpha = 0.4, color = "gray60", linewidth = 0.2) +
    # Occurrence points
    geom_point(data = sp_occ, aes(x = longitude, y = latitude, color = ecozone), alpha = 0.6, size = .5) +
    # Colors
    scale_fill_manual(values = sp_zone_colors, name = "Ecozone") +
    scale_color_manual(values = c(zone_colors, "Outside/Ocean" = "gray40"), 
                       name = "Occurrences") +
    # Labels
    labs(
      title = sp,
      subtitle = sprintf("n = %d occurrences across %d ecozones", 
                         nrow(sp_occ), length(sp_zones))
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1))
    )
  
  # Save individual map
  filename_safe <- str_replace_all(sp, "[^A-Za-z0-9_-]", "_")
  ggsave(sprintf("./ecozone_output/species_maps/%s.png", filename_safe),
         p_map, width = 14, height = 8, dpi = 150)
}

message(sprintf("Saved %d individual species maps", n_distinct(occ_results$species)))


# =============================================================================
# 9. EXPORT RESULTS
# =============================================================================

message("\nExporting results...")

export(
  list(
    "Species_Summary" = species_summary,
    "Ecozone_Summary" = ecozone_summary,
    "Species_Ecozone_Counts" = species_ecozone_summary,
    "Species_Ecozone_Matrix" = species_ecozone_wide,
    "All_Occurrences" = occ_results
  ),
  "./ecozone_output/ant_ecozone_analysis.xlsx"
)

message("Saved: ./ecozone_output/ant_ecozone_analysis.xlsx")

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("Total species analyzed: %d", n_distinct(occ_results$species)))
message(sprintf("Total occurrences: %d", nrow(occ_results)))
message(sprintf("Ecozones represented: %d", n_distinct(occ_results$ecozone)))
message(paste(rep("=", 60), collapse = ""))