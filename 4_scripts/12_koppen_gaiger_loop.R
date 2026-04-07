# -------------------------------------------------------------------------
# Load required libraries
# -------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(rgbif)
library(CoordinateCleaner)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(sf)
library(kgc)
library(tidyverse)
library(future.apply) # For parallel processing
library(parallel)

# Set up parallel processing using available cores (reserve 2 cores)
plan(multisession, workers = detectCores() - 2)

# -------------------------------------------------------------------------
# Set working directories and master file path
# -------------------------------------------------------------------------
base_dir <- "./Species"
master_file_path <- file.path(base_dir, "2_processed_data/master_species.xlsx")

# -------------------------------------------------------------------------
# Read the master species file (using columns: speciesKey and species)
# -------------------------------------------------------------------------
master_species <- rio::import(master_file) %>% 
  as_tibble() %>% 
  filter(!is.na(speciesKey) & !is.na(species))

# -------------------------------------------------------------------------
# Parallel loop through each species
# -------------------------------------------------------------------------
results <- future_lapply(seq_len(nrow(master_species)), function(i) {
  tryCatch({
    # Retrieve speciesKey and species name from master file
    speciesKey <- as.character(master_species$speciesKey[i])
    species_name <- master_species$species[i]
    message("Processing species: ", species_name, " (SpeciesKey: ", speciesKey, ")")
    
    # Create folder structure for species and Koppen-Geiger outputs
    species_folder <- file.path(base_dir, speciesKey)
    sdm_folder <- file.path(species_folder, "koppen_geiger")
    if (!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)
    if (!dir.exists(sdm_folder)) dir.create(sdm_folder, recursive = TRUE)
    
    # -----------------------------------------------------------------------
    # Retrieve GBIF Occurrence Data
    # -----------------------------------------------------------------------
    gbif_taxon_key <- tryCatch({
      result <- name_backbone(name = species_name)
      if (!is.null(result) && "usageKey" %in% names(result)) {
        result$usageKey
      } else {
        NA
      }
    }, error = function(e) {
      message("Error retrieving taxon key for species: ", species_name, " - ", e$message)
      NA
    })
    
    occ_data <- tryCatch({
      occ_data(taxonKey = gbif_taxon_key, limit = 100000)$data
    }, error = function(e) NULL)
    
    # Check if occurrence data is valid (at least 10 records)
    if (is.null(occ_data) || nrow(occ_data) < 10) {
      message("No occurrence data found for species: ", species_name)
      writeLines("No occurrence data found in GBIF for this species.", 
                 file.path(sdm_folder, "no_occurrence_data.txt"))
      next  # Proceed to the next species
    }
    
    # Process occurrence data if available
    occ_df <- tryCatch({
      occ_data %>%
        as_tibble() %>%
        dplyr::select(scientificName, decimalLongitude, decimalLatitude) %>%
        dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
        dplyr::mutate(species_name = species_name) %>%
        dplyr::select(species_name, decimalLongitude, decimalLatitude)
    }, error = function(e) {
      message("Error processing occurrence data for species: ", species_name)
      writeLines("Error processing occurrence data for this species.", 
                 file.path(sdm_folder, "error_in_data_processing.txt"))
      NULL
    })
    
    # Skip if no occurrence data was processed
    if (is.null(occ_df)) {
      next  # Proceed to the next species
    }
    
    names(occ_df) <- c("species", "decimalLongitude", "decimalLatitude")
    
    # Coordinate cleaning
    occ_df <- occ_df %>%
      dplyr::mutate(
        rndCoord.lat = RoundCoordinates(decimalLatitude),
        rndCoord.lon = RoundCoordinates(decimalLongitude)
      ) %>%
      dplyr::filter(!is.na(rndCoord.lat) & !is.na(rndCoord.lon))
    
    kg <- data.frame(occ_df, layer_climate = LookupCZ(occ_df))
    
    # Filter out specific unwanted climate types
    kg_data <- kg %>%
      dplyr::filter(!layer_climate %in% c("Ocean", "ET", "EF", "Climate Zone info missing"))
    
    # -----------------------------------------------------------------------
    # Read Koppen-Geiger Shapefile and filter by kg_data
    # -----------------------------------------------------------------------
    kg_path <- "C:/Users/dafl/Dropbox/Climat layers files/WC05_1975H_Koppen_Shapefile/WC05_1975H_Koppen.shp"
    shapefile_kg <- st_read(kg_path) %>%
      st_set_crs(st_crs("+proj=longlat +datum=WGS84 +no_defs"))
    
    subset_kg <- shapefile_kg %>% dplyr::filter(Koppen %in% kg_data$layer_climate)
    
    # -----------------------------------------------------------------------
    # Load World Boundaries
    # -----------------------------------------------------------------------
    world <- ne_countries(scale = "large", returnclass = "sf") %>%
      dplyr::filter(admin != "Antarctica") %>%
      st_transform("+proj=longlat +datum=WGS84 +no_defs")
    
    # -----------------------------------------------------------------------
    # Define manual colors for Koppen climate types
    # -----------------------------------------------------------------------
    codes <- data.frame(
      Koppen = c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 
                 'Cfa', 'Cfb', 'Cfc', 'Csa', 'Csb', 'Csc', 'Cwa', 'Cwb', 'Cwc', 
                 'Dfa', 'Dfb', 'Dfc', 'Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd', 
                 'Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF', 'ET', 'Ocean'),
      hex = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", 
              "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", 
              "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", 
              "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")
    )
    
    # Filter subset_kg to include only valid Koppen levels
    valid_koppen <- unique(codes$Koppen)
    subset_kg <- subset_kg %>%
      dplyr::filter(Koppen %in% valid_koppen)
    
    # Fix encoding and trim spaces
    species_name <- iconv(species_name, from = "UTF-8", to = "UTF-8")
    subset_kg$Koppen <- iconv(subset_kg$Koppen, from = "UTF-8", to = "UTF-8")
    subset_kg$Koppen <- trimws(subset_kg$Koppen)
    codes$Koppen <- trimws(codes$Koppen)
    
    # -----------------------------------------------------------------------
    # Summarize climate counts from occurrence data
    # -----------------------------------------------------------------------
    climate_counts <- kg_data %>%
      dplyr::group_by(layer_climate) %>%
      dplyr::summarise(Count = n()) %>%
      dplyr::ungroup()
    
    # Merge codes with climate_counts to ensure all Koppen climate types are represented
    merged_data <- merge(codes, climate_counts, by.x = "Koppen", by.y = "layer_climate", all.x = TRUE)
    merged_data[is.na(merged_data$Count), "Count"] <- 0  # Fill missing counts with zero
    
    # -----------------------------------------------------------------------
    # Define a custom theme for plots
    # -----------------------------------------------------------------------
    custom_theme <- theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold"),
      plot.subtitle = element_text(hjust = 0),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )
    
    # Define plot titles using bquote for italicized species name
    histogram_title <- bquote("Histogram of occurrences in Koppen-Geiger climate zones for " ~ italic(.(species_name)))
    plot_title <- bquote("Species occurrences and Koppen-Geiger climate zones for " ~ italic(.(species_name)))
    
    # -----------------------------------------------------------------------
    # Generate the climate map plot (plot1)
    # -----------------------------------------------------------------------
    plot1 <- ggplot() + 
      geom_sf(data = world, color = "#000000", fill = "gray100", size = 0.05, inherit.aes = FALSE) +
      geom_sf(data = subset_kg, aes(fill = Koppen, col = Koppen), size = 0.05) +
      scale_fill_manual(values = setNames(codes$hex, codes$Koppen)) +
      scale_color_manual(values = setNames(codes$hex, codes$Koppen)) +
      geom_point(data = occ_df, aes(x = decimalLongitude, y = decimalLatitude), size = 2, col = "black") +
      theme(legend.position = "bottom") +
      labs(
        title = paste("Koppen-Geiger Climate Zones and Species Occurrences:", species_name),
        x = "Longitude",
        y = "Latitude",
        fill = "Climate Zone",
        color = "Climate Zone"
      ) +
      custom_theme 
    
    # -----------------------------------------------------------------------
    # Generate the climate histogram plot (plot2)
    # -----------------------------------------------------------------------
    plot2 <- ggplot(data = merged_data, aes(x = Koppen, y = Count, fill = Koppen)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = merged_data$hex) +
      labs(title = histogram_title,
           subtitle = "Frequency distribution across different climate types",
           x = "Koppen Climate Type", 
           y = "Count",  
           caption = "VKM 2023\nData: GBIF, Beck et al. 2023, Natural Earth") +
      custom_theme +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))
    
    # -----------------------------------------------------------------------
    # Highlight specific Koppen groups in the histogram (if available)
    # -----------------------------------------------------------------------
    codes_to_highlight_c <- c("Cfb", "Cfc")
    codes_to_highlight_d <- c("Dfb", "Dfc")
    positions_c <- which(merged_data$Koppen %in% codes_to_highlight_c)
    positions_d <- which(merged_data$Koppen %in% codes_to_highlight_d)
    
    if(length(positions_c) > 0) {
      xmin <- min(positions_c) - 0.5
      xmax <- max(positions_c) + 0.5
      
      plot2 <- plot2 + 
        geom_rect(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text", x = mean(c(xmin, xmax)), y = Inf, label = "Norway", vjust = +1, size = 5, angle = 0)
    }
    
    if(length(positions_d) > 0) {
      xmin <- min(positions_d) - 0.5
      xmax <- max(positions_d) + 0.5
      
      plot2 <- plot2 + 
        geom_rect(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text", x = mean(c(xmin, xmax)), y = Inf, label = "Norway", vjust = +1, size = 5, angle = 0)
    }
    
    # -----------------------------------------------------------------------
    # Combine the two plots vertically using ggpubr
    # -----------------------------------------------------------------------
    combined_plot <- ggpubr::ggarrange(plot1, plot2, ncol = 1, nrow = 2, align = "v", heights = c(1, 1), widths  = c(1, 1))
    
    # Save the combined plot in the species-specific Koppen-Geiger folder
    output_file <- file.path(sdm_folder, paste0("koppen_classification_", speciesKey, ".png"))
    ggsave(plot = combined_plot, filename = output_file, dpi = 300, width = 30, height = 30, units = "cm", scale = 1, bg = "white")
    
    return(TRUE) # Return success indicator
  }, error = function(e) {
    message(paste("Error processing", i, ":", e$message))
    return(NULL) # Return NULL on error
  })
})

# Reset to default plan after processing is complete
plan(sequential)
