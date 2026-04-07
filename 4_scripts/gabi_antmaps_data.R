library(data.table)
library(tidyverse)

# Data downloaded from 
# https://github.com/nitishnarula/antmaps-app/blob/master/about.html
# https://antmaps.org/about.html

# fread is the fastest for big CSVs
gabi <- fread("./1_raw_data/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020.csv")
# Check size and structure
dim(gabi)
names(gabi)
head(gabi)

gabi <- as.data.frame(gabi)
gabi <- as_tibble(gabi)

# Fix valid_species_name: replace "." with " "
gabi$valid_species_name <- gsub("\\.", " ", gabi$valid_species_name)

head(gabi)

# 1. Check for dubious records
table(gabi$dubious)

dub <- c("Dubious", "Error", "Historical Dubious", "Redundant")

gabi_clean <- gabi %>% filter(!dubious %in% dub)

# -------------------------------------------------------------------------
library(rgbif)
library(parallel)
library(pbapply)

# Get unique species names
species_list <- unique(gabi_clean$valid_species_name)
cat("Unique species:", length(species_list), "\n")

# Function to get GBIF accepted key
get_gbif_key <- function(species_name) {
  tryCatch({
    result <- name_backbone(name = species_name, kingdom = "Animalia")
    data.frame(
      valid_species_name = species_name,
      acceptedKey = result$usageKey,
      status = result$status,
      matchType = result$matchType,
      gbif_name = result$canonicalName,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      valid_species_name = species_name,
      acceptedKey = NA,
      status = NA,
      matchType = NA,
      gbif_name = NA,
      stringsAsFactors = FALSE
    )
  })
}

# Setup multicore cluster (use n-1 cores)
n_cores <- detectCores() - 5
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(rgbif))

# Run with progress bar and parallel
gbif_keys <- pblapply(species_list, get_gbif_key, cl = cl)
gbif_keys <- do.call(rbind, gbif_keys)

# Stop cluster
stopCluster(cl)

# Check results
head(gbif_keys)
table(gbif_keys$status, useNA = "ifany")
cat("Species matched:", sum(!is.na(gbif_keys$acceptedKey)), "/", nrow(gbif_keys), "\n")

# Merge back to gabi_clean
gabi_clean <- left_join(gabi_clean, gbif_keys, by = "valid_species_name")

# -------------------------------------------------------------------------
rio::export(gabi_clean, "./2_processed_data/gabi_antmaps_data_clean.csv")
# check 
gabi_clean %>% filter(valid_species_name=="Camponotus cruentatus")

# plot check -------------------------------------------------------------------------
imparis <- gabi_clean %>% filter(valid_species_name=="Prenolepis imparis")
head(imparis)

library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "white") +
  geom_point(data = imparis, aes(x = dec_long, y = dec_lat, color = valid_species_name ), alpha = 0.6) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude")

imparis %>% count()
# END ---------------------------------------------------------------------

