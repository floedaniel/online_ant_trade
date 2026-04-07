##############################################################################
#   GLOBAL PCA (6 Bioclim Vars) + REAL SPECIES + GLOBAL + NORWAY ENVELOPES
#    Panel A is a fixed “Norway Envelope” (always the same) with species points
#    Panel B is each species' min–max envelope in Norway
#    Panel C & D are PCA-based global & Norway
#    Panel E & F are boxplot & scatter
##############################################################################

# -------------------------------------------------------------------------
# 0) Load Required Libraries
# -------------------------------------------------------------------------
library(geodata)       # WorldClim + GADM
library(terra)         # SpatRaster
library(raster)        # RasterStack for virtualspecies
library(sf)            # Spatial vector handling
library(rgbif)         # GBIF occurrences
library(dplyr)         # Data wrangling
library(ggplot2)       # Plotting
library(tidyr)         # Reshaping
library(cowplot)       # Combine plots
library(ade4)          # suprow() for PCA
library(virtualspecies)# PCA-based niche
library(tidyterra)     # Optional for geom_spatraster
library(rio)           # For importing master file

# -------------------------------------------------------------------------
# 1) Load Master Data and Set Base Directory for Species Outputs
# -------------------------------------------------------------------------
master_file <- "./2_processed_data/master_data.xlsx"
master_data <- rio::import(master_file) %>% 
  as_tibble() %>% 
  filter(!is.na(speciesKey) & !is.na(species))

# Define your species base folder (each species will be stored in its own speciesKey folder)
base_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"

# -------------------------------------------------------------------------
# 2) Norway Envelope (fixed) at 10m resolution: bio1,bio12,tmin,tmax + buffer
#    We'll show this map the same for every species, only points change
# -------------------------------------------------------------------------
bioclim_10m  <- geodata::worldclim_global(var="bio",  res=10, path=getwd())
tmin_10m_all <- geodata::worldclim_global(var="tmin", res=10, path=getwd())
tmax_10m_all <- geodata::worldclim_global(var="tmax", res=10, path=getwd())

norway_spat_10m <- geodata::gadm("GADM", country="NOR", level=0)
# Crop+mask for local stats
norway_bioclim <- mask(crop(bioclim_10m,  norway_spat_10m), norway_spat_10m)
norway_tmin    <- mask(crop(tmin_10m_all, norway_spat_10m), norway_spat_10m)
norway_tmax    <- mask(crop(tmax_10m_all, norway_spat_10m), norway_spat_10m)

# Extract single layers for reference
global_bio1  <- bioclim_10m[["wc2.1_10m_bio_1"]]
global_bio12 <- bioclim_10m[["wc2.1_10m_bio_12"]]
global_tmin  <- tmin_10m_all[["wc2.1_10m_tmin_01"]]
global_tmax  <- tmax_10m_all[["wc2.1_10m_tmax_07"]]

# Norway stats
bio1_norway  <- norway_bioclim[["wc2.1_10m_bio_1"]]
bio12_norway <- norway_bioclim[["wc2.1_10m_bio_12"]]
tmin_norway  <- norway_tmin[["wc2.1_10m_tmin_01"]]
tmax_norway  <- norway_tmax[["wc2.1_10m_tmax_07"]]

bio1_min  <- terra::minmax(bio1_norway)[1]
bio1_max  <- terra::minmax(bio1_norway)[2]
bio12_min <- terra::minmax(bio12_norway)[1]
bio12_max <- terra::minmax(bio12_norway)[2]
tmin_min  <- terra::minmax(tmin_norway)[1]
tmin_max  <- terra::minmax(tmin_norway)[2]
tmax_min  <- terra::minmax(tmax_norway)[1]
tmax_max  <- terra::minmax(tmax_norway)[2]

buffer <- 1
bio1_lower  <- bio1_min  - buffer
bio1_upper  <- bio1_max  + buffer
bio12_lower <- bio12_min - buffer
bio12_upper <- bio12_max + buffer
tmin_lower  <- tmin_min  - buffer
tmin_upper  <- tmin_max  + buffer
tmax_lower  <- tmax_min  - buffer
tmax_upper  <- tmax_max  + buffer

# Create a global mask for Norway's envelope
mask_norway <- (global_bio1 >= bio1_lower & global_bio1 <= bio1_upper) &
  (global_bio12>= bio12_lower & global_bio12<= bio12_upper) &
  (global_tmin>= tmin_lower   & global_tmin<= tmin_upper)
# optionally add tmax:
# & (global_tmax>= tmax_lower & global_tmax<= tmax_upper)

combined_norway_mask <- terra::ifel(mask_norway, 1, 0)

# We'll crop to lat> -60 if we want to remove Antarctica
norway_env_crop <- crop(combined_norway_mask, ext(-180,180,-60,90))
norway_env_crop <- stats::setNames(norway_env_crop, "value")

# Convert to df for plotting
norway_env_df <- as.data.frame(norway_env_crop, xy=TRUE, na.rm=FALSE)

base_norway_envelope_map <- ggplot() +
  geom_raster(data=norway_env_df, aes(x=x,y=y,fill=value)) +
  scale_fill_viridis_c(name="Norway Envelope") +
  coord_quickmap() +
  labs(
    title="Global Map of Norway Envelope (Fixed)",
    x="Longitude", y="Latitude"
  ) +
  theme_minimal()

# -------------------------------------------------------------------------
# 3) Build a Single Global PCA for 6 Vars (5m)
# -------------------------------------------------------------------------
all_bio_5m <- geodata::worldclim_global(var="bio", res=5, path=getwd())
# pick bio1,bio2,bio5,bio6,bio12,bio15
bio1_5m  <- all_bio_5m[["wc2.1_5m_bio_1"]]
bio2_5m  <- all_bio_5m[["wc2.1_5m_bio_2"]]
bio5_5m  <- all_bio_5m[["wc2.1_5m_bio_5"]]
bio6_5m  <- all_bio_5m[["wc2.1_5m_bio_6"]]
bio12_5m <- all_bio_5m[["wc2.1_5m_bio_12"]]
bio15_5m <- all_bio_5m[["wc2.1_5m_bio_15"]]

# Create a RasterStack (required by virtualspecies) using the raster package
bio1_r5  <- raster(bio1_5m)
bio2_r5  <- raster(bio2_5m)
bio5_r5  <- raster(bio5_5m)
bio6_r5  <- raster(bio6_5m)
bio12_r5 <- raster(bio12_5m)
bio15_r5 <- raster(bio15_5m)

env_stack_global <- stack(bio1_r5, bio2_r5, bio5_r5, bio6_r5, bio12_r5, bio15_r5)
names(env_stack_global) <- c("bio1","bio2","bio5","bio6","bio12","bio15")

global_pca_temp <- generateSpFromPCA(
  raster.stack  = env_stack_global,
  sample.points = TRUE,
  nb.points     = 10000,
  plot          = FALSE
)
my_pca_global <- global_pca_temp$details$pca

# -------------------------------------------------------------------------
# 4) For Norway PCA subsetting, convert boundary to sp if needed:
# -------------------------------------------------------------------------
norway_sp_10m <- as(norway_spat_10m, "Spatial")

# -------------------------------------------------------------------------
# 5) Loop Over Each Species from master_data
# -------------------------------------------------------------------------
for(i in seq_len(nrow(master_data))) {
  speciesKey <- as.character(master_data$speciesKey[i])
  species_name <- master_data$species[i]
  message("\n=== Processing species: ", species_name, " (SpeciesKey: ", speciesKey, ") ===")
  
  # Create species folder if it does not exist
  species_folder <- file.path(base_dir, speciesKey)
  if(!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)
  
  # --- (A) Retrieve Occurrence Data from GBIF
  gbif_out <- tryCatch(
    rgbif::occ_search(
      scientificName = species_name,
      limit         = 10000,
      fields        = c("decimalLatitude","decimalLongitude")
    ),
    error = function(e) NULL
  )
  if(is.null(gbif_out) || is.null(gbif_out$data) || nrow(gbif_out$data) == 0) {
    message("No occurrence data found for species: ", species_name, " => skipping.")
    next
  }
  occ_data <- gbif_out$data
  occ_data <- occ_data[!is.na(occ_data$decimalLongitude) & !is.na(occ_data$decimalLatitude), ]
  if(nrow(occ_data) == 0){
    message("No valid coordinates for species: ", species_name, " => skipping.")
    next
  }
  
  # (A1) Plot fixed Norway envelope with occurrence points
  plot_norway_env_fixed <- base_norway_envelope_map +
    geom_point(data=occ_data, aes(x=decimalLongitude, y=decimalLatitude),
               color="red", size=1.2, alpha=0.8) +
    labs(title=paste("Norway Envelope (Fixed) +", species_name, "Occurrences"))
  
  # -----------------------------------------------------------------------
  # (B) Species Envelope in Norway (using 10m layers)
  # -----------------------------------------------------------------------
  # Extract species climate values (bio1, bio12, tmin, tmax) at occurrence locations
  sp_4clim <- terra::extract(
    c(
      bioclim_10m[[c("wc2.1_10m_bio_1","wc2.1_10m_bio_12")]],
      tmin_10m_all[["wc2.1_10m_tmin_01"]],
      tmax_10m_all[["wc2.1_10m_tmax_07"]]
    ),
    occ_data[, c("decimalLongitude", "decimalLatitude")]
  )
  sp_4clim <- na.omit(sp_4clim)
  if(nrow(sp_4clim) == 0){
    message("No 10m climate extraction for species: ", species_name, " => skipping.")
    next
  }
  b1_min_sp  <- min(sp_4clim$wc2.1_10m_bio_1, na.rm=TRUE)
  b1_max_sp  <- max(sp_4clim$wc2.1_10m_bio_1, na.rm=TRUE)
  b12_min_sp <- min(sp_4clim$wc2.1_10m_bio_12, na.rm=TRUE)
  b12_max_sp <- max(sp_4clim$wc2.1_10m_bio_12, na.rm=TRUE)
  tmin_min_sp<- min(sp_4clim$wc2.1_10m_tmin_01, na.rm=TRUE)
  tmin_max_sp<- max(sp_4clim$wc2.1_10m_tmin_01, na.rm=TRUE)
  tmax_min_sp<- min(sp_4clim$wc2.1_10m_tmax_07, na.rm=TRUE)
  tmax_max_sp<- max(sp_4clim$wc2.1_10m_tmax_07, na.rm=TRUE)
  
  buffer2 <- 1
  b1lo   <- b1_min_sp - buffer2; b1hi   <- b1_max_sp + buffer2
  b12lo  <- b12_min_sp - buffer2; b12hi  <- b12_max_sp + buffer2
  tminlo <- tmin_min_sp - buffer2; tminhi <- tmin_max_sp + buffer2
  tmaxlo <- tmax_min_sp - buffer2; tmaxhi <- tmax_max_sp + buffer2
  
  mask_species <- (global_bio1  >= b1lo  & global_bio1  <= b1hi) &
    (global_bio12 >= b12lo & global_bio12 <= b12hi) &
    (global_tmin  >= tminlo & global_tmin  <= tminhi) &
    (global_tmax  >= tmaxlo & global_tmax  <= tmaxhi)
  
  combined_species_mask <- terra::ifel(mask_species, 1, 0)
  species_mask_norway <- mask(crop(combined_species_mask, norway_spat_10m), norway_spat_10m)
  
  species_df <- as.data.frame(species_mask_norway, xy=TRUE, na.rm=FALSE)
  names(species_df)[3] <- "suitable"
  
  plot_species_env_nor <- ggplot() +
    geom_raster(data=species_df, aes(x=x, y=y, fill=suitable)) +
    scale_fill_viridis_c(na.value="white", name="Suitability") +
    geom_sf(data=norway_spat_10m, fill=NA, color="black", linewidth=0.3) +
    labs(
      title=paste("Species Envelope in Norway:", species_name),
      x="Longitude", y="Latitude"
    ) +
    theme_minimal()
  
  # -----------------------------------------------------------------------
  # (C) Global PCA Map
  # -----------------------------------------------------------------------
  # Extract environmental data for occurrence points from the global raster stack
  occ_env <- raster::extract(env_stack_global, occ_data[, c("decimalLongitude", "decimalLatitude")])
  occ_env <- occ_env[complete.cases(occ_env), ]
  if(nrow(occ_env) < 2) {
    message("Not enough valid environmental data for PCA for species: ", species_name, " => skipping.")
    next
  }
  
  species_scores <- ade4::suprow(my_pca_global, occ_env)$lisup
  if(nrow(species_scores) < 2) {
    message("Not enough PCA scores for species: ", species_name, " => skipping.")
    next
  }
  sc_12 <- species_scores[, 1:2]
  my_mean <- apply(sc_12, 2, mean)
  my_sd   <- apply(sc_12, 2, sd)
  
  real_sp <- generateSpFromPCA(
    raster.stack = env_stack_global,
    pca          = my_pca_global,
    means        = my_mean,
    sds          = my_sd,
    rescale      = TRUE,
    plot         = FALSE
  )
  global_suit_df <- as.data.frame(real_sp$suitab.raster, xy=TRUE, na.rm=FALSE)
  names(global_suit_df)[3] <- "suit_global"
  
  plot_global_pca <- ggplot() +
    geom_raster(data=global_suit_df, aes(x=x, y=y, fill=suit_global)) +
    scale_fill_viridis_c(name="Suitability", na.value="white") +
    geom_point(data=occ_data, aes(x=decimalLongitude, y=decimalLatitude),
               color="red", size=1, alpha=0.5) +
    coord_quickmap() +
    labs(
      title=paste("Global PCA Suitability:", species_name),
      x="Longitude", y="Latitude"
    ) +
    theme_minimal()
  
  # -----------------------------------------------------------------------
  # (D) Norway PCA Map
  # -----------------------------------------------------------------------
  nor_pca_rast <- raster::mask(raster::crop(real_sp$suitab.raster, norway_spat_10m),
                               norway_spat_10m)
  nor_pca_df <- as.data.frame(nor_pca_rast, xy=TRUE, na.rm=FALSE)
  names(nor_pca_df)[3] <- "suit_nor_pca"
  
  plot_nor_pca <- ggplot() +
    geom_raster(data=nor_pca_df, aes(x=x, y=y, fill=suit_nor_pca)) +
    scale_fill_viridis_c(na.value="white", name="Suitability") +
    labs(
      title=paste("Norway PCA-based:", species_name),
      x="Longitude", y="Latitude"
    ) +
    geom_sf(data=st_as_sf(norway_spat_10m)) +
    theme_minimal()
  
  # -----------------------------------------------------------------------
  # (E) Boxplot & (F) Scatter in Norway
  # -----------------------------------------------------------------------
  # Prepare Norway climate data
  nor_box <- as.data.frame(norway_bioclim, xy=FALSE, na.rm=TRUE)[, c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]
  nor_box$tmin <- as.data.frame(norway_tmin, xy=FALSE, na.rm=TRUE)[, "wc2.1_10m_tmin_01"]
  nor_box$tmax <- as.data.frame(norway_tmax, xy=FALSE, na.rm=TRUE)[, "wc2.1_10m_tmax_07"]
  nor_box$type <- "Norwegian Climate"
  colnames(nor_box)[1:2] <- c("bio1", "bio12")
  
  sp_box <- data.frame(
    bio1  = sp_4clim$wc2.1_10m_bio_1,
    bio12 = sp_4clim$wc2.1_10m_bio_12,
    tmin  = sp_4clim$wc2.1_10m_tmin_01,
    tmax  = sp_4clim$wc2.1_10m_tmax_07,
    type  = paste(species_name, "Occurrences")
  )
  comb_box <- rbind(nor_box, sp_box)
  comb_long <- pivot_longer(
    comb_box,
    cols = c("bio1", "bio12", "tmin", "tmax"),
    names_to = "variable", values_to = "value"
  ) %>% filter(value > 0)
  
  plot_box <- ggplot(comb_long, aes(x=type, y=value, fill=type)) +
    geom_boxplot() +
    facet_wrap(~variable, scales="free") +
    labs(title=paste("Norway vs.", species_name, ": 4 Vars"),
         x="Data Source", y="Value") +
    theme_minimal() +
    coord_flip() +
    theme(legend.position="none")
  
  # Scatter plot with convex hulls
  nor_scatt <- nor_box %>% mutate(data_source="Norwegian Climate")
  sp_scatt  <- sp_box %>% mutate(data_source=paste(species_name, "Occurrences"))
  comb_scatter <- rbind(nor_scatt, sp_scatt) %>% filter(complete.cases(bio1, bio12))
  hulls <- comb_scatter %>%
    group_by(data_source) %>%
    group_modify(~ {
      if(nrow(.x) < 3) .x else .x[chull(.x$bio1, .x$bio12), ]
    }) %>% ungroup()
  
  plot_scatter <- ggplot(comb_scatter, aes(x=bio1, y=bio12, color=data_source)) +
    geom_polygon(data=hulls, aes(fill=data_source, group=data_source),
                 alpha=0.2, color=NA) +
    geom_point(alpha=0.5) +
    scale_y_sqrt() +
    labs(title=paste(species_name, ": bio1 vs. bio12"),
         subtitle="Scatter + Convex Hull") +
    theme_minimal()
  
  # -----------------------------------------------------------------------
  # (G) Combine 6 Panels:
  #    A) Fixed Norway Envelope + Occurrence Points
  #    B) Species Envelope in Norway
  #    C) Global PCA Suitability
  #    D) Norway PCA-based Suitability
  #    E) Boxplot
  #    F) Scatter
  # -----------------------------------------------------------------------
  combined_plots <- plot_grid(
    plot_norway_env_fixed, # A
    plot_species_env_nor,  # B
   # plot_global_pca,       # C
  #  plot_nor_pca,          # D
    plot_box,              # E
    plot_scatter,          # F
    labels = c("A", "B", "C", "D"),
    ncol = 2
  )
  
  # Save final plot in the species folder using speciesKey in the filename
  outfile <- file.path(species_folder, paste0(speciesKey, "_all_plots.png"))
  cowplot::save_plot(filename = outfile,
                     plot = combined_plots,
                     ncol = 1, base_height = 14)
  message("Saved plot for species: ", species_name, " -> ", outfile)
}
##############################################################################
# END SCRIPT
##############################################################################
