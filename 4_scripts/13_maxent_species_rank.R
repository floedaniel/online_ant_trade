# Set encoding to UTF-8 to handle Norwegian characters
Sys.setlocale("LC_ALL", "C")
options(encoding = "UTF-8")

library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(furrr)

# Setup parallel
plan(multisession, workers = 10)

# Last Norge
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, crs = 4326)
norway <- world %>% filter(name == "Norway")

# Finn alle species-mapper
species_dirs <- list.dirs(
  "./species",
  recursive = FALSE,
  full.names = TRUE
)

# Parallel loop over species
results <- future_map_dfr(species_dirs, function(species_dir) {
  
  tryCatch({
    species_name <- basename(species_dir)
    
    # Finn current og future filer
    current_file <- list.files(
      species_dir,
      pattern = "current_clamped_.*\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )[1]  # Ta første match
    
    future_file <- list.files(
      species_dir,
      pattern = "future_clamped_.*\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )[1]
    
    # Hvis begge filer eksisterer
    if(!is.na(current_file) && !is.na(future_file)) {
      
      # Current
      r_current <- rast(current_file)
      r_current_norway <- crop(r_current, norway)
      r_current_norway <- mask(r_current_norway, norway)
      vals_current <- values(r_current_norway, mat = FALSE)
      vals_current <- vals_current[!is.na(vals_current)]
      
      # Future
      r_future <- rast(future_file)
      r_future_norway <- crop(r_future, norway)
      r_future_norway <- mask(r_future_norway, norway)
      vals_future <- values(r_future_norway, mat = FALSE)
      vals_future <- vals_future[!is.na(vals_future)]
      
      # Beregn metrics
      current_mean <- mean(vals_current)
      future_mean <- mean(vals_future)
      current_max <- max(vals_current)
      future_max <- max(vals_future)
      
      # Continuing threat score: minimum av current og future mean
      # (arten må være trussel i begge perioder)
      continuing_threat <- min(current_mean, future_mean)
      
      # Alternatively: gjennomsnitt av begge
      avg_threat <- (current_mean + future_mean) / 2
      
      # Change in threat
      threat_change <- future_mean - current_mean
      threat_change_pct <- ((future_mean - current_mean) / current_mean) * 100
      
      tibble(
        species = species_name,
        current_mean = current_mean,
        current_max = current_max,
        future_mean = future_mean,
        future_max = future_max,
        continuing_threat = continuing_threat,
        avg_threat = avg_threat,
        threat_change = threat_change,
        threat_change_pct = threat_change_pct,
        n_cells = length(vals_current)
      )
    } else {
      # Missing files
      tibble(
        species = species_name,
        current_mean = NA,
        current_max = NA,
        future_mean = NA,
        future_max = NA,
        continuing_threat = NA,
        avg_threat = NA,
        threat_change = NA,
        threat_change_pct = NA,
        n_cells = NA,
        error = "Missing current or future file"
      )
    }
  }, error = function(e) {
    tibble(
      species = basename(species_dir),
      current_mean = NA,
      current_max = NA,
      future_mean = NA,
      future_max = NA,
      continuing_threat = NA,
      avg_threat = NA,
      threat_change = NA,
      threat_change_pct = NA,
      n_cells = NA,
      error = as.character(e)
    )
  })
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Sorter etter continuing threat (høyest først)
results <- results %>% 
  arrange(desc(continuing_threat)) %>%
  mutate(
    threat_category = case_when(
      continuing_threat >= 0.3 ~ "High",
      continuing_threat >= 0.1 ~ "Moderate",
      continuing_threat >= 0.01 ~ "Low",
      TRUE ~ "Minimal"
    )
  )

results <- results %>% filter(!current_mean=="NA")

# Lagre
write_csv(results, "species_maxent_threat_assessment.csv")

# Vis topp resultater
cat("\n=== TOPP 20 CONTINUING THREATS ===\n")

print(results %>% 
        select(species, continuing_threat, current_mean, future_mean, threat_change, threat_category) %>%
        head(20), n = 20)

# Vis arter med størst økning
cat("\n=== TOPP 20 INCREASING THREATS ===\n")

results_increasing <- results %>% 
  filter(!is.na(threat_change)) %>%
  arrange(desc(threat_change))

print(results_increasing %>% 
        select(species, threat_change, threat_change_pct, current_mean, future_mean, threat_category) %>%
        head(20), n = 20)

# Visualiser
library(ggplot2)

ggplot(results %>% filter(!is.na(continuing_threat)), 
       aes(x = current_mean, y = future_mean, color = threat_category)) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("High" = "red", "Moderate" = "orange", "Low" = "blue", "Minimal" = "green")) +
  labs(title = "Current vs Future Climate Suitability in Norway",
       subtitle = "Points above line = increasing threat",
       x = "Current Mean Probability",
       y = "Future Mean Probability",
       color = "Threat Level") +
  theme_minimal()

# ggsave("threat_comparison.png", width = 10, height = 8, dpi = 300)


# ployly ------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(plotly)

p <- ggplot(
  results %>% filter(!is.na(continuing_threat)),
  aes(
    x = current_mean,
    y = future_mean,
    color = threat_category,
    text = species
  )
) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    color = "black"
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    color = "black"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = c(
    "High" = "red",
    "Moderate" = "orange",
    "Low" = "blue",
    "Minimal" = "green"
  )) +
  labs(
    title = "Current vs Future Climate Suitability in Norway",
    subtitle = "Points above line = increasing threat",
    x = "Current Mean Probability",
    y = "Future Mean Probability",
    color = "Threat Level"
  ) +
  theme_minimal()

ggplotly(p, tooltip = c("text", "x", "y", "color"))

# -------------------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(furrr)
library(ggplot2)
library(plotly)

# Setup parallel
plan(multisession, workers = 10)

# Last Norge
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, crs = 4326)
norway <- world %>% filter(name == "Norway")

# Finn alle species-mapper
species_dirs <- list.dirs(
  
  "./species",
  recursive = FALSE,
  full.names = TRUE
)

# Parallel loop over species
results <- future_map_dfr(species_dirs, function(species_dir) {
  
  tryCatch({
    species_name <- basename(species_dir)
    
    # Finn current og future filer
    current_file <- list.files(
      species_dir,
      pattern = "current_clamped_.*\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )[1]
    
    future_file <- list.files(
      species_dir,
      pattern = "future_clamped_.*\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )[1]
    
    # Hvis begge filer eksisterer
    if (!is.na(current_file) && !is.na(future_file)) {
      
      # Current
      r_current <- rast(current_file)
      r_current_norway <- crop(r_current, norway)
      r_current_norway <- mask(r_current_norway, norway)
      vals_current <- values(r_current_norway, mat = FALSE)
      vals_current <- vals_current[!is.na(vals_current)]
      
      # Future
      r_future <- rast(future_file)
      r_future_norway <- crop(r_future, norway)
      r_future_norway <- mask(r_future_norway, norway)
      vals_future <- values(r_future_norway, mat = FALSE)
      vals_future <- vals_future[!is.na(vals_future)]
      
      # === BASIC METRICS ===
      current_mean <- mean(vals_current)
      future_mean <- mean(vals_future)
      current_max <- max(vals_current)
      future_max <- max(vals_future)
      
      # Continuing threat score
      continuing_threat <- min(current_mean, future_mean)
      avg_threat <- (current_mean + future_mean) / 2
      
      # Change in threat
      threat_change <- future_mean - current_mean
      threat_change_pct <- ((future_mean - current_mean) / current_mean) * 100
      
      # === HOTSPOT METRICS ===
      
      # Cell area (approximate for Norway's latitude ~62°N)
      cell_area_km2 <- prod(res(r_current_norway)) * 111^2 * cos(62 * pi / 180)
      
      # --- Current period ---
      # Percentiles (robust peak metrics)
      current_p95 <- as.numeric(quantile(vals_current, 0.95, na.rm = TRUE))
      current_p99 <- as.numeric(quantile(vals_current, 0.99, na.rm = TRUE))
      
      # Cells above thresholds
      current_n_above_0.5 <- sum(vals_current >= 0.5)
      current_n_above_0.7 <- sum(vals_current >= 0.7)
      current_n_above_0.9 <- sum(vals_current >= 0.9)
      current_prop_above_0.7 <- current_n_above_0.7 / length(vals_current)
      current_area_above_0.7_km2 <- current_n_above_0.7 * cell_area_km2
      current_area_above_0.9_km2 <- current_n_above_0.9 * cell_area_km2
      
      # --- Future period ---
      future_p95 <- as.numeric(quantile(vals_future, 0.95, na.rm = TRUE))
      future_p99 <- as.numeric(quantile(vals_future, 0.99, na.rm = TRUE))
      
      future_n_above_0.5 <- sum(vals_future >= 0.5)
      future_n_above_0.7 <- sum(vals_future >= 0.7)
      future_n_above_0.9 <- sum(vals_future >= 0.9)
      future_prop_above_0.7 <- future_n_above_0.7 / length(vals_future)
      future_area_above_0.7_km2 <- future_n_above_0.7 * cell_area_km2
      future_area_above_0.9_km2 <- future_n_above_0.9 * cell_area_km2
      
      # --- Combined hotspot scores ---
      hotspot_score <- (current_p95 + future_p95) / 2
      hotspot_area_change_km2 <- future_area_above_0.7_km2 - current_area_above_0.7_km2
      p95_change <- future_p95 - current_p95
      
      # Composite risk score: combines mean threat with peak suitability
      # Weights can be adjusted
      composite_risk <- (0.5 * avg_threat) + (0.5 * hotspot_score)
      
      tibble(
        species = species_name,
        # Basic metrics
        current_mean = current_mean,
        current_max = current_max,
        future_mean = future_mean,
        future_max = future_max,
        continuing_threat = continuing_threat,
        avg_threat = avg_threat,
        threat_change = threat_change,
        threat_change_pct = threat_change_pct,
        n_cells = length(vals_current),
        # Percentile metrics
        current_p95 = current_p95,
        current_p99 = current_p99,
        future_p95 = future_p95,
        future_p99 = future_p99,
        p95_change = p95_change,
        # Threshold-area metrics (current)
        current_n_above_0.5 = current_n_above_0.5,
        current_n_above_0.7 = current_n_above_0.7,
        current_n_above_0.9 = current_n_above_0.9,
        current_prop_above_0.7 = current_prop_above_0.7,
        current_area_above_0.7_km2 = current_area_above_0.7_km2,
        current_area_above_0.9_km2 = current_area_above_0.9_km2,
        # Threshold-area metrics (future)
        future_n_above_0.5 = future_n_above_0.5,
        future_n_above_0.7 = future_n_above_0.7,
        future_n_above_0.9 = future_n_above_0.9,
        future_prop_above_0.7 = future_prop_above_0.7,
        future_area_above_0.7_km2 = future_area_above_0.7_km2,
        future_area_above_0.9_km2 = future_area_above_0.9_km2,
        # Combined scores
        hotspot_score = hotspot_score,
        hotspot_area_change_km2 = hotspot_area_change_km2,
        composite_risk = composite_risk
      )
    } else {
      tibble(
        species = species_name,
        current_mean = NA, current_max = NA, future_mean = NA, future_max = NA,
        continuing_threat = NA, avg_threat = NA, threat_change = NA, threat_change_pct = NA,
        n_cells = NA, current_p95 = NA, current_p99 = NA, future_p95 = NA, future_p99 = NA,
        p95_change = NA, current_n_above_0.5 = NA, current_n_above_0.7 = NA,
        current_n_above_0.9 = NA, current_prop_above_0.7 = NA,
        current_area_above_0.7_km2 = NA, current_area_above_0.9_km2 = NA,
        future_n_above_0.5 = NA, future_n_above_0.7 = NA, future_n_above_0.9 = NA,
        future_prop_above_0.7 = NA, future_area_above_0.7_km2 = NA,
        future_area_above_0.9_km2 = NA, hotspot_score = NA,
        hotspot_area_change_km2 = NA, composite_risk = NA,
        error = "Missing current or future file"
      )
    }
  }, error = function(e) {
    tibble(
      species = basename(species_dir),
      current_mean = NA, current_max = NA, future_mean = NA, future_max = NA,
      continuing_threat = NA, avg_threat = NA, threat_change = NA, threat_change_pct = NA,
      n_cells = NA, current_p95 = NA, current_p99 = NA, future_p95 = NA, future_p99 = NA,
      p95_change = NA, current_n_above_0.5 = NA, current_n_above_0.7 = NA,
      current_n_above_0.9 = NA, current_prop_above_0.7 = NA,
      current_area_above_0.7_km2 = NA, current_area_above_0.9_km2 = NA,
      future_n_above_0.5 = NA, future_n_above_0.7 = NA, future_n_above_0.9 = NA,
      future_prop_above_0.7 = NA, future_area_above_0.7_km2 = NA,
      future_area_above_0.9_km2 = NA, hotspot_score = NA,
      hotspot_area_change_km2 = NA, composite_risk = NA,
      error = as.character(e)
    )
  })
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Filter out failed species
results <- results %>% filter(!is.na(current_mean))

# Add threat categories
results <- results %>%
  mutate(
    threat_category = case_when(
      continuing_threat >= 0.3 ~ "High",
      continuing_threat >= 0.1 ~ "Moderate",
      continuing_threat >= 0.01 ~ "Low",
      TRUE ~ "Minimal"
    ),
    hotspot_category = case_when(
      hotspot_score >= 0.7 ~ "Extreme",
      hotspot_score >= 0.5 ~ "High",
      hotspot_score >= 0.3 ~ "Moderate",
      TRUE ~ "Low"
    )
  )


# =============================================================================
# SORTED TABLES
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 1: TOP 20 BY CONTINUING THREAT (mean-based) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    arrange(desc(continuing_threat)) %>%
    select(species, continuing_threat, current_mean, future_mean, 
           threat_change, threat_category) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 2: TOP 20 BY HOTSPOT SCORE (p95-based peak suitability) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    arrange(desc(hotspot_score)) %>%
    select(species, hotspot_score, current_p95, future_p95, 
           current_mean, hotspot_category) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 3: TOP 20 BY COMPOSITE RISK (mean + hotspot combined) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    arrange(desc(composite_risk)) %>%
    select(species, composite_risk, avg_threat, hotspot_score, 
           threat_category, hotspot_category) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 4: TOP 20 BY AREA OF HIGH SUITABILITY (>=0.7, current) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    arrange(desc(current_area_above_0.7_km2)) %>%
    select(species, current_area_above_0.7_km2, future_area_above_0.7_km2,
           hotspot_area_change_km2, current_prop_above_0.7) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 5: TOP 20 BY EXTREME SUITABILITY (>=0.9, any cells) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    filter(current_n_above_0.9 > 0 | future_n_above_0.9 > 0) %>%
    arrange(desc(pmax(current_n_above_0.9, future_n_above_0.9))) %>%
    select(species, current_n_above_0.9, future_n_above_0.9,
           current_area_above_0.9_km2, future_area_above_0.9_km2, current_p99) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 6: TOP 20 INCREASING THREATS (future > current) ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    filter(!is.na(threat_change)) %>%
    arrange(desc(threat_change)) %>%
    select(species, threat_change, threat_change_pct, p95_change,
           hotspot_area_change_km2, current_mean, future_mean) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 7: HIDDEN THREATS (low mean, high hotspot score) ===\n")
cat("Species with localized but intense suitability\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    filter(current_mean < 0.1 & hotspot_score >= 0.3) %>%
    arrange(desc(hotspot_score)) %>%
    select(species, current_mean, hotspot_score, current_p95, 
           current_area_above_0.7_km2, hotspot_category) %>%
    head(20),
  n = 20
)

cat("\n", strrep("=", 70), "\n")
cat("=== TABLE 8: SUMMARY BY HOTSPOT CATEGORY ===\n")
cat(strrep("=", 70), "\n\n")

print(
  results %>%
    group_by(hotspot_category) %>%
    summarise(
      n_species = n(),
      mean_hotspot_score = mean(hotspot_score, na.rm = TRUE),
      mean_area_above_0.7 = mean(current_area_above_0.7_km2, na.rm = TRUE),
      total_species_with_extreme = sum(current_n_above_0.9 > 0, na.rm = TRUE)
    ) %>%
    arrange(desc(mean_hotspot_score))
)


# =============================================================================
# PLOTS
# =============================================================================

# --- PLOT 1: Current vs Future Mean (original) ---
p1 <- ggplot(
  results,
  aes(x = current_mean, y = future_mean, color = threat_category, text = species)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_color_manual(values = c(
    "High" = "#d62728", "Moderate" = "#ff7f0e", 
    "Low" = "#1f77b4", "Minimal" = "#2ca02c"
  )) +
  labs(
    title = "Current vs Future Climate Suitability (Mean)",
    subtitle = "Points above diagonal = increasing threat",
    x = "Current Mean Probability",
    y = "Future Mean Probability",
    color = "Threat Level"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
# ggsave("plot1_mean_comparison.png", p1, width = 10, height = 8, dpi = 300)


# --- PLOT 2: Current vs Future P95 (Hotspot comparison) ---
p2 <- ggplot(
  results,
  aes(x = current_p95, y = future_p95, color = hotspot_category, text = species)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_color_manual(values = c(
    "Extreme" = "#d62728", "High" = "#ff7f0e", 
    "Moderate" = "#1f77b4", "Low" = "#2ca02c"
  )) +
  labs(
    title = "Current vs Future Peak Suitability (95th Percentile)",
    subtitle = "Captures localized hotspots missed by mean values",
    x = "Current P95",
    y = "Future P95",
    color = "Hotspot Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)

# ggsave("plot2_p95_comparison.png", p2, width = 10, height = 8, dpi = 300)


# --- PLOT 3: Mean vs Hotspot Score (identify hidden threats) ---
p3 <- ggplot(
  results,
  aes(x = current_mean, y = hotspot_score, color = hotspot_category, text = species)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  annotate("rect", xmin = 0, xmax = 0.15, ymin = 0.3, ymax = 1, 
           alpha = 0.1, fill = "red") +
  annotate("text", x = 0.075, y = 0.95, label = "Hidden\nThreats", 
           color = "red", size = 3, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_color_manual(values = c(
    "Extreme" = "#d62728", "High" = "#ff7f0e", 
    "Moderate" = "#1f77b4", "Low" = "#2ca02c"
  )) +
  labs(
    title = "Mean Suitability vs Hotspot Score",
    subtitle = "Points above diagonal have localized high-suitability areas",
    x = "Current Mean Probability",
    y = "Hotspot Score (avg P95)",
    color = "Hotspot Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#print(p3)
# ggsave("plot3_mean_vs_hotspot.png", p3, width = 10, height = 8, dpi = 300)


# --- PLOT 4: Area of High Suitability (current vs future) ---
p4 <- ggplot(
  results %>% filter(current_area_above_0.7_km2 > 0 & future_area_above_0.7_km2 > 0),
  aes(x = current_area_above_0.7_km2, y = future_area_above_0.7_km2, 
      color = hotspot_category, text = species)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
 # scale_x_log10(labels = scales::comma) +
#  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = c(
    "Extreme" = "#d62728", "High" = "#ff7f0e", 
    "Moderate" = "#1f77b4", "Low" = "#2ca02c"
  )) +
  labs(
    title = "Area of High Suitability (>=0.7) in Norway",
    subtitle = "Log scale; points above diagonal = expanding suitable area",
    x = expression("Current Area (km"^2*")"),
    y = expression("Future Area (km"^2*")"),
    color = "Hotspot Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p4
# ggsave("plot4_area_comparison.png", p4, width = 10, height = 8, dpi = 300)


# --- PLOT 5: Composite Risk Distribution ---
p5 <- ggplot(results, aes(x = composite_risk, fill = hotspot_category)) +
  geom_histogram(bins = 30, alpha = 0.8, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = median(results$composite_risk, na.rm = TRUE), 
             linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(
    "Extreme" = "#d62728", "High" = "#ff7f0e", 
    "Moderate" = "#1f77b4", "Low" = "#2ca02c"
  )) +
  labs(
    title = "Distribution of Composite Risk Scores",
    subtitle = "Dashed line = median",
    x = "Composite Risk Score",
    y = "Number of Species",
    fill = "Hotspot Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p5)
# ggsave("plot5_risk_distribution.png", p5, width = 10, height = 6, dpi = 300)


# --- PLOT 6: Heatmap of Top Species by Multiple Metrics ---
top_species <- results %>%
  arrange(desc(composite_risk)) %>%
  head(30) %>%
  select(species, current_mean, future_mean, hotspot_score, 
         current_p95, future_p95, composite_risk) %>%
  pivot_longer(-species, names_to = "metric", values_to = "value")

p6 <- ggplot(top_species, aes(x = metric, y = reorder(species, value), fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", limits = c(0, 1)) +
  scale_x_discrete(labels = c(
    "current_mean" = "Curr Mean", "future_mean" = "Fut Mean",
    "hotspot_score" = "Hotspot", "current_p95" = "Curr P95",
    "future_p95" = "Fut P95", "composite_risk" = "Composite"
  )) +
  labs(
    title = "Top 30 Species by Composite Risk",
    subtitle = "Multiple metrics comparison",
    x = "", y = "", fill = "Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

print(p6)
# ggsave("plot6_heatmap_top30.png", p6, width = 10, height = 12, dpi = 300)

# -------------------------------------------------------------------------

results %>% 
  select(species, current_mean, future_mean, hotspot_score, 
         current_p95, future_p95, composite_risk) %>%
  pivot_longer(-species, names_to = "metric", values_to = "value") %>% 
  filter(metric=="composite_risk") %>% 
  arrange(desc(value)) %>% 
  write_csv(., "./5_outputs/tables/top_species_maxent_threat_composite_risk.csv")

# =============================================================================
# INTERACTIVE PLOTLY VERSIONS
# =============================================================================

# Interactive Plot 1: Mean comparison
ggplotly(p1, tooltip = c("text", "x", "y"))

# Interactive Plot 2: P95 comparison
ggplotly(p2, tooltip = c("text", "x", "y"))

# Interactive Plot 3: Mean vs Hotspot (best for finding hidden threats)
# ggplotly(p3, tooltip = c("text", "x", "y"))

# Interactive Plot 4: Area comparison
 ggplotly(p4, tooltip = c("text", "x", "y"))


# =============================================================================
# FUNCTION: Plot distribution for a single species
# =============================================================================

# Simple boxplot comparison for selected species

plot_species_boxplot <- function(species_names, species_dirs, norway) {
  
  all_data <- map_dfr(species_names, function(sp) {
    
    species_dir <- species_dirs[basename(species_dirs) == sp]
    if (length(species_dir) == 0) return(NULL)
    
    current_file <- list.files(species_dir, pattern = "current_clamped_.*\\.tif$", 
                               recursive = TRUE, full.names = TRUE)[1]
    future_file <- list.files(species_dir, pattern = "future_clamped_.*\\.tif$", 
                              recursive = TRUE, full.names = TRUE)[1]
    
    r_current <- rast(current_file)
    r_current_norway <- mask(crop(r_current, norway), norway)
    vals_current <- values(r_current_norway, mat = FALSE)
    vals_current <- vals_current[!is.na(vals_current)]
    
    r_future <- rast(future_file)
    r_future_norway <- mask(crop(r_future, norway), norway)
    vals_future <- values(r_future_norway, mat = FALSE)
    vals_future <- vals_future[!is.na(vals_future)]
    
    bind_rows(
      tibble(value = vals_current, period = "Current", species = sp),
      tibble(value = vals_future, period = "Future", species = sp)
    )
  })
  
  ggplot(all_data, aes(x = species, y = value, fill = period)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Current" = "#1f77b4", "Future" = "#d62728")) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_flip() +
    labs(
      title = "Climate Suitability Distribution",
      subtitle = "Box = IQR, whiskers = 1.5*IQR, top whisker approximates P95",
      x = NULL,
      y = "Suitability",
      fill = "Period"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Pick top 6 species by hotspot score
example_species <- results %>%
  arrange(desc(hotspot_score)) %>%
  slice(1:6) %>%
  pull(species)

p_box <- plot_species_boxplot(example_species, species_dirs, norway)

print(p_box)

# ggsave("boxplot_comparison.png", p_box, width = 10, height = 6, dpi = 300)


# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save full results
write_csv(results, "./5_outputs/tables/maxent_threat_category.csv")

# Save summary for high-risk species
# results %>%
  filter(hotspot_category %in% c("Extreme", "High") | threat_category == "High") %>%
  arrange(desc(composite_risk)) %>%
  write_csv("species_norway_high_risk_summary.csv")

cat("\n\nResults saved to:\n")
cat("- species_norway_threat_assessment_full.csv\n")
cat("- species_norway_high_risk_summary.csv\n")
