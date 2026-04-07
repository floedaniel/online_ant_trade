
# https://www.pnas.org/doi/10.1073/pnas.2016337118#supplementary-materials

df <- rio::import("C:/Users/dafl/Downloads/pnas.2016337118.sd03.csv") %>% as_tibble()

df <- df %>% filter(Sold=="Yes")

df

# -------------------------------------------------------------------------

df <- df %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  distinct()

# -------------------------------------------------------------------------

taxonomy_info <- rgbif::name_backbone_checklist(df$species) %>% as_tibble()

# -------------------------------------------------------------------------

df <- left_join(df, taxonomy_info)

 rio::export(df, "ants_sold_Jerome_etal_2021.xlsx")

df <-  rio::import("ants_sold_Jerome_etal_2021.xlsx")
# -------------------------------------------------------------------------
library(rgbif)
library(furrr)
library(dplyr)

# Function to check presence in a region (based on your working function)
check_presence <- function(species_name, region_type, region_code) {
  occ <- tryCatch({
    if(region_type == "country") {
      rgbif::occ_search(
        scientificName = species_name,
        country = region_code,
        limit = 1
      )
    } else {
      rgbif::occ_search(
        scientificName = species_name,
        continent = region_code,
        limit = 1
      )
    }
  }, error = function(e) NULL)
  
  if (is.null(occ) || is.null(occ$data) || nrow(occ$data) == 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

library(parallel)

# Total number of cores on your machine
total_cores <- detectCores()
print(paste("Total cores:", total_cores))

# Recommended: Use total cores minus 1-2 to keep system responsive
recommended_workers <- max(1, total_cores - 10)
print(paste("Recommended workers:", recommended_workers))

# Set up parallel processing
plan(multisession, workers = 10)

# Check each region one at a time with progress
cat("Checking Norway...\n")
df$Norway <- future_map_lgl(df$scientificName, ~check_presence(.x, "country", "NO"), .progress = TRUE)

cat("Checking Europe...\n")
df$Europe <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "Europe"), .progress = TRUE)

cat("Checking Asia...\n")
df$Asia <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "Asia"), .progress = TRUE)

cat("Checking Africa...\n")
df$Africa <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "Africa"), .progress = TRUE)

cat("Checking North America...\n")
df$North_America <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "North_America"), .progress = TRUE)

cat("Checking South America...\n")
df$South_America <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "South_America"), .progress = TRUE)

cat("Checking Oceania...\n")
df$Oceania <- future_map_lgl(df$scientificName, ~check_presence(.x, "continent", "Oceania"), .progress = TRUE)

# Reset to sequential processing
plan(sequential)

cat("Done!\n")

df


# -------------------------------------------------------------------------

df %>% filter(Norway==FALSE)
