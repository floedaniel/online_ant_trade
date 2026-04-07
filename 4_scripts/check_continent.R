

path <- "C:/Users/dafl/Downloads/0054907-251009101135966/0054907-251009101135966.csv"

griis_data <- rio::import(path)

griis_ants <- griis_data %>% filter(speciesKey %in% ant_data$usageKey) %>% as_tibble()

# -------------------------------------------------------------------------
library(rgbif)
library(dplyr)
library(purrr)

# species list
ant_species <- griis_ants$acceptedScientificName %>% unique() %>% na.omit()

# function to test presence on each continent
check_continent <- function(species, continent) {
  occ <- tryCatch(
    occ_search(
      scientificName = species,
      continent = continent,
      limit = 1
    ),
    error = function(e) NULL
  )
  !is.null(occ$data) && nrow(occ$data) > 0
}

# vector of continents as recognized by GBIF
continents <- c("Africa", "Asia", "Europe", "North America", "South America", "Oceania")

# check each species for each continent
ant_continent_presence <- map_dfr(ant_species, function(sp) {
  presence <- map_lgl(continents, ~ check_continent(sp, .x))
  tibble(species = sp, !!!setNames(as.list(presence), continents))
})

# join back to main data if needed
griis_ants <- griis_ants %>%
  left_join(ant_continent_presence, by = c("acceptedScientificName" = "species"))

# inspect
griis_ants %>% filter(Europe=="FALSE")


write_xlsx(griis_ants, "./2_processed_data/invasive_GRIIS_ants.xlsx")
