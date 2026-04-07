library(rgbif)
library(dplyr)
library(purrr)


clean <- rio::import("./2_processed_data/clean_ants.xlsx")
issue <- rio::import("./2_processed_data/ant_issues.xlsx")

df <- bind_rows(issue, clean)

# Simple function to check if species has ANY introduced occurrences
check_introduced <- function(species) {
  occ <- tryCatch(
    occ_search(
      scientificName = species,
      establishmentMeans = "introduced",
      limit = 1
    ),
    error = function(e) NULL
  )
  !is.null(occ$data) && nrow(occ$data) > 0
}

# Check each species for introduced status
ant_introduced_status <- map_dfr(df$species, function(sp) {
  tibble(
    species = sp,
    has_introduced_records = check_introduced(sp)
  )
})


introduced <- ant_introduced_status %>%  filter(has_introduced_records=="TRUE")     

has_introduced_records <- left_join(df, ant_introduced_status)


write_xlsx(has_introduced_records, "./2_processed_data/has_introduced_records.xlsx")
