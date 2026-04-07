# =========================================================================
# ANT SPECIES NAME VALIDATION AND SYNONYM LOOKUP
# =========================================================================

library(rgbif)
library(dplyr)
library(readxl)
library(purrr)
library(writexl)

# -------------------------------------------------------------------------
# 1. LOAD DATA FROM MULTIPLE SOURCES
# -------------------------------------------------------------------------

# Mdir species list
mdir_species <- rio::import("./1_raw_data/Artsliste maur 161025.xlsx") %>%
  mutate(source = "Mdir")

# Scraped data
path <- "./1_raw_data/scraped"
files <- list.files(path, pattern = "\\.xlsx$", full.names = TRUE)
scraped <- map_dfr(files, read_xlsx)

# Gippet data
gippet <- rio::import("./1_raw_data/pnas.2016337118.sd03.csv") %>%
  as_tibble() %>%
  filter(Sold == "Yes") %>%
  mutate(
    species = str_replace_all(species, "_", " "),
    source = "gippet"
  ) %>%
  select(species, source) %>% 
  distinct(species, .keep_all = TRUE)

antcheck <-  rio::import("./2_processed_data/AntCheck_data.xlsx") %>%
  as_tibble() %>%
  mutate(source = "antcheck") %>%
  select(source, "species"=species_name) %>% 
  distinct(species, .keep_all = TRUE)

Wang <- rio::import("./1_raw_data/Supplementary_Material_2_revised.csv") %>%
  as_tibble() %>%

  mutate(
    source = "Wang"
  ) %>%
  select(source, "species"=`Taxon code`) %>% 
  distinct(species, .keep_all = TRUE)

# -------------------------------------------------------------------------

# Combine all sources
data <- bind_rows(gippet, mdir_species, scraped, antcheck, Wang) %>%
  distinct(species, .keep_all = TRUE)

data

data %>% count(source)

# -------------------------------------------------------------------------
# 2. VALIDATE SPECIES NAMES WITH GBIF
# -------------------------------------------------------------------------

full <- name_backbone_checklist(data$species) %>% as_tibble()

species_key <- na.omit(full$speciesKey)

# -------------------------------------------------------------------------
# 3. GET FULL TAXONOMIC INFO FOR EACH SPECIES
# -------------------------------------------------------------------------

info_list <- list()

for (i in 1:length(species_key)) {
  key <- species_key[i]
  print(paste("Getting info:", i, "of", length(species_key)))
  info_list[[i]] <- name_usage(key = key)$data
}

# Combine results with different columns
bind_rows_base <- function(df_list) {
  all_cols <- unique(unlist(lapply(df_list, names)))
  df_list <- lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) df[[col]] <- NA
    return(df[, all_cols])
  })
  do.call(rbind, df_list)
}

results <- bind_rows_base(info_list)

results

# -------------------------------------------------------------------------
# 4. FILTER TO ANTS AND CLEAN DATA
# -------------------------------------------------------------------------

ants <- results %>% filter(family == "Formicidae")
cat("Removed", nrow(results) - nrow(ants), "non-ants\n")

# Add missing columns if they don't exist
if (!"accepted" %in% names(ants)) ants$accepted <- NA
if (!"acceptedKey" %in% names(ants)) ants$acceptedKey <- NA

# Create clean output
output <- ants %>%
  mutate(
    final_name = ifelse(is.na(accepted), scientificName, accepted),
    final_key = ifelse(is.na(acceptedKey), speciesKey, acceptedKey),
    is_synonym = !is.na(acceptedKey),
    has_issues = issues != ""
  ) %>%
  transmute(
    original_key = speciesKey,
    acceptedKey = final_key,
    original_name = scientificName,
    final_name,
    genus,
    species,
    rank,
    status = taxonomicStatus,
    is_synonym,
    has_issues,
    issues
  ) %>%
  filter(rank != "GENUS") %>%
  distinct(acceptedKey, .keep_all = TRUE)

# Replace empty strings with NA
output[output == ""] <- NA

# Split by data quality
clean_ants <- output %>% filter(is.na(issues))
clean_ants
ant_issues <- output %>% filter(!is.na(issues))
ant_issues
name_changes <- output %>% filter(original_key != acceptedKey | original_name != final_name)

# -------------------------------------------------------------------------
# 5. GET ALL SYNONYMS FOR EACH SPECIES
# -------------------------------------------------------------------------

all_synonyms <- data.frame()

for (i in 1:nrow(output)) {
  
  key <- output$acceptedKey[i]
  accepted_name <- output$final_name[i]
  
  if (is.na(key)) next
  
  print(paste("Getting synonyms:", i, "of", nrow(output)))
  
  syn <- tryCatch({
    name_usage(key = key, data = "synonyms")$data
  }, error = function(e) NULL)
  
  if (!is.null(syn) && nrow(syn) > 0) {
    syn$acceptedKey_searched <- key
    syn$accepted_name <- accepted_name
    all_synonyms <- bind_rows(all_synonyms, syn)
  }
}

all_synonyms <- as_tibble(all_synonyms)
all_synonyms
# -------------------------------------------------------------------------
# 6. CREATE COMPLETE NAME LOOKUP TABLE
# -------------------------------------------------------------------------

# Accepted names from synonyms
accepted_from_syn <- all_synonyms %>%
  distinct(acceptedKey_searched, accepted_name) %>%
  transmute(
    acceptedKey = acceptedKey_searched,
    name = accepted_name,
    genus = sub(" .*", "", accepted_name),
    status = "ACCEPTED",
    rank = "SPECIES"
  )

# Synonym names
synonym_names <- all_synonyms %>%
  transmute(
    acceptedKey = acceptedKey_searched,
    name = scientificName,
    genus,
    status = taxonomicStatus,
    rank
  )

# Combine
all_names <- bind_rows(accepted_from_syn, synonym_names) %>% distinct()
all_names

# Add species with NO synonyms
missing_species <- output %>%
  filter(!acceptedKey %in% all_names$acceptedKey) %>%
  transmute(
    acceptedKey,
    name = final_name,
    genus,
    status = "ACCEPTED",
    rank = "SPECIES"
  )

missing_species

# Final complete table
all_names_complete <- bind_rows(all_names, missing_species) %>%
  distinct() %>%
  arrange(acceptedKey, desc(status == "ACCEPTED"))

# -------------------------------------------------------------------------
# 7. SUMMARY STATISTICS
# -------------------------------------------------------------------------

cat("\n=== SUMMARY ===\n")
cat("Total species validated:", n_distinct(output$acceptedKey), "\n")
cat("Species with clean data:", n_distinct(clean_ants$acceptedKey), "\n")
cat("Species with issues:", n_distinct(ant_issues$acceptedKey), "\n")
cat("Species with name changes:", n_distinct(name_changes$acceptedKey), "\n")
cat("Species with synonyms:", n_distinct(all_synonyms$acceptedKey_searched), "\n")
cat("Total unique names (accepted + synonyms):", nrow(all_names_complete), "\n")

# -------------------------------------------------------------------------
# 8. EXPORT RESULTS
# -------------------------------------------------------------------------

write_xlsx(clean_ants, "./2_processed_data/clean_ants.xlsx")
write_xlsx(ant_issues, "./2_processed_data/ant_issues.xlsx")
write_xlsx(name_changes, "./2_processed_data/gbif_name_changes.xlsx")
write_xlsx(all_names_complete, "./2_processed_data/complete_ant_data.xlsx")

# END 
