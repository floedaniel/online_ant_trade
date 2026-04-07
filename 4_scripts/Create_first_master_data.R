# Load required libraries
library(europepmc)
library(rentrez)
library(rcrossref)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(rgbif)
library(rio)
library(stringr)
library(rgbif)
library(rio)
library(taxize)
library(tibble)
library(readxl)
library(tidyverse)

# mdir arter --------------------------------------------------------------

# Define the list of ant species from the raw Excel file
mdir_species <- rio::import("./1_raw_data/Artsliste maur 161025.xlsx") 

mdir_species  <- mdir_species %>%
  mutate(source = "Mdir")

# -------------------------------------------------------------------------

path <- "./1_raw_data/scraped"

# Get all xlsx files
files <- list.files(path, pattern = "\\.xlsx$", full.names = TRUE)

# Read and combine all files
scraped <- map_dfr(files, read_xlsx)

# -------------------------------------------------------------------------
gippet <- rio::import("./1_raw_data/pnas.2016337118.sd03.csv") %>% as_tibble()

gippet <- gippet %>% filter(Sold=="Yes")

gippet <- gippet %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  distinct()

gippet

gippet  <- gippet %>%
  mutate(source = "gippet")

data <- bind_rows(gippet, mdir_species, scraped)

data <- data %>% distinct(species, .keep_all = T)

sp <- unique(na.omit(data$species))

gna_df <- map_dfr(sp, function(x) {
  out <- tryCatch(gna_verifier(x), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) {
    tibble(query = x)
  } else {
    as_tibble(out) %>% mutate(query = x, .before = 1)
  }
})

# Join back to original data to restore the source column
df <- data %>%
  select(source, species) %>%
  left_join(gna_df, by = c("species" = "submittedName"))

# check 
unique(df$matchType)
df %>% filter(matchType== "PartialExact")
df %>% filter(matchType==  "Fuzzy" )
df %>% filter(matchType==  "NoMatch")
df %>% filter(is.na(matchType))

# filter 
df <- df %>% filter(!is.na(matchType))
df <- df %>% filter(!matchType=="NoMatch")
df <- df %>% filter(!taxonomicStatus=="N/A")

unique(df$taxonomicStatus)

synonyme <- df %>% dplyr::filter(taxonomicStatus=="Synonym")

synonyme_info <- rgbif::name_backbone_checklist(synonyme$species) %>% as_tibble()

splist <- df$matchedName %>% 
  rgbif::name_backbone_checklist() 

splist <- splist %>% 
  filter(!rank=="Genus") %>% 
  filter(status=="ACCEPTED") %>% 
  filter(matchType=="EXACT")

ant_data <- left_join(df, splist,  by = c("species" = "canonicalName"))

ant_data <- ant_data %>% distinct(species, .keep_all = T) %>% filter(family=="Formicidae")

ant_data %>%
  group_by(source) %>%
  summarise(
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  mutate(total_species = sum(n_species)) %>%
  arrange(desc(n_species))

# Extract unique species keys (assuming speciesKey is the identifier)
df_with_key <- ant_data %>% filter(!is.na(usageKey))  

df_no_key <- ant_data %>% filter(is.na(usageKey)) %>% distinct(species, .keep_all = T)

summary<-  df_with_key %>%
  group_by(family, genus) %>%
  summarise(
    n_species = n_distinct(species),
    sources = paste(unique(source), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

library(writexl)

write_xlsx(df_with_key, "./2_processed_data/df_with_key.xlsx")
write_xlsx(df_no_key,   "./2_processed_data/df_no_key.xlsx")
write_xlsx(summary,     "./2_processed_data/summary.xlsx")


# -------------------------------------------------------------------------

