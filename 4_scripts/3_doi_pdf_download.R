# =============================================================================
# Parallel DOI PDF download from multisource.xlsx files (per species folder)
# =============================================================================
library(metagear)
library(tidyverse)
library(furrr)
library(rio)

# --- 1. Find all multisource.xlsx files recursively ---
base_dir <- "./Species"

multisource_files <- list.files(
  path = base_dir,
  pattern = "_multisource\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)

message("Found ", length(multisource_files), " multisource.xlsx files")

if (length(multisource_files) == 0) {
  stop("No multisource.xlsx files found in ", base_dir)
}

# --- 2. Process each species file separately ---
for (species_file in multisource_files) {
  
  message("\n==============================================")
  message("Processing: ", species_file)
  
  # Read the file
  species_data <- tryCatch({
    rio::import(species_file)
  }, error = function(e) {
    message("Error reading file: ", e$message)
    next
  })
  
  if (is.null(species_data) || nrow(species_data) == 0) {
    message("No data in file, skipping...")
    next
  }
  
  # --- 3. Create PDF folder in same directory as xlsx file ---
  species_dir <- dirname(species_file)
  # Extract folder name above the file's directory (assuming structure Species/key/file)
  species_key <- basename(dirname(species_dir)) 
  pdf_folder <- file.path(species_dir, paste0(species_key, "_PDFs"))
  
  if (!dir.exists(pdf_folder)) {
    dir.create(pdf_folder, recursive = TRUE)
  }
  
  pdf_folder <- normalizePath(pdf_folder, winslash = "\\")
  # Append trailing slash for metagear internal logic consistency on Windows
  pdf_folder <- paste0(pdf_folder, "\\") 
  
  message("PDF folder: ", pdf_folder)
  
  # --- 4. Clean and prepare DOI data ---
  clean_df <- species_data %>%
    mutate(across(where(is.character), ~na_if(trimws(.x), ""))) %>%
    filter(!is.na(doi) & doi != "" & nchar(doi) > 5) %>%
    # Clean 'doi:' prefix if present, but retain 'https://doi.org/'
    mutate(doi = gsub("^doi:", "", doi, ignore.case = TRUE)) %>%
    mutate(doi = trimws(doi)) %>%
    distinct(doi, .keep_all = TRUE) %>%
    # Create safe filename by sanitizing the cleaned DOI string
    mutate(safe_filename = gsub("[/:\\*\\?\"<>\\|]", "_", doi))
  
  message("Found ", nrow(clean_df), " unique DOIs for this species")
  
  # --- 5. Skip existing files ---
  existing_files <- list.files(pdf_folder, pattern = "\\.pdf$", full.names = FALSE)
  existing_base <- gsub("\\.pdf$", "", existing_files)
  
  remaining_df <- clean_df %>%
    filter(!safe_filename %in% existing_base)
  
  message("Downloading ", nrow(remaining_df), " PDFs (", 
          nrow(clean_df) - nrow(remaining_df), " already exist)")
  
  if (nrow(remaining_df) == 0) {
    message("All PDFs already downloaded for this species!")
    next
  }
  
  # --- 6. Conservative parallel setup ---
  plan(multisession, workers = 10)
  options(timeout = 180)
  
  # --- 7. Robust download function - Accepts DOI and pre-sanitized filename ---
  safe_download <- function(doi, filename, out_dir) {
    tryCatch({
      result <- PDF_download(
        DOI = doi,
        directory = out_dir,
        theFileName = filename, # Use the pre-sanitized filename
        validatePDF = FALSE,
        WindowsProxy = TRUE
      )
      
      # Explicit check if file exists (metagear returns a list, not always T/F)
      expected_file <- file.path(gsub("\\\\$", "", out_dir), paste0(filename, ".pdf"))
      if (file.exists(expected_file)) {
        message("✅ Success: ", doi)
        return(TRUE)
      } else {
        message("❌ No file created: ", doi)
        return(FALSE)
      }
      
    }, error = function(e) {
      error_msg <- gsub("[\r\n\t]+", " ", as.character(e$message))
      message("❌ Failed: ", doi, " - ", substr(error_msg, 1, 100))
      return(FALSE)
    })
  }
  
  # --- 8. Process with parallel processing using pmap/future_pmap ---
  if (nrow(remaining_df) <= 10) {
    message("Processing sequentially...")
    results <- pmap_lgl(
      list(remaining_df$doi, remaining_df$safe_filename, list(pdf_folder)),
      safe_download
    )
  } else {
    message("Processing in parallel with 10 workers...")
    results <- future_pmap_lgl(
      list(remaining_df$doi, remaining_df$safe_filename, list(pdf_folder)),
      safe_download,
      .options = furrr_options(seed = NULL),
      .progress = TRUE
    )
  }
  
  plan(sequential)
  
  # --- 9. Summary for this species ---
  successful <- sum(results, na.rm = TRUE)
  message("Downloaded: ", successful, " out of ", nrow(remaining_df))
  
  final_count <- length(list.files(pdf_folder, pattern = "\\.pdf$"))
  message("Total PDFs for this species: ", final_count)
  
  # --- 10. Save download log for this species ---
  download_log <- remaining_df %>%
    mutate(
      download_success = results,
      download_date = Sys.time()
    )
  
  log_file <- file.path(species_dir, paste0(species_key, "_download_log.xlsx"))
  rio::export(download_log, log_file)
  message("Log saved: ", log_file)
}

message("\n==============================================")
message("All species processed!")
message("==============================================")