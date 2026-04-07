# =============================================================================
# Delete corrupt/invalid PDF files AND duplicates within the same folder
# =============================================================================
library(tidyverse)
library(pdftools)  # For PDF validation
library(digest)    # For file hashing to detect duplicates

# --- Set directory ---
pdf_dir <- "C:/Users/dafl/OneDrive - Folkehelseinstituttet/VKM Data/27.02.2025_maur_forprosjekt_biologisk_mangfold/data/species"

# --- Find all PDF files ---
pdf_files <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
message("Found ", length(pdf_files), " PDF files")

if (length(pdf_files) == 0) {
  message("No PDF files found in directory")
} else {
  
  # --- Function to check if PDF is valid ---
  is_pdf_valid <- function(file_path) {
    tryCatch({
      info <- pdf_info(file_path)
      if (is.null(info) || info$pages == 0) {
        return(FALSE)
      }
      return(TRUE)
    }, error = function(e) {
      return(FALSE)
    })
  }
  
  # --- Function to compute file hash ---
  get_file_hash <- function(file_path) {
    tryCatch({
      digest(file = file_path, algo = "md5")
    }, error = function(e) {
      NA_character_
    })
  }
  
  # --- Check all PDFs ---
  message("Checking PDF validity and computing hashes (this may take a while)...")
  
  file_status <- tibble(
    file_path = pdf_files,
    file_name = basename(pdf_files),
    folder = dirname(pdf_files),
    size_bytes = file.size(pdf_files),
    size_kb = round(size_bytes / 1024, 2)
  ) %>%
    filter(!is.na(size_bytes)) %>%
    mutate(
      is_valid = map_lgl(file_path, ~{
        message("Checking: ", basename(.x))
        is_pdf_valid(.x)
      }),
      file_hash = map_chr(file_path, get_file_hash)
    )
  
  # --- Find corrupt files ---
  corrupt_files <- file_status %>%
    filter(!is_valid)
  
  valid_files <- file_status %>%
    filter(is_valid)
  
  message("\n--- Validation Results ---")
  message("Valid PDFs: ", nrow(valid_files))
  message("Corrupt PDFs: ", nrow(corrupt_files))
  
  # --- Find duplicates WITHIN the same folder ---
  duplicates_to_delete <- valid_files %>%
    filter(!is.na(file_hash)) %>%
    group_by(folder, file_hash) %>%
    mutate(
      dup_count = n(),
      dup_rank = row_number()  # Keep first occurrence, mark rest as duplicates
    ) %>%
    ungroup() %>%
    filter(dup_count > 1, dup_rank > 1)  # Files to delete (not the first copy)
  
  # Show which files will be kept
  duplicates_kept <- valid_files %>%
    filter(!is.na(file_hash)) %>%
    group_by(folder, file_hash) %>%
    mutate(dup_count = n(), dup_rank = row_number()) %>%
    ungroup() %>%
    filter(dup_count > 1, dup_rank == 1)
  
  message("\n--- Duplicate Detection (within same folder) ---")
  message("Duplicate sets found: ", n_distinct(duplicates_to_delete$file_hash))
  message("Duplicate files to delete: ", nrow(duplicates_to_delete))
  
  # =========================================================================
  # STEP 1: Handle corrupt files
  # =========================================================================
  if (nrow(corrupt_files) > 0) {
    message("\n========================================")
    message("CORRUPT FILES TO DELETE:")
    message("========================================")
    corrupt_files %>%
      select(file_name, folder, size_kb) %>%
      arrange(folder, size_kb) %>%
      print(n = Inf)
    
    cat("\nDelete these", nrow(corrupt_files), "corrupt files? (y/n): ")
    response <- readline()
    
    if (tolower(response) %in% c("y", "yes")) {
      deleted_count <- 0
      for (i in 1:nrow(corrupt_files)) {
        tryCatch({
          file.remove(corrupt_files$file_path[i])
          deleted_count <- deleted_count + 1
          message("Deleted: ", corrupt_files$file_name[i])
        }, error = function(e) {
          message("Failed to delete: ", corrupt_files$file_name[i], " - ", e$message)
        })
      }
      message("Corrupt files deleted: ", deleted_count)
    } else {
      message("Corrupt file deletion cancelled")
    }
  } else {
    message("\nNo corrupt files found!")
  }
  
  # =========================================================================
  # STEP 2: Handle duplicates within same folder
  # =========================================================================
  if (nrow(duplicates_to_delete) > 0) {
    message("\n========================================")
    message("DUPLICATE FILES (within same folder):")
    message("========================================")
    
    # Show duplicates grouped by hash
    for (h in unique(duplicates_to_delete$file_hash)) {
      kept <- duplicates_kept %>% filter(file_hash == h)
      dupes <- duplicates_to_delete %>% filter(file_hash == h)
      
      message("\n--- Duplicate set in: ", kept$folder[1], " ---")
      message("  KEEPING: ", kept$file_name[1], " (", kept$size_kb[1], " KB)")
      for (j in 1:nrow(dupes)) {
        message("  DELETE:  ", dupes$file_name[j], " (", dupes$size_kb[j], " KB)")
      }
    }
    
    cat("\nDelete these", nrow(duplicates_to_delete), "duplicate files? (y/n): ")
    response <- readline()
    
    if (tolower(response) %in% c("y", "yes")) {
      deleted_count <- 0
      for (i in 1:nrow(duplicates_to_delete)) {
        tryCatch({
          file.remove(duplicates_to_delete$file_path[i])
          deleted_count <- deleted_count + 1
          message("Deleted duplicate: ", duplicates_to_delete$file_name[i])
        }, error = function(e) {
          message("Failed to delete: ", duplicates_to_delete$file_name[i], " - ", e$message)
        })
      }
      message("Duplicate files deleted: ", deleted_count)
    } else {
      message("Duplicate deletion cancelled")
    }
  } else {
    message("\nNo duplicates found within any folder!")
  }
  
  # =========================================================================
  # Final summary
  # =========================================================================
  remaining_valid <- nrow(valid_files) - nrow(duplicates_to_delete)
  
  message("\n========================================")
  message("FINAL SUMMARY")
  message("========================================")
  message("Original PDF count: ", length(pdf_files))
  message("Corrupt files: ", nrow(corrupt_files))
  message("Duplicate files: ", nrow(duplicates_to_delete))
  message("Unique valid files: ", remaining_valid)
  
  if (remaining_valid > 0) {
    # Recalculate for files that would remain
    final_valid <- valid_files %>%
      anti_join(duplicates_to_delete, by = "file_path")
    
    message("\nValid files size distribution:")
    message("Min size: ", round(min(final_valid$size_kb, na.rm = TRUE), 2), " KB")
    message("Max size: ", round(max(final_valid$size_kb, na.rm = TRUE), 2), " KB")
    message("Median size: ", round(median(final_valid$size_kb, na.rm = TRUE), 2), " KB")
  }
}