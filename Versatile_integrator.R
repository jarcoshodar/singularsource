library(Seurat)
library(here)

# Define the function to integrate Seurat objects
integrate_seurat_objects <- function(base_dir, choice) {
  rds_files <- list.files(path = file.path(base_dir, choice), pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  
  if (length(rds_files) == 0) {
    cat("No .rds files found in", file.path(base_dir, choice), "\n")
    return(NULL)
  }
  
  # Load Seurat objects

  seurat_objects <- lapply(rds_files, function(file) {
  seurat_obj <- readRDS(file)
  file_name <- basename(file)
  identifier <- tools::file_path_sans_ext(file_name)
  print(identifier)
  seurat_obj@meta.data$orig.ident <- identifier

  # Extract organ names
  organ <- gsub("^(.*?)-.*$", "\\1", identifier)

  # Extract sample numbers (if any)
  sample_number <- gsub("^.*-(\\d+)$", "\\1", identifier)
  sample_number[sample_number == identifier] <- NA  # Set to NA if no sample number found

  # Extract the remaining parts between "-"
  remaining_parts <- gsub("^.*?-(.*)$", "\\1", identifier)
  remaining_parts <- gsub("-?\\d+$", "", remaining_parts)
  remaining_parts_list <- strsplit(remaining_parts, "-")

  # Process remaining parts (sex and condition)
  sex <- ""
  condition <- ""
  
  for (part in remaining_parts_list[[1]]) {
    if (part %in% c("M", "F")) {
    # Check if there's both "M" and "F" in the identifier
    if (grepl("M", identifier) && grepl("F", identifier)) {
      condition <- paste(condition, part, sep = "-")
    } else {     
         sex <- part
    }
    } else {
      # Append part to the condition (in case there are multiple conditions)
      if (condition == "") {
        condition <- part
      } else {
        condition <- paste(condition, part, sep = "-")
      }
    }
  }
  
  seurat_obj@meta.data$organ <- organ
  seurat_obj@meta.data$sex <- sex
  seurat_obj@meta.data$condition <- condition
  seurat_obj@meta.data$sample_number <- sample_number

  return(seurat_obj)
})
 
  
  # Perform integration
  features <- SelectIntegrationFeatures(object.list = seurat_objects)
  seurat_list <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
  integrated_object <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features, reduction='rpca')
  integrated_object <- IntegrateData(anchorset = integrated_object, normalization.method = "SCT")
  
  # Save integrated object
  saveRDS(integrated_object, file = file.path(base_dir, paste0("Integrated_", choice, ".rds")))
  
  return(integrated_object)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Error: Please provide exactly two arguments: root directory and study name.")
}

root <- args[1]
study <- args[2]

# Update the results_dir based on the root and study arguments
results_dir <- file.path(root, study, "results")
directories <- list.dirs(path = results_dir, full.names = TRUE, recursive = TRUE)

# Remove .RData files from the list
directories <- directories[!grepl("\\.RData$", directories)]

# Get the immediate subdirectories of the results_dir
choices_directory <- directories

# Call the integrate_seurat_objects function in a loop with the choices_directory
for (choice in choices_directory) {
  subdir_name <- basename(choice)
  cat("Processing", subdir_name, "\n")
  integrated_object <- integrate_seurat_objects(results_dir, subdir_name)
}

