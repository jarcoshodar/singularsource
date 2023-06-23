library(WebGestaltR)

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
study <- args[2]

if (study == "CalorieRestriction") {
  organism <- "rnorvegicus"
} else {
  organism <- "mmusculus"
}

data_path <- file.path(root, study, "results", "Integrated")

enrichment_analysis <- function(interest_genes, background_genes, organism, database, csv_file) {
  enrichment_results <- list()
  
  # Extract the original filename (without extension) from the csv_file variable
  original_filename <- tools::file_path_sans_ext(basename(csv_file))
  
  # Create a label for the analysis being done
  if (database == "geneontology_Biological_Process") {
    analysis_label <- c("Positive_GO", "Negative_GO")
    # Read the original CSV file
    find_markers_output <- read.csv(csv_file, header = TRUE)
    
    # Subset upregulated and downregulated genes
    upregulated_genes <- find_markers_output[find_markers_output$avg_log2FC > 0, "Gene"]
    downregulated_genes <- find_markers_output[find_markers_output$avg_log2FC < 0, "Gene"]
    
    # Perform enrichment analysis for upregulated and downregulated genes separately
    interest_genes_list <- list(upregulated_genes, downregulated_genes)
  } else {
    analysis_label = c("KEGG")
    find_markers_output <- read.csv(csv_file, header = TRUE)
    interest_genes <- na.omit(find_markers_output['Gene'])
    interest_genes_list <- list(interest_genes)
  }
  
  for (i in 1:length(interest_genes_list)) {
    current_interest_genes <- interest_genes_list[[i]]
    # Remove the 'na.omit' attribute from the current_interest_genes object
    attr(current_interest_genes, "na.action") <- NULL
    
    # Remove the 'na.omit' attribute from the background_genes object
    attr(background_genes, "na.action") <- NULL
    
    tryCatch({
      enrichment_results[[paste0(original_filename, "_", analysis_label[i], "_Enrichment")]] <- WebGestaltR(
        organism = organism,
        enrichMethod = "ORA",
        enrichDatabase = database,
        interestGene = unname(unlist(current_interest_genes)),
        interestGeneType = "genesymbol",
        collapseMethod = "mean",
        referenceGene = unname(background_genes),
        referenceGeneType = "genesymbol",
        minNum = 10,
        maxNum = 500,
        sigMethod = "fdr",
        fdrMethod = "BH",
        fdrThr = 0.01,
        topThr = 10,
        isOutput = FALSE,
        hostName = "https://www.webgestalt.org/"
      )
      cat("Success:", paste0(original_filename, "_",  analysis_label[i], "_Enrichment"), "\n")
    }, error = function(e) {
      cat("Error in WebGestaltR call:", conditionMessage(e), "\n")
      enrichment_results[[paste0('FAILURE_', original_filename, '_', analysis_label[i], '_Enrichment')]] <- data.frame(
        Gene = character(),
        name = character(),
        log2FoldChange = character(),
        padj = character()
      )
    })
  }
  
  return(enrichment_results)
}



process_folder <- function(folder_path) {
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in csv_files) {
    find_markers_output <- read.csv(csv_file, header = TRUE)
    significant_genes <- find_markers_output[find_markers_output$p_val < 0.05,]
    background_genes <- find_markers_output$background # Get background_genes
    interest_genes <- na.omit(significant_genes$Gene) # Get interest_genes and remove NAs
    
    # Call the function for KEGG and GO_BP enrichment analysis
    enrichment_results_KEGG <- enrichment_analysis(interest_genes, background_genes, organism, "pathway_KEGG", csv_file)
    enrichment_results_GO_BP <- enrichment_analysis(interest_genes, background_genes, organism, "geneontology_Biological_Process", csv_file)
    
    # Combine the results from KEGG and GO_BP enrichment analyses
    enrichment_results <- c(enrichment_results_KEGG, enrichment_results_GO_BP)
    
    for (label in names(enrichment_results)) {
      result <- enrichment_results[[label]]
      if (nrow(result) == 0) {
        output_file <- paste0("FAILURE_", label, ".csv")
      } else {
        output_file <- paste0(label, ".csv")
      }
      
      # Create separate folders for Positive_GO, Negative_GO, and KEGG enrichments
      if (grepl("Positive_GO", label)) {
        folder <- file.path(folder_path, "pGO")
      } else if (grepl("Negative_GO", label)) {
        folder <- file.path(folder_path, "nGO")
      } else if (grepl("KEGG", label)) { # Updated condition for creating the KEGG folder
        folder <- file.path(folder_path, "KEGG")
      }
      
      if (!dir.exists(folder)) {
        dir.create(folder)
      }
      
      savepath = file.path(folder, basename(output_file))
      write.csv(result, savepath, row.names = FALSE)
    }
  }
}



if (length(args) == 3) {
  organ <- args[3]
  process_folder(file.path(data_path, organ))
} else {
  subfolders <- list.dirs(data_path, recursive = FALSE)
  
  for (subfolder in subfolders) {
    process_folder(subfolder)
  }
}
