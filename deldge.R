
library(Seurat)
library(DElegate)
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Please provide at least two arguments: root directory and study name.")
}
options(future.globals.maxSize = 89128960000)

root = args[1] # ~/Real_Data or #~/Documents/Real_Data
study <- args[2]

if (length(args) >= 3) {
  organs <- list(args[3])
} else {
  # If no organ is given, use all .rds files in the directory
  organ_files <- list.files(path = file.path(root, study, 'results'), pattern = "^Integrated_.*\\.rds$", full.names = FALSE)
  organs <- gsub("Integrated_(.*)\\.rds", "\\1", organ_files)
}

for (organ in organs) {
  basename_path = paste0('Integrated_',organ,'.rds')
  datapath = file.path(root,study,'results',basename_path)
  
  output_dir = file.path(root,study,'results3','Integrated',organ)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
  
  combined.sct = readRDS(datapath)
  
 # combined.sct <- PrepSCTFindMarkers(combined.sct,assay = "SCT",verbose = T)
  DefaultAssay(combined.sct) <- 'RNA'
  if (study == 'Parabiosis1') {
    library(stringr)
#combined.sct@meta.data$condition <- sub("^[0-9]", "", combined.sct@meta.data$condition)  # Remove the first character
    replacements <- c("HA" = "HetO", "IA" = "IsoO")
    combined.sct@meta.data$condition <- str_replace_all(combined.sct@meta.data$condition, replacements)   
   }
  if (study == 'Parabiosis3') {
  combined.sct@meta.data$condition[combined.sct@meta.data$condition %in% c("Isol-O", "Isor-O", "Non-O")] <- "Isonon-O"  	
  }  
  combined.sct@meta.data$cellcon <-paste0(combined.sct@active.ident, combined.sct@meta.data$condition)
  
  Idents(combined.sct) <- 'cellcon'
  
  unique_pairs <- unique(combined.sct@meta.data$cellcon)
  unique_pairs <- unique_pairs[!grepl("Y", unique_pairs)]
  background_df <- data.frame(background = rownames(combined.sct@assays$RNA)) #for later preparting for WebGestaltR
  
  # Initialize an empty list to store the results
# Initialize an empty list to store the results
marker_results <- list()

for (i in 1:(length(unique_pairs) - 1)) {
  for (j in (i + 1):length(unique_pairs)) {
    pair1 <- unique_pairs[i]
    pair2 <- unique_pairs[j]
    
    # Split the pairs by uppercase letters
    pair1_parts <- unlist(regmatches(pair1, gregexpr("[A-Z]+|[0-9]+|[a-z]+", pair1)))
    pair2_parts <- unlist(regmatches(pair2, gregexpr("[A-Z]+|[0-9]+|[a-z]+", pair2)))
    
    # Check if the pairs have the same cell identity and different condition
    if (pair1_parts[1] == pair2_parts[1] && pair1_parts[2] == pair2_parts[2] && pair1_parts[3] != pair2_parts[3]) {

      
      
      if ((grepl("Tg-O", pair1) & !grepl("Tg-Tre", pair1)) | (grepl("O", pair1) & !grepl("Het-O", pair1))) {
        ident1 <- pair2
        ident2 <- pair1
      } else {
        ident1 <- pair1
        ident2 <- pair2
      }
    # Perform FindMarkers on the subset data
    # if (grepl("Tg-O", pair1) & !grepl("Tg-Tre", pair1)) {
    #   ident1 <- pair2
    #   ident2 <- pair1
    # } else {
    #   ident1 <- pair1
    #   ident2 <- pair2
    # }
    #   
		# Perform FindMarkers on the subset data
		# if (grepl("O", pair1) & !grepl("Het-O", pair1)) {
		#   ident1 <- pair2
		#   ident2 <- pair1
		# } else {
		#   ident1 <- pair1
		#   ident2 <- pair2
		# }
		cat("Running SCT DGE for", ident1, "vs", ident2, "\n")
		
		result <- tryCatch({
			print("got to FindMarkers running")
		  res <- findDE(combined.sct, compare=c(ident1, ident2),method="deseq")
		  res <- res[,c("pvalue","log_fc","rate1","rate2","padj","feature")]
		  colnames(res) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","Gene")
		  res <- res[which(abs(res$avg_log2FC) >= 0.25 & ((res$pct.1 >= 0.1) | (res$pct.2 >= 0.1))),]
		  #res <- FindMarkers(combined.sct, assay = "SCT", ident.1 = ident1, ident.2 = ident2, verbose = T)
			print(head(res))
			res
		  },
		  error = function(e) {
			cat("Error in FindMarkers for", pair1, "vs", pair2, ":", toString(e), "\n")
			return(NULL)
		  }
		)
		
		# Add the background column for the WebGestaltR script
		if (!is.null(result)) {
		   print("got to after ! is null") 
		   print(class(result))
		   #result$Gene <- rownames(result)

		   # Create an empty data frame with the correct dimensions
		   merged_result <- data.frame(matrix(ncol = ncol(result) + 1, nrow = nrow(background_df)))
		   colnames(merged_result) <- c(colnames(result), "background")

		   # Fill the new data frame with the values from 'result'
		   merged_result[1:nrow(result), 1:ncol(result)] <- result

		   # Fill the 'background' column with the values from 'background_df'
		   merged_result$background <- background_df$background
		   
			# Store the result in the list with a combined name
			marker_results[[paste0(ident1, "_vs_", ident2)]] <- merged_result
		}
    }
  }
}



	for (name in names(marker_results)) {
	  output_file = paste0(name,".csv")
	    dir.create(output_dir)
      output_path = file.path(output_dir, output_file)
      cat("Saving to:", output_path, "\n")
	  write.csv(marker_results[[name]], file = file.path(output_dir, output_file), row.names = FALSE)
	}
}	

