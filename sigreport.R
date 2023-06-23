input_path = '~/Documents/realsig'

recursive_file_search <- function(path, pattern) {
  r_files <- c()
  
  # List all files in the current directory
  files <- list.files(path, full.names = TRUE)
  
  for (file_path in files) {
    # Get file info
    info <- file.info(file_path)
    
    # If this file is a directory, search it recursively
    if (!is.na(info$isdir) && info$isdir) {
      r_files <- c(r_files, recursive_file_search(file_path, pattern))
      
      # Otherwise, if it matches the pattern, add it to the list
    } else if (!is.na(info$isdir) && !info$isdir && grepl(pattern, file_path)) {
      r_files <- c(r_files, file_path)
    }
  }
  
  return(r_files)
}


rds_files <- recursive_file_search(input_path, "hotspots\\.rds")
for (hotspots in rds_files) {
  print(paste("Processing file:", hotspots))
  
  merge_scores <- function(sigresult_path) {
    print(paste("Merging scores for file:", sigresult_path))
    celltype_list <- readRDS(sigresult_path)
    
    # New list to store the successfully processed results
    processed_celltype_list <- list()
    
    # Existing lists for storing the merged results and changed genes
    merged_list <- list()
    changed_genes_list <- list()
    
    for(i in 1:length(celltype_list)) {
      condition_list <- celltype_list[[i]]
      
      # Avoid NA ruining loop
      if(is.null(condition_list) || length(condition_list) < 2) {
        warning(paste("Skipping cell type due to insufficient conditions:", names(celltype_list)[i]))
        next
      }
      
      condition1 <- condition_list[[1]]
      condition2 <- condition_list[[2]]
      
      # Do merge if there are no problems
      if(is.list(condition1) && !is.null(condition1$final_score) && is.list(condition2) && !is.null(condition2$final_score)) {
        merged_final_score <- merge(condition1$final_score, 
                                    condition2$final_score, 
                                    by = "Gene", 
                                    all = TRUE)
        merged_final_score = merged_final_score %>%  mutate_all(replace_na, 0)
        
        # Filter genes that changed state: right now, activated is above 0.7 and inhibited is below 0.3
        changed_genes <- merged_final_score[(merged_final_score$Activation_probability.x > 0.7 & merged_final_score$Activation_probability.y < 0.3) | 
                                              (merged_final_score$Activation_probability.x < 0.3 & merged_final_score$Activation_probability.y > 0.7), "Gene"]
        
        # Save the merged df and path in the merged_list
        merged_list[[names(celltype_list)[i]]] <- list(path = sigresult_path, merged_score = merged_final_score)
        
        # Save the genes that changed state in the changed_genes_list
        changed_genes_list[[names(celltype_list)[i]]] <- changed_genes
        
        # Save the processed cell type in the processed_celltype_list
        processed_celltype_list[[names(celltype_list)[i]]] <- condition_list
      } else {
        warning(paste("Skipping cell type due to missing final_score:", names(celltype_list)[i]))
      }
    }
    
    return(list(merged_data = merged_list, changed_genes = changed_genes_list, celltype_list = processed_celltype_list))
  }
  
  output <- merge_scores(hotspots)
  
  changed_genes_list = output$changed_genes
  celltype_list = output$celltype_list
  
  new_list <- list()
  
  names_to_iter <- names(celltype_list)
  
  for(name in names_to_iter) {
    first_condition <- celltype_list[[name]][[1]]
    new_list[[name]] <- first_condition
  }
  
  new_list <- lapply(celltype_list, function(x) x[[1]])
  names(new_list) <- names(celltype_list)
  
  # Loop over every name in new_list
  for(name in names(new_list)) {
    
    # Loop over both vis_net_A and vis_net_I
    for(net in c("vis_net_A", "vis_net_I")) {
      
      # Check if the name and network exist in new_list
      if(name %in% names(new_list) && net %in% names(new_list[[name]])) {
        
        # Loop over the elements in the vis_net_A or vis_net_I list
        for(i in seq_along(new_list[[name]][[net]])) {
          
          # Find the row where group is 'int' and get the corresponding label
          matching_rows <- new_list[[name]][[net]][[i]]$nodes[new_list[[name]][[net]][[i]]$nodes$group == "int", ]
          new_name <- na.omit(matching_rows$label)[1]
          
          # Check if the new name is not empty
          if(!is.na(new_name) && nchar(new_name) > 0) {
            
            # Set the name of the list element
            names(new_list[[name]][[net]])[i] <- new_name
          } else {
            warning(paste("Could not find 'int' group or label is NA for element", i, "in", net, "of", name))
          }
        }
        
        # Filtering Step:
        # Only keep the elements of vis_net_A and vis_net_I that match with the changed_genes_list
        new_list[[name]][[net]] <- new_list[[name]][[net]][names(new_list[[name]][[net]]) %in% changed_genes_list[[name]]]
        
      } else {
        warning(paste("Could not find", net, "in", name))
      }
    }
  }
  
celltype_list = new_list
#   # Loop over each cell type in the cell_type_list
#   for(name in names(cell_type_list)) {
#     
#     # Create empty lists to store all the nodes and edges dataframes
#     all_nodes <- list()
#     all_edges <- list()
#     
#     # Loop over both vis_net_A and vis_net_I
#     for(net in c("vis_net_A", "vis_net_I")) {
#       
#       # Check if the name and network exist in cell_type_list
#       if(name %in% names(cell_type_list) && net %in% names(cell_type_list[[name]])) {
#         
#         # Loop over the elements in the vis_net_A or vis_net_I list
#         for(i in seq_along(cell_type_list[[name]][[net]])) {
#           
#           # Add a 'signalling' column and 'int' column to the nodes and edges dataframes
#           cell_type_list[[name]][[net]][[i]]$nodes$signalling <- substr(net, 8, 8)  # 'A' or 'I'
#           cell_type_list[[name]][[net]][[i]]$nodes$int <- names(cell_type_list[[name]][[net]])[i]
#           
#           cell_type_list[[name]][[net]][[i]]$edges$signalling <- substr(net, 8, 8)  # 'A' or 'I'
#           cell_type_list[[name]][[net]][[i]]$edges$int <- names(cell_type_list[[name]][[net]])[i]
#           
#           # Store the nodes and edges dataframes in the all_nodes and all_edges lists
#           all_nodes[[length(all_nodes) + 1]] <- cell_type_list[[name]][[net]][[i]]$nodes
#           all_edges[[length(all_edges) + 1]] <- cell_type_list[[name]][[net]][[i]]$edges
#         }
#       }
#     }
#     
#     # Concatenate all the nodes and edges dataframes together
#     cell_type_list[[name]]$nodes <- do.call(rbind, all_nodes)
#     cell_type_list[[name]]$edges <- do.call(rbind, all_edges)
#     # Remove the vis_net_A and vis_net_I entries from the cell_type_list
#     cell_type_list[[name]] <- cell_type_list[[name]][!names(cell_type_list[[name]]) %in% c("vis_net_A", "vis_net_I")]
#   }
#                                             
# 
# celltype_list = new_list

  # Call the function to modify celltype_list
  # celltype_list <- modify_celltype_list(celltype_list, changed_genes_list)
  
# Extract the file path and name
file_info <- strsplit(hotspots, "/")[[1]]

# Find position of "Real_Data_2"
pos_real_data_2 <- which(file_info == "realsig")

# Build new name from the next position after "Real_Data_2"
path_info <- file_info[(pos_real_data_2 + 1):(length(file_info) - 1)]
name_info <- strsplit(file_info[length(file_info)], "\\.")[[1]][1]

# Build new name
new_name <- paste(path_info, collapse = "_")
new_path <- paste("/home/jarcos/Documents/Final_Mohammed/sigFINAL", paste(new_name, name_info, "filtered_elements_all.rds", sep = "_"), sep = "/")

if (length(celltype_list) > 0) {
  # Save the data
  saveRDS(celltype_list, new_path)
} else {
  print(paste("Skipped saving empty list for file:", hotspots))
}

}
