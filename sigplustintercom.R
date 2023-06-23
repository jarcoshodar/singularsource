library(InterCom)
library(SigHotSpotter)
library(visNetwork)
library(Seurat)
library(plyr)

run_intercom_analysis <- function(seur, species, out.path, ncores = 16, sig.cutoff = 0, z.score.cutoff = 0) {
  inpmat <- as.matrix(seur@assays$SCT@counts)

  anno.tbl <- data.frame(cell.name = rownames(seur@meta.data),
                         cell.type = unname(seur@active.ident),
                         stringsAsFactors = F)

  anno.tbl$cell.type <- gsub(" ", "\\.", anno.tbl$cell.type)
  anno.tbl$cell.type <- gsub("-", "\\.", anno.tbl$cell.type)
  anno.tbl$cell.type <- gsub("\\+", "\\.", anno.tbl$cell.type)
  anno.tbl$cell.name <- gsub("-", "\\.", anno.tbl$cell.name)

  colnames(inpmat) <- gsub("-", "\\.", colnames(inpmat))

  InterCom(data = inpmat,
           anno.tbl = anno.tbl,
           species = species,
           ncores = ncores,
           sig.cutoff = sig.cutoff,
           z.score.cutoff = z.score.cutoff,
           tissue.name = organ,
           out.path = out.path)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least two arguments (root folder and study directory) are required.")
}

root <- args[1]
study <- args[2]
organ <- ifelse(length(args) >= 3, args[3], NULL)
if (study == 'CalorieRestriction') {species = 'RAT'} else {species = 'MOUSE'}

run_sighotspotter_analysis <- function(combined.sct, species, out.path) {
  
  combined.sct@meta.data$cellcon <-paste0(combined.sct@active.ident, combined.sct@meta.data$condition)
  
  Idents(combined.sct) <- 'cellcon'
  
  # Read the csv files

  folder_path <- file.path(root, study, "results", "Integrated", organ)
  
  csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  final_list <- list()
  
  # Iterate through each csv file
  for (csv_file in csv_files) {
    
    # Extract the pair of cell conditions from the filename
    pair <- gsub("^(.*)\\.csv$", "\\1", basename(csv_file))
    origin_1 <- gsub("^(.*)_vs_(.*)$", "\\1", pair)
    origin_2 <- gsub("^(.*)_vs_(.*)$", "\\2", pair)
    
    # Read the differential SCT dge compiled earlier
    df <- read.csv(csv_file)
    
    df <- df[!is.na(df$Gene) & !is.na(df$avg_log2FC), ]
    
    # Generate the de dataframe
    de <- data.frame(Gene = df$Gene[df$p_val_adj < 0.01],
                     DE = df$avg_log2FC[df$p_val_adj < 0.01])
    
    de$DE <- sign(de$DE)
    
    de_path <- file.path(out.path, "de.txt")
    write.table(de, file = de_path, sep = "\t", row.names = F, col.names = T, quote = F)
    
    out_list <- list()
    
    # Iterate through the pair of cell conditions
    for (origin in c(origin_1, origin_2)) {
      
      print(origin)
      
      inp <- as.data.frame(combined.sct@assays$RNA@data[, rownames(combined.sct@meta.data)[combined.sct@meta.data$cellcon == origin]])
      inp <- cbind(rownames(inp), inp)
      colnames(inp)[1] <- "Gene"
      rownames(inp) <- NULL
      
      inp_path <- file.path(out.path, "inp.txt")
      write.table(inp, file = inp_path, sep = "\t", row.names = F, col.names = T, quote = F)
      
 
      invert_de <- if (origin == origin_1) FALSE else TRUE
  
      
      out_list[[origin]] <- tryCatch(SigHotSpotter_pipeline(species = species,
                                                        input_data = inp_path,
                                                        cutoff = 30,
                                                        DE_Genes_data = de_path,
                                                        percentile = 70,
                                                        invert_DE = invert_de), error = function(e) { NA })
    }
    final_list[[pair]] <- out_list
  }
  
  return(final_list)
}



input_path <- file.path(root, study, "results")
files <- list.files(input_path, pattern = "Integrated_.*\\.rds$", full.names = TRUE)

if (!is.null(organ)) {
  files <- grep(paste0("Integrated_", organ, "\\.rds$"), files, value = TRUE)
}

if (length(files) == 0) {
  stop("No matching RDS files found.")
}


for (file in files) {
  seur <- readRDS(file)
  organ <- gsub("Integrated_(.*?)\\.rds$", "\\1", basename(file))
  conditions <- unique(seur@meta.data$condition)
#  for (sit in conditions) {
#   subseur = subset(seur, subset = condition == sit)
#   out.path_intercom <- file.path(root, study, "results", organ, sit, "intercom")
#  dir.create(out.path_intercom, recursive = TRUE, showWarnings = FALSE)
#  run_intercom_analysis(subseur, species, out.path_intercom)
# }
  out.path_sighotspotter <- file.path(root, study, "results", organ, "sighotspotter")
  dir.create(out.path_sighotspotter, recursive = TRUE, showWarnings = FALSE)
  b = run_sighotspotter_analysis(seur, species, out.path_sighotspotter)
  saveRDS(b, paste0(out.path_sighotspotter,"/","hotspots.rds"))
}

