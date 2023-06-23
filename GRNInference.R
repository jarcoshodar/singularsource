library(gurobi) #slam is a dependency
library(GRNOpt)
#install with if necessary:
#devtools::install_git("https://gitlab.lcsb.uni.lu/andras.hartmann/GRNOptR.git", ref ="master")

library(stringr)
library(dplyr)
library(readr)

choicerule = "ID" 

species = "MOUSE"

if (species == "MOUSE") {
  pkn_path = '~/tools/trimicepkn.csv'
  tf_list = read.csv('/home/CICBIOGUNE/jarcos/tools/Mus_musculus_TF_LIST.txt', sep = '\t')
} else if (species == "RAT") {
  pkn_path = '~/tools/triratortmicepkn.csv'
  tf_list = read.csv('/home/CICBIOGUNE/jarcos/tools/Rattus_norvegicus_TF.txt', sep = '\t')
} else {
  stop("no human pkn yet")
}

invert_dge = FALSE # set this to FALSE if you do not want to invert the boolean values

#run the following only if you don't have already booleanized gene xpression
booleanize_csv_files <- function(folder_path) {
  # Target folder is booldata
  bool_folder <- file.path(folder_path, "booldata")
  if (!dir.exists(bool_folder)) {
    dir.create(bool_folder)
  }
  
  csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  
  for (csv_file in csv_files) {
    
    data <- read.csv(csv_file)
    # Remove NA in case any processing lead to different length of columns or similar
    #will lead erros somewhere if not
    data_clean <- data %>%
      filter(!is.na(Gene) & !is.na(avg_log2FC) & p_val_adj <= 0.05) 
    #assumes dge, at least elements of Gene already filtered by significance level!!!
    data_bool <- data_clean %>%
      dplyr::mutate(avg_log2FC_bool = ifelse(as.numeric(avg_log2FC) > 0, 1, 0)) %>%
      dplyr::select(Gene, avg_log2FC_bool)
    
    
    # save result on bool_(originalname)
    new_filename <- paste0("bool_", basename(csv_file))
    write_csv(data_bool, file.path(bool_folder, new_filename))
  }
}

#folder_path = "/home/CICBIOGUNE/jarcos/Documentos/Real_Data/CalorieRestriction/Integrated/WAT"
#again, skip if you already boolanized

#folder_path = "/home/CICBIOGUNE/jarcos/Documentos/Real_Data_2/Exercise/Integrated/Testis"


#pkn_path = '/home/CICBIOGUNE/jarcos/Documentos/TF interac/short_mm.csv'

pkn =  read.csv(pkn_path,sep = ',',header=TRUE,stringsAsFactors=FALSE,quote="") 
if (dim(pkn)[2] == 4){
  pkn = pkn[,-1]
}

colnames(pkn) <- c("TF", "effect", "target") 
tf_symbols <- tf_list$Symbol
pkn_filtered <- pkn[pkn$TF %in% tf_symbols, ]

pkn <- pkn_filtered
pkn <- pkn %>% distinct()
#pkn_path = '/home/CICBIOGUNE/jarcos/Documentos/TF interac/MOUSE_TF_TF_interaction.txt'
#pkn =  read.csv(pkn_path,sep = '\t',header=TRUE,stringsAsFactors=FALSE,quote="") #check if uppercase/sql conversion to symbol

#write.csv(pkn, '~/tools/may_rat_pkn.csv')
#write.csv(pkn, '~/tools/tf_only_may_rat_pkn.csv')


for (organ in list.dirs("/home/CICBIOGUNE/jarcos/Documentos/Real_Data_2/Reprogramming2/results3/Integrated", recursive=FALSE)) {
  {
  booleanize_csv_files(organ)
  
  bool_folder = file.path(organ,'booldata')
  CellTypes = list.files(bool_folder) #called because single cell data but can be whatever
  
  for (dir in CellTypes) {
    
    comp = substr(dir, 6, nchar(dir)) #remove start of bool versio name
    filename = paste0('GRNI5',comp) #outputfiles will be called GRNI(original dge name)
    bool_file = file.path(bool_folder,dir)
    bool_data = read.csv(bool_file,sep=',',header=TRUE,stringsAsFactors=FALSE,quote="")
    #again, prune gurobi is strict with column names or leads to uninformative errors
    colnames(bool_data) <- c("TF", "TF_value")   
    
    if (invert_dge) {
      bool_data$TF_value <- ifelse(bool_data$TF_value == 1, 0, 1)
    }
    
    #in our case, gene regulatory network inference, must be TFs
    #other matching maybe necessary based on your previous knowledge network
    #lack of matching can lead to errors
    bool_tfonly <- bool_data[bool_data$TF %in% unique(pkn$TF) | bool_data$TF %in% unique(pkn$target),]
    output_folder = file.path(organ, 'GRNI5')
    if(!dir.exists(file.path(output_folder))) {
      dir.create(file.path(output_folder))
    }
    #if this takes long, probably there was an error in filtering or column name steps:
    #(or your data is super big)
    tryCatch({
      result = prune_gurobi(pkn, bool_tfonly, rule = choicerule)
      print(dir)
      print(result)
      write.csv(result, file.path(output_folder,filename))
    },
    #errors to infer sometimes are expected and must be handlded
    error = function(e){
      print('failure')
      write.csv(bool_tfonly, file.path(output_folder, paste0('FAILGUROBI',dir)))
    })
  }
  }
}