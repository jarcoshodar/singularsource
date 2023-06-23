
library(Seurat)
library(DoubletFinder)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(cowplot) 
library(dplyr)
library(cluster, quietly = TRUE)
library(ggstatsplot)
library(SCINA)
library(plyr)
library(biomaRt)
library(densvis)
try(library(monocle3))
library(singleseqgset)
library(pheatmap)
library(SeuratWrappers)
library(GSA)
library(preprocessCore)
library(scCustomize)

### sascha functions start here 

## Functions ##
#seur: A Seurat object
#pass: An identifier of the number of times this function is running. Only for plotting.
#org: The organism the data is from. Only relevant if not all metadata info is there.
filterCells <- function(seur,mad.coeff = 3,pass = -1,org = "HUMAN", output_path = sample_dir){
  
  seur@meta.data$Barcodes <- rownames(seur@meta.data)
  #Calculate percent.mito and percent.ribo metadata columns if they are not there
  if(!any(colnames(seur@meta.data) == "percent.mito")){
    if(org == "HUMAN"){
      seur[["percent.mito"]] <- Seurat::PercentageFeatureSet(seur, pattern = "^MT-")
    } else if(org == "MOUSE"){
      seur[["percent.mito"]] <- Seurat::PercentageFeatureSet(seur, pattern = "^mt-")
    } else if(org == "RAT"){
      seur[["percent.mito"]] <- Seurat::PercentageFeatureSet(seur, pattern = "^[Mm]t-")
    }else {
      stop("The specified organism is not supported")
    }
  }
  if(!any(colnames(seur@meta.data) == "percent.ribo")){
    if(org == "HUMAN"){
      seur[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seur, pattern = "(^RPL|^RPS|^MRP)")
    } else if(org == "MOUSE"){
      seur[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seur, pattern = "(^Rpl|^Rps|^Mrp)")
    } else if(org == "RAT"){
      seur[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seur, pattern = "(^Rpl|^Rps|^[Mm]rp)")
    } else {
      stop("The specified organism is not supported")
    }
  }
  # Filtering cells based on percentage of mitochondrial transcripts
  
  non_na_cells <- !is.na(seur@meta.data$percent.mito)
  
  # Subset the Seurat object to keep only cells with non-NA percent.mito values
  seur <- subset(seur, cells = which(non_na_cells))
  
  cell.QC.stat <- seur@meta.data
  
  max.mito.thr <- median(cell.QC.stat$percent.mito) + mad.coeff*mad(cell.QC.stat$percent.mito)
  min.mito.thr <- median(cell.QC.stat$percent.mito) - mad.coeff*mad(cell.QC.stat$percent.mito)
  
  p1 <- ggplot(cell.QC.stat, aes(x=nFeature_RNA, y=percent.mito)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[2])," cells removed\n",
                                           as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[1])," cells remain"), x = 6000, y = -10)
  
  p <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)
  ggsave(paste0(output_path,"/Mitofilter_Marginal_Pass",pass,".png"),plot = p)
  
  cell.QC.stat <- cell.QC.stat %>% filter(percent.mito <= max.mito.thr) %>% filter(percent.mito >= min.mito.thr)
  
  # Filtering cells based on number of genes and transcripts detected
  # Set low and hight thresholds on the number of detected genes
  min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
  max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
  
  # Set hight threshold on the number of transcripts
  max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nCount_RNA))
  
  p1 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2)
  
  p <- ggMarginal(p1, type = "histogram", fill="lightgrey")
  ggsave(paste0(output_path,"/FeatureAndCountFilter_1_Pass",pass,".png"),plot = p)
  
  # Filter cells base on both metrics
  cell.QC.stat <- cell.QC.stat %>%
    filter(log10(nFeature_RNA) > min.features.thr) %>%
    filter(log10(nFeature_RNA) < max.features.thr) %>%
    filter(log10(nCount_RNA) < max.count.thr)
  
  lm.model <- lm(data = cell.QC.stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
  
  p2 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2) +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
    annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1], " QC passed cells"), x = 4, y = 3.8)
  
  p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
  ggsave(paste0(output_path,"/FeatureAndCountOutlier_2_Pass",pass,".png"),plot = p)
  
  # Cells to exclude lie below an intercept offset of -0.09
  cell.QC.stat$validCells <- log10(cell.QC.stat$nFeature_RNA) > (log10(cell.QC.stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
  
  p3 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(aes(colour = validCells), fill = "black",pch=21) +
    scale_color_manual(breaks = c("TRUE", "FALSE"),
                       values = c("black","firebrick1")) +
    geom_smooth(method="lm") +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
    theme(legend.position="none") +
    annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$validCells)[2]), " QC passed cells\n",
                                           as.numeric(table(cell.QC.stat$validCells)[1]), " QC filtered"), x = 4, y = 3.8)
  
  p <- ggMarginal(p3, type = "histogram", fill="lightgrey")
  
  ggsave(paste0(output_path,"/FeatureAndCountOutlier_3_Pass",pass,".png"),plot = p)
  
  # Remove invalid cells
  cell.QC.stat <- cell.QC.stat %>% filter(validCells)
  
  seur <- subset(seur, subset = Barcodes %in% cell.QC.stat$Barcodes)
  return(seur)
}


findOptimalResolution <- function(seur,pcs,output_path = sample_dir){
  #' Find Optimal Resolution for Clustering
  #'
  #' This function takes a Seurat object and pc comps number,
  #' estimates the Calinski-Harabasz index for a given range of resolutions,
  #' and ends up implementing the one with the highest score.
  #' This implementation thus rewards within and between cluster heterogeneity
  #' So far often several cell types are in a cluster though, so future
  #' functions implement ways of subsetting and merging clusters for labeling
  #' A plot of the Calinski-Harabasz index across resolutions is saved
  #' A modified Seurat object with the winning clustering implemented is returned
  #' 
  #' @param seur A Seurat object.
  #' @param pcs A integer range specifying the number of principal components to use. 
  #' Example use is 1:50 or 1:numPCs when numPCs was previously defined as 50 
  #' @param output_path (Optional) A string specifying the directory to save the plot. Default is `sample_dir`.
  #' Sample_dir is defined by the workflow in global environment for loops, but
  #' any other valid path can be given
  #' @return A modified Seurat object with the winning clustering implemented.
  #'
  #' @examples
  #' seur <- findOptimalResolution(seur, pcs = 1:50)
  
  cds <- as.cell_data_set(seur)
  emb <- Embeddings(object = seur@reductions$pca)[, pcs]
  resolutions <- c(rbind(10^seq(-7,0),5*10^seq(-7,0)))
  resolutions <- rev(rev(resolutions)[-1])
  idx_list <- list()
  for(i in 1:length(resolutions)){
    print(paste0("Resolution: ",resolutions[i]))
    cds <- cluster_cells(cds,
                         k = 10,
                         reduction_method = "UMAP",
                         cluster_method = "leiden",
                         resolution = resolutions[i],
                         num_iter = 5)
    part <- as.integer(clusters(cds))
    idx_list[[i]] <- fpc::calinhara(emb,part)
  }
  
  df <- data.frame(Resolution = resolutions,
                   Value = unlist(idx_list),
                   stringsAsFactors = F)
  
  p <- ggplot(data=df, aes(x=Resolution, y=Value)) +
    geom_line()+
    geom_point()+
    ylab("Calinski Harabasz Index")+
    xlab("Resolution")
  ggsave(paste0(output_path,"/Calinski_Harabasz_Index.png"),plot = p)
  
  res <- df$Resolution[which(df$Value == max(df$Value,na.rm = T))]
  cds <- cluster_cells(cds,
                       k = 10,
                       reduction_method = "UMAP",
                       cluster_method = "leiden",
                       resolution = res,
                       num_iter = 5)
  seur@meta.data$seurat_clusters <- clusters(cds)
  if (length(unique(clusters(cds))) > 30) {
    cds <- cluster_cells(cds, "UMAP")
    part <-partitions(cds, "UMAP")
    
    seur@meta.data$seurat_clusters <- part}
  
  Idents(seur) <- 'seurat_clusters'
  
  return(seur)
}

annotateCells <- function(seur,org = c("HUMAN","MOUSE","RAT"), tissue = "ALL",id.type = c("Symbol","Ensembl"),  output_path = paste0(sample_dir,"/UMAP_SCINA.png")){
  
  id.type <- match.arg(id.type, c("Symbol","Ensembl"))
  org <- match.arg(org, c("HUMAN","MOUSE","RAT"))
  
  inp_df <- as.matrix(seur@assays$RNA[,])
  inp_df <- log(inp_df+1)
  inp_df[] = normalize.quantiles(as.matrix(inp_df))
  
  if(org == "HUMAN"){
    markers <- cm_human
    mart_org <- "hsapiens_gene_ensembl"
    filter <- "hgnc_symbol"
  }else if(org == "MOUSE"){
    markers <- cm_mouse
    mart_org <- "mmusculus_gene_ensembl"
    filter <- "mgi_symbol"
  }else if(org == "RAT"){
    markers <- cm_rat
    mart_org <- "rnorvegicus_gene_ensembl"
    filter <- "rgd_symbol"
  }
  
  if(tissue != "ALL"){
    markers <- markers[which(grepl(tissue,markers$tissue_class)),]
  }
  markers <- markers[markers$cancer_type == "Normal",]
  markers <- markers[markers$cell_type != "Cancer cell",]
  markers <- markers[,c("cell_name","Symbol")]
  markers <- plyr::ddply(markers,.(cell_name,Symbol),nrow)
  markers <- markers[markers$Symbol != "",]
  markers <- markers[markers$V1 > 1,]
  markers <- markers[,c("cell_name","Symbol")]
  
  if(id.type == "Ensembl"){
    ##the end part of the loop is commented away. it's how the script originally ran
    ##to avoid problems with relying on Esembl being online and available for the computer cluster,
    ##it was turned into files to be read instead. you can comment away the load line and run it to get more recent data instead
    
    markers_file <- file.path(root, supplfiles, paste0("markers_", org, ".RData"))
    
    load(file = markers_file)
    # require(biomaRt)
    # mart <- useDataset(mart_org, useMart("ensembl"))
    # symbToEns <- getBM(filters= filter, attributes= c("ensembl_gene_id",filter),values=markers$Symbol,mart= mart)
    # markers <- merge(markers,symbToEns,by.x = 2, by.y = 2, all.x = F)
    # markers <- markers[,-1]
    # colnames(markers) <- c("cell_name","Symbol")
  }
  
  sigs <- lapply(unique(markers$cell_name),function(x){markers$Symbol[markers$cell_name == x]})
  names(sigs) <- unique(markers$cell_name)
  sigs <- lapply(sigs,function(x){intersect(x,rownames(inp_df))})
  sigs <- sigs[sapply(sigs,length) > 3]
  
  # Try SCINA with rm_overlap = F first, and if it fails, try with rm_overlap = T
  scina.results <- tryCatch({
    SCINA(inp_df,
          sigs,
          max_iter = 100,
          convergence_n = 10,
          convergence_rate = 0.99,
          sensitivity_cutoff = 0.9,
          rm_overlap = F,
          allow_unknown = T
    )
  }, error = function(e) {
    SCINA(inp_df,
          sigs,
          max_iter = 100,
          convergence_n = 10,
          convergence_rate = 0.99,
          sensitivity_cutoff = 0.9,
          rm_overlap = T,
          allow_unknown = T
    )
  })
  
  seur$SCINA_annot <- scina.results$cell_labels
  
  ggsave(output_path,DimPlot(seur, reduction = "umap", group.by = "SCINA_annot"))
  return(seur)
}
scina_counts <- function(cluster_number, seur = choice) {
  cluster_df <- data.frame(Cluster = seur@active.ident,
                           SCINA_annot = seur@meta.data$SCINA_annot)
  
  cluster_counts <- cluster_df %>%
    dplyr::filter(Cluster == as.character(cluster_number)) %>%
    dplyr::group_by(SCINA_annot) %>%
    dplyr::summarize(Count = n())
  
  return(cluster_counts)
}

scustom_feature_plots <- function(gene, seur = choice){
  DefaultAssay(seur) <-'RNA'
  plot <- try(FeaturePlot_scCustom(seur, gene, colors_use = viridis_plasma_dark_high, na_color = "lightgray"))
  DefaultAssay(seur) <-'integrated'
  return(plot)
}
bring_top_10 <- function(markers_df) {
  # Get the top 10 genes for each cluster
  top_10_genes_by_cluster <- markers_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(10, avg_log2FC) %>%
    dplyr::ungroup()
  
  # Convert the tibble to a list with one element per cluster
  genes_list <- split(top_10_genes_by_cluster$gene, top_10_genes_by_cluster$cluster)
  
  # Print the gene names for each cluster
  for (i in seq_along(genes_list)) {
    cat(paste("Cluster", i, ":", paste(genes_list[[i]], collapse = ", ")), "\n\n")
  }
  
  # Return the list of genes for further usage
  return(genes_list)
}


runDensUMAP <- function(seur, pcs, sample_dir = ".") {
  # Use tryCatch to handle errors within the function, since basilisk is capricious and might fail
  result <- tryCatch({
    require(densvis)
    emb <- Embeddings(object = seur@reductions$pca)[, pcs]
    fit <- densmap(emb,
                   dens_frac = 0.5,
                   dens_lambda = 0.5,
                   metric = "cosine",
                   n_epochs = 3000L,
                   random_state = 100L
    )
    colnames(fit) <- paste0("DensUMAP_", 1:2)
    rownames(fit) <- rownames(emb)
    seur[["densumap"]] <- CreateDimReducObject(embeddings = fit,
                                               key = "DensUMAP_",
                                               assay = DefaultAssay(seur)
    )
    return(seur)
  }, error = function(e) {
    # If an error occurs, create a file called "failed DensUMAP"
    write("DensUMAP failed", file = paste0(sample_dir, "/failed_DensUMAP.txt"))
    
    # Return the original Seurat object
    return(seur)
  })
  
  return(result)
}


findNumPCs <- function(seur){
  pct <- seur[["pca"]]@stdev / sum(seur[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
  
  return(min(co1,co2))
}



##jarcoshodar functions begin here

#sample_dir = '/home/jarcos/Documents/Real_Data/CalorieRestriction/Aorta'

ExploreCluster <- function(seur, cluster, runPCAUMAP = FALSE,
                           ucomps =50,comps=100, verb = T, k=20) {
  #' ExploreCluster
  #'
  #' This function takes a Seurat object and a cluster number or name, and
  #' optionally reruns PCA and UMAP on the selected cluster. Note that despite the
  #' option being offered, it is NOT recommended and it seems that using the
  #' clustering from the original data leads to better results.
  #' It returns a subsetted Seurat object with numbered subclusters
  #' branching from the original name separated by a dot.
  #' Follow up with RemapClusters to put them back in original object
  #'  recommended, before or after identifying and renaming them.
  #'
  #' @param seur A Seurat object of single-cell data with the 'integrated' assay.
  #' @param cluster An integer or name that can be found by the subset function.
  #' @param runPCAUMAP (Optional) A boolean indicating whether to rerun PCA and UMAP. Default is FALSE. Not recommended to set it to TRUE unless results with FALSE not satisfactory.
  #' @param comps (Optional) Number of components for PCA if ran. Default is 100.
  #' @param ucomps (Optional) Number of PCA components to use in UMAP. Default is 50.
  #' @param verb (Optional) A boolean to decide whether to print messages during PCA. Default is TRUE.
  #' @param k (Optional) Number of nearest neighbors for UMAP. Default is 20.
  #'
  #' @return A modified Seurat object of the selected cluster divided into subclusters (e.g., 2.1, 2.2, 2.3).
  #'
  #' @examples
  #' seur <- ExploreCluster(seur, 2)
  #' seur <- ExploreCluster(seur, 2), with no additional subclustering due to
  #' default runPCAUMAP to false
  #' returns subset seurat object with cells in cluster 2 divided in 2.1, 2.2, 2.3
  
  sub <- subset(seur, subset = seurat_clusters == cluster)
  
  if (runPCAUMAP) {
    numPCs = ucomps
    sub <- RunPCA(sub, assay = 'integrated', npcs = comps, verbose = verb)
    sub <- RunUMAP(sub, dims = 1:numPCs, n.neighbors = k)
  }
  
  sub <- findOptimalResolution(sub)
  sub$sub_backup <- sub$seurat_clusters
  return(sub)
}
#DimPlot(seur, label=TRUE)

#mcluster = 3 

#sub <- ExploreCluster(seur, mcluster, FALSE)

#DimPlot(sub, label=TRUE)

RemapClusters <- function(seur, sub, mcluster) {
  #' Remap Clusters
  #'
  #' Taking your original integrated Seurat object, a sub object from 
  #' ExploreCluster and a specific cluster that you already ran ExploreCLuster on,
  #' It creates a new metadata column, "RefinedClusters", putting them back in
  #' original object. In sub_backup you have all the original numeric clusters
  #' 
  #' If the cluster is a number RefinedClusters will have
  #' original seurat_cluster + . + subcluster number
  #' If it's already a cell type name it's just put into place
  #' Both cases matching by UMIs
  #'
  #' @param seur A Seurat object containing all the original data
  #' @param sub A Seurat object containing the subset clusters to be mapped back into the seur object.
  #' @param mcluster The cluster integer that was used to create the subset object.
  #'
  #' @return A modified Seurat object with the subset clusters remapped and added to the "RefinedClusters" metadata column.
  #'
  #' @examples
  #' seur <- RemapClusters(seur, sub, 3)
  #' 
  #' 
  
  is_convertible_to_numeric <-function(x) {
    !is.na(suppressWarnings(as.numeric(x)))
  }
  
  
  if ("sub_backup" %in% colnames(sub@meta.data)) {
    seur$sub_backup <- seur$seurat_clusters
    
    seur@meta.data[rownames(sub@meta.data), "sub_backup"] <- sub$sub_backup
    
  } else {
    seur$sub_backup <- seur$seurat_clusters
  }
  
  
  if (!"RefinedClusters" %in% colnames(seur@meta.data)) {
    seur$RefinedClusters <- as.character(seur$seurat_clusters)
  }
  
  # Find the indices of the relevant mcluster in RefinedClusters
  mcluster_indices <- which(seur$RefinedClusters == as.character(mcluster))
  
  # Check which rows in the sub@meta.data correspond to the mcluster_indices
  sub_rows <- intersect(rownames(sub@meta.data), rownames(seur@meta.data[mcluster_indices, ]))
  
  labs <- sapply(sub$seurat_clusters[sub_rows], function(cluster_label) {
    if (is_convertible_to_numeric(cluster_label)) { #only if numbers is the label 1 in 2 -> 2.1 etc. generated.
      paste0(mcluster, ".", cluster_label)
    } else {
      cluster_label
    }
  })
  
  # Update only the relevant cluster labels, avoid overwriting over repeated runs
  seur@meta.data[sub_rows, "RefinedClusters"] <- labs
  
  return(seur)
}

#seur <- RemapClusters(seur, sub, mcluster)

#DimPlot(seur, label=TRUE, group.by = "RefinedClusters")

#DimPlot(seur, label=TRUE)

#mcluster = 2

#sub <- ExploreCluster(seur, mcluster, FALSE)

#DimPlot(sub, label=TRUE)

#seur <- RemapClusters(seur, sub, mcluster)

#DimPlot(seur, label=TRUE, group.by = "RefinedClusters")

ProcessClusters <- function(seur, clusters, runPCAUMAP = rep(FALSE, length(clusters))) {
  #' Process Clusters
  #'
  #' This function takes an integrated Seurat object, a vector of clusters to process, and an optional
  #' boolean vector indicating whether to run PCA and UMAP for each cluster. It runs ExploreCluster and
  #' RemapClusters functions serially for each specified cluster, giving you more clusters only in the selected clusters list.
  #' 
  #' Intended workflow is to use this when initial exploration of gene markers 
  #' suggests that some clusters include too many different types of cell
  #' and you must break them down for proper identification.
  #'
  #' @param seur A Seurat object containing all the original data
  #' @param clusters A vector of integers specifying the clusters to process.
  #' @param runPCAUMAP (Optional) A boolean vector indicating whether to run PCA and UMAP for each cluster.
  #'                   Default is FALSE for each cluster. The length of runPCAUMAP must match the length of clusters.
  #'
  #' @return A modified Seurat object with refined clustering in the "RefinedClusters" metadata column.
  #'
  #' @examples
  #' clusters_to_process <- c(2, 3, 4)
  #' runPCAUMAP <- c(FALSE, TRUE, FALSE) #note false at all types recommended
  #' seur <- ProcessClusters(seur, clusters_to_process, runPCAUMAP)
  #' 
  if (length(runPCAUMAP) != length(clusters)) {
    stop("The length of runPCAUMAP must match the length of clusters.")
  }
  
  for (i in seq_along(clusters)) {
    mcluster <- clusters[i]
    sub <- ExploreCluster(seur, mcluster, runPCAUMAP[i])
    seur <- RemapClusters(seur, sub, mcluster)
  }
  return(seur)
}

#clusters_to_remerge <- list(c("2.1", "2.2", "2.4"))
#seur <- RmergeClusters(seur, clusters_to_remerge)


#DimPlot(seur, label=TRUE, group.by = "RefinedClusters")

#unique(combine)

ReorderRefClusters <- function(seur) {
  unique_labels <- unique(seur@meta.data$RefinedClusters)
  
  # Function to sort cluster labels in the desired order
  sort_clusters <- function(x) {
    as.numeric(unlist(lapply(strsplit(x, "\\."), function(y) paste0(y, collapse = "."))))
  }
  
  # Sort the unique cluster labels
  sorted_labels <- unique_labels[order(sort_clusters(unique_labels))]
  
  # Create a named vector to map the original cluster labels to new consecutive integers
  new_labels <- setNames(as.character(seq_along(sorted_labels)), sorted_labels)
  
  # Update the RefinedClusters with the new consecutive integer labels
  seur@meta.data$RefinedClusters <- new_labels[seur$RefinedClusters]
  #ensure ordered by number
  sorted_RefinedClusters <- sort(as.numeric(unique(seur@meta.data$RefinedClusters)), decreasing = FALSE)
  seur@meta.data$RefinedClusters <- factor(seur@meta.data$RefinedClusters, levels = as.character(sorted_RefinedClusters))

  return(seur)
} 
#seur <- ReorderRefinedClusters(seur)

truplot <- function (gene, seur = choice) {
  # Get gene names from the Seurat object
  gene_names <- rownames(seur@assays$RNA@counts)
  DefaultAssay(seur) <- 'RNA'
  # Define a function to create the FeaturePlot
  create_feature_plot <- function(gene, seur=choice) {
    plot <- FeaturePlot_scCustom(seur, features = c(gene), colors = viridis_plasma_dark_high, na_color = "lightgray")
    plot <- plot + ggtitle(gene)
    return(plot)
  }
  
  # Try to plot the gene, catch errors
  plot_result <- tryCatch({
    plot <- create_feature_plot(gene, seur)
    list(success = TRUE, plot = plot)
  }, error = function(e) {
    list(success = FALSE, message = e$message)
  })
  
  # If there was an error, search for the closest gene match and plot it
  if (!plot_result$success) {
    message("Error: ", plot_result$message)
    message("Attempting to find the closest gene match...")
    
    # Search for the closest match to the input gene
    matched_indices <- grep(gene, gene_names, ignore.case = TRUE)
    
    # Check if a match was found
    if (length(matched_indices) > 0) {
      # Update the gene variable with the matched gene name
      gene <- gene_names[matched_indices[1]]
      
      # Plot the matched gene
      message("Plotting the closest gene match: ", gene)
      plot <- create_feature_plot(gene, seur)
    } else {
      message("No close gene match found.")
      plot <- NULL
    }
  } else {
    plot <- plot_result$plot
  }
  
  return(plot)
}

find_top_markers <- function(seurat_obj, clusters) {
  # Initialize an empty list to store the top 10 genes for each cluster
  top_markers_list <- list()
  
  # Loop through the provided cluster names
  for (cluster in clusters) {
    # Run FindMarkers for the current cluster
    dge <- FindMarkers(seurat_obj, ident.1 = cluster)
    
    # Order the genes by absolute avg_log2FC and extract the top 10 genes
    top_10_genes <- dge[order(abs(dge$avg_log2FC), decreasing = TRUE), ][1:10, ]
    top_10_gene_names <- rownames(top_10_genes)
    
    # Store the top 10 genes in the list, using the cluster name as the key
    top_markers_list[[cluster]] <- top_10_gene_names
    
    # Print the top 10 genes for the current cluster
    cat("Top 10 genes for cluster", cluster, ":\n")
    print(top_10_gene_names)
    cat("\n")
    top_markers_list[[cluster]] <- dge
  }
  
  # Return the list with the top 10 genes for each cluster
  return(top_markers_list)
}


MergeClusters <- function(seur, clusters_to_merge, use_regex = FALSE) {
  #' Merge Clusters
  #'
  #' This function merges clusters in the RefinedClusters metadata column of a Seurat object.
  #' You can choose between two behaviors: merging based on regex patterns or merging specific clusters.
  #'
  #' @param seur A Seurat object containing the original data.
  #' @param clusters_to_merge A list of clusters to merge.
  #'   If use_regex is TRUE, clusters_to_merge should be a list of vectors, where each vector contains clusters to merge based on regex patterns.
  #'   If use_regex is FALSE, clusters_to_merge should be a list of pairs, where each pair contains two clusters to merge.
  #' @param use_regex A boolean indicating whether to use regex patterns for merging clusters (default is FALSE).
  #'   If use_regex is TRUE, clusters will be merged based on regex patterns.
  #'   If use_regex is FALSE, specific clusters will be merged as specified in the input list.
  #'   Note final behavior is different!!! FALSE and all subclusters will not give you equivalent to regex TRUE
  #'   
  #' @return A modified Seurat object with the specified clusters merged in the "RefinedClusters" metadata column.
  #'
  #' @examples
  #' # Merge clusters 2.1, 2.2, and 2.3 into a single cluster named 2 using the regex behavior
  #' seur <- MergeClusters(seur, list(c("2.1", "2.2", "2.3")), use_regex = TRUE)
  #'
  #' # Merge cluster 2.3 into cluster 2.1 using the specific cluster merging behavior
  #' seur <- MergeClusters(seur, list(c("2.1", "2.3")), use_regex = FALSE)
  if (use_regex) {
    for (cluster_group in clusters_to_merge) {
      # Create a regex pattern to match the specified clusters
      pattern <- paste0("^", paste(cluster_group, collapse = "|^"), "$")
      
      # Find the indices of the clusters to remerge
      remerge_indices <- grep(pattern, seur$RefinedClusters)
      
      # Extract the main cluster number and replace the RefinedClusters labels for the specified clusters
      main_cluster <- as.character(as.integer(strsplit(cluster_group[1], "\\.")[[1]][1]))
      seur@meta.data[remerge_indices, "RefinedClusters"] <- main_cluster
    }
  } else {
    for (pair in clusters_to_merge) {
      # Find the indices of the clusters to merge
      merge_indices <- which(seur$RefinedClusters == as.character(pair[2]))
      
      # Update the RefinedClusters labels for the specified clusters
      seur@meta.data[merge_indices, "RefinedClusters"] <- as.character(pair[1])
    }
    
    # Reorder the clusters after merging (cancelled)
    #seur <- ReorderRefClusters(seur)
  }
  
  return(seur)
}

#merge_pairs <- list(c(5, 7))
#seur <- MergeClusters(seur, merge_pairs)
#DimPlot(seur, label=TRUE, group.by = "RefinedClusters")


##Calorie Restriction functions from jarcoshodar start here
acronyms <- 'Fib, Fibroblast; M2, Macrophage II; ASC, Adipose-derived stem cell; M1, Macrophage I; Neu, Neutrophil; EC, Endothelial cell; DC1, Irf8- Dendritic cell; BC1, Krt15+ Basal cell; Neu (PT), Neutrophil (peripheral tissues); KC, Kupffer cell; T2, CD8+ T cell; T1, CD4-CD8- T cell; PT, Proximal tubule; SMC, Smooth muscle cell; NK, Natural killer cell; BC2, Krt14+ Basal cell; LOH, Thick ascending limb of the loop of Henle; HFC, Hair follicle cell; ProN, Pro-neutrophil; Pla, Plasmocyte; DC2, Irf8+ Dendritic cell; PB, Pro-B cell; Mon, Monocyte; NKT, Natural killer T cell; Mit2, Mitotic cell 2; LC, Langerhans cell; B, B cell; DLOH, Descending loop of Henle; CD2, Intercalated cell of the collecting duct, type A; Spi, Spinous cell; CD3, Intercalated cell of the collecting duct, type B; Mit1, Mitotic cell 1; Bas, Basophil; CC, Channel cell; CD1, Principal cell of the collecting duct; IB, Immature B cell; Hep, Hepatocyte; ProE, Pro-erythroblast; Ery, Erythroblast; Epi, Epithelial cell; LPB, Late pro-B cell; Cho, Cholangiocyte.'

enable_cross_check <- function(acronyms) {
  ref <- strsplit(acronyms, '; ')[[1]]
  todict <- list()
  for (i in ref) {
    mlem <- strsplit(i, ',')[[1]]
    todict[[mlem[1]]] <- mlem[2]
  }
  return(todict)
}

AcroDict = enable_cross_check(acronyms)

#their_excel = '/home/CICBIOGUNE/jarcos/Documentos/CR_current/supplfiles/inspiration.xlsx'

their_excel = '/home/CICBIOGUNE/jarcos/Documentos/inspiration.xlsx'

Expect_Cell_Types <- function(organ, pathtoattrbook) {
  if (organ == 'BM') {
    sheet <- 'Bone marrow'
  } else if (organ == 'Muscle') {
    sheet <- 'Skeletal muscle'
  } else {
    sheet <- organ
  }
  theirdata <- readxl::read_excel(pathtoattrbook, sheet = sheet, skip = 2)
  temp <- table(theirdata$'Cell Type')
  theyfound <- names(temp)
  print(theyfound)
  return(theyfound)
}

Markers_From_CR_excel <- function(pathtoattrbook) {
  sheet <- readxl::read_excel(pathtoattrbook, sheet = 'Marker genes', skip = 1)
  
  output <- list()
  
  cell <- sheet$'Cell Type'
  geneset <- sheet$'Marker Genes'
  for (i in seq_along(cell)) {
    j <- geneset[i]
    j <- gsub(' ', '', j)
    j <- gsub('+', '', j)
    j <- tolower(j)
    j <- strsplit(j, ',')[[1]]
    j <- j[!grepl('d/e/g', j)]
    j <- unique(j)
    output[[cell[i]]] <- j
  }
  return(output)
}



FilterMarkers <- function(all_markers, expected, acronym_dict) {
  filtered_markers <- list()
  
  for (short_name in expected) {
    full_name <- acronym_dict[[short_name]]
    if (!is.null(full_name)) {
      # Calculate string distances between full_name and all_markers names
      string_dists <- stringdist::stringdist(tolower(full_name), tolower(names(all_markers)), method = "jw")
      
      # Find the index of the closest match
      closest_index <- which.min(string_dists)
      
      # Add the closest match to the filtered_markers
      closest_full_name <- names(all_markers)[closest_index]
      filtered_markers[[closest_full_name]] <- all_markers[[closest_full_name]]
    }
  }
  
  return(filtered_markers)
}

CreatePlots <- function(seur, filtered_markers) {
  DefaultAssay(seur) <- 'RNA'
  plots <- list()
  
  for (cell_type in names(filtered_markers)) {
    markers <- filtered_markers[[cell_type]]
    title <- cell_type
    
    cleaned_markers <- c()
    for (marker in markers) {
      if (grepl("\\+$", marker)) {
        marker <- gsub("\\+$", "", marker)
      } else if (grepl("\\-$", marker)) {
        title <- paste(title, paste(toupper(substr(marker, 1, 1)), tolower(substr(marker, 2, nchar(marker) - 1)), "SHOULD NOT", sep = ""))
        marker <- gsub("\\-$", "", marker)
      }
      marker <- paste0(toupper(substr(marker, 1, 1)), tolower(substr(marker, 2, nchar(marker))))
      cleaned_markers <- append(cleaned_markers, marker)
    }
    
    tryCatch({
      p <- FeaturePlot_scCustom(seur, features = cleaned_markers,colors_use = viridis_plasma_dark_high, na_color = "lightgray") + ggtitle(title)
      plots[[cell_type]] <- p
    }, error = function(e) {
      message(paste("Error creating plot for", cell_type, ":", e$message))
    })
  }
  DefaultAssay(seur) <- 'integrated'
  return(plots)
}

format_iden_list <- function(iden_list) {
  # Split the string into lines
  lines <- strsplit(iden_list, "\n")[[1]]
  
  # Extract the second part of each line (after the space)
  strings <- sapply(lines, function(line) {
    parts <- strsplit(line, " ")[[1]]
    return(parts[2])
  })
  
  # Create a character vector
  result <- c(strings)
  
  return(result)
}
update_cell_ids <- function(seurat_obj, new_idens = NULL) {
  # If new_idens is provided, create a named vector from it
  if (!is.null(new_idens)) {
    new_idens_vector <- as.character(new_idens)
    names(new_idens_vector) <- as.character(levels(seurat_obj@active.ident))
  }
  
  # Check if RefinedClusters exists in the metadata
  if ("RefinedClusters" %in% colnames(seurat_obj@meta.data)) {
    column_to_use <- "RefinedClusters"
  } else {
    column_to_use <- "active.ident"
  }
  
  # Update cell_ids based on either RefinedClusters or active.ident
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(cell_ids = if (is.null(new_idens)) {
      mapvalues(as.character(get(column_to_use)), from = names(new_idens_vector), to = new_idens_vector)
    } else {
      as.character(get(column_to_use))
    })
  
  # Set the active identity to the new cell_ids column
  Idents(seurat_obj) <- seurat_obj@meta.data$cell_ids
  
  return(seurat_obj)
}
umap <- function(seur, genes = NULL) {
  if (is.null(genes)) {
    if ("RefinedClusters" %in% colnames(seur@meta.data)) {
      group_by_col <- "RefinedClusters"
      plot <- DimPlot(seur, label = TRUE, group.by = group_by_col)
    } else {
      DimPlot(seur, label = TRUE)
    }

  } else {
    viridis_plasma_dark_high <- viridis(length(genes), end = 0.9, begin = 0.1, direction = 1)
    plot <- FeaturePlot(seur, features = genes, colors_use = viridis_plasma_dark_high, na_color = "lightgray")
  }
  return(plot)
}
