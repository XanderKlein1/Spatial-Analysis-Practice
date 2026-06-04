#custom_heatmaps.R
#Functions to generate heatmaps of seurat objects.

library(Seurat)
library(pheatmap)
library(ggplot2)
library(patchwork)
#install.packages('devtools')
devtools::install_github('immunogenomics/presto')


#Generates a heatmap of clusters x genes for the top 5 differentially expressed genes in each cluster.
#
#Parameters:
#object --> the seurat object containing all spatial data
#n_genes --> the number of genes per cluster that should be included in the heatmap
create_cluster_heatmap <- function(object, n_genes) {
  marker_list <- lapply(levels(object), function(i) FindMarkers(object, 
                                                                   features = VariableFeatures(object), 
                                                                   ident.1 = i, 
                                                                   min.pct = 0.01, 
                                                                   logfc.threshold = 0.25, test.use = "LR"))
  #Filter to show only the top 5 genes per cluster and remove any duplicates:
  filtered_list <- lapply(marker_list, function(cluster) {
    rownames(cluster[1:n_genes, ])
  })
  gene_list <- unique(do.call(c, filtered_list))
  
  #Create the cluster x gene matrix that we will plot:
  df <- FetchData(object, vars = c(gene_list, 'ident'))
  #Calculate the average expression across each of the genes in gene_list, for each cluster.
  means <- lapply(levels(intestine), function(clust) {
      lapply(df[df$ident==clust, -ncol(df)], mean) 
    })
  #unwind the list (of lists) into a matrix:
  #row = cluster (first row is cluster 1)
  #col = gene
  means_mtx <- do.call(rbind, lapply(means, unlist))
  rownames(means_mtx) <- lapply(c(1:22), function(i) { paste("Cluster", i)})
  
  #generate the heatmap of clusters x genes:
  pheatmap(means_mtx,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Cluster gene expression",
           color = colorRampPalette(c("deepskyblue", "white", "red"))(100),
           breaks = seq(-3,3,length.out = 101))
}

