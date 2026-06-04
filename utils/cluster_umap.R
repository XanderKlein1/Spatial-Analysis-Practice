#cluster_umap.R -------------------
#Provides clustering functionality for spatial transcriptomics data.
#
#Includes functionality for:
# - finding neighborhoods
# - defining clusters
# - projecting clusters into 2 dimensions (UMAP)
#
#Notes:
# - defines reusable function(s) only -- no top level execution.

library(Seurat)

cluster_umap <- function(object, num_pca, cluster_resolution) {
  object <- FindNeighbors(object, dims = 1:num_pca)
  object <- FindClusters(object, resolution = cluster_resolution)
  object <- RunUMAP(object, dims = 1:num_pca)
  object
}