
## calculate UMAP
calcUMAP <- function(s, UMAP.resolution = 0.06) {
    s <- SCTransform(s)
    s <- RunPCA(s, verbose = FALSE)
    s <- FindNeighbors(s, dims = 1:10)
    s <- FindClusters(s, resolution = UMAP.resolution)
    s <- RunUMAP(s, dims = 1:10)
}
