### Author: Matthew Aaron Loberg
### Date: June 25, 2024
### Script: 24-0625_DoubletFinder_Function.R
### Adapted from: 23-1025_DoubletFinder_Function.R
# Only difference is the addition of a dir.create(outputdir)

# Goal: make a function for detecting doublets using scDblFinder

Doublet_Detection <- function(SO, outputdir, resolution = 0.4){
  dir.create(outputdir)

  SO_Processed <- SO %>% SCTransform(vst.flavor = "v2") %>%
    RunPCA(npcs = 50)

  ggsave(file.path(outputdir, "scDblFinder_elbow.png"),
         ElbowPlot(SO_Processed, ndims = 50),
         height = 5, width = 5, dpi = 600)

  SO_Processed <- SO_Processed %>% RunUMAP(dims = 1:30, reduction = "pca") %>%
    FindNeighbors(dims = 1:30, reduction = "pca") %>%
    FindClusters(resolution = resolution)

  ggsave(file.path(outputdir, "scDblFinder_UMAP_Clusters.png"),
         UMAPPlot(SO_Processed), height = 5, width = 5, dpi = 600)

  # Set the default assay to RNA
  DefaultAssay(SO_Processed) <- "RNA"

  # Run scDblFinder on raw RNA acounts with the identified clusters
  sce <- scDblFinder(GetAssayData(SO_Processed, slot = "counts"), clusters = Idents(SO_Processed)) # Should isolate RAW counts - make sure that DefaultAssay is set to RNA to isolate raw counts

  SO_Processed$scDblFinder.score <- sce$scDblFinder.score
  SO_Processed$scDblFinder.class <- sce$scDblFinder.class

  ggsave(file.path(outputdir, "scDblFinder_UMAP_scDblFinder.png"),
         FeaturePlot(SO_Processed, features = c("scDblFinder.score")) + UMAPPlot(SO_Processed, group.by = "scDblFinder.class"),
         width = 5, height = 10, dpi = 600)

  # Add doublet classification to original Seurat Object
  SO$scDblFinder.class <- SO_Processed$scDblFinder.class
  #SO <- SO %>% subset(scDblFinder.class == "singlet") COULD subset down to just singlets before returning ... will NOT do this
  #SO$scDblFinder.class <- NULL # if subsetting, should remove class variable (no utility)
  return(SO) # Return original Seurat Object
}
