### Author: Matthew Aaron Loberg
### Date: October 26, 2023
### Script: 23-1026_SingleR_Prediction_Function.R

### Goal:
# Run SingleR Prediction for human primary cell atlas (HPCA), Myeloid reference, CD4 reference, and CD8 reference

SingleR_Predictions <- function(SO_List, outputdir, savedir){
  # First run single-cell T and Myeloid References
  # Need to establish references
  # Establish CD4 Reference as SCE
  CD4_Reference <- readRDS("~/Chu_2023_Nature_CD4_Reference.RDS") %>% Seurat::as.SingleCellExperiment()
  # Establish CD8 reference as SCE
  CD8_Reference <- readRDS("~/Chu_2023_Nature_CD8_Reference.RDS") %>% Seurat::as.SingleCellExperiment()
  # Establish Myeloid Reference (saved as SCE so no need to convert)
  Myeloid_THCA_Reference <- readRDS("~/GSE154763_THCA_Myeloid_Normalized_Reference.RDS")

  # Now we will run through the HPCA data
  hpca.se <- celldex::HumanPrimaryCellAtlasData()

  # Generate a list object to store myeloid and T cell SingleR annotations
  Myeloid_T_Annotations <- list()
  Myeloid_Annotations <- list()
  CD4_Annotations <- list()
  CD8_Annotations <- list()

  # Generate a list object to store HPCA SingleR annotations
  HPCA_Annotations_Main <- list()
  HPCA_Annotations_Fine <- list()


  # Run a for loop across entire list for determining cell type annotations
  for(i in 1:length(SO_List)){
    Myeloid_T_Annotations[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                                   ref = list(CD4 = CD4_Reference, CD8 = CD8_Reference, Myeloid = Myeloid_THCA_Reference),
                                                   labels = list(CD4_Reference$cell.type, CD8_Reference$cell.type, Myeloid_THCA_Reference$MajorCluster),
                                                   de.method = "wilcox",
                                                   de.n = 25, # Note, default here is 10 - not sure what is optimal (in example vignette there was only a small difference between 10 and 50; I will split the difference w/ 25)
                                                   aggr.ref = TRUE)
    Myeloid_Annotations[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                                 ref = Myeloid_THCA_Reference,
                                                 labels = Myeloid_THCA_Reference$MajorCluster,
                                                 de.method = "wilcox",
                                                 de.n = 25,
                                                 aggr.ref = TRUE)
    CD4_Annotations[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                             ref = CD4_Reference,
                                             labels = CD4_Reference$cell.type,
                                             de.method = "wilcox",
                                             de.n = 25,
                                             aggr.ref = TRUE)
    CD8_Annotations[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                             ref = CD8_Reference,
                                             labels = CD8_Reference$cell.type,
                                             de.method = "wilcox",
                                             de.n = 25,
                                             aggr.ref = TRUE)
    HPCA_Annotations_Main[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                                   ref = hpca.se,
                                                   labels = hpca.se$label.main)
    HPCA_Annotations_Fine[[i]] <- SingleR::SingleR(test = SO_List[[i]] %>% Seurat::as.SingleCellExperiment(assay = "SCT"),
                                                   ref = hpca.se,
                                                   labels = hpca.se$label.fine)
    SO_List[[i]]$SingleR_Myeloid_T <- Myeloid_T_Annotations[[i]]$labels
    SO_List[[i]]$SingleR_Myeloid <- Myeloid_Annotations[[i]]$labels
    SO_List[[i]]$SingleR_CD4 <- CD4_Annotations[[i]]$labels
    SO_List[[i]]$SingleR_CD8 <- CD8_Annotations[[i]]$labels
    SO_List[[i]]$SingleR_HPCA_Main <- HPCA_Annotations_Main[[i]]$labels
    SO_List[[i]]$SingleR_HPCA_Fine <- HPCA_Annotations_Fine[[i]]$labels
    TEMP_SO <- SO_List[[i]]
    TEMP_SO <- TEMP_SO %>% RunPCA(npcs = 50) %>% RunUMAP(dims = 1:30, reduction = "pca")
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "Myeloid_T_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_Myeloid_T"),
                    width = 7, height = 5, dpi = 600)
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "Myeloid_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_Myeloid"),
                    width = 7, height = 5, dpi = 600)
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "CD4_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_CD4"),
                    width = 7, height = 5, dpi = 600)
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "CD8_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_CD8"),
                    width = 7, height = 5, dpi = 600)
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "HPCA_Main_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_HPCA_Main"),
                     width = 7, height = 5, dpi = 600)
    ggplot2::ggsave(file.path(outputdir, SO_List[[i]]$orig.ident[1], "HPCA_Fine_Annotations_UMAP.png"),
                    UMAPPlot(TEMP_SO, group.by = "SingleR_HPCA_Fine"),
                    width = 7, height = 5, dpi = 600)
  }
  rm(CD4_Reference, CD8_Reference, Myeloid_THCA_Reference, hpca.se, TEMP_SO)
  saveRDS(Myeloid_T_Annotations, file = file.path(savedir, "Myeloid_T_Annotations.RDS"))
  saveRDS(Myeloid_Annotations, file = file.path(savedir, "Myeloid_Annotations.RDS"))
  saveRDS(CD4_Annotations, file = file.path(savedir, "CD4_Annotations.RDS"))
  saveRDS(CD8_Annotations, file = file.path(savedir, "CD8_Annotations.RDS"))
  saveRDS(HPCA_Annotations_Main, file = file.path(savedir, "HPCA_Annotations_Main.RDS"))
  saveRDS(HPCA_Annotations_Fine, file = file.path(savedir, "HPCA_Annotations_Fine.RDS"))
  return(SO_List)
}
