### Author: Matthew Aaron Loberg
### Date: August 21, 2024
### Script: FastMNN_Fibroblast_Subclustering_3000_Features.R
### Source Script Name: "24-0821_FastMNN_Fibroblast_Subclustering_3000_Features.R"

# # Goal: Subcluster fibroblast cells (use 24-0819 Fibroblast IDs from script 24-0819_FastMNN_3000.R)
# Here, I will subcluster fibroblasts

##### Load packages #####
library(Seurat)
library(SeuratWrappers)
library(tidyverse) # for ggsave
library(RColorBrewer) # For plot customization

# Read in merged data
Merged_SO <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_for_FastMNN.RDS")

# Read in Fibroblast IDs for subsetting prior to integration
Fibroblast_IDs <- readRDS(file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_Fibroblast_Labels.RDS")

# Subset by Fibroblast IDs
# Label cell names
Merged_SO$Cell_Names <- colnames(Merged_SO)
Merged_SO <- Merged_SO %>% subset(Cell_Names %in% Fibroblast_IDs)

# Test that the new subsetted object is the appropriate length
ncol(Merged_SO)
# ncols is 24,463 which is the right number

# Run FastMNN
Merged_SO <- NormalizeData(Merged_SO)
Merged_SO <- FindVariableFeatures(Merged_SO, nfeatures = 3000)
# There is some confusion on my end on whether NormalizeData should be run prior to or after merge
# In the Seurat CCA integration, it would be normal to do an individual SCTransform prior to integration
# It could make more sense to do a library based multi-sample normalization rather than normalizing after merging
# However, the tutorials that I have seen perform normalization AFTER merging the objects
# See tutorial here:
# https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.md
# Another point of confusion for me is how to select variable features (e.g., individual, etc.)
# There seem to be conflicting pieces of adivce on this. For instance, see the following GitHub issue thread:
# https://github.com/satijalab/seurat-wrappers/issues/15
# For now, I will force FastMNN to use 5000 features and go from there
Merged_SO_FastMNN <- RunFastMNN(object.list = SplitObject(Merged_SO, split.by = "Identifier"), features = 3000)
# Note: Matrix package version needed to be 1.6-1 for RunUMAP with this old version of Seurat
# Ran the following command (after deleting the Matrix package):
# remotes::install_version("Matrix", version = "1.6-1")
Merged_SO_FastMNN <- RunUMAP(Merged_SO_FastMNN, reduction = "mnn", dims = 1:30)
Merged_SO_FastMNN <- FindNeighbors(Merged_SO_FastMNN, reduction = "mnn", dims = 1:30)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.1)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.2)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.3)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.4)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.5)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 1)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.6)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.7)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.9)

# Make a Histology simplified label
# I am doing this here because I want to print out the clusters by their histology
# Histology simplified label
Merged_SO_FastMNN$Histology_Simplified <- Merged_SO_FastMNN$Histology
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$Histology_Simplified
Merged_SO_FastMNN$Histology_Simplified <- Idents(Merged_SO_FastMNN) # Workaround to force levels
levels(Merged_SO_FastMNN$Histology_Simplified)
levels(Merged_SO_FastMNN$Histology_Simplified) <- c("ATC", "PTC", "Paratumor/Normal", "PTC", "PTC")
levels(Merged_SO_FastMNN$Histology_Simplified)

# I have chosen 6 colors from the alphabet palette (see below)
#cols <- c("#0075DC", "#FE8F42", "#766C95", "#990000", "#426600", "#808080")
cols <- c("#0075DC", "#FE8F42", "#783FC1", "#990000", "#426600", "#808080")
#cols <- c("#0075DC", "#FF5005", "#783FC1", "#990000", "#426600", "#808080")
savedir = "outputs/Fibroblast_Subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/Clustering"
# Plot 0.1 Resolution
ggsave(file.path(savedir, "UMAP_0.1_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.1", cols = DiscretePalette(n = 3), label = TRUE),
       height = 5, width = 6, dpi = 600)
# Plot 0.2 Resolution
ggsave(file.path(savedir, "UMAP_0.2_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.2", cols = DiscretePalette(n = 6), label = TRUE),
       height = 5, width = 6, dpi = 600)
# Plot 0.2 Resolution with labels turned OFF (for paper figures)
ggsave(file.path(savedir, "UMAP_0.2_Resolution_NOT_Labeled.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.2", cols = cols, label = FALSE),
       height = 5, width = 6, dpi = 600)
# Plot 0.2 Resolution for cells that come from ATC samples ONLY
ggsave(file.path(savedir, "UMAP_0.2_Resolution_NOT_Labeled_ATC.png"),
       DimPlot(Merged_SO_FastMNN %>% subset(Histology_Simplified == "ATC"),
               group.by = "RNA_snn_res.0.2", cols = cols, label = FALSE),
       height = 5, width = 6, dpi = 600)
# Plot 0.2 Resolution for cells that come from PTC samples ONLY
ggsave(file.path(savedir, "UMAP_0.2_Resolution_NOT_Labeled_PTC.png"),
       DimPlot(Merged_SO_FastMNN %>% subset(Histology_Simplified == "PTC"),
               group.by = "RNA_snn_res.0.2", cols = cols, label = FALSE),
       height = 5, width = 6, dpi = 600)
# Plot 0.2 Resolution for cells that come from Paratumor/normal samples ONLY
ggsave(file.path(savedir, "UMAP_0.2_Resolution_NOT_Labeled_Normal.png"),
       DimPlot(Merged_SO_FastMNN %>% subset(Histology_Simplified == "Paratumor/Normal"),
               group.by = "RNA_snn_res.0.2", cols = cols, label = FALSE),
       height = 5, width = 6, dpi = 600)
# Plot higher resolution clustering
ggsave(file.path(savedir, "UMAP_0.3_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.3", cols = DiscretePalette(n = 8), label = TRUE),
       height = 5, width = 6, dpi = 600)
ggsave(file.path(savedir, "UMAP_0.4_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.4", cols = DiscretePalette(n = 10), label = TRUE),
       height = 5, width = 6, dpi = 600)
ggsave(file.path(savedir, "UMAP_0.5_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.5", cols = DiscretePalette(n = 11), label = TRUE),
       height = 5, width = 6, dpi = 600)
ggsave(file.path(savedir, "UMAP_0.6_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.6", cols = DiscretePalette(n = 11), label = TRUE),
       height = 5, width = 6, dpi = 600)
ggsave(file.path(savedir, "UMAP_0.7_Resolution.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.7", cols = DiscretePalette(n = 13), label = TRUE),
       height = 5, width = 6, dpi = 600)


##### GOING WITH 0.2 RESOLUTION #####
### Set Labels of 0.2 resolution clustering
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.2
levels(Merged_SO_FastMNN)
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% RenameIdents("0" = "iPVCAF",
                                                        "1" = "iCAF",
                                                        "2" = "myCAF",
                                                        "3" = "dPVCAF",
                                                        "4" = "APOE+_CAF",
                                                        "5" = "iCAF2")
levels(Merged_SO_FastMNN)
Merged_SO_FastMNN$CAF_Labels <- Idents(Merged_SO_FastMNN)

# Save new annotations
savedir = "outputs/Fibroblast_Subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/DimPlots"
ggsave(file.path(savedir, "CAF_Labels.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "CAF_Labels", cols = DiscretePalette(n = 6), label = FALSE),
       height = 5, width = 5, dpi = 600)


###### PULL OUT myCAF IDs to SUBCLUSTER myCAFs ######
DimPlot(Merged_SO_FastMNN, group.by = "CAF_Labels")
myCAF_IDs <- colnames(Merged_SO_FastMNN %>% subset(CAF_Labels == "myCAF"))
saveRDS(myCAF_IDs, file = "data_in_use/Integrated_Data/24-0821_3000_FastMNN_Fibroblast_Subclustering_myCAF_Labels.RDS")

##### PULL OUT FIBROBLAST SUBCLUSTERING IDs FOR MAPPING BACK ONTO INTEGRATED ATLAS #####
CAF_Labels <- Merged_SO_FastMNN$CAF_Labels
saveRDS(CAF_Labels, file = "data_in_use/Integrated_Data/24-0821_3000_FastMNN_Fibroblast_Subclustering_CAF_Labels.RDS")

##### saveRDS
saveRDS(Merged_SO_FastMNN, file = "~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

##### SESSION INFO #####
sessionInfo()
# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
#
# Matrix products: default
#
#
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8
#
# time zone: America/Chicago
# tzcode source: internal
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] RColorBrewer_1.1-3    SeuratWrappers_0.2.0  SeuratDisk_0.0.0.9020 SeuratObject_4.1.4    Seurat_4.4.0
# [6] lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2
# [11] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.15.0      jsonlite_1.8.7         magrittr_2.0.3         spatstat.utils_3.0-4   farver_2.1.1
# [6] ragg_1.2.6             vctrs_0.6.4            ROCR_1.0-11            spatstat.explore_3.2-5 htmltools_0.5.7
# [11] sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-21     htmlwidgets_1.6.2      ica_1.0-3
# [16] plyr_1.8.9             plotly_4.10.3          zoo_1.8-12             igraph_2.0.3           mime_0.12
# [21] lifecycle_1.0.4        pkgconfig_2.0.3        rsvd_1.0.5             Matrix_1.6-1           R6_2.5.1
# [26] fastmap_1.1.1          fitdistrplus_1.1-11    future_1.33.0          shiny_1.8.0            digest_0.6.33
# [31] colorspace_2.1-0       patchwork_1.2.0        tensor_1.5             irlba_2.3.5.1          textshaping_0.3.7
# [36] labeling_0.4.3         progressr_0.14.0       fansi_1.0.5            spatstat.sparse_3.0-3  timechange_0.2.0
# [41] httr_1.4.7             polyclip_1.10-6        abind_1.4-5            compiler_4.3.1         remotes_2.4.2.1
# [46] bit64_4.0.5            withr_2.5.2            MASS_7.3-60            tools_4.3.1            lmtest_0.9-40
# [51] httpuv_1.6.12          future.apply_1.11.0    goftest_1.2-3          glue_1.6.2             nlme_3.1-162
# [56] promises_1.2.1         grid_4.3.1             Rtsne_0.16             cluster_2.1.4          reshape2_1.4.4
# [61] generics_0.1.3         hdf5r_1.3.8            gtable_0.3.4           spatstat.data_3.0-3    tzdb_0.4.0
# [66] data.table_1.14.8      hms_1.1.3              sp_2.1-1               utf8_1.2.4             spatstat.geom_3.2-7
# [71] RcppAnnoy_0.0.21       ggrepel_0.9.4          RANN_2.6.1             pillar_1.9.0           later_1.3.1
# [76] splines_4.3.1          lattice_0.21-8         survival_3.5-5         bit_4.0.5              deldir_1.0-9
# [81] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3          scattermore_1.2
# [86] matrixStats_1.1.0      stringi_1.8.1          lazyeval_0.2.2         codetools_0.2-19       BiocManager_1.30.22
# [91] cli_3.6.1              uwot_0.1.16            xtable_1.8-4           reticulate_1.34.0      systemfonts_1.0.5
# [96] munsell_0.5.0          Rcpp_1.0.11            globals_0.16.2         spatstat.random_3.2-1  png_0.1-8
# [101] parallel_4.3.1         ellipsis_0.3.2         listenv_0.9.0          viridisLite_0.4.2      scales_1.3.0
# [106] ggridges_0.5.4         leiden_0.4.3.1         crayon_1.5.2           rlang_1.1.2            cowplot_1.1.1
