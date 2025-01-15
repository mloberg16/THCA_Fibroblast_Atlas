### Author: Matthew Aaron Loberg
### Date: November 4, 2024
### Script: 24-1104_FastMNN_3000_RCTD_Sparse_Matrix_Cell_Labels.R

###### Goal: #####
# Pull out sparse matrix of cell counts + cell type labels with fibroblast subclustering included and pEMT PTC designation included
# This is for RCTD deconvolution of spatial data
# Will EXCLUDE iCAF2 from these counts for the purpose of deconvolution

###### Reference: #####
# Used 24-0212_FastMNN.R as a reference for doing this
# In 24-0212_FastMNN.R, I had a portion of the script where I pulled out sparse counts for RCTD deconvolution of spatial data

#### 24-1104 update
# Adding in PTC pEMT designation as well
# This will help me decode PTC

##### Load packages ####
library(Seurat)

##### Load Data #####
Merged_SO_FastMNN <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_AFTER_FastMNN_3000.RDS")

##### Extract sparse matrix and cell IDs for RCTD deconvolution analysis of spatial data with iCAF2 excluded #####
# Exclude iCAF2
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% subset(Broad_Labels_CAFs_pEMT_Included != "iCAF2")
table(Merged_SO_FastMNN$Broad_Labels_CAFs_pEMT_Included)

##### Extract sparse matrix and cell IDs for RCTD deconvolution analysis of spatial data #####
# extract as a sparse matrix
raw_count_sparse_matrix <- Merged_SO_FastMNN@assays$RNA@counts

# save the sparse matrix
saveRDS(raw_count_sparse_matrix, "~/24-1104_FastMNN_scRNA_Sparse_RNA_Counts_iCAF2_Excluded.RDS")
rm(raw_count_sparse_matrix)

# set the Idents of Merged_SO_FastMNN to "Broad_Labels_CAFS_Included" variable
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$Broad_Labels_CAFs_pEMT_Included
head(Idents(Merged_SO_FastMNN))
levels(Idents(Merged_SO_FastMNN))

cell_types <- Idents(Merged_SO_FastMNN)
levels(cell_types)
saveRDS(cell_types, "~/24-1104_FastMNN_scRNA_cell_type_labels_iCAF2_Excluded.RDS")

##### Cleaning up #####
rm(list = ls())

######## Session Info #############
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
#   [1] RColorBrewer_1.1-3   lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4
# [6] purrr_1.0.2          readr_2.1.4          tidyr_1.3.0          tibble_3.2.1         ggplot2_3.5.0
# [11] tidyverse_2.0.0      SeuratWrappers_0.2.0 SeuratObject_4.1.4   Seurat_4.4.0
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
# [46] withr_2.5.2            MASS_7.3-60            tools_4.3.1            lmtest_0.9-40          httpuv_1.6.12
# [51] future.apply_1.11.0    goftest_1.2-3          glue_1.6.2             nlme_3.1-162           promises_1.2.1
# [56] grid_4.3.1             Rtsne_0.16             cluster_2.1.4          reshape2_1.4.4         generics_0.1.3
# [61] gtable_0.3.4           spatstat.data_3.0-3    tzdb_0.4.0             data.table_1.14.8      hms_1.1.3
# [66] sp_2.1-1               utf8_1.2.4             spatstat.geom_3.2-7    RcppAnnoy_0.0.21       ggrepel_0.9.4
# [71] RANN_2.6.1             pillar_1.9.0           later_1.3.1            splines_4.3.1          lattice_0.21-8
# [76] survival_3.5-5         deldir_1.0-9           tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2
# [81] gridExtra_2.3          scattermore_1.2        matrixStats_1.1.0      stringi_1.8.1          lazyeval_0.2.2
# [86] codetools_0.2-19       BiocManager_1.30.22    cli_3.6.1              uwot_0.1.16            xtable_1.8-4
# [91] reticulate_1.34.0      systemfonts_1.0.5      munsell_0.5.0          Rcpp_1.0.11            globals_0.16.2
# [96] spatstat.random_3.2-1  png_0.1-8              parallel_4.3.1         ellipsis_0.3.2         listenv_0.9.0
# [101] viridisLite_0.4.2      scales_1.3.0           ggridges_0.5.4         leiden_0.4.3.1         rlang_1.1.2
# [106] cowplot_1.1.1


