#### Author: Matthew Aaron Loberg
#### Date: September 17, 2024
#### Script: "FastMNN_3000_Plot_Marker_Genes_TEMPLATE.R"

# This is adapted from 24-0819_FastMNN.R

# Goal: Script just for plotting feature plots amongst the large UMAP
# Any features of interest
# TEMPLATE is for GitHub upload; replace "INSERT_GENE_NAME" with desired genes to plot

##### Load Packages #####
library(Seurat)
library(SeuratWrappers)
library(tidyverse) # for ggsave
library(RColorBrewer) # for custom plot coloring

##### Load Data #####
Merged_SO_FastMNN <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_After_FastMNN_3000.RDS")

##### FeaturePlots #####
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/FeaturePlots"
dir.create(savedir)
gene <- list("INSERT_GENE_1_NAME",
             "INSERT_GENE_2_NAME",
             "...",
             "INSERT_GENE_N_NAME")

for(i in 1:length(gene)){
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot.png")), # Standard blue plot
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]])),
         height = 5, width = 5, dpi = 600)
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot_Gradient.png")), # Blue -> Red gradient plot
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]]))  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
         height = 5, width = 5, dpi = 600)
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot_Gradient_Red.png")), # Red plot with q90 cutoff
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]]), cols = c("gray", "red"), max.cutoff = 'q90', order = FALSE),
         height = 5, width = 5, dpi = 600)
}

##### Cleaning Up #####
rm(list = ls())

######## Session Info #############
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] beepr_1.3                   SeuratWrappers_0.2.0        RColorBrewer_1.1-3          scDblFinder_1.14.0
# [5] SingleCellExperiment_1.22.0 SingleR_2.2.0               celldex_1.10.1              SummarizedExperiment_1.30.2
# [9] Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1
# [13] S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.1.0
# [17] SoupX_1.6.2                 SeuratDisk_0.0.0.9020       SeuratObject_4.1.4          Seurat_4.4.0
# [21] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4
# [25] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1
# [29] ggplot2_3.5.0               tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.0-3         bitops_1.0-7                  httr_1.4.7                    tools_4.3.1
# [5] sctransform_0.4.1             ResidualMatrix_1.10.0         utf8_1.2.4                    R6_2.5.1
# [9] lazyeval_0.2.2                uwot_0.1.16                   withr_2.5.2                   sp_2.1-1
# [13] gridExtra_2.3                 progressr_0.14.0              cli_3.6.1                     textshaping_0.3.7
# [17] spatstat.explore_3.2-5        labeling_0.4.3                spatstat.data_3.0-3           ggridges_0.5.4
# [21] pbapply_1.7-2                 Rsamtools_2.16.0              systemfonts_1.0.5             scater_1.28.0
# [25] parallelly_1.36.0             limma_3.56.2                  rstudioapi_0.15.0             RSQLite_2.3.3
# [29] generics_0.1.3                BiocIO_1.10.0                 ica_1.0-3                     spatstat.random_3.2-1
# [33] Matrix_1.6-1                  ggbeeswarm_0.7.2              fansi_1.0.5                   abind_1.4-5
# [37] lifecycle_1.0.4               yaml_2.3.7                    edgeR_3.42.4                  BiocFileCache_2.8.0
# [41] Rtsne_0.16                    grid_4.3.1                    blob_1.2.4                    promises_1.2.1
# [45] dqrng_0.3.1                   ExperimentHub_2.8.1           crayon_1.5.2                  miniUI_0.1.1.1
# [49] lattice_0.21-8                beachmat_2.16.0               cowplot_1.1.1                 KEGGREST_1.40.1
# [53] pillar_1.9.0                  metapod_1.8.0                 rjson_0.2.21                  xgboost_1.7.5.1
# [57] future.apply_1.11.0           codetools_0.2-19              leiden_0.4.3.1                glue_1.6.2
# [61] data.table_1.14.8             remotes_2.4.2.1               vctrs_0.6.4                   png_0.1-8
# [65] gtable_0.3.4                  cachem_1.0.8                  S4Arrays_1.0.6                mime_0.12
# [69] survival_3.5-5                audio_0.1-11                  statmod_1.5.0                 bluster_1.10.0
# [73] interactiveDisplayBase_1.38.0 ellipsis_0.3.2                fitdistrplus_1.1-11           ROCR_1.0-11
# [77] nlme_3.1-162                  bit64_4.0.5                   filelock_1.0.2                RcppAnnoy_0.0.21
# [81] irlba_2.3.5.1                 vipor_0.4.5                   KernSmooth_2.23-21            colorspace_2.1-0
# [85] DBI_1.1.3                     tidyselect_1.2.0              bit_4.0.5                     compiler_4.3.1
# [89] curl_5.1.0                    BiocNeighbors_1.18.0          hdf5r_1.3.8                   DelayedArray_0.26.7
# [93] plotly_4.10.3                 rtracklayer_1.60.1            scales_1.3.0                  lmtest_0.9-40
# [97] rappdirs_0.3.3                digest_0.6.33                 goftest_1.2-3                 spatstat.utils_3.0-4
# [101] XVector_0.40.0                htmltools_0.5.7               pkgconfig_2.0.3               sparseMatrixStats_1.12.2
# [105] dbplyr_2.3.4                  fastmap_1.1.1                 rlang_1.1.2                   htmlwidgets_1.6.2
# [109] shiny_1.8.0                   DelayedMatrixStats_1.22.6     farver_2.1.1                  zoo_1.8-12
# [113] jsonlite_1.8.7                BiocParallel_1.34.2           BiocSingular_1.16.0           RCurl_1.98-1.13
# [117] magrittr_2.0.3                scuttle_1.10.3                GenomeInfoDbData_1.2.10       patchwork_1.2.0
# [121] munsell_0.5.0                 Rcpp_1.0.11                   viridis_0.6.4                 reticulate_1.34.0
# [125] stringi_1.8.1                 zlibbioc_1.46.0               MASS_7.3-60                   AnnotationHub_3.8.0
# [129] plyr_1.8.9                    parallel_4.3.1                listenv_0.9.0                 ggrepel_0.9.4
# [133] deldir_1.0-9                  Biostrings_2.68.1             splines_4.3.1                 tensor_1.5
# [137] hms_1.1.3                     locfit_1.5-9.8                igraph_2.0.3                  spatstat.geom_3.2-7
# [141] reshape2_1.4.4                ScaledMatrix_1.8.1            BiocVersion_3.17.1            XML_3.99-0.15
# [145] scran_1.28.2                  BiocManager_1.30.22           batchelor_1.16.0              tzdb_0.4.0
# [149] httpuv_1.6.12                 RANN_2.6.1                    polyclip_1.10-6               future_1.33.0
# [153] scattermore_1.2               rsvd_1.0.5                    xtable_1.8-4                  restfulr_0.0.15
# [157] later_1.3.1                   viridisLite_0.4.2             ragg_1.2.6                    memoise_2.0.1
# [161] beeswarm_0.4.0                AnnotationDbi_1.62.2          GenomicAlignments_1.36.0      cluster_2.1.4
# [165] timechange_0.2.0              globals_0.16.2
