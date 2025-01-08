### Author: Matthew Aaron Loberg
### Date: November 4, 2024
### Script: 24-1104_RCTD_Reference_Builder.R

##### Goal: #####
# Build RCTD Reference from newly created annotations From ~420k integrated cells with FastMNN
# Include CAF subclustering labels and broad cell type IDs
# This version will also include a split of PTC into PTC and pEMT-PTC groups

#### Reference/Updates from previous #####
# This is an update of 24-0930, 24-0212, and 23-1207 RCTD scripts
# I am running these with an updated atlas from the prior reference builder scripts
# This is with the 24-0819 atlas integration
# 24-0930 was run with the most recent atlas but without stromal subclusters
# In this 24-1104 update, I am including stromal subclusters (myCAF, iCAF, APOE+ PVL, pericyte, vSMC)
# I am running this on the ATLAS EXCLUDING iCAF2 (largely from one sample)
# I think that the exclusion of iCAF2 will give the best marker genes for iCAF

#### 24-1104 Update ####
# Building a NEW reference with pEMT included

#### Load required packages ####
library(spacexr)
library(Matrix)
library(Seurat)
library(tidyverse)

#### Load data for reference ####
reference_count_data <- readRDS(file = "~/24-1104_FastMNN_scRNA_Sparse_RNA_Counts_iCAF2_Excluded.RDS")
reference_cell_types <- readRDS(file = "~/24-1104_FastMNN_scRNA_cell_type_labels_iCAF2_Excluded.RDS")

# problems with the factor levels of reference_cell_types containing "NK/T"
# Here I will re-code the factor levels
levels(reference_cell_types)
reference_cell_types <- recode(reference_cell_types,
                               "NK/T" = "NKT")
levels(reference_cell_types)

#####Generate reference from count data + cell types #####
# Original ran into an error due to "NK/T" -> error went away when I changed it to "NKT"
reference <- spacexr::Reference(reference_count_data, reference_cell_types) 
# Note that NKT cells were downsampled from 167,528 to 10,000 cells

###### Save Reference #####
saveRDS(reference, file = "~/24-1104_FastMNN_RCTD_Reference_iCAF2_Excluded.RDS")

###### Cleaning up #####
rm(list = ls())

##### Session Info #####
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
#   [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.4       
# [7] tidyr_1.3.0        tibble_3.2.1       ggplot2_3.5.0      tidyverse_2.0.0    SeuratObject_4.1.4 Seurat_4.4.0      
# [13] Matrix_1.6-1       spacexr_2.2.1     
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.2            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.21       spatstat.geom_3.2-7    matrixStats_1.1.0      ggridges_0.5.4         compiler_4.3.1        
# [11] png_0.1-8              vctrs_0.6.4            reshape2_1.4.4         pkgconfig_2.0.3        fastmap_1.1.1         
# [16] ellipsis_0.3.2         utf8_1.2.4             promises_1.2.1         tzdb_0.4.0             jsonlite_1.8.7        
# [21] goftest_1.2-3          later_1.3.1            spatstat.utils_3.0-4   irlba_2.3.5.1          parallel_4.3.1        
# [26] cluster_2.1.4          R6_2.5.1               ica_1.0-3              stringi_1.8.1          RColorBrewer_1.1-3    
# [31] spatstat.data_3.0-3    reticulate_1.34.0      parallelly_1.36.0      lmtest_0.9-40          scattermore_1.2       
# [36] Rcpp_1.0.11            iterators_1.0.14       tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
# [41] sctransform_0.4.1      timechange_0.2.0       httpuv_1.6.12          splines_4.3.1          igraph_2.0.3          
# [46] tidyselect_1.2.0       rstudioapi_0.15.0      abind_1.4-5            spatstat.random_3.2-1  doParallel_1.0.17     
# [51] codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.2-5 listenv_0.9.0          lattice_0.21-8        
# [56] plyr_1.8.9             withr_2.5.2            shiny_1.8.0            ROCR_1.0-11            Rtsne_0.16            
# [61] future_1.33.0          survival_3.5-5         polyclip_1.10-6        fitdistrplus_1.1-11    pillar_1.9.0          
# [66] KernSmooth_2.23-21     foreach_1.5.2          plotly_4.10.3          generics_0.1.3         sp_2.1-1              
# [71] hms_1.1.3              munsell_0.5.0          scales_1.3.0           globals_0.16.2         xtable_1.8-4          
# [76] glue_1.6.2             lazyeval_0.2.2         tools_4.3.1            data.table_1.14.8      RANN_2.6.1            
# [81] leiden_0.4.3.1         cowplot_1.1.1          grid_4.3.1             colorspace_2.1-0       nlme_3.1-162          
# [86] patchwork_1.2.0        cli_3.6.1              spatstat.sparse_3.0-3  fansi_1.0.5            viridisLite_0.4.2     
# [91] uwot_0.1.16            gtable_0.3.4           digest_0.6.33          progressr_0.14.0       ggrepel_0.9.4         
# [96] htmlwidgets_1.6.2      htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7             mime_0.12             
# [101] MASS_7.3-60  
