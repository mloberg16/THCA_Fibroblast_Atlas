### Author: Matthew A. Loberg
### Date: November 4, 2024
### Script: Peds_Samples_Seurat_PreProcessing.R
### Source Script Name: 24-1104_Peds_Samples_Seurat_PreProcessing.R

# Purpose: New pediatric Visium sequencing data just obtained from Vantage
# Here, I will read the data into R studio and begin basic processing of the data + save as .RDS

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Peds samples and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in a list of pediatric sample data directories for processing as a list
data_dir <- c('Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0001_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0002_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0003_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0004_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0005_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0006_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0007_v2',
              'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/count_11133-CP-0008_v2')# Set directory to load from

# Initiate Spatial_SOs as a list
# Will store the seurat objects in this list
Spatial_SOs <- list()

# Read in Seurat objects from data directory
for(i in 1:length(data_dir)){
  Spatial_SOs[[i]] <- Seurat::Load10X_Spatial(data.dir = data_dir[i], slice = "slice1") 
  Spatial_SOs[[i]]$orig.ident <- paste0("Peds0", i)
}

rm(data_dir)

# QC plots
for(i in 1:length(Spatial_SOs)){
  
  # Violin plot (raw counts)
  ViolinPlot <- VlnPlot(Spatial_SOs[[i]], features = "nCount_Spatial", raster = FALSE) + NoLegend() + theme(axis.text = element_text(face = "bold", size = 15))
  ggsave(paste0("outputs/Peds0", i, "/24-1104_preprocessing/QC/24-1104_Raw_Counts_Violin.png"),
         ViolinPlot,
         width = 4, height = 5, dpi = 600, create.dir = TRUE)
  
  # Spatial raw counts plot
  SpatialPlot <- SpatialFeaturePlot(Spatial_SOs[[i]], features = "nCount_Spatial") + theme(legend.position = "right",
                                                                                           legend.text = element_text(face = "bold", size = 15),
                                                                                           legend.title = element_text(face = "bold", size = 15))
  ggsave(paste0("outputs/Peds0", i, "/24-1104_preprocessing/QC/24-1104_Raw_Counts_Spatial.png"),
         SpatialPlot,
         width = 7, height = 5, dpi = 600, create.dir = TRUE)
}

# Cleaning up
rm(ViolinPlot, SpatialPlot)

# Save .RDS objects
for(i in 1:length(Spatial_SOs)){
  saveRDS(Spatial_SOs[[i]], paste0("Data_in_Use/Belcher_Peds_Visium/Processed/Peds0", i, "v2/24-1104_Peds0", i, "_Raw_PreProcessed.rds"))
}

### Chapter 3: Cleaning up 
rm(list = ls())

### Chapter 4: Session Info 
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
#   [1] hdf5r_1.3.8        patchwork_1.2.0    lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
# [7] purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.5.0      tidyverse_2.0.0   
# [13] SeuratObject_4.1.4 Seurat_4.4.0       Matrix_1.6-1       spacexr_2.2.1     
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.7         magrittr_2.0.3         ggbeeswarm_0.7.2      
# [6] spatstat.utils_3.0-4   farver_2.1.1           ragg_1.2.6             vctrs_0.6.4            ROCR_1.0-11           
# [11] spatstat.explore_3.2-5 htmltools_0.5.7        sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-21    
# [16] htmlwidgets_1.6.2      ica_1.0-3              plyr_1.8.9             plotly_4.10.3          zoo_1.8-12            
# [21] igraph_2.0.3           mime_0.12              lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
# [26] R6_2.5.1               fastmap_1.1.1          fitdistrplus_1.1-11    future_1.33.0          shiny_1.8.0           
# [31] digest_0.6.33          colorspace_2.1-0       tensor_1.5             irlba_2.3.5.1          textshaping_0.3.7     
# [36] labeling_0.4.3         progressr_0.14.0       fansi_1.0.5            spatstat.sparse_3.0-3  timechange_0.2.0      
# [41] httr_1.4.7             polyclip_1.10-6        abind_1.4-5            compiler_4.3.1         bit64_4.0.5           
# [46] withr_2.5.2            doParallel_1.0.17      maps_3.4.1.1           MASS_7.3-60            tools_4.3.1           
# [51] vipor_0.4.5            lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.12          future.apply_1.11.0   
# [56] goftest_1.2-3          quadprog_1.5-8         glue_1.6.2             nlme_3.1-162           promises_1.2.1        
# [61] grid_4.3.1             Rtsne_0.16             cluster_2.1.4          reshape2_1.4.4         generics_0.1.3        
# [66] gtable_0.3.4           spatstat.data_3.0-3    tzdb_0.4.0             data.table_1.14.8      hms_1.1.3             
# [71] sp_2.1-1               utf8_1.2.4             spatstat.geom_3.2-7    RcppAnnoy_0.0.21       ggrepel_0.9.4         
# [76] RANN_2.6.1             foreach_1.5.2          pillar_1.9.0           pals_1.8               later_1.3.1           
# [81] splines_4.3.1          lattice_0.21-8         bit_4.0.5              survival_3.5-5         deldir_1.0-9          
# [86] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3          scattermore_1.2       
# [91] matrixStats_1.1.0      stringi_1.8.1          lazyeval_0.2.2         codetools_0.2-19       cli_3.6.1             
# [96] uwot_0.1.16            systemfonts_1.0.5      xtable_1.8-4           reticulate_1.34.0      munsell_0.5.0         
# [101] dichromat_2.0-0.1      Rcpp_1.0.11            globals_0.16.2         spatstat.random_3.2-1  mapproj_1.2.11        
# [106] png_0.1-8              ggrastr_1.0.2          parallel_4.3.1         ellipsis_0.3.2         listenv_0.9.0         
# [111] viridisLite_0.4.2      scales_1.3.0           ggridges_0.5.4         leiden_0.4.3.1         rlang_1.1.2           
# [116] cowplot_1.1.1  
