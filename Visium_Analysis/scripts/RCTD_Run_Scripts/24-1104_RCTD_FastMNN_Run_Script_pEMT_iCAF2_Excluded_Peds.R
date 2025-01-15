### Author: Matthew Aaron Loberg
### Date: November 4, 2024
### Script: 24-1104_RCTD_FastMNN_Run_Script.R

#### Goal: ####
# Run The RCTD_Processing function for ALL Thy objects
# Peds objects included

#### 24-1001 Updates ####
# Including pt.size.factor variable that is a list of 20 values; different value for each sample that I previously determined
# Calling updated function scripts for processing and output plots
# Adding in a function call to "24-1001_SpatialFeaturePlotBlend.R", which will allow me to plot multiple features at the same time

#### 24-1104 Updates ####
# I have segregated my PTC to include a "pEMT" group that is enriched for the pEMT gene signature
# I will re-run RCTD analysis with the NEW pEMT object
# Include plots with the NEW pEMT object for each sample

#### Load required packages ####
library(spacexr)
library(Matrix)
library(Seurat)
library(tidyverse)
library(patchwork) # required to wrap plots in SpatialFeaturePlotBlend

#### Load required function scripts ####
source("scripts/RCTD_Functions/24-1001_RCTD_Analysis_Processing_Function.R") # No update from prior
source("scripts/RCTD_Functions/RCTD_Editable_Plotting_Functions.R") # No update from prior
source("scripts/RCTD_Functions/24-1104_iCAF2_Excluded_RCTD_Output_Plots.R") # Updated for new objects
source('scripts/Function_Scripts/24-1001_SpatialFeaturePlotBlend.R') # No update from prior

#### Load Reference Data ####
reference <- readRDS(file = "~/24-1104_FastMNN_RCTD_Reference_iCAF2_Excluded.RDS")

#### Load in readdirs ####
readdirs <- list(
                 #"Data_in_Use/2021_JHU_Data/Thy1_Processed/22-0913_Thy1_Raw_PreProcessed.rds",
                 #"Data_in_Use/2021_JHU_Data/Thy2_Processed/22-1129_Thy2_Raw_PreProcessed.rds",
                 #"Data_in_Use/2021_JHU_Data/Thy3_Processed/22-1129_Thy3_Raw_PreProcessed.rds",
                 #"Data_in_Use/2021_JHU_Data/Thy4_Processed/22-1129_Thy4_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy5_Processed/22-0825_Thy5_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy6_Processed/22-0915_Thy6_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy7_Processed/22-1122_Thy7_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy8_Processed/22-1122_Thy8_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy9_ManualAlign_Processed/22-1003_Thy9_ManualAlign_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy10_Processed/22-0915_Thy10_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy11_Processed/22-0905_Thy11_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy12_Processed/22-0915_Thy12_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy13_Processed/22-0915_Thy13_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy14_Processed/22-1129_Thy14_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy15_Processed/22-1129_Thy15_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy16_Processed/22-1129_Thy16_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy17_Processed/22-1129_Thy17_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy18_Processed/22-1129_Thy18_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy19_Processed/22-1129_Thy19_Raw_PreProcessed.rds",
                 #"Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy20_Processed/22-1129_Thy20_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds01v2/24-1104_Peds01_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds02v2/24-1104_Peds02_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds03v2/24-1104_Peds03_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds04v2/24-1104_Peds04_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds05v2/24-1104_Peds05_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds06v2/24-1104_Peds06_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds07v2/24-1104_Peds07_Raw_PreProcessed.rds",
                 "Data_in_Use/Belcher_Peds_Visium/Processed/Peds08v2/24-1104_Peds08_Raw_PreProcessed.rds"
                 )

# Create a list length 20 of the pt.size.factors to use for the objects
# This is a 24-1001 addition after the pt.size.factors being WAY off
pt.size.factor <- list(1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65)

# Create save directory
dir.create("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded")

##### Loop through all the readdirs and do the analysis on all of the readdirs ####
for(i in 1:length(readdirs)){
  # Create resultsdir
  resultsdir <- paste0("outputs/Peds0", i, "_RCTD/24-1104_iCAF2_Excluded/")
  dir.create(resultsdir)
  
  # Next, load spatial data as "Thy_Spatial_SO"
  Thy_Spatial_SO <- readRDS(file = readdirs[[i]])
  
  # Now generate Thy_RCTD_Obj from Thy_Spatial_SO
  Thy_RCTD_Obj <- RCTD_Processing(reference = reference,
                                  Raw_SO = Thy_Spatial_SO,
                                  resultsdir = resultsdir)
  
  # Save the RCTD object
  saveRDS(Thy_RCTD_Obj, file = paste0( "~/24-1104_Peds0", i, "_FastMNN_RCTD_Obj_iCAF2_Excluded.RDS")) # Saving to the "documents" folder of the lab PC to save space on OneDrive; will back up on s drive, SharePoint
  
  # Make plots
  Thy_Spatial_SO <- RCTD_Plotting(RCTD_Obj = Thy_RCTD_Obj, SO = Thy_Spatial_SO, resultsdir = resultsdir, pt.size.factor = pt.size.factor[[i]], alpha = c(0.01, 1.0))
  
  # Save Thy_Spatial_SO
  saveRDS(Thy_Spatial_SO, file = paste0("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds0", i, ".RDS"))
  
}

##### Cleaning Up #####
rm(list = ls())

# Adding additional plots
Thy_Spatial_SOs <- list()
placeholders <- c(1:20, 1:8)
for(i in 1:length(placeholders)){
  if(i <= 20){
    Thy_Spatial_SOs[[i]] <- readRDS(file = paste0("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy", i, ".RDS"))
  }
  else if(i > 20){
    Thy_Spatial_SOs[[i]] <- readRDS(file = paste0("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds0", placeholders[i], ".RDS"))
  }
}

pt.size.factor <- c(3, 3, 3, 3, 1.6, 1.7, 1.7, 1.7, 1.7, 1.6, 
                       1.7, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.9, 1.6,
                       1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65)
for(i in 1:length(Thy_Spatial_SOs)){
  if(i <= 20){
    resultsdir <- paste0("outputs/Thy", i, "_RCTD/24-1104_iCAF2_Excluded/")
  }
  else if(i > 20){
    resultsdir <- paste0("outputs/Peds0", placeholders[i], "_RCTD/24-1104_iCAF2_Excluded/")
  }
  Spatial_Feature <- SpatialFeaturePlot(Thy_Spatial_SOs[[i]], features = c("dPVCAF_RCTD"), pt.size.factor = pt.size.factor[i], alpha = c(1.0))
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/dPVCAF_RCTD_Spatial_Feature_Plot_alpha_1.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
}

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
#   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.7         magrittr_2.0.3         spatstat.utils_3.0-4  
# [6] farver_2.1.1           ragg_1.2.6             vctrs_0.6.4            ROCR_1.0-11            spatstat.explore_3.2-5
# [11] htmltools_0.5.7        sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-21     htmlwidgets_1.6.2     
# [16] ica_1.0-3              plyr_1.8.9             plotly_4.10.3          zoo_1.8-12             igraph_2.0.3          
# [21] mime_0.12              lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3        R6_2.5.1              
# [26] fastmap_1.1.1          fitdistrplus_1.1-11    future_1.33.0          shiny_1.8.0            digest_0.6.33         
# [31] colorspace_2.1-0       patchwork_1.2.0        tensor_1.5             irlba_2.3.5.1          textshaping_0.3.7     
# [36] labeling_0.4.3         progressr_0.14.0       fansi_1.0.5            spatstat.sparse_3.0-3  timechange_0.2.0      
# [41] httr_1.4.7             polyclip_1.10-6        abind_1.4-5            compiler_4.3.1         withr_2.5.2           
# [46] doParallel_1.0.17      maps_3.4.1.1           MASS_7.3-60            tools_4.3.1            lmtest_0.9-40         
# [51] httpuv_1.6.12          future.apply_1.11.0    goftest_1.2-3          quadprog_1.5-8         glue_1.6.2            
# [56] nlme_3.1-162           promises_1.2.1         grid_4.3.1             Rtsne_0.16             cluster_2.1.4         
# [61] reshape2_1.4.4         generics_0.1.3         gtable_0.3.4           spatstat.data_3.0-3    tzdb_0.4.0            
# [66] data.table_1.14.8      hms_1.1.3              sp_2.1-1               utf8_1.2.4             spatstat.geom_3.2-7   
# [71] RcppAnnoy_0.0.21       ggrepel_0.9.4          RANN_2.6.1             foreach_1.5.2          pillar_1.9.0          
# [76] pals_1.8               later_1.3.1            splines_4.3.1          lattice_0.21-8         survival_3.5-5        
# [81] deldir_1.0-9           tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3         
# [86] scattermore_1.2        matrixStats_1.1.0      stringi_1.8.1          lazyeval_0.2.2         codetools_0.2-19      
# [91] cli_3.6.1              uwot_0.1.16            systemfonts_1.0.5      xtable_1.8-4           reticulate_1.34.0     
# [96] munsell_0.5.0          dichromat_2.0-0.1      Rcpp_1.0.11            globals_0.16.2         spatstat.random_3.2-1 
# [101] mapproj_1.2.11         png_0.1-8              parallel_4.3.1         ellipsis_0.3.2         listenv_0.9.0         
# [106] viridisLite_0.4.2      scales_1.3.0           ggridges_0.5.4         leiden_0.4.3.1         rlang_1.1.2           
# [111] cowplot_1.1.1  
