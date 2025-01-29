### Author: Matthew Aaron Loberg
### Date: 24-1118
### Script: BayesSpace_Enhanced_Feature_Plots.R

### Goal: 
# I will read in SCE and SCE Enhanced objects for all Adult Visium images
# From there I will look at enhanced features of interest

### 24-1118 update
# Adding in custom high color for RGS5 as "#0075DC" and POSTN as "#783FC1" to match their respective fibroblast clusters
# Add in markers color with argument "high" in featurePlot function

### Load required packages:
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

### Read in objects
# Set read directory list for BayesSpace sce objects (1 for each adult thyroid Visium sample)
sce_readdirs <- list("Data_in_Use/Processed_Outputs/Thy1_Processed/Thy1_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy2_Processed/Thy2_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy3_Processed/Thy3_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy4_Processed/Thy4_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy5_Processed/Thy5_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy6_Processed/Thy6_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy7_Processed/Thy7_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy8_Processed/Thy8_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy9_ManualAlign_Processed/Thy9_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy10_Processed/Thy10_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy11_Processed/Thy11_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy12_Processed/Thy12_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy13_Processed/Thy13_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy14_Processed/Thy14_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy15_Processed/Thy15_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy16_Processed/Thy16_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy17_Processed/Thy17_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy18_Processed/Thy18_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy19_Processed/Thy19_BayesSpace.RDS",
                     "Data_in_Use/Processed_Outputs/Thy20_Processed/Thy20_BayesSpace.RDS")

# Set read directory list for BayesSpace sce.enhanced objects (1 for each adult thyroid Visium sample)
sce.enhanced_readdirs <- list("Data_in_Use/Processed_Outputs/Thy1_Processed/Thy1_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy2_Processed/Thy2_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy3_Processed/Thy3_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy4_Processed/Thy4_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy5_Processed/Thy5_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy6_Processed/Thy6_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy7_Processed/Thy7_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy8_Processed/Thy8_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy9_ManualAlign_Processed/Thy9_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy10_Processed/Thy10_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy11_Processed/Thy11_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy12_Processed/Thy12_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy13_Processed/Thy13_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy14_Processed/Thy14_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy15_Processed/Thy15_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy16_Processed/Thy16_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy17_Processed/Thy17_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy18_Processed/Thy18_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy19_Processed/Thy19_BayesSpace_Enhanced.RDS",
                              "Data_in_Use/Processed_Outputs/Thy20_Processed/Thy20_BayesSpace_Enhanced.RDS")

# Set markers to plot and colors to plot them in 
markers <- c("RGS5", "POSTN") # INSERT MARKERS OF INTEREST
colors <- c("#0075DC", "#783FC1") # INSERT DESIRED COLORS FOR MARKERS OF INTEREST

# Set savedir for output plots
savedir <- "outputs/BayesSpace_Feature_Plots/Thy"

# For loop to make plots for each sample
for(i in 1:length(sce_readdirs)){

  # Generate save directory
  dir.create(paste0(savedir,i))

  # Read in sce and sce.enhanced object
  sce <- readRDS(sce_readdirs[[i]])
  sce.enhanced <- readRDS(sce.enhanced_readdirs[[i]])

  # Enhance features of interest
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                       feature_names = markers,
                                       nrounds = 0)

  # Plot each marker
  for(n in 1:length(markers)){
    tryCatch({
      ggsave(paste0(paste0(paste0(savedir,i), paste0("/", markers[n])), "_Enhanced_FeaturePlot.png"),
             featurePlot(sce.enhanced, markers[n], high = colors[n]),
             width = 5, height = 5, dpi = 600)
      
    }, error = function(e){
      
      print("Error plotting gene:")
      print(paste0(sce_readdirs[i], ": ", markers[n]))  
      print(e)
      
    })
  }
}

### Cleaning up
rm(list = ls())

### SESSION INFO
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BayesSpace_1.10.1           ggplot2_3.5.0               SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2
# [5] Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1             
# [9] S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.1.0          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.1.3                 bitops_1.0-7              gridExtra_2.3             sandwich_3.0-2            DirichletReg_0.7-1       
# [6] rlang_1.1.2               magrittr_2.0.3            scater_1.28.0             compiler_4.3.1            RSQLite_2.3.3            
# [11] DelayedMatrixStats_1.22.6 vctrs_0.6.4               pkgconfig_2.0.3           crayon_1.5.2              fastmap_1.1.1            
# [16] dbplyr_2.3.4              XVector_0.40.0            scuttle_1.10.3            utf8_1.2.4                ggbeeswarm_0.7.2         
# [21] miscTools_0.6-28          purrr_1.0.2               bit_4.0.5                 bluster_1.10.0            zlibbioc_1.46.0          
# [26] cachem_1.0.8              beachmat_2.16.0           jsonlite_1.8.7            blob_1.2.4                rhdf5filters_1.12.1      
# [31] DelayedArray_0.26.7       Rhdf5lib_1.22.1           BiocParallel_1.34.2       irlba_2.3.5.1             parallel_4.3.1           
# [36] cluster_2.1.4             R6_2.5.1                  limma_3.56.2              xgboost_1.7.5.1           Rcpp_1.0.11              
# [41] assertthat_0.2.1          zoo_1.8-12                Matrix_1.6-1              igraph_2.0.3              tidyselect_1.2.0         
# [46] rstudioapi_0.15.0         abind_1.4-5               viridis_0.6.4             maxLik_1.5-2              codetools_0.2-19         
# [51] curl_5.1.0                lattice_0.21-8            tibble_3.2.1              withr_2.5.2               coda_0.19-4              
# [56] BiocFileCache_2.8.0       mclust_6.0.1              pillar_1.9.0              filelock_1.0.2            generics_0.1.3           
# [61] RCurl_1.98-1.13           sparseMatrixStats_1.12.2  munsell_0.5.0             scales_1.3.0              glue_1.6.2               
# [66] metapod_1.8.0             tools_4.3.1               data.table_1.14.8         BiocNeighbors_1.18.0      ScaledMatrix_1.8.1       
# [71] locfit_1.5-9.8            scran_1.28.2              rhdf5_2.44.0              grid_4.3.1                edgeR_3.42.4             
# [76] colorspace_2.1-0          GenomeInfoDbData_1.2.10   beeswarm_0.4.0            BiocSingular_1.16.0       vipor_0.4.5              
# [81] Formula_1.2-5             cli_3.6.1                 rsvd_1.0.5                fansi_1.0.5               S4Arrays_1.0.6           
# [86] viridisLite_0.4.2         dplyr_1.1.4               gtable_0.3.4              digest_0.6.33             ggrepel_0.9.4            
# [91] dqrng_0.3.1               memoise_2.0.1             lifecycle_1.0.4           httr_1.4.7                statmod_1.5.0            
# [96] bit64_4.0.5  
