# Author: Matthew Aaron Loberg
# Date: November 30, 2024
# Script: 24-1130_iCAF2_Broad_Lu_ATC09_pEMT_RCTD_Output_Plots_V2.R

### Script info 
# Generates RCTD output plots for the iCAF2_Excluded object w/ broad labels outside of CAF subclusters (myCAF; iCAF; APOE_CAF; dPVCAF; iPVCAF)

RCTD_Plotting <- function(RCTD_Obj, SO, resultsdir, pt.size.factor = 3, alpha = 1){

  results <- RCTD_Obj@results
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- RCTD_Obj@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- RCTD_Obj@spatialRNA
  dir.create(resultsdir) # Create results directory
  dir.create(paste0(resultsdir, "RCTD_SpatialFeature")) # Prevents the need to create.dir while the code runs
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  plot_weights_unthreshold_edited(cell_type_names, spatialRNA, resultsdir, norm_weights)
  #plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
  #plot_all_cell_types is not working -> I am not sure why, tbh -_-

  # Pull out NKT from column 1 of norm_weights
  SO$NKT_RCTD <- norm_weights[,1]
  SO$NKT_RCTD[is.na(SO$NKT_RCTD)] <- 0
  # Save PNG of NKT RCTD deconvolution
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("NKT_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/NKT_RCTD_Spatial_Feature_Plot.png"),
       Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out B_Cell from column 2 of norm_weights
  SO$B_Cell_RCTD <- norm_weights[,2]
  SO$B_Cell_RCTD[is.na(SO$B_Cell_RCTD)] <- 0
  # Save PNG of B_Cell
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("B_Cell_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/B_Cell_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Myeloid from column 3 of norm_weights
  SO$Myeloid_RCTD <- norm_weights[,3]
  SO$Myeloid_RCTD[is.na(SO$Myeloid_RCTD)] <- 0
  # Save PNG of Myeloid
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("Myeloid_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/Myeloid_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Plasma from column 4 of norm_weights
  SO$Plasma_RCTD <- norm_weights[,4]
  SO$Plasma_RCTD[is.na(SO$Plasma_RCTD)] <- 0
  # Save PNG of Plasma
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("Plasma_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/Plasma_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Endothelial from column 5 of norm_weights
  SO$Endothelial_RCTD <- norm_weights[,5]
  SO$Endothelial_RCTD[is.na(SO$Endothelial_RCTD)] <- 0
  # Save PNG of Endothelial
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("Endothelial_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/Endothelial_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Thyrocyte from column 6 of norm_weights
  SO$Thyrocyte_RCTD <- norm_weights[,6]
  SO$Thyrocyte_RCTD[is.na(SO$Thyrocyte_RCTD)] <- 0
  # Save PNG of Thyrocyte
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("Thyrocyte_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/Thyrocyte_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out pDC from column 7 of norm_weights
  SO$pDC_RCTD <- norm_weights[,7]
  SO$pDC_RCTD[is.na(SO$pDC_RCTD)] <- 0
  # Save PNG of pDC
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("pDC_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pDC_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out PTC from column 8 of norm_weights
  SO$PTC_RCTD <- norm_weights[,8]
  SO$PTC_RCTD[is.na(SO$PTC_RCTD)] <- 0
  # Save PNG of PTC
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("PTC_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/PTC_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out iPVCAF from column 9 of norm_weights
  SO$iPVCAF_RCTD <- norm_weights[,9]
  SO$iPVCAF_RCTD[is.na(SO$iPVCAF_RCTD)] <- 0
  # Save PNG of iPVCAF
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("iPVCAF_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/iPVCAF_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out iCAF from column 10 of norm_weights
  SO$iCAF_RCTD <- norm_weights[,10]
  SO$iCAF_RCTD[is.na(SO$iCAF_RCTD)] <- 0
  # Save PNG of CDK6_ATC
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("iCAF_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/iCAF_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out myCAF from column 11 of norm_weights
  SO$myCAF_RCTD <- norm_weights[,11]
  SO$myCAF_RCTD[is.na(SO$myCAF_RCTD)] <- 0
  # Save PNG of myCAF
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("myCAF_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/myCAF_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out dPVCAF from column 12 of norm_weights
  SO$dPVCAF_RCTD <- norm_weights[,12]
  SO$dPVCAF_RCTD[is.na(SO$dPVCAF_RCTD)] <- 0
  # Save PNG of dPVCAF
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("dPVCAF_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/dPVCAF_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out APOE_CAF from column 13 of norm_weights
  SO$APOE_CAF_RCTD <- norm_weights[,13]
  SO$APOE_CAF_RCTD[is.na(SO$APOE_CAF_RCTD)] <- 0
  # Save PNG of APOE_CAF
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("APOE_CAF_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/APOE_CAF_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Lu_ATC09_ATC from column 14 of norm_weights
  SO$ATC_RCTD <- norm_weights[,14]
  SO$ATC_RCTD[is.na(SO$ATC_RCTD)] <- 0
  # Save PNG of ATC
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("ATC_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/ATC_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  # Pull out Lu_ATC09_ATC_pEMT from column 15 of norm_weights
  SO$pEMT_ATC_RCTD <- norm_weights[,15]
  SO$pEMT_ATC_RCTD[is.na(SO$pEMT_ATC_RCTD)] <- 0
  # Save PNG of pEMT_ATC
  Spatial_Feature <- SpatialFeaturePlot(SO, features = c("pEMT_ATC_RCTD"), pt.size.factor = pt.size.factor, alpha = alpha)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pEMT_ATC_RCTD_Spatial_Feature_Plot.png"),
         Spatial_Feature, width = 5, height = 5, dpi = 600)
  
  ## pEMT spatial feature plot blends with ATC and myCAF
  # With ATC
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("pEMT_ATC_RCTD", "ATC_RCTD"), pt.size.factor = pt.size.factor)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pEMT_ATC_ATC_RCTD_Spatial_Feature_Blend_Default_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
  
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("pEMT_ATC_RCTD", "ATC_RCTD"), pt.size.factor = pt.size.factor,
                                                   top_left = "blue", bottom_right = "orange",
                                                   bottom_left = "white", top_right = "#FF0000")
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pEMT_ATC_ATC_RCTD_Spatial_Feature_Blend_Alt_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
  
  # With myCAF
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("pEMT_ATC_RCTD", "myCAF_RCTD"), pt.size.factor = pt.size.factor)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pEMT_ATC_myCAF_RCTD_Spatial_Feature_Blend_Default_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
  
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("pEMT_ATC_RCTD", "myCAF_RCTD"), pt.size.factor = pt.size.factor,
                                                   top_left = "blue", bottom_right = "orange",
                                                   bottom_left = "white", top_right = "#FF0000")
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/pEMT_ATC_myCAF_RCTD_Spatial_Feature_Blend_Alt_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
  
  # myCAF Spatial Feature Blend Plots
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("myCAF_RCTD", "ATC_RCTD"), pt.size.factor = pt.size.factor)
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/myCAF_ATC_RCTD_Spatial_Feature_Blend_Default_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
  
  Spatial_Feature_Blend <- SpatialFeaturePlotBlend(object = SO, features = c("myCAF_RCTD", "ATC_RCTD"), pt.size.factor = pt.size.factor,
                           top_left = "blue", bottom_right = "orange",
                           bottom_left = "white", top_right = "#FF0000")
  ggsave(file.path(resultsdir, "RCTD_SpatialFeature/myCAF_ATC_RCTD_Spatial_Feature_Blend_Alt_Color.png"),
         Spatial_Feature_Blend, width =18, height = 5, dpi = 600)
                                                   
return(SO)
}
