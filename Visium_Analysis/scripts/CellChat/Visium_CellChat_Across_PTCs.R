# Author: Matthew Aaron Loberg
# Date: December 8, 2024
# Script: CellChat_Across_PTCs.R
# Source Script Name: 24-1208_CellChat_Across_PTCs.R

# Goal: Run CellChat across PTC Visium samples that have PTC, pEMT-PTC, and myCAF populations to infer L-R interactions with spatial info

# More info: 
# CellChat requires each individual visium capture area (barcode) to be labeled as one cell type
# This breaks the assumption/the reality that each capture area actually contains a mix of cell types
# However, it can still be used to infer spatial L-R interactions
# The way I will do this is to label each individual barcode based on the cell population (by RCTD deconvolution) that contributes the most to the composition of the spot 
# I will do this for all of the PTCs
# I will only proceed with CellChat analysis on the PTCs that have at least 10 spots with PTC AND pEMT-PTC

# Sample info: 
# The PTCs are Peds01 - Peds08 (pediatric samples) and Thy7, Thy15 - Thy18 (adult samples)

##### Load Packages For Pre-Processing #####
library(Seurat)
library(tidyverse)
library(patchwork) # required to wrap plots in SpatialFeaturePlotBlend

# Peds01
Peds01 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds01.RDS") # Read in Seurat object containing RCTD deconvolved percents as meta data
pt.size.factor <- 1.65
outputdir <- "outputs/Peds01_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds01@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max cell population
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column to the seurat object
Peds01$max_meta_column <- max_column

# View a table of the cell populations that compose the max
table(Peds01$max_meta_column)

# Custom colors for plotting
cols <- c("APOE_CAF_RCTD" = "#426600", 
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")

# Save a spatial plot with the barcodes colored based on which cell type is contributing the most to the composition of the barcode
ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds01, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds01, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds01 (not enough barcodes)
rm(list = ls())


# Peds02
Peds02 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds02.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds02_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds02@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds02$max_meta_column <- max_column

table(Peds02$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600", 
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          #"pDC_RCTD" = "#94FFB5", # no pDC max area
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds02, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds02, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds02 (not enough barcodes)
rm(list = ls())



# Peds03
Peds03 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds03.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds03_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds03@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds03$max_meta_column <- max_column

table(Peds03$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          #"dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          #"iCAF_RCTD" = "#FE8F42",
          #"iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          #"pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds03, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds03, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds03 (not enough barcodes)
rm(list = ls())


# Peds04
Peds04 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds04.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds04_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds04@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds04$max_meta_column <- max_column

table(Peds04$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds04, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds04_CellChat_Subset, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD"),
               pt.size = 0),
       width = 6, height = 3, dpi = 600)

# Will proceed with L-R interaction analysis; get rid of not using but KEEP Peds04
rm(meta_data, max_column, outputdir, pt.size.factor)

# Run SCTransform on Peds04
Peds04 <- Peds04 %>% SCTransform(assay = "Spatial", 
                                 return.only.var.genes = FALSE, 
                                 verbose = TRUE,
                                 vst.flavor = "v2")
Peds04_CellChat_Subset <- Peds04
rm(Peds04)



# Peds05
Peds05 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds05.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds05_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds05@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds05$max_meta_column <- max_column

table(Peds05$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds05, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds05, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds05 (not enough barcodes)
rm(meta_data, max_column, outputdir, pt.size.factor, Peds05)


# Peds06
Peds06 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds06.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds06_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds06@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds06$max_meta_column <- max_column

table(Peds06$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds06, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds06, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds06 (not enough barcodes)
rm(meta_data, max_column, outputdir, pt.size.factor, Peds06)




# Peds07
Peds07 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds07.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds07_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"
meta_data <- Peds07@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# for each barcode, define a max
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds07$max_meta_column <- max_column

table(Peds07$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds07, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)


ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds07, 
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# remove - not proceeding with Peds07 (not enough barcodes)
rm(meta_data, max_column, outputdir, pt.size.factor, Peds07)


Peds08 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds08.RDS")
pt.size.factor <- 1.65
outputdir <- "outputs/Peds08_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"

meta_data <- Peds08@meta.data
meta_data <- Peds08@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Peds08$max_meta_column <- max_column

table(Peds08$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600", 
          "B_Cell_RCTD" = "#993F00", 
          "Endothelial_RCTD" = "#2BCE48",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Peds08, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds08, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# will proceed with Cell Chat - remove things NOT using
rm(meta_data, max_column, pt.size.factor)

# Run SCTransform 
Peds08 <- Peds08 %>% SCTransform(assay = "Spatial", 
                                 return.only.var.genes = FALSE, 
                                 verbose = FALSE,
                                 vst.flavor = "v2")

# Need only things with a minimum count of 10
table(Peds08$max_meta_column)
Peds08_CellChat_Subset <- Peds08 %>% subset(max_meta_column != "APOE_CAF_RCTD" &
                                              max_meta_column != "B_Cell_RCTD" &
                                              max_meta_column != "NKT_RCTD" &
                                              max_meta_column != "Plasma_RCTD")
table(Peds08_CellChat_Subset$max_meta_column) # shows that this is appropriately down to groups with 10+ barcodes
rm(Peds08)

# repeat violin with CellChat subset
ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Peds08_CellChat_Subset, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)










########## Thy7 ##############
Thy7 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy7.RDS")
pt.size.factor <- 1.7
outputdir <- "outputs/Thy7_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"

meta_data <- Thy7@meta.data
meta_data <- Thy7@meta.data[,c("APOE_CAF_RCTD",
                                 "ATC_RCTD", 
                                 "B_Cell_RCTD", 
                                 "dPVCAF_RCTD",
                                 "Endothelial_RCTD", 
                                 "iCAF_RCTD",
                                 "iPVCAF_RCTD",
                                 "NKT_RCTD", 
                                 "myCAF_RCTD",
                                 "Myeloid_RCTD", 
                                 "pDC_RCTD", 
                                 "pEMT_PTC_RCTD",
                                 "Plasma_RCTD", 
                                 "PTC_RCTD",
                                 "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Thy7$max_meta_column <- max_column

table(Thy7$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Thy7, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Thy7, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# will proceed with Cell Chat - remove things NOT using
rm(meta_data, max_column, pt.size.factor)

# Run SCTransform 
Thy7 <- Thy7 %>% SCTransform(assay = "Spatial", 
                                 return.only.var.genes = FALSE, 
                                 verbose = FALSE,
                                 vst.flavor = "v2")

# Need only things with a minimum count of 10
table(Thy7$max_meta_column)
Thy7_CellChat_Subset <- Thy7 %>% subset(max_meta_column != "iPVCAF_RCTD" &
                                              max_meta_column != "NKT_RCTD")
table(Thy7_CellChat_Subset$max_meta_column) # shows that this is appropriately down to groups with 10+ barcodes
rm(Thy7)

# repeat violin with CellChat subset
ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln_CellChat_Subset.png"),
       VlnPlot(Thy7_CellChat_Subset, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)





########## Thy15 ##############
Thy15 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy15.RDS")
pt.size.factor <- 1.6
outputdir <- "outputs/Thy15_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"

meta_data <- Thy15@meta.data
meta_data <- Thy15@meta.data[,c("APOE_CAF_RCTD",
                               "ATC_RCTD", 
                               "B_Cell_RCTD", 
                               "dPVCAF_RCTD",
                               "Endothelial_RCTD", 
                               "iCAF_RCTD",
                               "iPVCAF_RCTD",
                               "NKT_RCTD", 
                               "myCAF_RCTD",
                               "Myeloid_RCTD", 
                               "pDC_RCTD", 
                               "pEMT_PTC_RCTD",
                               "Plasma_RCTD", 
                               "PTC_RCTD",
                               "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Thy15$max_meta_column <- max_column

table(Thy15$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Thy15, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Thy15, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# will proceed with Cell Chat - remove things NOT using
rm(meta_data, max_column, pt.size.factor)

# Run SCTransform 
Thy15 <- Thy15 %>% SCTransform(assay = "Spatial", 
                             return.only.var.genes = FALSE, 
                             verbose = FALSE,
                             vst.flavor = "v2")

# Need only things with a minimum count of 10
table(Thy15$max_meta_column)
Thy15_CellChat_Subset <- Thy15 %>% subset(max_meta_column != "APOE_CAF_RCTD")
table(Thy15_CellChat_Subset$max_meta_column) # shows that this is appropriately down to groups with 10+ barcodes
rm(Thy15)

# repeat violin with CellChat subset
ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln_CellChat_Subset.png"),
       VlnPlot(Thy15_CellChat_Subset, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)




########## Thy16 ##############
Thy16 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy16_V2.RDS")
pt.size.factor <- 1.6
outputdir <- "outputs/Thy16_RCTD/24-1104_iCAF2_Excluded_V2/RCTD_SpatialFeature/"

meta_data <- Thy16@meta.data
meta_data <- Thy16@meta.data[,c("APOE_CAF_RCTD",
                               "ATC_RCTD", 
                               "B_Cell_RCTD", 
                               "dPVCAF_RCTD",
                               "Endothelial_RCTD", 
                               "iCAF_RCTD",
                               "iPVCAF_RCTD",
                               "NKT_RCTD", 
                               "myCAF_RCTD",
                               "Myeloid_RCTD", 
                               "pDC_RCTD", 
                               "pEMT_PTC_RCTD",
                               "Plasma_RCTD", 
                               "PTC_RCTD",
                               "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Thy16$max_meta_column <- max_column

table(Thy16$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Thy16, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Thy16, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# will proceed with Cell Chat - remove things NOT using
rm(meta_data, max_column, pt.size.factor)

rm(Thy16) # NOT ENOUGH pEMT or myeloid to continue



########## Thy17 ##############
Thy17 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy17.RDS")
pt.size.factor <- 1.6
outputdir <- "outputs/Thy17_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"

meta_data <- Thy17@meta.data
meta_data <- Thy17@meta.data[,c("APOE_CAF_RCTD",
                                "ATC_RCTD", 
                                "B_Cell_RCTD", 
                                "dPVCAF_RCTD",
                                "Endothelial_RCTD", 
                                "iCAF_RCTD",
                                "iPVCAF_RCTD",
                                "NKT_RCTD", 
                                "myCAF_RCTD",
                                "Myeloid_RCTD", 
                                "pDC_RCTD", 
                                "pEMT_PTC_RCTD",
                                "Plasma_RCTD", 
                                "PTC_RCTD",
                                "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Thy17$max_meta_column <- max_column

table(Thy17$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Thy17, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Thy17, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# will proceed with Cell Chat - remove things NOT using
rm(meta_data, max_column, pt.size.factor)

# Run SCTransform 
Thy17 <- Thy17 %>% SCTransform(assay = "Spatial", 
                               return.only.var.genes = FALSE, 
                               verbose = FALSE,
                               vst.flavor = "v2")

# Need only things with a minimum count of 10
table(Thy17$max_meta_column)
Thy17_CellChat_Subset <- Thy17 %>% subset(max_meta_column != "APOE_CAF_RCTD" &
                                          max_meta_column != "Endothelial_RCTD" &
                                          max_meta_column != "NKT_RCTD")
table(Thy17_CellChat_Subset$max_meta_column) # shows that this is appropriately down to groups with 10+ barcodes
rm(Thy17)

# repeat violin with CellChat subset
ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln_CellChat_Subset.png"),
       VlnPlot(Thy17_CellChat_Subset, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)




########## Thy18 ##############
Thy18 <- readRDS(file = "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy18.RDS")
pt.size.factor <- 1.6
outputdir <- "outputs/Thy18_RCTD/24-1104_iCAF2_Excluded/RCTD_SpatialFeature/"

meta_data <- Thy18@meta.data
meta_data <- Thy18@meta.data[,c("APOE_CAF_RCTD",
                                "ATC_RCTD", 
                                "B_Cell_RCTD", 
                                "dPVCAF_RCTD",
                                "Endothelial_RCTD", 
                                "iCAF_RCTD",
                                "iPVCAF_RCTD",
                                "NKT_RCTD", 
                                "myCAF_RCTD",
                                "Myeloid_RCTD", 
                                "pDC_RCTD", 
                                "pEMT_PTC_RCTD",
                                "Plasma_RCTD", 
                                "PTC_RCTD",
                                "Thyrocyte_RCTD")]

# Apply function to get the column name with the maximum value for each cell
max_column <- apply(meta_data, 1, function(x) colnames(meta_data)[which.max(x)])

# Add this as a new metadata column
Thy18$max_meta_column <- max_column

table(Thy18$max_meta_column)

cols <- c("APOE_CAF_RCTD" = "#426600",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          "dPVCAF_RCTD" = "#990000",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42",
          "iPVCAF_RCTD" = "#0075DC",
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "NKT_RCTD" = "#F0A0FF", 
          "pDC_RCTD" = "#94FFB5",
          "pEMT_PTC_RCTD" = "black", 
          "Plasma_RCTD" = "#005C31", 
          "PTC_RCTD" = "#FFCC99",
          "Thyrocyte_RCTD" = "lightgrey")


ggsave(file.path(outputdir, "24-1208_max_meta_column.png"),
       SpatialDimPlot(Thy18, 
                      group.by = "max_meta_column", 
                      cols = cols,
                      pt.size.factor = pt.size.factor),
       width = 6, height = 5, dpi = 600)

ggsave(file.path(outputdir, "24-1208_myCAF_PTC_pEMT_Vln.png"),
       VlnPlot(Thy18, pt.size = 0,
               group.by = "max_meta_column",
               cols = cols, 
               features = c("myCAF_RCTD", "PTC_RCTD", "pEMT_PTC_RCTD")),
       width = 6, height = 3, dpi = 600)

# no cell chat for Thy18 (almost exclusively PTC) - proceed to CellChat w/ out
rm(meta_data, max_column, pt.size.factor, Thy18, max_column, cols, outputdir)






################### Cell Chat analysis ########################
# See CellChat Tutorial - was used to write this code
### Load Required Packages
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Create a list of Cell_Chat_SOs (These are the Seurat Objects that contained PTC and pEMT-PTC popluations)
Cell_Chat_SOs <- list(Peds04_CellChat_Subset, 
                      Peds08_CellChat_Subset,
                      Thy7_CellChat_Subset,
                      Thy15_CellChat_Subset,
                      Thy17_CellChat_Subset)

# Clean up individual objecs
rm(Peds04_CellChat_Subset, 
   Peds08_CellChat_Subset,
   Thy7_CellChat_Subset,
   Thy15_CellChat_Subset,
   Thy17_CellChat_Subset)



# Define expression data.input for CellChat analysis for each SO in the list
data.input <- list()

for(i in 1:length(Cell_Chat_SOs)){
  data.input[[i]] <- Seurat::GetAssayData(Cell_Chat_SOs[[i]], slot = "data", assay = "SCT")
}

# Define the meta data for each SO
meta <- list()

for(i in 1:length(Cell_Chat_SOs)){
  Idents(Cell_Chat_SOs[[i]]) <- Cell_Chat_SOs[[i]]$max_meta_column
  meta[[i]] <- data.frame(labels = Seurat::Idents(Cell_Chat_SOs[[i]]),
                          samples = Cell_Chat_SOs[[i]]$orig.ident[1],
                          row.names = names(Seurat::Idents(Cell_Chat_SOs[[i]])))
  meta[[i]]$samples <- factor(meta[[i]]$samples)
}


# Load spatial transcriptomics information
spatial.locs <- list()
scalefactors <- c(235.6634359013668,  # Peds04 scale factor
                  235.6667708600977,  # Peds08 scale factor
                  239.13277820640744, # Thy7 scale factor
                  235.46194252361263, # Thy15 scale factor
                  235.47463457815473) # Thy17 scale factor
spot.size <- 65
conversion.factor <- c()
spatial.factors <- list()

for(i in 1:length(Cell_Chat_SOs)){
  spatial.locs[[i]] <- Seurat::GetTissueCoordinates(Cell_Chat_SOs[[i]], scale = NULL, cols = c("imagerow", "imagecol"))
  conversion.factor[i] <- spot.size/scalefactors[i]
  spatial.factors[[i]] <- data.frame(ratio = conversion.factor[i],
                                     tol = spot.size/2)
}


# compute d.spatial
d.spatial <- list()
for(i in 1:length(Cell_Chat_SOs)){
  d.spatial[[i]] <- computeCellDistance(coordinates = spatial.locs[[i]], ratio = spatial.factors[[i]]$ratio, tol = spatial.factors[[i]]$tol)
  print(min(d.spatial[[i]][d.spatial[[i]]!=0])) # this value should approximately equal 100um for 10X Visium data - and it is for me exactly 100 for Peds and close for adults!
}


# Create Cell Chat object
cellchat <- list()

for(i in 1:length(Cell_Chat_SOs)){
  cellchat[[i]] <- createCellChat(object = data.input[[i]],
                                  meta = meta[[i]],
                                  group.by = "labels",
                                  datatype = "spatial",
                                  coordinates = spatial.locs[[i]],
                                  spatial.factors = spatial.factors[[i]])
}

# checking that cell chat objects appropriately generated
cellchat[[1]]
cellchat[[2]]
cellchat[[3]]
cellchat[[4]]
cellchat[[5]]

##### RUNNING CELL CHAT

##### CHOOSING DATABASES
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# There are a number of options to choose from; I am detailing them below

# Option:
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Option: 
# Only uses the Secreted Signaling from CellChatDB v1
# CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# Option: 
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# Option: 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling) that can be only estimated from gene expression data.

# I am going to go with USE ALL EXCEPT FOR NON-PROTEIN SIGNALING
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
for(i in 1:length(cellchat)){
  cellchat[[i]]@DB <- CellChatDB.use
}



##### PREPROCESSING 
# subset the expression data of signaling genes for saving computation cost

for(i in 1:length(cellchat)){
  cellchat[[i]] <- subsetData(cellchat[[i]]) # This step is necessary even if using the whole database
}

future::plan("multisession", workers = 4) 

for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% identifyOverExpressedGenes(do.fast = FALSE)
}

for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% identifyOverExpressedInteractions(variable.both = F)
}


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


###### COMPUTE COMMUNICATION PROBABILITY + INFER CELLULAR COMMUNICATION NETWORK

for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% computeCommunProb(type = "triMean",
                                                       distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                                                       contact.dependent = TRUE, contact.range = 100) 
}


# default of 10 cells required in each group for cell-cell communication
for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% filterCommunication(min.cells = 10)
}



##### Extract the inferred cellular communication network as a data frame

# Background info from tutorial (https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html): 
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

###### INFERENCE AT SIGNALING PATHWAY LEVEL #####
for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% computeCommunProbPathway()
}


###### CALCULATE AGGREGATED cell-cell communication network #####
for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% aggregateNet()
}


outputdir <- "outputs/CellChat/"
dir.create(outputdir)

col1 <- c(          "Plasma_RCTD" = "#005C31", 
                    "iCAF_RCTD" = "#FE8F42",
                    "Myeloid_RCTD" = "#808080", 
                    "dPVCAF_RCTD" = "#990000",
                    "NKT_RCTD" = "#F0A0FF",
                    "iPVCAF_RCTD" = "#0075DC",
                    "Thyrocyte_RCTD" = "lightgrey",
                    "Endothelial_RCTD" = "#2BCE48",
          "PTC_RCTD" = "#FFCC99",
          "pEMT_PTC_RCTD" = "black",
          "ATC_RCTD" = "lightblue",
          "B_Cell_RCTD" = "#993F00",
          
          "pDC_RCTD" = "#94FFB5",

          "myCAF_RCTD" = "#783FC1",


          "APOE_CAF_RCTD" = "#426600"
            )

col2 <- c("PTC_RCTD" = "#FFCC99",
          "pEMT_PTC_RCTD" = "black", 
          "myCAF_RCTD" = "#783FC1",
          "Myeloid_RCTD" = "#808080", 
          "Endothelial_RCTD" = "#2BCE48"
)

col3 <- c("PTC_RCTD" = "#FFCC99",
          "myCAF_RCTD" = "#783FC1",
          "Thyrocyte_RCTD" = "lightgrey",
          "pEMT_PTC_RCTD" = "black", 
          "Myeloid_RCTD" = "#808080",
          "APOE_CAF_RCTD" = "#426600",
          "Endothelial_RCTD" = "#2BCE48",
          "iCAF_RCTD" = "#FE8F42", 
          "dPVCAF_RCTD" = "#990000" 
          )

col4 <- c("PTC_RCTD" = "#FFCC99",
          "pEMT_PTC_RCTD" = "black", 
          "Endothelial_RCTD" = "#2BCE48",
          "dPVCAF_RCTD" = "#990000",
          "B_Cell_RCTD" = "#993F00",
          "iCAF_RCTD" = "#FE8F42",
          "Thyrocyte_RCTD" = "lightgrey",
          "NKT_RCTD" = "#F0A0FF",
          "Myeloid_RCTD" = "#808080", 
          "myCAF_RCTD" = "#783FC1")

col5 <- c("PTC_RCTD" = "#FFCC99",
          "Myeloid_RCTD" = "#808080", 
          "myCAF_RCTD" = "#783FC1",
          "pEMT_PTC_RCTD" = "black")
  


cols <- list(col1, col2, col3, col4, col5)

objects <- c("Peds04", "Peds08", "Thy7", "Thy15", "Thy17")
for(i in 1:length(cellchat)){
  
  groupSize <- as.numeric(table(cellchat[[i]]@idents))
  
  dir.create(paste0(outputdir, objects[i]))
  
  png(file = paste0(outputdir, objects[i], "/Interaction_Number.png"), width = 5, height = 5, units = "in", res = 2000)
  print(netVisual_circle(cellchat[[i]]@net$count, 
                         vertex.weight = groupSize, 
                         weight.scale = T, 
                         label.edge= F, 
                         title.name = "Number of interactions",
                         color.use = cols[[i]],
                         vertex.label.cex = 0.00000000001))
  dev.off()
  
  png(file = paste0(outputdir, objects[i], "/Interaction_Weight.png"), width = 5, height = 5, units = "in", res = 2000)
  print(netVisual_circle(cellchat[[i]]@net$weight, 
                         vertex.weight = groupSize, 
                         weight.scale = T, 
                         label.edge= F, 
                         title.name = "Interaction weights/strength",
                         color.use = cols[[i]],
                         vertex.label.cex = 0.00000000001))
  dev.off()
}

png("outputs/Lu_etal_2024_Analysis_Outputs/24-0924_Individual_Analysis/ATC09_Lu/CellChat/Interaction_Number.png", width = 5, height = 5, units = "in", res = 600)
groupSize <- as.numeric(table(cellChat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
png("outputs/Lu_etal_2024_Analysis_Outputs/24-0924_Individual_Analysis/ATC09_Lu/CellChat/Interaction_Weight.png", width = 5, height = 5, units = "in", res = 600)
groupSize <- as.numeric(table(cellChat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")



df.net <- subsetCommunication(cellchat[[1]])

df.net.p <- subsetCommunication(cellchat[[1]], slot.name = "netP")


Interaction_Weights <- data.frame("Interaction" = c("myCAF_PTC", "myCAF_PTC", "myCAF_PTC", "myCAF_PTC",
                                                    "myCAF_pEMT", "myCAF_pEMT", "myCAF_pEMT", "myCAF_pEMT"),
                                  "Weight" = c(cellchat[[1]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                               cellchat[[2]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                               cellchat[[3]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                               cellchat[[5]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                               cellchat[[1]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                               cellchat[[2]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                               cellchat[[3]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                               cellchat[[5]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")]),
                                  "Sample" = c("Peds04", "Peds08", "Thy7", "Thy17",
                                               "Peds04", "Peds08", "Thy7", "Thy17"))



Interaction_Weights_stats <- data.frame("PTC_myCAF" = c(cellchat[[1]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                                        cellchat[[2]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                                        cellchat[[3]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")],
                                                        cellchat[[5]]@net$weight[c("myCAF_RCTD"), c("PTC_RCTD")]),
                                       "pEMT_myCAF" = c(cellchat[[1]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                                        cellchat[[2]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                                        cellchat[[3]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")],
                                                        cellchat[[5]]@net$weight[c("myCAF_RCTD"), c("pEMT_PTC_RCTD")]))


plot <- ggplot(Interaction_Weights, aes(Interaction, Weight)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Interaction),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_point(aes(),
             #position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  geom_line(aes(group = Sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("black", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Distance") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("myCAF_PTC", "myCAF_pEMT")) +
  scale_y_continuous(breaks = c(0, 4, 8, 12),
                     limits = c(0, 12)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/CellChat/24-1208_PTCs/24-1208_PTC_pEMT_myCAF_Interaction_Weights.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

Interaction_Weights_Stats <- pivot_wider(Interaction_Weights, names_from = Interaction, values_from = Weight)
wilcox.test(Interaction_Weights_stats$PTC_myCAF, Interaction_Weights_stats$pEMT_myCAF, paired = TRUE) # two-sided t-test: 0.006445




##### VISUALIZATION OF CELL-CELL COMMUNICATION NETWORK
pathways.show <- c("TGFb") 
# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute the network centrality scores

for(i in 1:length(cellchat)){
  cellchat[[i]] <- cellchat[[i]] %>% netAnalysis_computeCentrality(slot.name = "netP")
}

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "COLLAGEN", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "LAMININ", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "TENASCIN", width = 8, height = 2.5, font.size = 10)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2), remove.isolate = FALSE)
#> Comparing communications on a single object
#> 

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2), signaling = c("TGFb", "CCL", "CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2), signaling = c("COLLAGEN", "LAMININ", "TENASCIN"), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = c(1:3), signaling = c("COLLAGEN", "LAMININ", "TENASCIN"), remove.isolate = FALSE)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(3), targets.use = c(1:2), slot.name = "netP", legend.pos.x = 10)

runCellChatApp(cellchat)

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2), signaling = c("COLLAGEN"), remove.isolate = FALSE)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat[[3]], pattern = "outgoing", cluster.cols = FALSE, height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat[[3]], pattern = "incoming", cluster.cols = FALSE, height = 16)
ht1 + ht2

# repeat with just TOP 10 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat[[3]], 
                                         signaling = c("COLLAGEN", "FN1", "LAMININ", "THBS", "TENASCIN",
                                                       "APP", "MK", "PERIOSTIN", "CD99", "ADGRG"),
                                         pattern = "outgoing", 
                                         cluster.cols = FALSE, 
                                         height = 8)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat[[3]],
                                         signaling = c("COLLAGEN", "FN1", "LAMININ", "THBS", "TENASCIN",
                                                       "APP", "MK", "PERIOSTIN", "CD99", "ADGRG"),
                                         pattern = "incoming", 
                                         cluster.cols = FALSE, 
                                         height = 8)
ht1 + ht2

# extract netp probablilities
cellchat@netP$prob
cellchat@netP$prob[,,c("TENASCIN")] # Extracts just TENASCIN pathway

cellchat@netP$centr
