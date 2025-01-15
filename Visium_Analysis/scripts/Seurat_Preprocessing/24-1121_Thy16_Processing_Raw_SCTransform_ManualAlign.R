# Author: Matthew A. Loberg
# Date: November 21th, 2024
# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will read the data into R studio and begin basic processing of the data

# Thy16

# 24-1121 Update
# Re-running with manual align file from Lana Olson from 24-1120

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy16 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy16 data
data_dir <- 'Data_in_Use/August_2022_VANTAGE_Visium_Run/Raw_SpaceRanger_Outputs/Thy16_manAlign' # Set directory to load from
Thy16 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy16
Thy16$orig.ident <- "Thy16"
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers"
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Thy16@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Thy16@images[["slice1"]]@coordinates[["tissue"]])
Thy16@images[["slice1"]]@coordinates[["row"]] <- as.integer(Thy16@images[["slice1"]]@coordinates[["row"]])
Thy16@images[["slice1"]]@coordinates[["col"]] <- as.integer(Thy16@images[["slice1"]]@coordinates[["col"]])
Thy16@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Thy16@images[["slice1"]]@coordinates[["imagerow"]])
Thy16@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Thy16@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy16, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy16, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy16_QC/24-1121_Thy16_Processing_Raw_SCTransform/24-1121_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy16_QC/24-1121_Thy16_Processing_Raw_SCTransform/24-1121_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy16 as an RDS
saveRDS(Thy16, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy16_Processed/24-1121_Thy16_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands
# The scTransform is returning an error because one of my cells has no UMIs ... I need to filter out cells (spots) that have no UMIs before proceeding
min(Thy16$nCount_Spatial) # Shows that there are spots with 0 counts
Thy16 <- subset(Thy16, nCount_Spatial > 0) # Subset out only spots that have > 0 counts
Thy16 <- SCTransform(Thy16, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE) # Now I can perform SCTransform without error

# Save SCTransformed Thy16 as an RDS
saveRDS(Thy16, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy16_Processed/24-1121_Thy16_SCTransformed_All_Genes.rds")

SpatialFeaturePlot(Thy16, features = c("FN1"), pt.size.factor = 1.6)
SpatialFeaturePlot(Thy16, features = c("TG"), pt.size.factor = 1.7)
SpatialFeaturePlot(Thy16, features = c("FAP"))
SpatialFeaturePlot(Thy16, features = c("ACTA2"))
SpatialFeaturePlot(Thy16, features = c("RGS5"))
SpatialFeaturePlot(Thy16, features = c("KRT8"), alpha = 0)

# Cleaning up
rm(Thy16)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy16

