### Author: Matthew A. Loberg
### Date: November 29th, 2022
### Script: Thy18_Processing_Raw_SCTransform.R
### Source Script Name: 22-1129_Thy18_Processing_Raw_SCTransform.R

# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will read the data into R studio and begin basic processing of the data

# Goal: 
# Here, I will read the data into R studio and begin basic processing of the data
# I will save a seurat object as a .RDS, which I will use for future analysis

# Thy18

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy18 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy18 data from SpaceRanger output folder
data_dir <- 'Data_in_Use/Raw_SpaceRanger_Outputs/Thy18' # Set directory to load from
Thy18 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy18
Thy18$orig.ident <- "Thy18"
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers"
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Thy18@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Thy18@images[["slice1"]]@coordinates[["tissue"]])
Thy18@images[["slice1"]]@coordinates[["row"]] <- as.integer(Thy18@images[["slice1"]]@coordinates[["row"]])
Thy18@images[["slice1"]]@coordinates[["col"]] <- as.integer(Thy18@images[["slice1"]]@coordinates[["col"]])
Thy18@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Thy18@images[["slice1"]]@coordinates[["imagerow"]])
Thy18@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Thy18@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy18, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy18, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy18_QC/22-1129_Thy18_Processing_Raw_SCTransform/22-1129_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy18_QC/22-1129_Thy18_Processing_Raw_SCTransform/22-1129_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy18 Seurat Object as a .RDS
saveRDS(Thy18, "Data_in_Use/Processed_Outputs/Thy18_Processed/22-1129_Thy18_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands


Thy18 <- SCTransform(Thy18, 
                     vst.flavor = "v2",
                     assay = "Spatial", 
                     return.only.var.genes = FALSE, 
                     verbose = FALSE)

# Save SCTransformed Thy18 Seurat Object as a .RDS
saveRDS(Thy18, "Data_in_Use/Processed_Outputs/Thy18_Processed/22-1129_Thy18_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Thy18)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy18

