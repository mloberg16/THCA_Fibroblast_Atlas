# Author: Matthew A. Loberg
# Date: September 13th, 2022
# Script: Thy1_Processing_Raw_SCTransform.R
# Source Script Name: 22-0913_Thy1_Processing_Raw_SCTransform.R

# Goal: 
# Here, I will read the data into R studio and begin basic processing of the data
# I will save a seurat object as a .RDS, which I will use for future analysis

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy1 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy1 data
data_dir <- 'Data_in_Use/Raw_SpaceRanger_Outputs/Thy1' # Set directory to load from
Thy1 <- Seurat::Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy1
Thy1$orig.ident <- "Thy1"
# Cleaning up
rm(data_dir)


# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy1, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy1, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy1_QC/22-0913_Thy1_Processing_Raw_SCTransform/22-0913_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy1_QC/22-0913_Thy1_Processing_Raw_SCTransform/22-0913_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy1 Seurat Object as an RDS (prior to SCTransform - just contains raw count data)
saveRDS(Thy1, "Data_in_Use/Processed_Outputs/Thy1_Processed/22-0913_Thy1_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands
Thy1 <- SCTransform(Thy1, 
                    vst.flavor = "v2",
                    assay = "Spatial", 
                    return.only.var.genes = FALSE, 
                    verbose = FALSE)

# Save SCTransformed Thy1 Seurat Object as an RDS for future use
saveRDS(Thy1, "Data_in_Use/Processed_Outputs/Thy1_Processed/22-0913_Thy1_SCTransformed_All_Genes.rds")

# Just testing spatial feature plot
SpatialFeaturePlot(Thy1, features = "TGFB1")
SpatialFeaturePlot(Thy1, features = "GREM1")

# Cleaning up
rm(Thy1)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy1
# Can use this SCTransformed assay to addmodule score, perform other analyses

