# Author: Matthew A. Loberg
# Date: November 29th, 2022
# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will read the data into R studio and begin basic processing of the data

# Thy20

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy20 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy20 data
data_dir <- 'Data_in_Use/August_2022_VANTAGE_Visium_Run/Raw_SpaceRanger_Outputs/8405-CP-0016_S14-26091_3_Thy20' # Set directory to load from
Thy20 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy20
Thy20$orig.ident <- "Thy20"
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers"
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Thy20@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Thy20@images[["slice1"]]@coordinates[["tissue"]])
Thy20@images[["slice1"]]@coordinates[["row"]] <- as.integer(Thy20@images[["slice1"]]@coordinates[["row"]])
Thy20@images[["slice1"]]@coordinates[["col"]] <- as.integer(Thy20@images[["slice1"]]@coordinates[["col"]])
Thy20@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Thy20@images[["slice1"]]@coordinates[["imagerow"]])
Thy20@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Thy20@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy20, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy20, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy20_QC/22-1129_Thy20_Processing_Raw_SCTransform/22-1129_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy20_QC/22-1129_Thy20_Processing_Raw_SCTransform/22-1129_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy20 as an RDS
saveRDS(Thy20, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy20_Processed/22-1129_Thy20_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands


Thy20 <- SCTransform(Thy20, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Save SCTransformed Thy20 as an RDS
saveRDS(Thy20, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy20_Processed/22-1129_Thy20_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Thy20)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy20

