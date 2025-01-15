# Author: Matthew A. Loberg
# Date: September 15th, 2022
# Purpose: New Visium sequencing data just obtained from Vantage
# Here, I will read the data into R studio and begin basic processing of the data

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy10 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy10 data
data_dir <- 'Data_in_Use/August_2022_VANTAGE_Visium_Run/8405-CP-0006_S10-12119_1P_Thy10' # Set directory to load from
Thy10 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy10
Thy10$orig.ident <- "Thy10"
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers"
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Thy10@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Thy10@images[["slice1"]]@coordinates[["tissue"]])
Thy10@images[["slice1"]]@coordinates[["row"]] <- as.integer(Thy10@images[["slice1"]]@coordinates[["row"]])
Thy10@images[["slice1"]]@coordinates[["col"]] <- as.integer(Thy10@images[["slice1"]]@coordinates[["col"]])
Thy10@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Thy10@images[["slice1"]]@coordinates[["imagerow"]])
Thy10@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Thy10@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy10, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy10, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy10_QC/22-0915_Thy10_Processing_Raw_SCTransform/22-0915_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy10_QC/22-0915_Thy10_Processing_Raw_SCTransform/22-0915_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy10 as an RDS
saveRDS(Thy10, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Thy10_Processed/22-0915_Thy10_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands
Thy10 <- SCTransform(Thy10, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Save SCTransformed Thy10 as an RDS
saveRDS(Thy10, "Data_in_Use/August_2022_VANTAGE_Visium_Run/Thy10_Processed/22-0915_Thy10_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Thy10)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy10

