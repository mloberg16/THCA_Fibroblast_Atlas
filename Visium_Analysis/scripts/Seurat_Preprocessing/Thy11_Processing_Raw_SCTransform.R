### Author: Matthew A. Loberg
### Date: September 5th, 2022
### Script: Thy11_Processing_Raw_SCTransform.R
### Source Script Name: 22-0905_Thy11_Processing_Raw_SCTransform.R

### Goal: 
# Here, I will read the data into R studio and begin basic processing of the data
# I will save a seurat object as a .RDS, which I will use for future analysis

# Thy11

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy11 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy11 data
data_dir <- 'Data_in_Use/Raw_SpaceRanger_Outputs/Thy11' # Set directory to load from
Thy11 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy11
Thy11$orig.ident <- "Thy11"
# Cleaning up
rm(data_dir)

# I was having a problem with spatial feature plot coordinates being "characters" instead of "integers"
# When I ran spatial feature plot I was getting the following error: 
# "Error in FUN(left, right) : non-numeric argument to binary operator
# According to stack overflow, the following code should fix the issue 
# See line: https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
Thy11@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(Thy11@images[["slice1"]]@coordinates[["tissue"]])
Thy11@images[["slice1"]]@coordinates[["row"]] <- as.integer(Thy11@images[["slice1"]]@coordinates[["row"]])
Thy11@images[["slice1"]]@coordinates[["col"]] <- as.integer(Thy11@images[["slice1"]]@coordinates[["col"]])
Thy11@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(Thy11@images[["slice1"]]@coordinates[["imagerow"]])
Thy11@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(Thy11@images[["slice1"]]@coordinates[["imagecol"]])

# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy11, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy11, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy11_QC/22-0905_Thy11_Processing_Raw_SCTransform/22-0905_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy11_QC/22-0905_Thy11_Processing_Raw_SCTransform/22-0905_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy11 Seurat Object as a .RDS
saveRDS(Thy11, "Data_in_Use/Processed_Outputs/Thy11_Processed/22-0905_Thy11_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands
Thy11 <- SCTransform(Thy11, 
                     vst.flavor = "v2",
                     assay = "Spatial", 
                     return.only.var.genes = FALSE, 
                     verbose = FALSE)

# Save SCTransformed Thy11 Seurat Object as a .RDS
saveRDS(Thy11, "Data_in_Use/Processed_Outputs/Thy11_Processed/22-0905_Thy11_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Thy11)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy11

