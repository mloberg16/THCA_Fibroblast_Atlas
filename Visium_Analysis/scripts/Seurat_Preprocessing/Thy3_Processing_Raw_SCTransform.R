### Author: Matthew A. Loberg
### Date: November 29th, 2022
### Script: Thy3_Processing_Raw_SCTransform.R
### Source Script Name: 22-1129_Thy3_Processing_Raw_SCTransform.R

### Goal: 
# Here, I will read the data into R studio and begin basic processing of the data
# I will save a seurat object as a .RDS, which I will use for future analysis

# Thy3

#### Chapter 1: Loading Packages ####
# Load required packages
library(Seurat)
library(hdf5r) # required to read in data file
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)

#### Chapter 2: Reading in Thy3 and looking at raw count data by violin and SpatialFeaturePlot ####

# Load in Thy3 data
data_dir <- 'Data_in_Use/Thy3' # Set directory to load from
Thy3 <- Load10X_Spatial(data.dir = data_dir, slice = "slice1") # Load Thy3
Thy3$orig.ident <- "Thy3"
# Cleaning up
rm(data_dir)


# Visualize raw count data as a violin plot. Note, 'raster = FALSE' for image quality
plot1 <- VlnPlot(Thy3, features = "nCount_Spatial", raster = FALSE) + NoLegend()

# Visiualize spatial location of raw count data
plot2 <- SpatialFeaturePlot(Thy3, features = "nCount_Spatial") + theme(legend.position = "right")

# Visiualize spatial location and violin plot on the same graph 
test <- wrap_plots(plot1, plot2)

# Format plot1 and plot2
plot1 <- plot1 + 
  theme(
    axis.text = element_text(face = "bold", size = 15)
  )
ggsave("outputs/Thy3_QC/22-1129_Thy3_Processing_Raw_SCTransform/22-1129_Raw_Counts_Violin.png",
       plot1,
       width = 4, height = 5, dpi = 600)

plot2 <- plot2 + theme(
  legend.text = element_text(face = "bold", size = 15),
  legend.title = element_text(face = "bold", size = 15)
)
ggsave("outputs/Thy3_QC/22-1129_Thy3_Processing_Raw_SCTransform/22-1129_Raw_Counts_Spatial.png",
       plot2,
       width = 7, height = 5, dpi = 600)

# Cleaning up
rm(plot1, plot2, test)

# Save raw Thy3 Seurat Object as an RDS
saveRDS(Thy3, "Data_in_Use/Thy3_Processed/22-1129_Thy3_Raw_PreProcessed.rds")

#### Chapter 3: Data Transformation ####
# I will perform data transformation with SCTransform
# See SCTransform vignette here: https://satijalab.org/seurat/articles/sctransform_vignette.html
# See SCTransform paper here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# There are many discussions on whether return.only.var.genes should be set to TRUE/FALSE.
# In the vignette from satijalab, they recommend FALSE for best performance. 
# I need to do more reading to see how this affects addModuleScore and other commands
Thy3 <- SCTransform(Thy3, 
                    vst.flavor = "v2",
                    assay = "Spatial", 
                    return.only.var.genes = FALSE, 
                    verbose = FALSE)

# Save SCTransformed Thy3 Seurat Object as an RDS
saveRDS(Thy3, "Data_in_Use/Thy3_Processed/22-1129_Thy3_SCTransformed_All_Genes.rds")

# Cleaning up
rm(Thy3)

# I will start subsequent analysis by loading the SCTransformed/processed version of Thy3

