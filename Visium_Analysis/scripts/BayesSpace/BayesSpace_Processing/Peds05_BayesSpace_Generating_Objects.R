### Author: Matthew Aaron Loberg
### Date: 24-0524
### Script: Peds05_BayesSpace_Generating_Object.R
### Source Script Name: 24-0524_Peds05_BayesSpace_Generating_Object.R

### Goal:
# Today I will generate SCE + Clustering for Peds05
# I will also generate SCE.enhanced and enhanced clustering

# I will be using the BayesSpace vignette/guide: 
# https://edward130603.github.io/BayesSpace/articles/BayesSpace.html


### Load required packages:
library(SingleCellExperiment); 
library(ggplot2)
library(BayesSpace)
library(Seurat)
#library(SeuratData) # think this might be necessary for coordinate transfer?
library(hdf5r) 
library(beepr) # for making beep noise

### Load data
# Set directory for reading data
Dir <- "Data_in_Use/Raw_SpaceRanger_Outputs/Peds_05"

# Note that in contrast to Seurat, SCE requires both a filtered_feature_bc_matrix folder and a spatial folder
# The filtered_feature_bc_matrix folder should contain the following files: 
# barcodes.tsv.gz
# features.tsv.gz
# matrix.mtx.gz

# Load Peds05 as a single cell experiment
# Note: originally got the following error: 
# cannot open file 'Data_in_Use/Belcher_Peds_Visium/Raw_SpaceRanger_Outputs/Peds_01/spatial/tissue_positions_list.csv': No such file or directory
# My "tissue_positions_list.csv" file was just "tissue_positions.csv" -> I added the "_list" to it and it worked
# In order to get it to work, I also had to DELETE THE FIRST ROW
# I only knew to do this because I have had this problem prior
sce <- readVisium(Dir)

#### Pre-processing data ####
set.seed(102)
# Note, platform was previously set to "ST" and that was WRONG
# Platform needs to be set to "Visium" -> this is the reason that I was not finding any neighbors
sce <- spatialPreprocess(sce, platform = "Visium",
                         n.PCs = 10, n.HVGs = 2000, log.normalize = TRUE)

# Clustering
# Selection of number of clusters
# qTune
# Testing qs 2 through 10
# I believe that d = the dimensionality of the data set?? (e.g., number of PCs)
sce <- qTune(sce, qs = seq(2,10), platform = "Visium", d = 10) # set the seq paramater to the minimum and maximum number of expected clusters
qPlot(sce); beep()

# First elbow at 4; going with 4 for now

set.seed(149)
sce <- spatialCluster(sce, q=4, platform="Visium", d=10,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)

Peds05_Clustering <- clusterPlot(sce)
ggsave("outputs/Peds05_BayesSpace/24-0524_Peds05_BayesSpace_Clustering/24-0524_Peds05_BayesSpace_Clustering.png",
       Peds05_Clustering,
       width = 4, height = 5, dpi = 600)
rm(Peds05_Clustering)


sce.enhanced <- spatialEnhance(sce, q=4, platform="Visium", d=10,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

Peds05_Enhanced_Clustering <- clusterPlot(sce.enhanced)
ggsave("outputs/Peds05_BayesSpace/24-0524_Peds05_BayesSpace_Clustering/14-0524_Peds05_BayesSpace_EnhancedClustering.png",
       Peds05_Enhanced_Clustering,
       width = 4, height = 5, dpi = 600)
rm(Thy1_Enhanced_Clustering)

# Save sce and sce.enhanced as R objects to prevent need to rerun script
saveRDS(sce, "Data_in_Use/Processed_Outputs/Peds05/Peds05_BayesSpace.RDS")
saveRDS(sce.enhanced, "Data_in_Use//Processed_Outputs/Peds05/Peds05_BayesSpace_Enhanced.RDS")

# Cleaning up
rm(list = ls()); beep()

