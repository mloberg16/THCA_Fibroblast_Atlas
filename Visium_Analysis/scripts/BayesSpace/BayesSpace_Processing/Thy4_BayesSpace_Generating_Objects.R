### Author: Matthew Aaron Loberg
### Date: 23-0113
### Script: Thy4_BayesSpace_Generating_Objects.R
### Source Script Name: 23-0113_Thy4_BayesSpace_Generating_Object.R

### Goal:
# Today I will generate SCE + Clustering for Thy4
# I will also generate SCE.enhanced and enhanced clustering

# I will be using the BayesSpace vignette/guide: 
# https://edward130603.github.io/BayesSpace/articles/BayesSpace.html


### Load required packages:
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Seurat)
#library(SeuratData) # think this might be necessary for coordinate transfer?
library(hdf5r) 
library(beepr)

### Load date
# Set directory for reading data
Dir <- "Data_in_Use/Thy4"

# Note that in contrast to Seurat, SCE requires both a filtered_feature_bc_matrix folder and a spatial folder
# The filtered_feature_bc_matrix folder should contain the following files: 
# barcodes.tsv.gz
# features.tsv.gz
# matrix.mtx.gz

# Load Thy4 as a single cell experiment
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

# Note: the qPlot leveled off at 7. Will go w/ 7 clusters. 

#
set.seed(149)
sce <- spatialCluster(sce, q=7, platform="Visium", d=10,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)

Thy4_Clustering <- clusterPlot(sce)
ggsave("outputs/Thy4_BayesSpace/23-0113_Thy4_BayesSpace_Clustering/23-0113_Thy4_BayesSpace_Clustering.png",
       Thy4_Clustering,
       width = 4, height = 5, dpi = 600)
rm(Thy4_Clustering)


sce.enhanced <- spatialEnhance(sce, q=7, platform="Visium", d=10,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

Thy4_Enhanced_Clustering <- clusterPlot(sce.enhanced)
ggsave("outputs/Thy4_BayesSpace/23-0113_Thy4_BayesSpace_Clustering/23-0113_Thy4_BayesSpace_EnhancedClustering.png",
       Thy4_Enhanced_Clustering,
       width = 4, height = 5, dpi = 600)
rm(Thy4_Enhanced_Clustering)

# Save sce and sce.enhanced as R objects to prevent need to rerun script
saveRDS(sce, "Data_in_Use/Processed_Outputs/Thy4_Processed/Thy4_BayesSpace.RDS")
saveRDS(sce.enhanced, "Data_in_Use/Processed_Outputs/Thy4_Processed/Thy4_BayesSpace_Enhanced.RDS"); beep()

# Cleaning up
rm(list = ls())
