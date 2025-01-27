# Author: Matthew Aaron Loberg
# Date: June 30, 2023
# Script: 23-0630_GSM5814583_ATC_08_Transformation_and_Clustering.R

#### Load Packages ####
library(tidyverse)
library(Seurat)

### Load data ####
# Load the raw SO generated from "23-0630_GSM5814583_ATC_08_Intital_Attempt.R"
SO <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814583_ATC08.RDS")

#### SCTransform ####
# Here, I will be trying out the V2 of SCTransform, recently released
# SCTransform V2 guide:
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
# SCTransform guide:
# https://satijalab.org/seurat/articles/sctransform_vignette.html
# Original SCTransform paper:
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
# SCTransform V2 Paper:
# https://www.biorxiv.org/content/10.1101/2021.07.07.451498v1.full.pdf
# Note: SCTransform v2 requires glmGamPoi (install from BiocManager)
# I had some trouble installing glmGamPoi
# Updating to R version 4.3.1 (from 4.2.2) solved the errors I was having

SO <- SCTransform(SO, vst.flavor = "v2", verbose = TRUE)

#### PCA, Dimensionality Reduction, and Clustering ####
# Run PCA
SO <- RunPCA(SO)
# Vizualize PCA results
DimPlot(SO)
DimHeatmap(SO)
ElbowPlot(SO, ndims = 50)
# Based on elbow plot will proceed with 30 dimensions
SO <- RunUMAP(SO, dims = 1:30, verbose = FALSE)
SO <- FindNeighbors(SO, dims = 1:30, verbose = FALSE)
SO <- FindClusters(SO, verbose = FALSE, resolution = 0.5)
DimPlot(SO, label = TRUE) + NoLegend()
FeaturePlot(SO, features = c("KRT8", "KRT5", "FAP"))

# Read in cell type annotation
Annotation <- read.table("data_in_use/Lu_etal_2023_ATC_scRNA/GSE193581_celltype_annotation.txt", sep = "\t")
Annotation <- Annotation %>% subset(sample.ID == "ATC08")
SO$External_Cell_Type <- Annotation$celltype

DimPlot(SO, group.by = "External_Cell_Type")
SO_Myeloid <- SO %>% subset(External_Cell_Type == "Myeloid cell") # 28 "samples" (cells) (of 8994)
SO_Fibroblast <- SO %>% subset(External_Cell_Type == "Fibroblast") # 55 "samples"
SO_Malignant <- SO %>% subset(External_Cell_Type == "Malignant cell") # 235 "samples"
SO_B <- SO %>% subset(External_Cell_Type == "B cell") # 171 "samples"
SO_NK <- SO %>% subset(External_Cell_Type == "NK cell") # 248 "samples"
SO_T <- SO %>% subset(External_Cell_Type == "T cell")
SO_T # Print out info -> 8,257 "samples"
rm(SO_Myeloid, SO_Fibroblast, SO_Malignant, SO_B, SO_NK, SO_T)
