### Author: Matthew Loberg
### Script: 22-1101_Luo_Create_SOs.R
### Source Script Name: 22-1101_Create_SOs_Cleaned.R
### Date: November 1, 2022
### Goal: Read in raw data files for each sample and create Seurat objects

library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)

### Now that I have separate count matrices for each sample, I will start working with them as Seurat Objects. here I will turn them into Seurat objects.

# First, I will make a function for making Seurat Objects out of raw tables and add in mitochondrial percentage
CreateSO <- function(scRawData){
  rownames(scRawData) <- scRawData$GENE # set the GENE name column equal to the rownames
  scRawData <- scRawData[ , 2:ncol(scRawData)] # remove the GENE column (first column) from the data
  SO <- CreateSeuratObject(counts = scRawData, min.cells = 3, min.genes = 200, project = "Luo_scRNAseq") # create Seurat object
  SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")
  return(SO)
}

### Create Seurat objects
# See supplementary Table 3 of Luo et al. paper for sample index
## Create Seurat Objects for ATCs
# Note: ATCs from core needle biopsies
# Create Seurat object for ATC1 (ATC_WYF_raw_counts)
ATC_WYF_raw_counts <- readRDS("data_in_use/ATC_WYF_raw_counts.RDS")
ATC_WYF_ATC1_SO <- CreateSO(ATC_WYF_raw_counts)
saveRDS(ATC_WYF_ATC1_SO, file = "data_in_use/ATC_WYF_ATC1_SO.RDS") # Save Seurat Object of ATC_WYF

# Create Seurat object for ATC2 (ATC_MSQ_raw_counts)
ATC_MSQ_raw_counts <- readRDS("data_in_use/ATC_MSQ_raw_counts.RDS")
ATC_MSQ_ATC2_SO <- CreateSO(ATC_MSQ_raw_counts)
saveRDS(ATC_MSQ_ATC2_SO, file = "data_in_use/ATC_MSQ_ATC2_SO.RDS") # Save Seurat Object of ATC_MSQ

# Create Seurat object for ATC3 (ATC_LJ_raw_counts)
ATC_LJ_raw_counts <- readRDS("data_in_use/ATC_LJ_raw_counts.RDS")
ATC_LJ_ATC3_SO <- CreateSO(ATC_LJ_raw_counts)
saveRDS(ATC_LJ_ATC3_SO, file = "data_in_use/ATC_LJ_ATC3_SO.RDS") # Save Seurat Object of ATC_LJ

## Create Suerat Objects for PTCs
# XYH1 is one sample w/ two specimens
PTC_XHY1_raw_counts <- readRDS("data_in_use/PTC_XHY1_raw_counts.RDS")
PTC_XHY1_SO <- CreateSO(PTC_XHY1_raw_counts)
saveRDS(PTC_XHY1_SO, file = "data_in_use/PTC_XHY1_SO.RDS")
PTC_XHY2_raw_counts <- readRDS("data_in_use/PTC_XHY2_raw_counts.RDS")
PTC_XHY2_SO <- CreateSO(PTC_XHY2_raw_counts)
saveRDS(PTC_XHY2_SO, file = "data_in_use/PTC_XHY2_SO.RDS")

# WJL1 (listed as YJL on supplemental table 3 of the paper) is one sample with two specimens
PTC_WJL1_raw_counts <- readRDS("data_in_use/PTC_WJL1_raw_counts.RDS")
PTC_WJL1_SO <- CreateSO(PTC_WJL1_raw_counts)
saveRDS(PTC_WJL1_SO, file = "data_in_use/PTC_WJL1_SO.RDS")
PTC_WJL2_raw_counts <- readRDS("data_in_use/PTC_WJL2_raw_counts.RDS")
PTC_WJL2_SO <- CreateSO(PTC_WJL2_raw_counts)
saveRDS(PTC_WJL2_SO, file = "data_in_use/PTC_WJL2_SO.RDS")

PTC_XTZ_raw_counts <- readRDS("data_in_use/PTC_XTZ_raw_counts.RDS")
PTC_XTZ_SO <- CreateSO(PTC_XTZ_raw_counts)
saveRDS(PTC_XTZ_SO, file = "data_in_use/PTC_XTZ_SO.RDS")

# Create Seurat Object for normal
NOM_XTZ_raw_counts <- readRDS("data_in_use/NOM_XTZ_raw_counts.RDS")
NOM_XTZ_SO <- CreateSO(NOM_XTZ_raw_counts)
saveRDS(NOM_XTZ_SO, file = "data_in_use/NOM_XTZ_SO.RDS")

# Cleaning up:
rm(list = ls())
