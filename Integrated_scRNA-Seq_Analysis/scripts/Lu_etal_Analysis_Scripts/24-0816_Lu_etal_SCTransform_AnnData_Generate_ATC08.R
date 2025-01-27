### Author: Matthew Aaron Loberg
### Date: August 16, 2024
### Script: 24-0816_Lu_etal_SCTransform_AnnData_Generate.R
### Goal: Incorporate normal, PTC, and ATC scRNA-sequencing data from Lu et al. into seurat objects and save them as h5ad for use as AnnData in scanorama integration

# In brief:
# Creating SCTransformed Seurat Objects and AnnData objects for Lu et al.

#### 23-1031 Update ####
# Today I am re-running this script but including the following: scDblFinder, SoupX, and SingleR annotations
# Modifying to include these and re-naming everything to 23-1031
# Save location will be the PC

#### 23-1108 Update ####
# Today I am re-running WITHOUT SoupX

#### 24-0814 Update ####
#### I am re-running for JUST ATC08 with updated QC requiring > 500 nCount_RNA
# Also using several updated source scripts purely for output directory management
# ATC08 was the only object from Lu et al. that did NOT already meet the QC requirement of > 500 nCount_RNA

#### 24-0816 Update ####
# Running with new Basic QC that has nCount >= 500 (instead of > 500)
# Was having issues with exclusion of some cells with exactly 500 nCount

#### Load packages ####
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SoupX)
library(celldex) # Provides access to several reference datasets ("Pokedex for Cell Types) - celldex source: https://bioconductor.org/packages/3.17/data/experiment/vignettes/celldex/inst/doc/userguide.html
library(SingleR) # Load in SingleR - make sure to install the Bioconductor version as there is an old version on Github
# Packages required for scDblFinder
library(SingleCellExperiment)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
# scDblFinder
library(scDblFinder)

##### Run source function scripts ####
source("function_scripts/24-0816_Seurat_Basic_QC.R")
source("function_scripts/24-0625_DoubletFinder_Function.R")
source("function_scripts/24-0814_SingleR_Prediction_Function.R")
# source("function_scripts/23-1031_AmbientRNA_Function.R") # Not runnning SoupX

#### Load previously generated raw seurat objects ####
# Load data and specify metadata (e.g., Paper, Histology, orig.ident)
# Adding paper specific metadata as well
ATC_08_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814583_ATC08.RDS")
ATC_08_Lu$orig.ident <- "ATC08_Lu"
ATC_08_Lu$Paper <- "Lu"
ATC_08_Lu$Histology <- "ATC"
ATC_08_Lu$Histology_Subtype <- "Squamoid_ATC"
ATC_08_Lu$Treatment_Hx <- "Atezolizumab/Vermuafenib/Cobimetinib_3_Months"
ATC_08_Lu$BRAFV600E <- "YES"
ATC_08_Lu$TERT_MUT <- "NO"
ATC_08_Lu$TP53_MUT <- "YES"
ATC_08_Lu$RAS_MUT <- "NO"
ATC_08_Lu$Age <- 70
ATC_08_Lu$Sex <- "F"

# Make the SOs containing sequencing data from Lu et al. into a list
Lu_SOs <- list(ATC_08_Lu)

# Cleaning up
rm(ATC_08_Lu)

#### Doublet Detection ####
outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/24-0816_scDblFinder/"
dir.create(outputdir)
# Run doublet detection
for(i in 1:length(Lu_SOs)){
  Lu_SOs[[i]] <- Doublet_Detection(Lu_SOs[[i]], outputdir = paste0(outputdir, Lu_SOs[[i]]$orig.ident[1]))
}
rm(outputdir, i)

# Statistics for doublets + subsetting
Lu_Doublet_Info <- list()
for(i in 1:length(Lu_SOs)){
  Lu_Doublet_Info[[i]] <- table(Lu_SOs[[i]]$scDblFinder.class)
  Lu_SOs[[i]] <- Lu_SOs[[i]] %>% subset(scDblFinder.class == "singlet")
  Lu_SOs[[i]]$scDblFinder.class <- NULL
}
saveRDS(Lu_Doublet_Info, file = "data_in_use/Lu_etal_2023_ATC_scRNA/24-0816_Lu_Doublet_Info.rds")
rm(Lu_Doublet_Info, i)

#### SoupX Ambient RNA Detection ####
#outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/23-1031_SoupX/"
#for(i in 1:length(Lu_SOs)){
#  Lu_SOs[[i]] <- AmbientRNA_Processing(Lu_SOs[[i]], outputdir = paste0(outputdir, Lu_SOs[[i]]$orig.ident[1]))
#}
#rm(outputdir, i)
# I AM NOT RUNNING SoupX THIS TIME

# QC plots + mitochondrial subsetting
QC_OutputDir <- "outputs/Lu_etal_2023_Analysis_Outputs/24-0816_All_Sample_QC/"
dir.create(QC_OutputDir)
for(i in 1:length(Lu_SOs)){
  Lu_SOs[[i]] <- Seurat_Basic_QC(Lu_SOs[[i]], outputdir = paste0(QC_OutputDir,Lu_SOs[[i]]$orig.ident[1]))
}

# Run SCTransform for all Lu_Sos using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Lu_SOs <- lapply(X = Lu_SOs,
                 FUN = SCTransform,
                 variable.features.n = 7000,
                 return.only.var.genes = FALSE,
                 vst.flavor = "v2")

#### SingleR Annotations ####
### Run SingleR on the Lu_SOs list
### But First ... I would like to RUN SingleR to do some deconvolution to add to this data!
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Lu_etal_2023_ATC_scRNA/24-0816_SingleR_Predictions/" # Set savedir for the tables with prediction info
dir.create(savedir)
outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/24-0816_SingleR_UMAPs/" # set the outputdir for sample umaps
dir.create(outputdir)
Lu_SOs <- SingleR_Predictions(SO_List = Lu_SOs, outputdir = outputdir, savedir = savedir)

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Lu_SOs)){
  print(length(rownames(GetAssayData(Lu_SOs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Lu_SOs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Lu_SOs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Lu_SOs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Lu_SOs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad

# Save SCTransformed data for each of the Lu_SOs
savedir <- "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/"
for(i in 1:length(Lu_SOs)){
  savedir_temp <- paste0(savedir, Lu_SOs[[i]]$orig.ident[1])
  SaveH5Seurat(Lu_SOs[[i]], filename = paste0(savedir_temp, "/24-0816_SCTransformed.h5Seurat"), overwrite = TRUE)
  #Convert(paste0(savedir_temp, "/24-0816_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE) # do not need h5ad right now
}

rm(list = ls())
