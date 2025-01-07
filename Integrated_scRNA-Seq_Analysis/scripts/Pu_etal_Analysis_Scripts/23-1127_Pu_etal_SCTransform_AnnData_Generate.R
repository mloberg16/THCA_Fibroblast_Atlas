# Author: Matthew Aaron Loberg
# Date: November 27th, 2023
# Script: 23-1127_Pu_etal_SCTransform_AnnData_Generate.R

# Creating SCTransformed Seurat Objects and AnnData objects for Pu et al.
# I will use these to merge with other samples + run FastMNN
# NOTE: NO LONGER DOING AnnData objects - just the Seurat object for future integration

# 23-1108 Update
# NOT RUNNING SoupX
# Seurat version 4 is REQUIRED for use of SeuratDisk

# 23-1121 Update
# Running ALL Pu Samples this time ... including lymph nodes, mets, etc.
# Note, script finished + run on 23-1122, so output files are from 23-1122

# 23-1127 Update
# I am re-running this with a new Pu Create SO function that will name Cell IDs without overlapping identities
# This way when I merge all of the objects together to create an integrated dataset, none of them will need to have cell IDs get changed

### Load packages
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

# Run source function scripts
source("function_scripts/23-1127_Pu_Create_SO.R") # New Create SO function that fixes the error of overlapping names
source("function_scripts/23-1009_Seurat_Basic_QC.R")
source("function_scripts/23-1025_DoubletFinder_Function.R")
source("function_scripts/23-1026_SingleR_Prediction_Function.R")
#source("function_scripts/23-1031_AmbientRNA_Function.R") # No longer running SoupX, so will NOT run this (for now)

# Read Directories in as a list
# Read Directories
Pu_readdirs <- list("data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585102_PTC01_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585103_PTC01_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585104_PTC02_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585105_PTC02_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585106_PTC02_LeftLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585107_PTC03_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585108_PTC03_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585109_PTC03_LeftLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585110_PTC03_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585111_PTC04_SC",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585112_PTC05_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585113_PTC05_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585114_PTC05_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585115_PTC06_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585116_PTC07_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585117_PTC08_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585118_PTC08_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585119_PTC09_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585120_PTC09_P",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585121_PTC10_T",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585122_PTC10_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585123_PTC11_RightLN",
                    "data_in_use/Pu_etal_2021_PTC_scRNA/GSE184362_RAW/GSM5585124_PTC11_SC")

# Create Seurat Object in a list for each directory
Pu_PTCs <- list()
for(i in 1:length(Pu_readdirs)){
  Pu_PTCs[[i]] <- Pu_Create_SO(readdir = Pu_readdirs[[i]])
  # Note sure why orig.ident is not working in the function but I will run it myself here
  Pu_PTCs[[i]]$orig.ident <- paste("Pu", substr(Pu_readdirs[[i]], 61, nchar(Pu_readdirs[[i]])), sep = '_')
  # Pu_PTCs[[i]]$Histology <- "PTC" # Including paratumors (Normal), LNs, and SCs, so removing this for now (included in old version where I was only including PTCs)
  Pu_PTCs[[i]]$Paper <- "Pu"
}
# Will get the following warning: "Feature names cannot have underscores ('_'), replacing with dashes ('-')
# This warning is NOT a problem
rm(Pu_readdirs)

# Explore Pu_PTCs to make sure it is what I want - confirmed, commenting out for now
#Pu_PTCs[[1]]
#Pu_PTCs[[2]]$orig.ident
#Pu_PTCs[[3]]$orig.ident

#### Doublet Detection ####
outputdir <- "outputs/Pu_etal_Analysis_Outputs/23-1127_scDblFinder/"
# Run doublet detection
for(i in 1:length(Pu_PTCs)){
  Pu_PTCs[[i]] <- Doublet_Detection(Pu_PTCs[[i]], outputdir = paste0(outputdir, Pu_PTCs[[i]]$orig.ident[1]))
}
rm(outputdir, i)

# Statistics for doublets + subsetting
Pu_Doublet_Info <- list()
for(i in 1:length(Pu_PTCs)){
  Pu_Doublet_Info[[i]] <- table(Pu_PTCs[[i]]$scDblFinder.class)
  Pu_PTCs[[i]] <- Pu_PTCs[[i]] %>% subset(scDblFinder.class == "singlet")
  Pu_PTCs[[i]]$scDblFinder.class <- NULL
}
saveRDS(Pu_Doublet_Info, file = "data_in_use/Pu_etal_2021_PTC_scRNA/23-1127_Pu_Doublet_Info.rds")
rm(Pu_Doublet_Info, i)

#### COMMENTING OUT - NOT RUNNING SoupX ####
#### SoupX Ambient RNA Detection ####
#outputdir <- "outputs/Pu_etal_Analysis_Outputs/23-1031_SoupX/"
#for(i in 1:length(Pu_PTCs)){
#  Pu_PTCs[[i]] <- AmbientRNA_Processing(Pu_PTCs[[i]], outputdir = paste0(outputdir, Pu_PTCs[[i]]$orig.ident[1]))
#}
#rm(outputdir, i)
# A note on the Rho values:
# PTC01_T: 0.04;
# PTC_02_T: 0.20;
# PTC_03_T: 0.09;
# PTC_05_T: 0.10;
# PTC_08_T: 0.01;
# PTC_09_T: 0.01;
# PTC-10_T: 0.02;
# The following warning: In sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w], :  'giveCsparse' is deprecated; setting repr="T" for you
# Should NOT need to worry about the warning

# Can check and see corrected counts
# Pu_PTCs[[1]]@assays$RNA@counts # Commenting out for now - do not want this to run every time

#### Basic QC ####
# QC plots + mitochondrial subsetting
QC_OutputDir <- "outputs/Pu_etal_Analysis_Outputs/23-1127_All_Sample_QC/"
for(i in 1:length(Pu_PTCs)){
  Pu_PTCs[[i]] <- Seurat_Basic_QC(Pu_PTCs[[i]], outputdir = paste0(QC_OutputDir,Pu_PTCs[[i]]$orig.ident[1]))
}
rm(QC_OutputDir, i)

#### Normalization with SCTransform ####
# Run SCTransform for all Pu PTCs using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Pu_PTCs <- lapply(X = Pu_PTCs,
                  FUN = SCTransform,
                  variable.features.n = 7000,
                  return.only.var.genes = FALSE,
                  vst.flavor = "v2")

### Run SingleR on the Pu_PTCs list
### But First ... I would like to RUN SingleR to do some deconvolution to add to this data!
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Pu_etal_2021_PTC_scRNA/23-1127_SingleR_Predictions/" # Set savedir for the tables with prediction info
outputdir <- "outputs/Pu_etal_Analysis_Outputs/23-1127_SingleR_UMAPs/" # set the outputdir for sample umaps
Pu_PTCs <- SingleR_Predictions(SO_List = Pu_PTCs, outputdir = outputdir, savedir = savedir)

### Add in MetaData about samples from the paper
# BRAF Mutation Status
Pu_PTCs[[1]]$BRAFV600E <- "YES" # This is PTC01_T
Pu_PTCs[[2]]$BRAFV600E <- "YES_Para" # This is PTC01_P
Pu_PTCs[[3]]$BRAFV600E <- "YES" # This is PTC02_T
Pu_PTCs[[4]]$BRAFV600E <- "YES_Para" # This is PTC02_P
Pu_PTCs[[5]]$BRAFV600E <- "YES" # This is PTC02_LeftLN
Pu_PTCs[[6]]$BRAFV600E <- "NO"  # This is PTC03_T
Pu_PTCs[[7]]$BRAFV600E <- "NO_Para" # This is Pu_PTC03_P
Pu_PTCs[[8]]$BRAFV600E <- "NO" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$BRAFV600E <- "NO" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$BRAFV600E <- "YES" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$BRAFV600E <- "NO"  # This is PTC05_T
Pu_PTCs[[12]]$BRAFV600E <- "NO_Para" # This is Pu_PTC05_P
Pu_PTCs[[13]]$BRAFV600E <- "NO" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$BRAFV600E <- "YES" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$BRAFV600E <- "NO" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$BRAFV600E <- "NO"  # This is PTC08_T
Pu_PTCs[[17]]$BRAFV600E <- "NO_Para" # This is Pu_PTC08_P
Pu_PTCs[[18]]$BRAFV600E <- "YES" # This is PTC09_T
Pu_PTCs[[19]]$BRAFV600E <- "YES_Para" # This is Pu_PTC09_P
Pu_PTCs[[20]]$BRAFV600E <- "YES" # This is PTC10_T
Pu_PTCs[[21]]$BRAFV600E <- "YES" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$BRAFV600E <- "YES" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$BRAFV600E <- "YES" # This is Pu_PTC11_SC

# TERT Mutation Status (note: TERT mutations in patients 4 (subcutaneous met) and 6 (lymph node))
Pu_PTCs[[1]]$TERT_MUT <- "NO" # This is PTC01_T
Pu_PTCs[[2]]$TERT_MUT <- "NO_Para" # This is PTC01_P
Pu_PTCs[[3]]$TERT_MUT <- "NO" # This is PTC02_T
Pu_PTCs[[4]]$TERT_MUT <- "NO_Para" # This is PTC02_P
Pu_PTCs[[5]]$TERT_MUT <- "NO" # This is PTC02_LeftLN
Pu_PTCs[[6]]$TERT_MUT <- "NO"  # This is PTC03_T
Pu_PTCs[[7]]$TERT_MUT <- "NO_Para" # This is Pu_PTC03_P
Pu_PTCs[[8]]$TERT_MUT <- "NO" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$TERT_MUT <- "NO" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$TERT_MUT <- "YES" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$TERT_MUT <- "NO"  # This is PTC05_T
Pu_PTCs[[12]]$TERT_MUT <- "NO_Para" # This is Pu_PTC05_P
Pu_PTCs[[13]]$TERT_MUT <- "NO" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$TERT_MUT <- "YES" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$TERT_MUT <- "NO" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$TERT_MUT <- "NO"  # This is PTC08_T
Pu_PTCs[[17]]$TERT_MUT <- "NO_Para" # This is Pu_PTC08_P
Pu_PTCs[[18]]$TERT_MUT <- "NO" # This is PTC09_T
Pu_PTCs[[19]]$TERT_MUT <- "NO_Para" # This is Pu_PTC09_P
Pu_PTCs[[20]]$TERT_MUT <- "NO" # This is PTC10_T
Pu_PTCs[[21]]$TERT_MUT <- "NO" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$TERT_MUT <- "NO" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$TERT_MUT <- "NO" # This is Pu_PTC11_SC

# RAS Mutation Status (note: Sequenced NRAS, KRAS, HRAS - no mutations in any of them)
Pu_PTCs[[1]]$RAS_MUT <- "NO" # This is PTC01_T
Pu_PTCs[[2]]$RAS_MUT <- "NO_Para" # This is PTC01_P
Pu_PTCs[[3]]$RAS_MUT <- "NO" # This is PTC02_T
Pu_PTCs[[4]]$RAS_MUT <- "NO_Para" # This is PTC02_P
Pu_PTCs[[5]]$RAS_MUT <- "NO" # This is PTC02_LeftLN
Pu_PTCs[[6]]$RAS_MUT <- "NO"  # This is PTC03_T
Pu_PTCs[[7]]$RAS_MUT <- "NO_Para" # This is Pu_PTC03_P
Pu_PTCs[[8]]$RAS_MUT <- "NO" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$RAS_MUT <- "NO" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$RAS_MUT <- "NO" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$RAS_MUT <- "NO"  # This is PTC05_T
Pu_PTCs[[12]]$RAS_MUT <- "NO_Para" # This is Pu_PTC05_P
Pu_PTCs[[13]]$RAS_MUT <- "NO" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$RAS_MUT <- "NO" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$RAS_MUT <- "NO" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$RAS_MUT <- "NO"  # This is PTC08_T
Pu_PTCs[[17]]$RAS_MUT <- "NO_Para" # This is Pu_PTC08_P
Pu_PTCs[[18]]$RAS_MUT <- "NO" # This is PTC09_T
Pu_PTCs[[19]]$RAS_MUT <- "NO_Para" # This is Pu_PTC09_P
Pu_PTCs[[20]]$RAS_MUT <- "NO" # This is PTC10_T
Pu_PTCs[[21]]$RAS_MUT <- "NO" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$RAS_MUT <- "NO" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$RAS_MUT <- "NO" # This is Pu_PTC11_SC

# Location
Pu_PTCs[[1]]$Location <- "Thyroid" # This is PTC01_T
Pu_PTCs[[2]]$Location <- "Paratumor" # This is PTC01_P
Pu_PTCs[[3]]$Location <- "Thyroid" # This is PTC02_T
Pu_PTCs[[4]]$Location <- "Paratumor" # This is PTC02_P
Pu_PTCs[[5]]$Location <- "Lymph Node" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Location <- "Thyroid"  # This is PTC03_T
Pu_PTCs[[7]]$Location <- "Paratumor" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Location <- "Lymph Node" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Location <- "Lymph Node" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Location <- "Skin Cutaneous" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Location <- "Thyroid"  # This is PTC05_T
Pu_PTCs[[12]]$Location <- "Paratumor" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Location <- "Lymph Node" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Location <- "Lymph Node" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Location <- "Lymph Node" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Location <- "Thyroid"  # This is PTC08_T
Pu_PTCs[[17]]$Location <- "Paratumor" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Location <- "Thyroid" # This is PTC09_T
Pu_PTCs[[19]]$Location <- "Paratumor" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Location <- "Thyroid" # This is PTC10_T
Pu_PTCs[[21]]$Location <- "Lymph Node" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Location <- "Lymph Node" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Location <- "Skin Cutaneous" # This is Pu_PTC11_SC

# Histology
Pu_PTCs[[1]]$Histology <- "PTC" # This is PTC01_T
Pu_PTCs[[2]]$Histology <- "Normal" # This is PTC01_P
Pu_PTCs[[3]]$Histology <- "PTC" # This is PTC02_T
Pu_PTCs[[4]]$Histology <- "Normal" # This is PTC02_P
Pu_PTCs[[5]]$Histology <- "PTC" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Histology <- "PTC"  # This is PTC03_T
Pu_PTCs[[7]]$Histology <- "Normal" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Histology <- "PTC" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Histology <- "PTC" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Histology <- "PTC" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Histology <- "FVPTC"  # This is PTC05_T
Pu_PTCs[[12]]$Histology <- "Normal" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Histology <- "FVPTC" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Histology <- "tcPTC" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Histology <- "FVPTC" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Histology <- "PTC"  # This is PTC08_T
Pu_PTCs[[17]]$Histology <- "Normal" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Histology <- "PTC" # This is PTC09_T
Pu_PTCs[[19]]$Histology <- "Normal" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Histology <- "PTC" # This is PTC10_T
Pu_PTCs[[21]]$Histology <- "PTC" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Histology <- "PTC" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Histology <- "PTC" # This is Pu_PTC11_SC

# Histology_Subtype
Pu_PTCs[[1]]$Histology_Subtype <- "cPTC" # This is PTC01_T
Pu_PTCs[[2]]$Histology_Subtype <- "Normal" # This is PTC01_P
Pu_PTCs[[3]]$Histology_Subtype <- "cPTC" # This is PTC02_T
Pu_PTCs[[4]]$Histology_Subtype <- "Normal" # This is PTC02_P
Pu_PTCs[[5]]$Histology_Subtype <- "cPTC" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Histology_Subtype <- "cPTC"  # This is PTC03_T
Pu_PTCs[[7]]$Histology_Subtype <- "Normal" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Histology_Subtype <- "cPTC" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Histology_Subtype <- "cPTC" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Histology_Subtype <- "cPTC" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Histology_Subtype <- "FVPTC"  # This is PTC05_T
Pu_PTCs[[12]]$Histology_Subtype <- "Normal" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Histology_Subtype <- "FVPTC" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Histology_Subtype <- "tcPTC" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Histology_Subtype <- "FVPTC" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Histology_Subtype <- "cPTC"  # This is PTC08_T
Pu_PTCs[[17]]$Histology_Subtype <- "Normal" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Histology_Subtype <- "cPTC" # This is PTC09_T
Pu_PTCs[[19]]$Histology_Subtype <- "Normal" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Histology_Subtype <- "cPTC" # This is PTC10_T
Pu_PTCs[[21]]$Histology_Subtype <- "cPTC" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Histology_Subtype <- "cPTC" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Histology_Subtype <- "cPTC" # This is Pu_PTC11_SC

# TNM Stage
Pu_PTCs[[1]]$TNM_Stage <- "T4aN0M0" # This is PTC01_T
Pu_PTCs[[2]]$TNM_Stage <- "T4aN0M0" # This is PTC01_P
Pu_PTCs[[3]]$TNM_Stage <- "T4aN1bM0" # This is PTC02_T
Pu_PTCs[[4]]$TNM_Stage <- "T4aN1bM0" # This is PTC02_P
Pu_PTCs[[5]]$TNM_Stage <- "T4aN1bM0" # This is PTC02_LeftLN
Pu_PTCs[[6]]$TNM_Stage <- "T1bN1bM0"  # This is PTC03_T
Pu_PTCs[[7]]$TNM_Stage <- "T1bN1bM0" # This is Pu_PTC03_P
Pu_PTCs[[8]]$TNM_Stage <- "T1bN1bM0" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$TNM_Stage <- "T1bN1bM0" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$TNM_Stage <- "rT0N0M1" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$TNM_Stage <- "T4aN1bM1"  # This is PTC05_T
Pu_PTCs[[12]]$TNM_Stage <- "T4aN1bM1" # This is Pu_PTC05_P
Pu_PTCs[[13]]$TNM_Stage <- "T4aN1bM1" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$TNM_Stage <- "rT0an1bM0" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$TNM_Stage <- "rT0N1bM0" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$TNM_Stage <- "T4aN1bM0"  # This is PTC08_T
Pu_PTCs[[17]]$TNM_Stage <- "T4aN1bM0" # This is Pu_PTC08_P
Pu_PTCs[[18]]$TNM_Stage <- "T1bN1aM0" # This is PTC09_T
Pu_PTCs[[19]]$TNM_Stage <- "T1bN1aM0" # This is Pu_PTC09_P
Pu_PTCs[[20]]$TNM_Stage <- "T4aN1bM1" # This is PTC10_T
Pu_PTCs[[21]]$TNM_Stage <- "T4aN1bM1" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$TNM_Stage <- "rT0N1bM1" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$TNM_Stage <- "rT0N1bM1" # This is Pu_PTC11_SC

# Concomitant Hashimoto's
Pu_PTCs[[1]]$Hashimotos <- "YES" # This is PTC01_T
Pu_PTCs[[2]]$Hashimotos <- "YES" # This is PTC01_P
Pu_PTCs[[3]]$Hashimotos <- "YES" # This is PTC02_T
Pu_PTCs[[4]]$Hashimotos <- "YES" # This is PTC02_P
Pu_PTCs[[5]]$Hashimotos <- "YES" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Hashimotos <- "NO"  # This is PTC03_T
Pu_PTCs[[7]]$Hashimotos <- "NO" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Hashimotos <- "NO" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Hashimotos <- "NO" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Hashimotos <- "YES" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Hashimotos <- "NO"  # This is PTC05_T
Pu_PTCs[[12]]$Hashimotos <- "NO" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Hashimotos <- "NO" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Hashimotos <- "NO" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Hashimotos <- "NO" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Hashimotos <- "YES" # This is Pu_PTC08_T
Pu_PTCs[[17]]$Hashimotos <- "YES" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Hashimotos <- "NO" # This is PTC09_T
Pu_PTCs[[19]]$Hashimotos <- "NO" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Hashimotos <- "NO" # This is PTC10_T
Pu_PTCs[[21]]$Hashimotos <- "NO" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Hashimotos <- "NO" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Hashimotos <- "NO" # This is Pu_PTC11_SC

# Treatment History? (Binary: YES or NO)
Pu_PTCs[[1]]$Tx_Hx_Binary <- "NO" # This is PTC01_T
Pu_PTCs[[2]]$Tx_Hx_Binary <- "NO" # This is PTC01_P
Pu_PTCs[[3]]$Tx_Hx_Binary <- "NO" # This is PTC02_T
Pu_PTCs[[4]]$Tx_Hx_Binary <- "NO" # This is PTC02_P
Pu_PTCs[[5]]$Tx_Hx_Binary <- "NO" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Tx_Hx_Binary <- "NO"  # This is PTC03_T
Pu_PTCs[[7]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Tx_Hx_Binary <- "YES" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Tx_Hx_Binary <- "NO"  # This is PTC05_T
Pu_PTCs[[12]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Tx_Hx_Binary <- "YES" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Tx_Hx_Binary <- "YES" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC08_T
Pu_PTCs[[17]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Tx_Hx_Binary <- "NO" # This is PTC09_T
Pu_PTCs[[19]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Tx_Hx_Binary <- "NO" # This is PTC10_T
Pu_PTCs[[21]]$Tx_Hx_Binary <- "NO" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Tx_Hx_Binary <- "YES" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Tx_Hx_Binary <- "YES" # This is Pu_PTC11_SC

# Treatment History? (Descriptive)
Pu_PTCs[[1]]$Tx_Hx_Descriptive <- "None" # This is PTC01_T
Pu_PTCs[[2]]$Tx_Hx_Descriptive <- "None" # This is PTC01_P
Pu_PTCs[[3]]$Tx_Hx_Descriptive <- "None" # This is PTC02_T
Pu_PTCs[[4]]$Tx_Hx_Descriptive <- "None" # This is PTC02_P
Pu_PTCs[[5]]$Tx_Hx_Descriptive <- "None" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Tx_Hx_Descriptive <- "None"  # This is PTC03_T
Pu_PTCs[[7]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Tx_Hx_Descriptive <- "2x operations; 1x iodine ablation; TSH suppression" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Tx_Hx_Descriptive <- "None"  # This is PTC05_T
Pu_PTCs[[12]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Tx_Hx_Descriptive <- "Total thyroidectomy; TSH suppression" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Tx_Hx_Descriptive <- "Hemithyroidectomy; TSH suppression" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC08_T
Pu_PTCs[[17]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Tx_Hx_Descriptive <- "None" # This is PTC09_T
Pu_PTCs[[19]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Tx_Hx_Descriptive <- "None" # This is PTC10_T
Pu_PTCs[[21]]$Tx_Hx_Descriptive <- "None" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Tx_Hx_Descriptive <- "Total thyroidectomy; 3x iodine ablation; TSH suppression" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Tx_Hx_Descriptive <- "Total thyroidectomy; 3x iodine ablation; TSH suppression" # This is Pu_PTC11_SC

# Disease Characteristics
Pu_PTCs[[1]]$Disease_Characteristics <- "Tumor invading trachea, muscle, and skin" # This is PTC01_T
Pu_PTCs[[2]]$Disease_Characteristics <- "Adjacent normal of tumor invading trachea, muscle, and skin" # This is PTC01_P
Pu_PTCs[[3]]$Disease_Characteristics <- "Primary intrathyroidal tumor; lateral neck mets" # This is PTC02_T
Pu_PTCs[[4]]$Disease_Characteristics <- "Adjacent normal of intrathyroidal tumor with lateral neck mets" # This is PTC02_P
Pu_PTCs[[5]]$Disease_Characteristics <- "Lateral neck mets of primary intrathyroidal tumor" # This is PTC02_LeftLN
Pu_PTCs[[6]]$Disease_Characteristics <- "Primary intrathyroidal tumor; laternal neck mets"  # This is PTC03_T
Pu_PTCs[[7]]$Disease_Characteristics <- "Adjacent normal of intrathyroidal tumor with lateral neck mets" # This is Pu_PTC03_P
Pu_PTCs[[8]]$Disease_Characteristics <- "Lateral neck mets of primary intrathyroidal tumor" # This is Pu_PTC03_LeftLN
Pu_PTCs[[9]]$Disease_Characteristics <- "Lateral neck mets of primary intrathyroidal tumor" # This is Pu_PTC03_RightLN
Pu_PTCs[[10]]$Disease_Characteristics <- "Subcutaneous recurrence" # This is Pu_PTC04_SC
Pu_PTCs[[11]]$Disease_Characteristics <- "Tumor invading muscle, esophagus, and recurrent laryngeal nerve; lateral neck mets"  # This is PTC05_T
Pu_PTCs[[12]]$Disease_Characteristics <- "Adjacent normal of tumor invading muscle, esophagus, and recurrent laryngeal with lateral neck matastases" # This is Pu_PTC05_P
Pu_PTCs[[13]]$Disease_Characteristics <- "Lateral neck mets of tumor invading muscle, esophagus, and recurrent laryngeal nerve" # This is Pu_PTC05_RightLN
Pu_PTCs[[14]]$Disease_Characteristics <- "Recurrent right neck mets" # This is Pu_PTC06_RightLN
Pu_PTCs[[15]]$Disease_Characteristics <- "Recurrent bilateral neck mets" # This is Pu_PTC07_RightLN
Pu_PTCs[[16]]$Disease_Characteristics <- "Tumor with extrathyroidal extension; lateral neck mets" # This is Pu_PTC08_T
Pu_PTCs[[17]]$Disease_Characteristics <- "Adjacent normal of tumor with extrathyroidal extension and lateral neck matastases" # This is Pu_PTC08_P
Pu_PTCs[[18]]$Disease_Characteristics <- "Intrathyroidal tumor" # This is PTC09_T
Pu_PTCs[[19]]$Disease_Characteristics <- "Adjacent normal of intrathyroidal tumor" # This is Pu_PTC09_P
Pu_PTCs[[20]]$Disease_Characteristics <- "Tumor invading hypopharynx and thyroid cartilage; lateral neck mets; lung mets" # This is PTC10_T
Pu_PTCs[[21]]$Disease_Characteristics <- "Lateral neck mets of tumor invading hypopharynx and thyroid cartilage also with lung mets" # This is Pu_PTC10_RightLN
Pu_PTCs[[22]]$Disease_Characteristics <- "Lymph noed of tumor with subcutaneous recurrence and lung mets" # This is Pu_PTC11_RightLN
Pu_PTCs[[23]]$Disease_Characteristics <- "Subcutaneous recurrence also with lung mets" # This is Pu_PTC11_SC

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Pu_PTCs)){
  print(length(rownames(GetAssayData(Pu_PTCs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Pu_PTCs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Pu_PTCs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Pu_PTCs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Pu_PTCs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad
# I am going to save to the PC to save space in my OneDrive. I will format this the same way that the 2023_Integrated_scRNA-Seq_Aanalysis project is formatted

# Save SCTransformed data for each of the Pu PTCs
savedir <- "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/"
for(i in 1:length(Pu_PTCs)){
  savedir_temp <- paste0(savedir, Pu_PTCs[[i]]$orig.ident[1])
  SaveH5Seurat(Pu_PTCs[[i]], filename = paste0(savedir_temp, "/23-1127_SCTransformed.h5Seurat"), overwrite = TRUE)
  #Convert(paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE) # To save space, I will NOT be saving the AnnData object here ... can always come back and run in future
}
rm(list = ls())


#### SESSION INFO ####
# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
#
# Matrix products: default
#
#
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C
# [5] LC_TIME=English_United States.utf8
#
# time zone: America/Chicago
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] scDblFinder_1.14.0          SingleCellExperiment_1.22.0
# [3] SingleR_2.2.0               celldex_1.10.1
# [5] SummarizedExperiment_1.30.2 Biobase_2.60.0
# [7] GenomicRanges_1.52.1        GenomeInfoDb_1.36.4
# [9] IRanges_2.34.1              S4Vectors_0.38.2
# [11] BiocGenerics_0.46.0         MatrixGenerics_1.12.3
# [13] matrixStats_1.1.0           SoupX_1.6.2
# [15] SeuratDisk_0.0.0.9020       beepr_1.3
# [17] RColorBrewer_1.1-3          lubridate_1.9.3
# [19] forcats_1.0.0               stringr_1.5.1
# [21] dplyr_1.1.4                 purrr_1.0.2
# [23] readr_2.1.4                 tidyr_1.3.0
# [25] tibble_3.2.1                ggplot2_3.5.0
# [27] tidyverse_2.0.0             SeuratWrappers_0.2.0
# [29] SeuratObject_4.1.4          Seurat_4.4.0
#
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.0-3         bitops_1.0-7
# [3] httr_1.4.7                    tools_4.3.1
# [5] sctransform_0.4.1             utf8_1.2.4
# [7] R6_2.5.1                      ResidualMatrix_1.10.0
# [9] lazyeval_0.2.2                uwot_0.1.16
# [11] withr_2.5.2                   sp_2.1-1
# [13] prettyunits_1.2.0             gridExtra_2.3
# [15] progressr_0.14.0              cli_3.6.1
# [17] textshaping_0.3.7             spatstat.explore_3.2-5
# [19] labeling_0.4.3                spatstat.data_3.0-3
# [21] ggridges_0.5.4                pbapply_1.7-2
# [23] Rsamtools_2.16.0              systemfonts_1.0.5
# [25] scater_1.28.0                 parallelly_1.36.0
# [27] limma_3.56.2                  rstudioapi_0.15.0
# [29] RSQLite_2.3.3                 generics_0.1.3
# [31] BiocIO_1.10.0                 ica_1.0-3
# [33] spatstat.random_3.2-1         Matrix_1.6-1
# [35] ggbeeswarm_0.7.2              fansi_1.0.5
# [37] abind_1.4-5                   lifecycle_1.0.4
# [39] edgeR_3.42.4                  yaml_2.3.7
# [41] BiocFileCache_2.8.0           Rtsne_0.16
# [43] grid_4.3.1                    blob_1.2.4
# [45] dqrng_0.3.1                   promises_1.2.1
# [47] ExperimentHub_2.8.1           crayon_1.5.2
# [49] miniUI_0.1.1.1                lattice_0.21-8
# [51] beachmat_2.16.0               cowplot_1.1.1
# [53] KEGGREST_1.40.1               metapod_1.8.0
# [55] pillar_1.9.0                  rjson_0.2.21
# [57] xgboost_1.7.5.1               future.apply_1.11.0
# [59] codetools_0.2-19              leiden_0.4.3.1
# [61] glue_1.6.2                    data.table_1.14.8
# [63] remotes_2.4.2.1               vctrs_0.6.4
# [65] png_0.1-8                     gtable_0.3.4
# [67] cachem_1.0.8                  S4Arrays_1.0.6
# [69] mime_0.12                     survival_3.5-5
# [71] audio_0.1-11                  statmod_1.5.0
# [73] bluster_1.10.0                interactiveDisplayBase_1.38.0
# [75] ellipsis_0.3.2                fitdistrplus_1.1-11
# [77] ROCR_1.0-11                   nlme_3.1-162
# [79] bit64_4.0.5                   progress_1.2.2
# [81] filelock_1.0.2                RcppAnnoy_0.0.21
# [83] irlba_2.3.5.1                 vipor_0.4.5
# [85] KernSmooth_2.23-21            colorspace_2.1-0
# [87] DBI_1.1.3                     tidyselect_1.2.0
# [89] bit_4.0.5                     compiler_4.3.1
# [91] curl_5.1.0                    BiocNeighbors_1.18.0
# [93] hdf5r_1.3.8                   DelayedArray_0.26.7
# [95] plotly_4.10.3                 rtracklayer_1.60.1
# [97] scales_1.3.0                  lmtest_0.9-40
# [99] rappdirs_0.3.3                digest_0.6.33
# [101] goftest_1.2-3                 spatstat.utils_3.0-4
# [103] XVector_0.40.0                htmltools_0.5.7
# [105] pkgconfig_2.0.3               sparseMatrixStats_1.12.2
# [107] dbplyr_2.3.4                  fastmap_1.1.1
# [109] rlang_1.1.2                   htmlwidgets_1.6.2
# [111] shiny_1.8.0                   DelayedMatrixStats_1.22.6
# [113] farver_2.1.1                  zoo_1.8-12
# [115] jsonlite_1.8.7                BiocParallel_1.34.2
# [117] BiocSingular_1.16.0           RCurl_1.98-1.13
# [119] magrittr_2.0.3                scuttle_1.10.3
# [121] GenomeInfoDbData_1.2.10       patchwork_1.2.0
# [123] munsell_0.5.0                 Rcpp_1.0.11
# [125] viridis_0.6.4                 reticulate_1.34.0
# [127] stringi_1.8.1                 zlibbioc_1.46.0
# [129] MASS_7.3-60                   MAST_1.26.0
# [131] AnnotationHub_3.8.0           plyr_1.8.9
# [133] parallel_4.3.1                listenv_0.9.0
# [135] ggrepel_0.9.4                 deldir_1.0-9
# [137] Biostrings_2.68.1             splines_4.3.1
# [139] tensor_1.5                    hms_1.1.3
# [141] locfit_1.5-9.8                igraph_2.0.3
# [143] spatstat.geom_3.2-7           reshape2_1.4.4
# [145] ScaledMatrix_1.8.1            BiocVersion_3.17.1
# [147] XML_3.99-0.15                 scran_1.28.2
# [149] BiocManager_1.30.22           tzdb_0.4.0
# [151] httpuv_1.6.12                 batchelor_1.16.0
# [153] RANN_2.6.1                    polyclip_1.10-6
# [155] future_1.33.0                 scattermore_1.2
# [157] rsvd_1.0.5                    xtable_1.8-4
# [159] restfulr_0.0.15               later_1.3.1
# [161] viridisLite_0.4.2             ragg_1.2.6
# [163] memoise_2.0.1                 beeswarm_0.4.0
# [165] AnnotationDbi_1.62.2          GenomicAlignments_1.36.0
# [167] cluster_2.1.4                 timechange_0.2.0
# [169] globals_0.16.2
