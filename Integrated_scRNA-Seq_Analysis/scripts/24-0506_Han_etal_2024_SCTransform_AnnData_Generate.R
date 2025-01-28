### Author: Matthew Aaron Loberg
### Date: May 6, 2024
### Script: 24-0506_Han_etal_2024_SCTransform_AnnData_Generate.R

# Creating SCTransformed Seurat Objects and AnnData objects for Han et al.
# I will use these to merge with other samples + run FastMNN
# This script is adapted from 23-1127_Pu_etal_SCTransform_AnnData_Generate.R
# Updates (below) through 23-1127 refer to the development of the Pu et al. script
# Updates (below) after 23-1127 (e.g., 24-0506 update) refer to the development of this Han et al. script. 

# 23-1108 Update
# NOT RUNNING SoupX
# Seurat version 4 is REQUIRED for use of SeuratDisk

# 23-1121 Update
# Running ALL Pu Samples this time ... including lymph nodes, mets, etc.
# Note, script finished + run on 23-1122, so output files are from 23-1122

# 23-1127 Update
# I am re-running this with a new Pu Create SO function that will name Cell IDs without overlapping identities

# 24-0506 Update
# I am adapting the 23-1127 Pu et al. script to the new Han et al. 2024 ATC (BRAF mutant) samples
# Paper link: https://insight.jci.org/articles/view/173712

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
source("function_scripts/24-0506_Han24_Create_SO.R") # New Create SO function that fixes the error of overlapping names
source("function_scripts/23-1009_Seurat_Basic_QC.R")
source("function_scripts/23-1025_DoubletFinder_Function.R")
source("function_scripts/23-1026_SingleR_Prediction_Function.R")
#source("function_scripts/23-1031_AmbientRNA_Function.R") # No longer running SoupX, so will NOT run this (for now)

# Read Directories in as a list
# Read Directories
Han_readdirs <- list("data_in_use/Han_etal_2024_ATC_scRNA/RAW_DATA/ATC34",
                     "data_in_use/Han_etal_2024_ATC_scRNA/RAW_DATA/ATC35",
                     "data_in_use/Han_etal_2024_ATC_scRNA/RAW_DATA/ATC36",
                     "data_in_use/Han_etal_2024_ATC_scRNA/RAW_DATA/ATC37")

# Create Seurat Object in a list for each directory
Han_ATCs <- list()
for(i in 1:length(Han_readdirs)){
  Han_ATCs[[i]] <- Han_Create_SO(readdir = Han_readdirs[[i]])
  # Note sure why orig.ident is not working in the function but I will run it myself here
  Han_ATCs[[i]]$orig.ident <- paste("Han24", substr(Han_readdirs[[i]], 46, nchar(Han_readdirs[[i]])), sep = '_')
  Han_ATCs[[i]]$Histology <- "ATC"
  Han_ATCs[[i]]$BRAFV600E <- "YES"
  Han_ATCs[[i]]$Paper <- "Han24"
}
# Will get the following warning: "Feature names cannot have underscores ('_'), replacing with dashes ('-')
# This warning is NOT a problem

# Cleaning up
rm(Han_readdirs)

# Explore Han_ATCs to make sure it is what I want - confirmed, commenting out for now
#Han_ATCs[[1]]
#Han_ATCs[[2]]$orig.ident
#Han_ATCs[[3]]$orig.ident
#Han_ATCs[[4]]$orig.ident

#### Doublet Detection ####
outputdir <- "outputs/Han_etal_2024_Analysis_Outputs/24-0506_scDblFinder/"
# Run doublet detection
for(i in 1:length(Han_ATCs)){
  Han_ATCs[[i]] <- Doublet_Detection(Han_ATCs[[i]], outputdir = paste0(outputdir, Han_ATCs[[i]]$orig.ident[1]))
}
rm(outputdir, i)

# Statistics for doublets + subsetting
Han_Doublet_Info <- list()
for(i in 1:length(Han_ATCs)){
  Han_Doublet_Info[[i]] <- table(Han_ATCs[[i]]$scDblFinder.class)
  Han_ATCs[[i]] <- Han_ATCs[[i]] %>% subset(scDblFinder.class == "singlet")
  Han_ATCs[[i]]$scDblFinder.class <- NULL
}
Han_Doublet_Info
saveRDS(Han_Doublet_Info, file = "data_in_use/Han_etal_2024_ATC_scRNA/24-0506_Han24_Doublet_Info.rds")
rm(Han_Doublet_Info, i)

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
QC_OutputDir <- "outputs/Han_etal_2024_Analysis_Outputs/24-0506_All_Sample_QC/"
for(i in 1:length(Han_ATCs)){
  Han_ATCs[[i]] <- Seurat_Basic_QC(Han_ATCs[[i]], outputdir = paste0(QC_OutputDir,Han_ATCs[[i]]$orig.ident[1]))
}
rm(QC_OutputDir, i)

#### Normalization with SCTransform ####
# Run SCTransform for all Pu PTCs using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Han_ATCs <- lapply(X = Han_ATCs,
                   FUN = SCTransform,
                   variable.features.n = 7000,
                   return.only.var.genes = FALSE,
                   vst.flavor = "v2")

### Run SingleR on the Pu_PTCs list
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Han_etal_2024_ATC_scRNA/24-0506_SingleR_Predictions/" # Set savedir for the tables with prediction info
outputdir <- "outputs/Han_etal_2024_Analysis_Outputs/24-0506_SingleR_UMAPs/" # set the outputdir for sample umaps
Han_ATCs <- SingleR_Predictions(SO_List = Han_ATCs, outputdir = outputdir, savedir = savedir)

### Add in MetaData about samples from the paper
# TERT Mutation Status
Han_ATCs[[1]]$TERT_MUT <- "YES" # This is Han24_ATC34
Han_ATCs[[2]]$TERT_MUT <- "YES" # This is Han24_ATC35
Han_ATCs[[3]]$TERT_MUT <- "YES" # This is Han24_ATC36
Han_ATCs[[4]]$TERT_MUT <- "YES" # This is Han24_ATC37

# RAS Mutation Status (note: assuming all have NO RAS Mutation as all BRAFV600E+)
Han_ATCs[[1]]$RAS_MUT <- "NO" # This is Han24_ATC34
Han_ATCs[[2]]$RAS_MUT <- "NO" # This is Han24_ATC35
Han_ATCs[[3]]$RAS_MUT <- "NO" # This is Han24_ATC36
Han_ATCs[[4]]$RAS_MUT <- "NO" # This is Han24_ATC37

# Location
Han_ATCs[[1]]$Location <- "Thyroid" # This is Han24_ATC34 (Primary tumor)
Han_ATCs[[2]]$Location <- "Thyroid" # This is Han24_ATC35 (Primary tumor)
Han_ATCs[[3]]$Location <- "Thyroid" # This is Han24_ATC36 (Primary tumor)
Han_ATCs[[4]]$Location <- "Thyroid" # This is Han24_ATC37 (Primary tumor)

# Histology_Subtype
Han_ATCs[[1]]$Histology_Subtype <- "ATC" # This is Han24_ATC34
Han_ATCs[[2]]$Histology_Subtype <- "Spindled_ATC" # This is Han24_ATC35
Han_ATCs[[3]]$Histology_Subtype <- "ATC_Focal_PTC" # This is Han24_ATC36
Han_ATCs[[4]]$Histology_Subtype <- "ATC" # This is Han24_ATC37

# TNM Stage
Han_ATCs[[1]]$TNM_Stage <- "Unknown" # This is Han24_ATC34
Han_ATCs[[2]]$TNM_Stage <- "Unknown" # This is Han24_ATC35
Han_ATCs[[3]]$TNM_Stage <- "Unknown" # This is Han24_ATC36
Han_ATCs[[4]]$TNM_Stage <- "Unknown" # This is Han24_ATC37

# Concomitant Hashimoto's
Han_ATCs[[1]]$Hashimotos <- "NO" # This is Han24_ATC34
Han_ATCs[[2]]$Hashimotos <- "NO" # This is Han24_ATC35
Han_ATCs[[3]]$Hashimotos <- "NO" # This is Han24_ATC36
Han_ATCs[[4]]$Hashimotos <- "NO" # This is Han24_ATC37

# Treatment History? (Binary: YES or NO)
Han_ATCs[[1]]$Tx_Hx_Binary <- "NO" # This is Han24_ATC34
Han_ATCs[[2]]$Tx_Hx_Binary <- "NO" # This is Han24_ATC35
Han_ATCs[[3]]$Tx_Hx_Binary <- "NO" # This is Han24_ATC36
Han_ATCs[[4]]$Tx_Hx_Binary <- "YES" # This is Han24_ATC37

# Treatment History? (Descriptive)
Han_ATCs[[1]]$Tx_Hx_Descriptive <- "None" # Han24_ATC34
Han_ATCs[[2]]$Tx_Hx_Descriptive <- "None" # Han24_ATC35
Han_ATCs[[3]]$Tx_Hx_Descriptive <- "None" # Han24_ATC36
Han_ATCs[[4]]$Tx_Hx_Descriptive <- "Famitinib+Camrelizumab" # Han24_ATC37

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Han_ATCs)){
  print(length(rownames(GetAssayData(Han_ATCs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Han_ATCs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Han_ATCs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Han_ATCs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Han_ATCs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad
# I am going to save to the PC to save space in my OneDrive. I will format this the same way that the 2023_Integrated_scRNA-Seq_Aanalysis project is formatted

# Save SCTransformed data for each of the Han 2024 ATCs
savedir <- "data_in_use/Han_etal_2024_ATC_scRNA/Processed_Data/Individual_Samples/"
for(i in 1:length(Han_ATCs)){
  savedir_temp <- paste0(savedir, Han_ATCs[[i]]$orig.ident[1])
  SaveH5Seurat(Han_ATCs[[i]], filename = paste0(savedir_temp, "/24-0506_SCTransformed.h5Seurat"), overwrite = TRUE)
  #Convert(paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE) # To save space, I will NOT be saving the AnnData object here ... can always come back and run in future
}
rm(list = ls())

sessionInfo()
