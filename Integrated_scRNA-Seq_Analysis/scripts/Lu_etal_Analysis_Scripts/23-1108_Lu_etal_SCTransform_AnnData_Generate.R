### Author: Matthew Aaron Loberg
### Date: November 8, 2023
### Script: 23-1108_Lu_etal_SCTransform_AnnData_Generate.R

### Goal: Incorporate normal, PTC, and ATC scRNA-sequencing data from Lu et al. into seurat objects and save them as h5Seurat for future use and h5ad for use as AnnData objects

# In brief:
# Creating SCTransformed Seurat Objects and AnnData objects for Lu et al. normal (paratumors), PTCs, and ATCs

### Adapted from original script: 23-1016_Lu_etal_SCTransform_AnnData_Generate.R

#### 23-1031 Update ####
# Today I am re-running this script but including the following: scDblFinder, SoupX, and SingleR annotations
# Modifying to include these and re-naming everything to 23-1031
# Save location will be the PC

#### 23-1108 Update ####
# Today I am re-running WITHOUT SoupX

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
source("function_scripts/23-1009_Seurat_Basic_QC.R")
source("function_scripts/23-1025_DoubletFinder_Function.R")
source("function_scripts/23-1026_SingleR_Prediction_Function.R")
# source("function_scripts/23-1031_AmbientRNA_Function.R") # Not runnning SoupX

#### Load previously generated raw seurat objects ####
# Load data and specify metadata (e.g., Paper, Histology, orig.ident)
# Adding paper specific metadata as well
PTC_01_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814574_PTC01.RDS")
PTC_01_Lu$orig.ident <- "PTC01_Lu"
PTC_01_Lu$Paper <- "Lu"
PTC_01_Lu$Histology <- "PTC"
PTC_01_Lu$Histology_Subtype <- "Classical_PTC"
PTC_01_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
PTC_01_Lu$Age <- 33
PTC_01_Lu$Sex <- "F"

PTC_02_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814575_PTC02.RDS")
PTC_02_Lu$orig.ident <- "PTC02_Lu"
PTC_02_Lu$Paper <- "Lu"
PTC_02_Lu$Histology <- "PTC"
PTC_02_Lu$Histology_Subtype <- "Follicular_Variant_PTC"
PTC_02_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
PTC_02_Lu$Age <- 69
PTC_02_Lu$Sex <- "M"

PTC_03_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814576_PTC03.RDS")
PTC_03_Lu$orig.ident <- "PTC03_Lu"
PTC_03_Lu$Paper <- "Lu"
PTC_03_Lu$Histology <- "PTC"
PTC_03_Lu$Histology_Subtype <- "Tall_Cell_PTC"
PTC_03_Lu$Treatment_Hx <- "Dabrafenib/Trametenib_5_Months"
PTC_03_Lu$BRAFV600E <- "YES"
PTC_03_Lu$TERT_MUT <- "YES"
PTC_03_Lu$TP53_MUT <- "NO"
PTC_03_Lu$RAS_MUT <- "NO"
PTC_03_Lu$Age <- 67
PTC_03_Lu$Sex <- "M"

NORM_03_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814577_NORM03.RDS")
NORM_03_Lu$orig.ident <- "NORM03_Lu"
NORM_03_Lu$Paper <- "Lu"
NORM_03_Lu$Histology <- "Normal"
NORM_03_Lu$Histology_Subtype <- "Adjacent_Normal_Tall_Cell_PTC"
NORM_03_Lu$Treatment_Hx <- "Dabrafenib/Trametenib_5_Months"
# Normal so NOT including mutations
NORM_03_Lu$Age <- 67
NORM_03_Lu$Sex <- "M"

PTC_04_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814578_PTC04.RDS")
PTC_04_Lu$orig.ident <- "PTC04_Lu"
PTC_04_Lu$Paper <- "Lu"
PTC_04_Lu$Histology <- "PTC"
PTC_04_Lu$Histology_Subtype <- "Classical_PTC"
PTC_04_Lu$Treatment_Hx <- "Untreated"
PTC_04_Lu$BRAFV600E <- "NO"
PTC_04_Lu$TERT_MUT <- "NO"
PTC_04_Lu$TP53_MUT <- "NO"
PTC_04_Lu$RAS_MUT <- "NO"
PTC_04_Lu$Age <- 15
PTC_04_Lu$Sex <- "F"
PTC_04_Lu$Pediatric <- "YES"
PTC_04_Lu$Fusion <- "Probably"

PTC_05_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814579_PTC05.RDS")
PTC_05_Lu$orig.ident <- "PTC05_Lu"
PTC_05_Lu$Paper <- "Lu"
PTC_05_Lu$Histology <- "PTC"
PTC_05_Lu$Histology_Subtype <- "Classical_PTC"
PTC_05_Lu$Treatment_Hx <- "Untreated"
PTC_05_Lu$BRAFV600E <- "YES"
PTC_05_Lu$TERT_MUT <- "NO"
PTC_05_Lu$TP53_MUT <- "NO"
PTC_05_Lu$RAS_MUT <- "NO"
PTC_05_Lu$Age <- 23
PTC_05_Lu$Sex <- "F"

PTC_06_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814580_PTC06.RDS")
PTC_06_Lu$orig.ident <- "PTC06_Lu"
PTC_06_Lu$Paper <- "Lu"
PTC_06_Lu$Histology <- "PTC"
PTC_06_Lu$Histology_Subtype <- "Classical_PTC"
PTC_06_Lu$Treatment_Hx <- "Untreated"
PTC_06_Lu$BRAFV600E <- "YES"
PTC_06_Lu$TERT_MUT <- "NO"
PTC_06_Lu$TP53_MUT <- "NO"
PTC_06_Lu$RAS_MUT <- "NO"
PTC_06_Lu$Age <- 40
PTC_06_Lu$Sex <- "M"

PTC_07_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814581_PTC07.RDS")
PTC_07_Lu$orig.ident <- "PTC07_Lu"
PTC_07_Lu$Paper <- "Lu"
PTC_07_Lu$Histology <- "PTC"
PTC_07_Lu$Histology_Subtype <- "Classical_PTC"
PTC_07_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
PTC_07_Lu$Age <- 32
PTC_07_Lu$Sex <- "F"

NORM_07_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814582_NORM07.RDS")
NORM_07_Lu$orig.ident <- "NORM07_Lu"
NORM_07_Lu$Paper <- "Lu"
NORM_07_Lu$Histology <- "Normal"
NORM_07_Lu$Histology_Subtype <- "Adjacent_Normal_Classical_PTC"
NORM_07_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
NORM_07_Lu$Age <- 32
NORM_07_Lu$Sex <- "F"

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

ATC_09_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814584_ATC09.RDS")
ATC_09_Lu$orig.ident <- "ATC09_Lu"
ATC_09_Lu$Paper <- "Lu"
ATC_09_Lu$Histology <- "ATC"
ATC_09_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_09_Lu$Treatment_Hx <- "Untreated"
ATC_09_Lu$BRAFV600E <- "YES"
ATC_09_Lu$TERT_MUT <- "YES"
ATC_09_Lu$TP53_MUT <- "NO"
ATC_09_Lu$RAS_MUT <- "NO"
ATC_09_Lu$Age <- 51
ATC_09_Lu$Sex <- "F"

ATC_10_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814585_ATC10.RDS")
ATC_10_Lu$orig.ident <- "ATC10_Lu"
ATC_10_Lu$Paper <- "Lu"
ATC_10_Lu$Histology <- "ATC"
ATC_10_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_10_Lu$Treatment_Hx <- "Untreated"
ATC_10_Lu$BRAFV600E <- "NO"
ATC_10_Lu$TERT_MUT <- "NO"
ATC_10_Lu$TP53_MUT <- "YES"
ATC_10_Lu$RAS_MUT <- "NO"
ATC_10_Lu$Age <- 68
ATC_10_Lu$Sex <- "M"

ATC_11_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814586_ATC11.RDS")
ATC_11_Lu$orig.ident <- "ATC11_Lu"
ATC_11_Lu$Paper <- "Lu"
ATC_11_Lu$Histology <- "ATC"
ATC_11_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_11_Lu$Treatment_Hx <- "Untreated"
ATC_11_Lu$BRAFV600E <- "NO"
ATC_11_Lu$TERT_MUT <- "YES"
ATC_11_Lu$TP53_MUT <- "NO"
ATC_11_Lu$RAS_MUT <- "YES"
ATC_11_Lu$Age <- 59
ATC_11_Lu$Sex <- "M"

ATC_12_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814587_ATC12.RDS")
ATC_12_Lu$orig.ident <- "ATC12_Lu"
ATC_12_Lu$Paper <- "Lu"
ATC_12_Lu$Histology <- "ATC"
ATC_12_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_12_Lu$Treatment_Hx <- "Untreated"
ATC_12_Lu$BRAFV600E <- "NO"
ATC_12_Lu$TERT_MUT <- "YES"
ATC_12_Lu$TP53_MUT <- "YES"
ATC_12_Lu$RAS_MUT <- "YES"
ATC_12_Lu$Age <- 51
ATC_12_Lu$Sex <- "F"

ATC_13_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814588_ATC13.RDS")
ATC_13_Lu$orig.ident <- "ATC13_Lu"
ATC_13_Lu$Paper <- "Lu"
ATC_13_Lu$Histology <- "ATC"
ATC_13_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_13_Lu$Treatment_Hx <- "Untreated"
ATC_13_Lu$BRAFV600E <- "NO"
ATC_13_Lu$TERT_MUT <- "NO"
ATC_13_Lu$TP53_MUT <- "NO"
ATC_13_Lu$RAS_MUT <- "YES"
ATC_13_Lu$Age <- 69
ATC_13_Lu$Sex <- "F"

ATC_14_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814589_ATC14.RDS")
ATC_14_Lu$orig.ident <- "ATC14_Lu"
ATC_14_Lu$Paper <- "Lu"
ATC_14_Lu$Histology <- "ATC"
ATC_14_Lu$Histology_Subtype <- "Squamoid_ATC"
ATC_14_Lu$Treatment_Hx <- "Untreated"
ATC_14_Lu$BRAFV600E <- "NO"
ATC_14_Lu$TERT_MUT <- "YES"
ATC_14_Lu$TP53_MUT <- "YES"
ATC_14_Lu$RAS_MUT <- "NO"
ATC_14_Lu$Age <- 59
ATC_14_Lu$Sex <- "M"

ATC_15_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814590_ATC15.RDS")
ATC_15_Lu$orig.ident <- "ATC15_Lu"
ATC_15_Lu$Paper <- "Lu"
ATC_15_Lu$Histology <- "ATC"
ATC_15_Lu$Histology_Subtype <- "Squamoid_ATC"
ATC_15_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
ATC_15_Lu$Age <- 48
ATC_15_Lu$Sex <- "F"

ATC_17_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814591_ATC17.RDS")
ATC_17_Lu$orig.ident <- "ATC17_Lu"
ATC_17_Lu$Paper <- "Lu"
ATC_17_Lu$Histology <- "ATC"
ATC_17_Lu$Histology_Subtype <- "Spindled_ATC"
ATC_17_Lu$Treatment_Hx <- "Untreated"
ATC_17_Lu$BRAFV600E <- "NO"
ATC_17_Lu$TERT_MUT <- "YES"
ATC_17_Lu$TP53_MUT <- "YES"
ATC_17_Lu$RAS_MUT <- "YES"
ATC_17_Lu$Age <- 57
ATC_17_Lu$Sex <- "F"

ATC_18_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814592_ATC18.RDS")
ATC_18_Lu$orig.ident <- "ATC18_Lu"
ATC_18_Lu$Paper <- "Lu"
ATC_18_Lu$Histology <- "ATC"
ATC_18_Lu$Histology_Subtype <- "Unknown_ATC"
ATC_18_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
ATC_18_Lu$Age <- 65
ATC_18_Lu$Sex <- "F"

NORM_18_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814593_NORM18.RDS")
NORM_18_Lu$orig.ident <- "NORM18_Lu"
NORM_18_Lu$Paper <- "Lu"
NORM_18_Lu$Histology <- "Normal"
NORM_18_Lu$Histology_Subtype <- "Adjacent_Normal_Unknown_ATC"
NORM_18_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
NORM_18_Lu$Age <- 65
NORM_18_Lu$Sex <- "F"

# Norm 19 collected from patient with thymic cancer
# I have excluded it in the past for that reason
# I am INCLUDING it here
NORM_19_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814594_NORM19.RDS")
NORM_19_Lu$orig.ident <- "NORM19_Lu"
NORM_19_Lu$Paper <- "Lu"
NORM_19_Lu$Histology <- "Normal"
NORM_19_Lu$Histology_Subtype <- "Normal_Thyroid_Thymic_Cancer"
NORM_19_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
NORM_19_Lu$Age <- 53
NORM_19_Lu$Sex <- "F"

NORM_20_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814595_NORM20.RDS")
NORM_20_Lu$orig.ident <- "NORM20_Lu"
NORM_20_Lu$Paper <- "Lu"
NORM_20_Lu$Histology <- "Normal"
NORM_20_Lu$Histology_Subtype <- "Adjacent_Normal_Unknown_PTC"
NORM_20_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
NORM_20_Lu$Age <- 37
NORM_20_Lu$Sex <- "F"

NORM_21_Lu <- readRDS(file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814596_NORM21.RDS")
NORM_21_Lu$orig.ident <- "NORM21_Lu"
NORM_21_Lu$Paper <- "Lu"
NORM_21_Lu$Histology <- "Normal"
NORM_21_Lu$Histology_Subtype <- "Adjacent_Normal_Unknown_PTC"
NORM_21_Lu$Treatment_Hx <- "Untreated"
# Mutation status unknown so not inputing
NORM_21_Lu$Age <- 47
NORM_21_Lu$Sex <- "F"

# Make the SOs containing sequencing data from Lu et al. into a list
Lu_SOs <- list(
  PTC_01_Lu,
  PTC_02_Lu,
  PTC_03_Lu,
  NORM_03_Lu,
  PTC_04_Lu,
  PTC_05_Lu,
  PTC_06_Lu,
  PTC_07_Lu,
  NORM_07_Lu,
  ATC_08_Lu,
  ATC_09_Lu,
  ATC_10_Lu,
  ATC_11_Lu,
  ATC_12_Lu,
  ATC_13_Lu,
  ATC_14_Lu,
  ATC_15_Lu,
  ATC_17_Lu,
  ATC_18_Lu,
  NORM_18_Lu,
  NORM_19_Lu,
  NORM_20_Lu,
  NORM_21_Lu)

# Cleaning up
rm(PTC_01_Lu,
   PTC_02_Lu,
   PTC_03_Lu,
   NORM_03_Lu,
   PTC_04_Lu,
   PTC_05_Lu,
   PTC_06_Lu,
   PTC_07_Lu,
   NORM_07_Lu,
   ATC_08_Lu,
   ATC_09_Lu,
   ATC_10_Lu,
   ATC_11_Lu,
   ATC_12_Lu,
   ATC_13_Lu,
   ATC_14_Lu,
   ATC_15_Lu,
   ATC_17_Lu,
   ATC_18_Lu,
   NORM_18_Lu,
   NORM_19_Lu,
   NORM_20_Lu,
   NORM_21_Lu)

#### Doublet Detection ####
outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/23-1108_scDblFinder/"
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
saveRDS(Lu_Doublet_Info, file = "data_in_use/Lu_etal_2023_ATC_scRNA/23-1108_Lu_Doublet_Info.rds")
rm(Lu_Doublet_Info, i)

#### SoupX Ambient RNA Detection ####
#outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/23-1031_SoupX/"
#for(i in 1:length(Lu_SOs)){
#  Lu_SOs[[i]] <- AmbientRNA_Processing(Lu_SOs[[i]], outputdir = paste0(outputdir, Lu_SOs[[i]]$orig.ident[1]))
#}
#rm(outputdir, i)
# I AM NOT RUNNING SoupX THIS TIME

# QC plots + mitochondrial subsetting
QC_OutputDir <- "outputs/Lu_etal_2023_Analysis_Outputs/23-1108_All_Sample_QC/"
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
savedir <- "data_in_use/Lu_etal_2023_ATC_scRNA/23-1108_SingleR_Predictions/" # Set savedir for the tables with prediction info
outputdir <- "outputs/Lu_etal_2023_Analysis_Outputs/23-1108_SingleR_UMAPs/" # set the outputdir for sample umaps
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
  SaveH5Seurat(Lu_SOs[[i]], filename = paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), overwrite = TRUE)
  Convert(paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE)
}

rm(list = ls())
