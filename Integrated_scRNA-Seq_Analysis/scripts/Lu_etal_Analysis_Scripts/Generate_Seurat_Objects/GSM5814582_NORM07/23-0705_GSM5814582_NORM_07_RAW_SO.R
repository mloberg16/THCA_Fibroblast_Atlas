# Author: Matthew Aaron Loberg
# Date: July 5, 2023
# Script: 23-0705_GSM5814582_NORM_07_RAW_SO.R

# Goals:
# Read in and create a Seurat object out of NORM_07
# Perform initial clustering attempt of this object
# Read in cell type labels and see how those map on

#### Load Packages ####
library(Seurat)
library(tidyverse)

#### Read in Raw Data nad Convert to Seurat Object ####
# Read data from the .txt.gz file
NORM_07_Counts <- read.table(gzfile("data_in_use/Lu_etal_2023_ATC_scRNA/GSE193581_RAW/GSM5814582_NORM07_UMI.txt.gz"))
# Notes on the data: there are 38,224 rows (genes) and 3,226 (cells)
# See that here:
nrow(NORM_07_Counts)
ncol(NORM_07_Counts)

# Turn the raw count data into a Seurat object
SO <- CreateSeuratObject(counts = NORM_07_Counts, project = "Lu et al. 2023", orig.ident = "NORM_07", min.cells = 3, min.features = 200)

# Print out SO features
SO

# Check orig.ident
head(SO$orig.ident)

# Now that SO is created, can remove raw counts
rm(NORM_07_Counts)

#### QC: selecting cells for further analysis ####
# Label percent mitochondrial
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")
# Look at the variability of this
quantile(SO$percent.mt)

# A note from Lue et al. 2023 on filtering out based on mitochondrial reads:
# Given the common problem of low viability when single cells are isolated from ATC tumor tissues,
# cells that had 30% or greater of fractions of mitochondrial expression were filtered out to remove dying cells.
# This percentage seems somewhat high compared to what is standard in the field
# Cells with < 200 or > 6000 genes detected were also removed
# In total, ~70% of cells were retained for downstream analysis

# visualize QC
QC_Plots <- VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("outputs/Lu_etal_2023_Analysis_Outputs/Individual_Sample_Analysis_Outputs/GSM5814582_NORM07/NORM07_QC_VLN.png",
       QC_Plots,
       width = 10, height = 5, dpi = 600)

# For now I will not remove anything more than what they already filtered
# I may come back and consider a more restrictive % mitochondrial cutoff

# Look at association between mitochondrial percent and nCount_RNA
plot1 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Note that there is a negative correlation here (more mitochondrial counts = fewer overall RNA counts)

# Add in Cell Type Annotation
Annotation <- read.table("data_in_use/Lu_etal_2023_ATC_scRNA/GSE193581_celltype_annotation.txt", sep = "\t")
Annotation <- Annotation %>% subset(sample.ID == "NORM07") # Not sure why there is a T in the annotation file
SO$External_Cell_Type <- Annotation$celltype

# Save the raw Seurat Object
saveRDS(SO, file = "data_in_use/Lu_etal_2023_ATC_scRNA/RAW_SOs/GSM5814582_NORM07.RDS")

# cleaning up
rm(list = ls())


