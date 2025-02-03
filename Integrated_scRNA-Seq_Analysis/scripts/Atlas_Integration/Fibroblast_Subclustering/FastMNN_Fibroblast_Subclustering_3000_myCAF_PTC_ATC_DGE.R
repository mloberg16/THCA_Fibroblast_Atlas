### Author: Matthew Aaron Loberg
### Date: December 12, 2024
### Script: FastMNN_Fibroblast_Subclustering_3000_myCAF_PTC_ATC_DGE.R
### Source Script Name: 24-1212_FastMNN_myCAF_PTC_ATC_DGE.R

# From November 16, 2024:
# Goal: Fibroblast Subclustering Differential Gene Expression Analysis

### Adapted from 24-1126_FastMNN_Fibroblast_Subclustering_3000_DGE_Analysis.R
# Goal for that script was differential gene expression analysis

### 24-1212 Update
# DGE within myCAF between PTC and ATC

### Load packages
library(Seurat)
library(SeuratWrappers)
library(tidyverse) # for ggsave
library(RColorBrewer) # For plot customization

### Load data
Merged_SO_FastMNN <- readRDS(file = "~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

### DGE Analysis within myCAFs
# Subset to myCAF
table(Merged_SO_FastMNN$CAF_Labels)
myCAF <- Merged_SO_FastMNN %>% subset(CAF_Labels == "myCAF")

# Check how many myCAFs are from paratumor/normal
table(myCAF$Histology_Simplified) # 43

# remove the 43 myCAFs from paratumor/normal
myCAF <- myCAF %>% subset(Histology_Simplified != "Paratumor/Normal")
table(myCAF$Histology_Simplified)

# Now, perform DGE analysis between PTC/ATC for myCAF - first for ATC
Idents(myCAF) <- myCAF$Histology_Simplified
myCAF_ATC_Markers <- myCAF %>% FindMarkers(ident.1 = c("ATC"), test = "MAST")
saveRDS(myCAF_ATC_Markers, file = "data_in_use/Integrated_Data/24-1212_myCAF_ATC_PTC_Markers.RDS")

myCAF_ATC_Markers_pct1 <- myCAF_ATC_Markers %>% subset(pct.1 > 0.5)



### Volcano Plot
# Load required packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

# Set thresholds and colors for for plotting
keyvals <- ifelse(
  myCAF_ATC_Markers$avg_log2FC < -1 & myCAF_ATC_Markers$p_val_adj < 0.05, '#FFCC99',
  ifelse(myCAF_ATC_Markers$avg_log2FC > 1 & myCAF_ATC_Markers$p_val_adj < 0.05, '#0075DC',
         'grey'))

keyvals[is.na(keyvals)] <- 'grey'

names(keyvals)[keyvals == '#0075DC'] <- 'Up'
names(keyvals)[keyvals == '#FFCC99'] <- 'Down'
names(keyvals)[keyvals == 'grey'] <- 'NS'

# Make volcano plot
myCAF_ATC_PTC_Volcano <- EnhancedVolcano(myCAF_ATC_Markers,
                              lab = rownames(myCAF_ATC_Markers),
                              legendPosition = "none",
                              x = 'avg_log2FC',
                              xlim = c(-4, 4),
                              y = 'p_val_adj',
                              pCutoff = 0.05,
                              FCcutoff = 1.0,
                              pointSize = 1.8,
                              colCustom = keyvals,
                              title = "myCAF PTC vs ATC",
                              subtitle = NULL,
                              titleLabSize = 30,
                              axisLabSize = 25,
                              # subtitleLabSize = 24,
                              # labsize = 10,
                              selectLab = c("FAP",
                                            "LRRC15",
                                            "RGS5",
                                            "ACTA2",
                                            "APOE",
                                            "THBS2",
                                            "CD274",
                                            "SPARCL",
                                            "ASPN",
                                            "WNT5A",
                                            "CXCL8",
                                            "CXCL6",
                                            "CXCL3",
                                            "HSPA6",
                                            "TAGLN",
                                            "CXCL1",
                                            #"HLA-A",
                                            "HLA-B",
                                            "HLA-C",
                                            "COL8A1",
                                            "COL10A1",
                                            "COL11A1",
                                            "ADIRF",
                                            "TPM2",
                                            "CD82",
                                            "ADM"),
                              boxedLabels = TRUE,
                              drawConnectors = TRUE,
                              colConnectors = 'black',
                              widthConnectors = 1,
                              labSize = 4,
                              lengthConnectors = unit(0.20, "npc"),
                              arrowheads = FALSE)
myCAF_ATC_PTC_Volcano <-
  myCAF_ATC_PTC_Volcano + theme(axis.title = element_text(face = "bold", size = 40),
                     axis.text = element_text(face = "bold", size = 30))
# ggsave the volcano plot
ggsave("outputs/Fibroblast_Subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/Volcano/myCAF_ATC_vs_PTC.png",
       width = 7,
       height = 7,
       myCAF_ATC_PTC_Volcano,
       dpi = 600)
