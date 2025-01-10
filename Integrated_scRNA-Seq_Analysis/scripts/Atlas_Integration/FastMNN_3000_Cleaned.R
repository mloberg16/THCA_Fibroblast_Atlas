#### Author: Matthew Aaron Loberg
#### Date: August 19, 2024
#### Script: "FastMNN_3000_Cleaned.R"

# Goal: Integrate scRNA-Seq data sets using Fast MNN

# Notes on FastMNN
# Object integration may be dependent on order of integration (i.e., whether some subpopulations are lost)
# here is a thread on the issue:
# https://support.bioconductor.org/p/127635/
# Essentially, you would want to first merge samples that contain the highest diversity of populations ... that would prevent merging of disparate populations
# This feels highly problematic to me due to the inherent diversity of mATC, iATC, etc.

# 23-1120 Update
# See previous attempt on 23-1117
# ON 23-1117, I was attempting to do this using the most recent version of Seurat (V5 Seurat/Seurat object) and SeuratWrappers (V0.3.1)
# Today, I have deleted these versions, re-installed old versions, and am attempting using the V4 Seurat/SeuratObject and V0.2.0 SeuratWrappers


# Updates process
# Attempting with Seurat V4 and OLD Seurat wrappers
# Prior to running this, I deleted Seruat V5 ... here I will reinstall Seurat V4
# install.packages('Seurat', repos = c('https://satijalab.r-universe.dev', 'https://cloud.r-project.org'))
# The above did not work -> it installed 5.0.1, which is NOT what I want
# Here I will force install an older version:
#remotes::install_version(package = 'Seurat', version = package_version('4.4.0'))
#remotes::install_version(package = 'SeuratObject', version = package_version('4.1.4'))

# Attempt to install OLD version of Seurat wrappers
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
# This DID indeed force SeuratWrappers to install the older 0.2.0 version ... so this MAY WORK!! :)

# The following downloads the most recent version of Seurat Wrappers:
# remotes::install_github('satijalab/seurat-wrappers')

##### LOAD REQUIED PACKAGES #####
library(Seurat)
library(SeuratWrappers)
library(tidyverse) # for ggsave
library(RColorBrewer) # for custom plot coloring
library(beepr)

##### LOAD DATA FOR INTEGRATION #####
# Load previously generated Merged single-cell object
# See script: FastMNN_Object_Generation.R for merging of single-cell objects
Merged_SO <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_for_FastMNN.RDS")

##### Run FastMNN #####
Merged_SO <- NormalizeData(Merged_SO)
Merged_SO <- FindVariableFeatures(Merged_SO, nfeatures = 3000)
# There is some confusion on my end on whether NormalizeData should be run prior to or after merge
# In the Seurat CCA integration, it would be normal to do an individual SCTransform prior to integration
# It could make more sense to do a library based multi-sample normalization rather than normalizing after merging
# However, the tutorials that I have seen perform normalization AFTER merging the objects
# See tutorial here:
# https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.md
# Another point of confusion for me is how to select variable features (e.g., individual, etc.)
# There seem to be conflicting pieces of advice on this. For instance, see the following GitHub issue thread:
# https://github.com/satijalab/seurat-wrappers/issues/15
# For now, I will force FastMNN to use 3000 features and go from there
# I should try additional versions of his script with alternate # of variable feautres
# 24-0818 update - I have tried different # of features in previous atlas building iterations
# The 3000 feature iteration seemed to do the best at balancing having enough biological signal and not too much noise to isolate major cell types
Merged_SO_FastMNN <- RunFastMNN(object.list = SplitObject(Merged_SO, split.by = "Identifier"), features = 3000)
saveRDS(Merged_SO_FastMNN, file = "~/24-0819_Merged_SOs_scRNA_AFTER_FastMNN_3000.RDS")
# Note: Matrix package version needed to be 1.6-1 for RunUMAP with this old version of Seurat
# Ran the following command (after deleting the Matrix package):
# remotes::install_version("Matrix", version = "1.6-1")

##### Dimensionality Reduction Analysis #####
# Run UMAP using 30 of the mnn dimensions
Merged_SO_FastMNN <- RunUMAP(Merged_SO_FastMNN, reduction = "mnn", dims = 1:30)
# Find nearest neighbors using 30 of the mnn dimensions
Merged_SO_FastMNN <- FindNeighbors(Merged_SO_FastMNN, reduction = "mnn", dims = 1:30)
# Look at a number of DIFFERENT clustering resolutions
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN) # Defaults to 0.8 resolution
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.4)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.5)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 0.6)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 1.5)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 1.7)
Merged_SO_FastMNN <- FindClusters(Merged_SO_FastMNN, resolution = 2.0)
# Just looking at a few DimPlots to see how the integration looks
DimPlot(Merged_SO_FastMNN)
DimPlot(Merged_SO_FastMNN, group.by = c("SingleR_Myeloid", "External_Cell_Type"), ncol = 2, cols = DiscretePalette(n = 13, palette = "alphabet2"))
DimPlot(Merged_SO_FastMNN, group.by = c("SingleR_Myeloid", "External_Cell_Type", "SingleR_CD4", "SingleR_CD8"), ncol = 4)
DimPlot(Merged_SO_FastMNN, group.by = c("orig.ident", "Paper"))

##### Save Clustering Results #####
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/Clustering"
dir.create(savedir)
ggsave(file.path(savedir, "UMAP_0.4_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.4", cols = DiscretePalette(n = 25), label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_0.5_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.5", cols = DiscretePalette(n = 27), label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_0.6_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.6", cols = DiscretePalette(n = 30), label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_0.8_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.8", cols = DiscretePalette(n = 33), label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_1.5_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.1.5", cols = DiscretePalette(n = 44), label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_1.7_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.1.7", label = TRUE),
       height = 5, width = 7)
ggsave(file.path(savedir, "UMAP_2.0_Res.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.2", label = TRUE),
       height = 5, width = 7)

##### DimPlots #####
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/DimPlots"
dir.create(savedir)
# First is looking at Lu External Cell Type DimPlot and other DimPots (SingleR Dimplots, Paper, Histology, etc.)
ggsave(file.path(savedir, "Lu_External_Cell_Type.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "External_Cell_Type"),
       height = 5, width = 5)
ggsave(file.path(savedir, "SingleR_Myeloid.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "SingleR_Myeloid"),
       height = 5, width = 5)
ggsave(file.path(savedir, "SingleR_CD4.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "SingleR_CD4"),
       height = 5, width = 7)
ggsave(file.path(savedir, "SingleR_CD8.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "SingleR_CD8"),
       height = 5, width = 7)
ggsave(file.path(savedir, "SingleR_Myeloid_T.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "SingleR_Myeloid_T", cols = DiscretePalette(n = 39)),
       height = 5, width = 8)
ggsave(file.path(savedir, "Paper.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "Paper"),
       height = 5, width = 5)
ggsave(file.path(savedir, "Histology.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "Histology"),
       height = 5, width = 5)

### Make Paper Simplified for DimPlots
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$Paper
levels(Idents(Merged_SO_FastMNN))
Merged_SO_FastMNN$Paper_Simplified <- Idents(Merged_SO_FastMNN)
levels(Merged_SO_FastMNN$Paper_Simplified)
levels(Merged_SO_FastMNN$Paper_Simplified) <- c("Hong", "Han", "Lu", "Luo", "Pu", "Lee", "Wang")
levels(Merged_SO_FastMNN$Paper_Simplified)
# Save "Paper_Simplified" DimPlot
ggsave(file.path(savedir, "Paper_Simplified.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "Paper_Simplified", cols = DiscretePalette(n = 6)),
       height = 5, width = 6)


# Second, I will load in IA Broad and add it
# I will also make a simplified version of IA_Broad
# Then I will make DimPlots

# Add IA_Broad
IA_Broad_Combined <- readRDS(file = "data_in_use/Integrated_Data/24-0818_IA_Broad_Combined_V4.RDS")
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% AddMetaData(IA_Broad_Combined, col.name = "IA_Broad")
rm(IA_Broad_Combined)

# Add IA_Broad Simplified
Merged_SO_FastMNN$IA_Broad_Simplified <- Merged_SO_FastMNN$IA_Broad
levels(Merged_SO_FastMNN$IA_Broad_Simplified)
levels(Merged_SO_FastMNN$IA_Broad_Simplified) <- c("NKT", "Thyrocyte", "B_Cell", "Stromal", "Myeloid", "B_Plasma", "Endothelial", "Stromal",
                                                   "NM_Thyrocyte", "Stromal", "Endothelial", "C-Cell", "Myeloid", "B_Plasma", "Myeloid", "Unknown_AT17",
                                                   "Tumor/Stroma", "Unknown", "Immune_Prolif", "B/T_Mix")
table(Merged_SO_FastMNN$IA_Broad_Simplified)
ggsave(file.path(savedir, "IA_Broad.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "IA_Broad", cols = DiscretePalette(n = 16)),
       height = 5, width = 5)
ggsave(file.path(savedir, "IA_Broad_Simplified.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "IA_Broad_Simplified", cols = DiscretePalette(n = 10)),
       height = 5, width = 5)



##### FeaturePlots #####
# here, i will plot features of interest within the Merged_SO_FastMNN object
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/FeaturePlots"
dir.create(savedir)
gene <- list("CD3E", "CD3G", "CD4", "FOXP3", "IL2RA", "NKG7", "CD8A", # NKT markers
             "MS4A1", "JCHAIN", "IGHG1", "IGHG4", "DERL3", "CD79A", "CD19", "MZB1", "IGHD", "TCL1A", "CD27", "TNFRSF13B", "BCL6", "AICDA", "MKI67", # B cell and plasma cell markers
             "SPP1", "MRC1", "CD68", "LYZ", "ITGAM", "FOLR2", "SELENOP", "MARCO", "S100A8", "S100A9", "KIT", "LILRA4", "LAMP3", "CD14", #
             "KRT8", "TG", "TPO", "TSHR", "KRT5", "TNC", "KRT19", "FN1", "PAX8", "EPCAM", "CDK6", "CREB3L1", # Epithelial markers
             "KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8", "KRT16", "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3", "SFN", "EPCAM", # Epithelial markers from HNSCC Puram et al. 2017 (some redundancy with other epithelial markers)
             "VWF", "PECAM1", # Endothelial markers
             "FAP", "LUM", "SFRP2", "ACTA2", "RGS5", "CXCL12", "DCN", "COL1A1") # Fibroblast/Pericyte markers

for(i in 1:length(gene)){
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot.png")), # default blue color
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]])),
         height = 5, width = 5, dpi = 600)
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot_Gradient.png")), # blue/red color scale
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]]))  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
         height = 5, width = 5, dpi = 600)
  ggsave(file.path(savedir, paste0(gene[[i]], "_Integrated_Feature_Plot_Gradient_Red.png")), # red color with a 90th percentile cutoff for max
         FeaturePlot(Merged_SO_FastMNN, features = c(gene[[i]]), cols = c("gray", "red"), max.cutoff = 'q90', order = FALSE),
         height = 5, width = 5, dpi = 600)
}

##### Split to ATC Tumor/Stroma #####
# Here, I am going to look at IA Broad Simplified just within ATC/Fibroblast clusters and within NON-ATC/Fibroblast clusters
# I I am just looking at how well separated these groups are

# Plot JUST ATC/Fibroblast group
ATC_Fibroblast <- Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6 == 12 |
                                               RNA_snn_res.0.6 == 15 |
                                               RNA_snn_res.0.6 == 17 |
                                               RNA_snn_res.0.6 == 20 |
                                               RNA_snn_res.0.6 == 25)

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/DimPlots"
ggsave(file.path(savedir, "ATC_Fibroblast_IA_Broad_Simplified.png"),
       DimPlot(ATC_Fibroblast, group.by = "IA_Broad_Simplified", cols = DiscretePalette(n = 9)),
       height = 5, width = 5)
table(ATC_Fibroblast$IA_Broad_Simplified)

# Plot everything that is NOT ATC/Fibroblast
NOT_ATC_Fibroblast <- Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6 != 12 &
                                                     RNA_snn_res.0.6 != 15 &
                                                     RNA_snn_res.0.6 != 17 &
                                                     RNA_snn_res.0.6 != 20 &
                                                     RNA_snn_res.0.6 != 25)

ggsave(file.path(savedir, "NOT_ATC_Fibroblast_IA_Broad_Simplified.png"),
       DimPlot(NOT_ATC_Fibroblast, group.by = "IA_Broad_Simplified", cols = DiscretePalette(n = 10)),
       height = 5, width = 5)

# Cleaning up
rm(ATC_Fibroblast, NOT_ATC_Fibroblast, Fibroblast_IDs)

##### ASSIGN IDENTS BASED ON 0.6 RESOLUTION CLUSTERING #####
# After review with VW, I decided that 0.6 resolution was the most appropriate
# Here, I am assigning BROAD Cell Type IDs to this group
# The "ATC_Fibroblast_Mix" group will be separated by the individual sample analysis
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6
Merged_SO_FastMNN <- RenameIdents(object = Merged_SO_FastMNN,
                                  "0"  = "NK/T",
                                  "1"  = "NK/T",
                                  "3"  = "NK/T",
                                  "8"  = "NK/T",
                                  "11" = "NK/T",
                                  "22" = "NK/T",
                                  "18" = "Plasma",
                                  "23" = "Plasma",
                                  "24" = "pDC",
                                  "7"  = "Myeloid",
                                  "10" = "Myeloid",
                                  "13" = "Myeloid",
                                  "16" = "Myeloid",
                                  "9"  = "Endothelial",
                                  "21" = "Endothelial",
                                  "6"  = "B_Cell",
                                  "14" = "B_Cell",
                                  "28" = "B_Cell", # This is GC_B_Cell and Cycling B_Cell mix
                                  "5"  = "Thyrocyte",
                                  "27" = "Thyrocyte",
                                  "2"  = "PTC",
                                  "4"  = "PTC",
                                  "19" = "PTC",
                                  "26" = "PTC",
                                  "29" = "PTC",
                                  "12" = "ATC_Fibroblast_Mix",
                                  "15" = "ATC_Fibroblast_Mix",
                                  "17" = "ATC_Fibroblast_Mix",
                                  "20" = "ATC_Fibroblast_Mix",
                                  "25" = "ATC_Fibroblast_Mix")
# Cluster labels Intermediate will be the cluster labels PRIOR to adding in the ATC_Fibroblast labels
Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Intermediate <- Idents(Merged_SO_FastMNN)

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/DimPlots"
ggsave(file.path(savedir, "RNA_snn_res.0.6_ClusterLabels_Intermediate.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.6_ClusterLabels_Intermediate", cols = DiscretePalette(n = 9)),
       height = 5, width = 5)

## Now make combined labels for Fibroblasts + ATC as well (using IA_Broad_Simplified)
ATC_Fibroblast <- Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Intermediate == "ATC_Fibroblast_Mix")
NOT_ATC_Fibroblast <- Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Intermediate != "ATC_Fibroblast_Mix")

# Update ATC_Fibroblast IDs
ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final <- ATC_Fibroblast$IA_Broad_Simplified
table(ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final)
DimPlot(ATC_Fibroblast, group.by = "RNA_snn_res.0.6_ClusterLabels_Final", cols = DiscretePalette(n = 10))
Idents(ATC_Fibroblast) <- ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final
levels(ATC_Fibroblast)
new_names <- c("NK/T", "ATC", "B_Cell", "Fibroblast", "Myeloid", "Plasma",
               "Endothelial", "Thyrocyte", "NK/T")
names(new_names) <- levels(ATC_Fibroblast)
ATC_Fibroblast <- ATC_Fibroblast %>% RenameIdents(new_names)
DimPlot(ATC_Fibroblast)
ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final <- Idents(ATC_Fibroblast)
RNA_snn_res.0.6_ClusterLabels_Final <- rbind(as.data.frame(ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final) %>% dplyr::rename("RNA_snn_res.0.6_ClusterLabels_Final" = "ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Final"),
                                             as.data.frame(NOT_ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Intermediate) %>% dplyr::rename("RNA_snn_res.0.6_ClusterLabels_Final" = "NOT_ATC_Fibroblast$RNA_snn_res.0.6_ClusterLabels_Intermediate"))

Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final <- RNA_snn_res.0.6_ClusterLabels_Final

#DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.6_ClusterLabels_Final", cols = DiscretePalette(n = 10)) # Commented out as repeating in plot below

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/DimPlots"

ggsave(file.path(savedir, "RNA_snn_res.0.6_ClusterLabels_Final.png"),
       DimPlot(Merged_SO_FastMNN, group.by = "RNA_snn_res.0.6_ClusterLabels_Final", cols = DiscretePalette(n = 10), label = FALSE),
       height = 5, width = 5)

rm(ATC_Fibroblast, NOT_ATC_Fibroblast, RNA_snn_res.0.6_ClusterLabels_Final)

### Save a Fibroblast label list for subclustering
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final
table(Idents(Merged_SO_FastMNN))
Fibroblast_IDs <- colnames(Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Final == "Fibroblast"))
saveRDS(Fibroblast_IDs, file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_Fibroblast_Labels.RDS")

### Save Myeloid label list for subclustering
Myeloid_IDs <- colnames(Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Final == "Myeloid"))
saveRDS(Myeloid_IDs, file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_Myeloid_Labels.RDS")

### Save NK/T label list for subclustering
NKT_IDs <- colnames(Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Final == "NK/T"))
saveRDS(NKT_IDs, file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_NKT_Labels.RDS")

### Save ATC label list for subclustering
ATC_IDs <- colnames(Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Final == "ATC"))
saveRDS(ATC_IDs, file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_ATC_Labels.RDS")

### Save B/Plasma labels
B_Plasma_IDs <- colnames(Merged_SO_FastMNN %>% subset(RNA_snn_res.0.6_ClusterLabels_Final == "B_Cell" |
                                                      RNA_snn_res.0.6_ClusterLabels_Final == "Plasma"))
saveRDS(B_Plasma_IDs, file = "data_in_use/Integrated_Data/24-0819_3000_FastMNN_B_Plasma_Labels.RDS")

# here are ATC Features from Han et al. :
# gene <- list("CDKN2A", "COL6A1", "COL6A2", "COL6A3", "FABP5", "GNG11", "NNMT", "SRGN", "STMN1", "TPM2", "TPM4", "TUBA1A") # These are the "12 ATC gene markers"
# NOT using those right here

##### Feature Dot Plot #####
Bubble_Plot_Broad_Markers <- c("TPO", "TG", "EPCAM", #Epithelial/PTC
                               "KRT19", "FN1", "CDKN2A", "CREB3L1",# PTC/ATC
                               #"TNC", "KRT5", "CDKN2A", # ATC
                               "DCN", "COL1A1", # Fibroblast
                               "VWF", "PECAM1", # Endothelial
                               "NKG7", "CD3E", # NKT
                               "MS4A1",  "CD79A", # B-Cells
                               "MZB1", "IGHG1", # Plasma
                               "LILRA4", "CLEC4C", # pDC
                               "LYZ", "CD14") # Myeloid

# Make a feature dotplot (using prior Lu et al. myeloid markers DotPlot from 23-0726 as a reference)
# Need to reorder factor levels for Idents for this to work
# See example here: https://github.com/satijalab/seurat/issues/711
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final
Merged_SO_FastMNN@active.ident <- factor(Merged_SO_FastMNN@active.ident,
                                         levels = c("Myeloid",
                                                    "pDC",
                                                    "Plasma",
                                                    "B_Cell",
                                                    "NK/T",
                                                    "Endothelial",
                                                    "Fibroblast",
                                                    "ATC",
                                                    "PTC",
                                                    "Thyrocyte"))
Feature_DotPlot <- DotPlot(Merged_SO_FastMNN, features = Bubble_Plot_Broad_Markers, assay = "RNA", dot.scale = 5.05) +
  scale_colour_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file = "outputs/24-0819_Atlas_Integration/FastMNN_3000/DotPlots/24-0819_Broad_Cell_Type_Markers.png",
       Feature_DotPlot, width = 9, height = 3, dpi = 600, create.dir = TRUE)
rm(Feature_DotPlot)

##### Add a TDS Module Score #####
# TDS is Thyroid Differentiation Score
# TDS was first defined in the 2014 Cell TCGA paper
# See supplemental table S5F (they are the 16 genes defined here)
TDS_Genes_List <- list(c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
                         "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"))
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% Seurat::AddModuleScore(features = TDS_Genes_List,
                                                                  name = "TDS_Mod_Score")

# Save a TDS module score plot
savedir <- "outputs/24-0819_Atlas_Integration/FastMNN_3000/ModuleScorePlots"
dir.create(savedir)
ggsave(file.path(savedir, "TDS_Module_Score.png"),
       FeaturePlot(Merged_SO_FastMNN, features = c("TDS_Mod_Score1"), reduction = "umap",
                   cols = c("gray", "red"), min.cutoff = 'q10', max.cutoff = 'q90', order = FALSE),
       height = 5, width = 6, dpi = 600)

##### Add Han et al. ATC Module Score #####
# See Han et al. 2024 JCI Insight paper for the 12 genes that they indicate label ATC
ATC_Gene_List <- list(c("CDKN2A", "COL6A1", "COL6A2", "COL6A3", "FABP5", "GNG11", "NNMT", "SRGN", "STMN1", "TPM2", "TPM4", "TUBA1A"))
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% Seurat::AddModuleScore(features = ATC_Gene_List,
                                                                  name = "Han_ATC_Mod_Score")

# Save a Han ATC module score plot
savedir <- "outputs/24-0819_Atlas_Integration/FastMNN_3000/ModuleScorePlots"
dir.create(savedir)
ggsave(file.path(savedir, "Han_ATC_Module_Score.png"),
       FeaturePlot(Merged_SO_FastMNN, features = c("Han_ATC_Mod_Score1"), reduction = "umap",
                   cols = c("gray", "red"), min.cutoff = 'q10', max.cutoff = 'q90', order = FALSE),
       height = 5, width = 6, dpi = 600)

##### Module Score Violin Plots #####
### Plot fibroblast scores
module_scores <- list("TDS", "Han_ATC")

# Set save directory
savedir <- "outputs/24-0819_Atlas_Integration/FastMNN_3000/ModuleScorePlots"
dir.create(savedir)

for(i in 1:length(module_scores)){

  # labeled with CAF_Labels
  Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final
  ggsave(file.path(savedir, paste0(module_scores[[i]], "_Module_Score_Vln.png")),
         Seurat::VlnPlot(Merged_SO_FastMNN,
                         features = c(paste0(module_scores[[i]], "_Mod_Score1")),
                         cols = DiscretePalette(n = 10),
                         pt.size = 0),
         height = 5, width = 5, dpi = 600)
}


##### SAVE #####
saveRDS(Merged_SO_FastMNN, file = "~/24-0819_Merged_SOs_scRNA_AFTER_FastMNN_3000.RDS")

######## Session Info #############
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
#
# Matrix products: default
#
#
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8
#
# time zone: America/Chicago
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] beepr_1.3                   SeuratWrappers_0.2.0        RColorBrewer_1.1-3          scDblFinder_1.14.0
# [5] SingleCellExperiment_1.22.0 SingleR_2.2.0               celldex_1.10.1              SummarizedExperiment_1.30.2
# [9] Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1
# [13] S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.1.0
# [17] SoupX_1.6.2                 SeuratDisk_0.0.0.9020       SeuratObject_4.1.4          Seurat_4.4.0
# [21] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4
# [25] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1
# [29] ggplot2_3.5.0               tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.0-3         bitops_1.0-7                  httr_1.4.7                    tools_4.3.1
# [5] sctransform_0.4.1             ResidualMatrix_1.10.0         utf8_1.2.4                    R6_2.5.1
# [9] lazyeval_0.2.2                uwot_0.1.16                   withr_2.5.2                   sp_2.1-1
# [13] gridExtra_2.3                 progressr_0.14.0              cli_3.6.1                     textshaping_0.3.7
# [17] spatstat.explore_3.2-5        labeling_0.4.3                spatstat.data_3.0-3           ggridges_0.5.4
# [21] pbapply_1.7-2                 Rsamtools_2.16.0              systemfonts_1.0.5             scater_1.28.0
# [25] parallelly_1.36.0             limma_3.56.2                  rstudioapi_0.15.0             RSQLite_2.3.3
# [29] generics_0.1.3                BiocIO_1.10.0                 ica_1.0-3                     spatstat.random_3.2-1
# [33] Matrix_1.6-1                  ggbeeswarm_0.7.2              fansi_1.0.5                   abind_1.4-5
# [37] lifecycle_1.0.4               yaml_2.3.7                    edgeR_3.42.4                  BiocFileCache_2.8.0
# [41] Rtsne_0.16                    grid_4.3.1                    blob_1.2.4                    promises_1.2.1
# [45] dqrng_0.3.1                   ExperimentHub_2.8.1           crayon_1.5.2                  miniUI_0.1.1.1
# [49] lattice_0.21-8                beachmat_2.16.0               cowplot_1.1.1                 KEGGREST_1.40.1
# [53] pillar_1.9.0                  metapod_1.8.0                 rjson_0.2.21                  xgboost_1.7.5.1
# [57] future.apply_1.11.0           codetools_0.2-19              leiden_0.4.3.1                glue_1.6.2
# [61] data.table_1.14.8             remotes_2.4.2.1               vctrs_0.6.4                   png_0.1-8
# [65] gtable_0.3.4                  cachem_1.0.8                  S4Arrays_1.0.6                mime_0.12
# [69] survival_3.5-5                audio_0.1-11                  statmod_1.5.0                 bluster_1.10.0
# [73] interactiveDisplayBase_1.38.0 ellipsis_0.3.2                fitdistrplus_1.1-11           ROCR_1.0-11
# [77] nlme_3.1-162                  bit64_4.0.5                   filelock_1.0.2                RcppAnnoy_0.0.21
# [81] irlba_2.3.5.1                 vipor_0.4.5                   KernSmooth_2.23-21            colorspace_2.1-0
# [85] DBI_1.1.3                     tidyselect_1.2.0              bit_4.0.5                     compiler_4.3.1
# [89] curl_5.1.0                    BiocNeighbors_1.18.0          hdf5r_1.3.8                   DelayedArray_0.26.7
# [93] plotly_4.10.3                 rtracklayer_1.60.1            scales_1.3.0                  lmtest_0.9-40
# [97] rappdirs_0.3.3                digest_0.6.33                 goftest_1.2-3                 spatstat.utils_3.0-4
# [101] XVector_0.40.0                htmltools_0.5.7               pkgconfig_2.0.3               sparseMatrixStats_1.12.2
# [105] dbplyr_2.3.4                  fastmap_1.1.1                 rlang_1.1.2                   htmlwidgets_1.6.2
# [109] shiny_1.8.0                   DelayedMatrixStats_1.22.6     farver_2.1.1                  zoo_1.8-12
# [113] jsonlite_1.8.7                BiocParallel_1.34.2           BiocSingular_1.16.0           RCurl_1.98-1.13
# [117] magrittr_2.0.3                scuttle_1.10.3                GenomeInfoDbData_1.2.10       patchwork_1.2.0
# [121] munsell_0.5.0                 Rcpp_1.0.11                   viridis_0.6.4                 reticulate_1.34.0
# [125] stringi_1.8.1                 zlibbioc_1.46.0               MASS_7.3-60                   AnnotationHub_3.8.0
# [129] plyr_1.8.9                    parallel_4.3.1                listenv_0.9.0                 ggrepel_0.9.4
# [133] deldir_1.0-9                  Biostrings_2.68.1             splines_4.3.1                 tensor_1.5
# [137] hms_1.1.3                     locfit_1.5-9.8                igraph_2.0.3                  spatstat.geom_3.2-7
# [141] reshape2_1.4.4                ScaledMatrix_1.8.1            BiocVersion_3.17.1            XML_3.99-0.15
# [145] scran_1.28.2                  BiocManager_1.30.22           batchelor_1.16.0              tzdb_0.4.0
# [149] httpuv_1.6.12                 RANN_2.6.1                    polyclip_1.10-6               future_1.33.0
# [153] scattermore_1.2               rsvd_1.0.5                    xtable_1.8-4                  restfulr_0.0.15
# [157] later_1.3.1                   viridisLite_0.4.2             ragg_1.2.6                    memoise_2.0.1
# [161] beeswarm_0.4.0                AnnotationDbi_1.62.2          GenomicAlignments_1.36.0      cluster_2.1.4
# [165] timechange_0.2.0              globals_0.16.2
