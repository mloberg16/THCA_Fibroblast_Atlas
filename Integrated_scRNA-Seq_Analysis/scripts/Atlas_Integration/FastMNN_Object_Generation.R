### Author: Matthew Aaron Loberg
### Date: August 19, 2024
### Script: FastMNN_Object_Generation.R

# Goal: Merge all data sets with annotations (but excluding SCTransformed data) in preparation for running FastMNN on the data

# FastMNN guide: https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html

### 24-0721 Update ###
# Attempting with a new object order:
# Lee_AT9
# Lee_AT20
# Han35
# Lu ATC09
# Lu ATC08
# Han36
# Han37
# Lee_AT13 (lots of fibroblasts, small # of tumor)
# Lee_AT16 (lots of fibroblasts)
# Lu ATC11
# Han34
# Lee AT17 (weird)
# Lu ATC15
# Lu ATC18
# Lu ATC10
# Lu ATC17
# Lu ATC12
# Lu ATC13
# Lee_PT3
# Lee_PT7
# Lee_PT8
# Lee_PT9
# Lee_PT10
# Lee_PT12
# Lee_PT5
# PTC01
# ATC14
# PTC07
# PTC05
# PTC04
# PTC06
# PTC03
# PTC02
# Luo ATCs first
# Pu PTCs

# 24-0818 Update
# Merging with NEW samples after doing an nCount_RNA >= 500 filtering step

# 24-0819 Update
# Excluding C-Cells from ATC-LJ PIOR to Atlas Integration

##### Load packages #####
library(tidyverse)
library(Seurat)
library(SeuratDisk)

##### Load readdirs for SCTransformed Seurat Objects #####
Readdirs <- list(

  # Ordered Samples
  # To replicate this analysis, it is IMPERATIVE that the sample order remains the same
  # FastMNN is sensitive to sample order
  # I tried to put the samples with the highest diversity of cell populations first to do the best job of preserving biological differences
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/AT9/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/AT20/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Han_etal_2024_ATC_scRNA/Individual_Samples/ATC35/24-0506_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC09_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC08_Lu/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Han_etal_2024_ATC_scRNA/Individual_Samples/ATC36/24-0506_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Han_etal_2024_ATC_scRNA/Individual_Samples/ATC37/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/AT13/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/AT16/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC11_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Han_etal_2024_ATC_scRNA/Individual_Samples/ATC34/24-0506_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/AT17/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC15_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC18_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC10_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC17_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC12_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC13_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT3/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT7/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT8/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT9/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT10/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT12/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/PT5/24-0625_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC01_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/ATC14_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC07_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC05_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC04_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC06_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC03_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/PTC02_Lu/23-1108_SCTransformed.h5Seurat",

  # Luo Samples
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_ATC_LJ/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_ATC_MSQ/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_ATC_WYF/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_PTC_WJL1/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_PTC_WJL2/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_PTC_XHY1/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_PTC_XHY2/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_PTC_XTZ/24-0816_SCTransformed.h5Seurat",


  # Pu Samples
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC01_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC01_P/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC02_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC02_P/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC02_LeftLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC03_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC03_P/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC03_LeftLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC03_RightLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC04_SC/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC05_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC05_P/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC05_RightLN/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC06_RightLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC07_RightLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC08_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC08_P/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC09_T/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC09_P/23-1127_SCTransformed.h5Seurat",
  #"~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC10_T/23-1127_SCTransformed.h5Seurat", # Excluding PTC10T from Pu TODAY due to quality issues
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC10_RightLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC11_RightLN/23-1127_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Pu_etal_2021_PTC_scRNA/Processed_Data/Individual_Samples/Pu_PTC11_SC/24-0816_SCTransformed.h5Seurat",

  # Luo normal
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/Luo_NOM_XTZ/24-0816_SCTransformed.h5Seurat",

  # Lu normal Samples
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM03_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM07_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM18_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM19_Lu/23-1108_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM20_Lu/23-1108_SCTransformed.h5Seurat",
  #"~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lu_etal_2023_ATC_scRNA/Processed_Data/Individual_Samples/NORM21_Lu/23-1108_SCTransformed.h5Seurat", # Excluded due to quality issues

  # Lee et al. 2024 Normal samples
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/N3-GEX/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy01/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy04/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy05/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy06/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy10/24-0816_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Lee_etal_2024_Normal_scRNA/Processed_Data/Individual_Samples/Thy15/24-0816_SCTransformed.h5Seurat",

  # Wang last because I trust it the least
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T1L/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T1R/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T2L/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T2R/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T3L/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_T3R/24-0201_SCTransformed.h5Seurat",
  "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/Wang_NT/24-0201_SCTransformed.h5Seurat"

)

##### Create a list of Seurat objects to merge from the 83 readdirs #####
merge_SOs <- list()
# Will also be adding in metadata ... read that in here
for(i in 1:length(Readdirs)){
  merge_SOs[[i]] <- LoadH5Seurat(Readdirs[[i]])
  # Recreate Seurat object without the SCTransformed info ... just raw counts and meta data
  merge_SOs[[i]] <- CreateSeuratObject(counts = merge_SOs[[i]]@assays$RNA@counts,
                                       min.cells = 5,
                                       min.features = 200,
                                       meta.data = merge_SOs[[i]]@meta.data)

}
rm(Readdirs, i)



##### Add Identifier here for eventual FastMNN #####
# Identifier ONLY exists for Luo paper samples ... I will add it to all of the other paper objects
# Identifier exists to deal with the problem of some of the Luo tumors being split into multiple samples
# "Identifier" is essentially orig.ident with an additional identifier added to tumors that had multiple samples
# The orig.ident for these is set to be identical
# Here, I will make Identifier <- orig.ident for the rest of the samples
for(i in 1:length(merge_SOs)){
  if(merge_SOs[[i]]$Paper[1] != "Luo"){ # Checks to see if first cell is in Luo Paper
    merge_SOs[[i]]$Identifier <- merge_SOs[[i]]$orig.ident
  }
  print(merge_SOs[[i]]$Identifier[1])
  print(merge_SOs[[i]]$orig.ident[1])
}

##### Removing C-Cells from Luo_ATC_LJ #####
# C-cells are a parafollicular cell population that we are NOT interested in
# C-cells are marked clearly by CALCA and CALCB
# Luo_ATC_LJ is the only sample with C-Cells
# Because of this, they integrate poorly (basically in the middle of NK/T)
# I will just exclude C-Cells from the atlas in Luo_ATC_LJ
# Luo_ATC_LJ is samples 34 of merge_SOs (see readdirs above)

# Reading in C-Cell meta data and isolating Luo_ATC_LJ meta data specifically
Luo_IA_Broad <- readRDS("data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/24-0817_Luo_IA_Broad_Meta_Data_V4.RDS")
Luo_ATC_LJ <- Luo_IA_Broad[[1]]
rm(Luo_IA_Broad)
Luo_ATC_LJ <- Luo_ATC_LJ %>% subset(Luo_ATC_LJ != "C-Cell") # Subset down to just IDs that are NOT C-Cells

# Isolate Luo_ATC_LJ from merge_SOs as Luo_ATC_LJ_SO
Luo_ATC_LJ_SO <- merge_SOs[[34]]
ncol(Luo_ATC_LJ_SO) # 682 cells prior to subsetting
length(Luo_ATC_LJ) # 653 cells that are NOT C-Cells

# Subset Luo_ATC_LJ_SO down to columns whose meta data are NOT C-Cells
Luo_ATC_LJ_SO$colnames <- colnames(Luo_ATC_LJ_SO)
Luo_ATC_LJ_SO <- Luo_ATC_LJ_SO %>% subset(colnames %in% names(Luo_ATC_LJ))
ncol(Luo_ATC_LJ_SO) # new length of Luo_ATC_LJ_SO is 653
Luo_ATC_LJ_SO$colnames <- NULL # reset colnames back to NULL

# Replace merge_SOs slot 34 with the NEW Luo_ATC_LJ_SO and clean up
merge_SOs[[34]] <- Luo_ATC_LJ_SO
ncol(merge_SOs[[34]]) # 653
rm(Luo_ATC_LJ_SO, i, Luo_ATC_LJ)

##### merge Seruat Objects for FastMNN #####
# Set the merged seurat object to the first Seurat object of merge_SOs
merged_SOs <- merge_SOs[[1]]
# Merge the first Seurat object with the remaining 83 Seurat objects
for(i in 2:length(merge_SOs)){
  merged_SOs <- merge(x = merged_SOs, y = merge_SOs[[i]])
}

# Note: this process will work different in Seurat V4 vs Seruat v5
# I am currently using Seurat v4
# In Seurat V5, the normalized data will be kept as separate "layers" that can then be used for integration

# Just printing out colnames for my own verification purposes
for(i in 1:length(merge_SOs)){
  print(colnames(merge_SOs[[i]])[1])
}

# clean up merge_SOs
rm(merge_SOs)

# This was my PREVIOUS approach to the identifier problem
# Look above prior to the merge for my NEW/EFFICIENT approach to the identifier problem
# The Identifier column contains NAs
# Fill out the NAs with the orig.ident
# I will use Identifier instead of orig.ident to do the FastMNN correction
# for(i in 1:ncol(merged_SOs)){
#   if(is.na(merged_SOs$Identifier[i])){
#     merged_SOs$Identifier[i] <- merged_SOs$orig.ident[i]
#   }
# }

# Save merged RDS
saveRDS(merged_SOs, file = "~/24-0819_Merged_SOs_scRNA_for_FastMNN.RDS")

##### SESSION INFO #####
sessionInfo()
# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] SeuratDisk_0.0.0.9020 SeuratObject_4.1.4    Seurat_4.4.0          lubridate_1.9.3       forcats_1.0.0
# [6] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.4           tidyr_1.3.0
# [11] tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.2            magrittr_2.0.3
# [6] RcppAnnoy_0.0.21       spatstat.geom_3.2-7    matrixStats_1.1.0      ggridges_0.5.4         compiler_4.3.1
# [11] png_0.1-8              vctrs_0.6.4            reshape2_1.4.4         hdf5r_1.3.8            crayon_1.5.2
# [16] pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2         utf8_1.2.4             promises_1.2.1
# [21] tzdb_0.4.0             bit_4.0.5              jsonlite_1.8.7         goftest_1.2-3          later_1.3.1
# [26] spatstat.utils_3.0-4   irlba_2.3.5.1          parallel_4.3.1         cluster_2.1.4          R6_2.5.1
# [31] ica_1.0-3              stringi_1.8.1          RColorBrewer_1.1-3     spatstat.data_3.0-3    reticulate_1.34.0
# [36] parallelly_1.36.0      lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.11            tensor_1.5
# [41] future.apply_1.11.0    zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.12          Matrix_1.6-1
# [46] splines_4.3.1          igraph_2.0.3           timechange_0.2.0       tidyselect_1.2.0       rstudioapi_0.15.0
# [51] abind_1.4-5            spatstat.random_3.2-1  codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.2-5
# [56] listenv_0.9.0          lattice_0.21-8         plyr_1.8.9             shiny_1.8.0            withr_2.5.2
# [61] ROCR_1.0-11            Rtsne_0.16             future_1.33.0          survival_3.5-5         polyclip_1.10-6
# [66] fitdistrplus_1.1-11    pillar_1.9.0           KernSmooth_2.23-21     plotly_4.10.3          generics_0.1.3
# [71] sp_2.1-1               hms_1.1.3              munsell_0.5.0          scales_1.3.0           globals_0.16.2
# [76] xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2         tools_4.3.1            data.table_1.14.8
# [81] RANN_2.6.1             leiden_0.4.3.1         cowplot_1.1.1          grid_4.3.1             colorspace_2.1-0
# [86] nlme_3.1-162           patchwork_1.2.0        cli_3.6.1              spatstat.sparse_3.0-3  fansi_1.0.5
# [91] viridisLite_0.4.2      uwot_0.1.16            gtable_0.3.4           digest_0.6.33          progressr_0.14.0
# [96] ggrepel_0.9.4          htmlwidgets_1.6.2      htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7
# [101] mime_0.12              bit64_4.0.5            MASS_7.3-60
