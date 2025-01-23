### Author: Matthew Aaron Loberg
### Date: November 25, 2024
### Script: FastMNN_Fibroblast_Subclustering_3000_Add_Module_Scores.R
### Source Script Name: 24-1125_FastMNN_Fibroblast_Subclustering_3000_Updated_Module_Scores_Cords_Hornburg.R

# Goal: Call AddFibroblastScores() function to add stromal subpopulation scores from published papers to my stromal subclustering seurat object

### Load packages
library(Seurat)
library(tidyverse)

### Load source script with module score addition function
# Load source
source("function_scripts/24-1125_AddFibroblastModuleScores.R")

### Read in data
# Read RDS of Stromal Subclustering Seurat Object
Merged_SO_FastMNN <- readRDS("~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

### Add scores
# Add fibroblast scores by calling the Add Fibroblast Scores function script
Merged_SO_FastMNN <- Merged_SO_FastMNN %>% AddFibroblastScores()

### Overwrite stromal subclustering object .RDS with new .RDS containing fibroblast module scores
saveRDS(Merged_SO_FastMNN, file = "~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

### Cleaning up
rm(list = ls())
