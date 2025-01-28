### Author: Matthew Loberg
### Script: 22-1031_Luo_ReadData.R
### Source Script Name: 22-1031_ReadData.R
### Goal: Read in scRNA-sequencing data from Luo at all and work on formatting it for analysis with Seurat and eventual integration

## Data organization
# Data is formatted as a raw count matrix for all samples: "thyroid_count_matrix.csv"
# I will take this matrix and just save it as a .RDS here for future splitting into individual objects, etc.

## Using the following tutorial for loading data
# https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/loading-your-own-data-in-seurat-reanalyze-a-different-dataset/
# The biggest challenge is that the samples are merged in this data file
# Luckily, the cell IDs are prefixed with sample IDs, so I can track individual cells to individual samples

# Load packages
library(tidyverse) # contains dplyr
library(Matrix)

# Load count table
raw_counts <- read.table(file = "data_in_use/thyroid_count_matrix.csv", sep = ",", header = TRUE)

# Save this file as an RDS for faster load times
saveRDS(raw_counts, file = "data_in_use/thyroid_raw_count_matrix.rds")
rm(raw_counts)
