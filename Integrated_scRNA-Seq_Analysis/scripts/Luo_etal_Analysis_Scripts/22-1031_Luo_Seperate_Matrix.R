### Author: Matthew Loberg
### Script: 22-1031_Luo_Seperate_Matrix.R
### Source Script Name: 22-1031_Seperate_Matrix.R
### Date: October 31, 2022
### Goal: Read in the matrix (saved as an RDS in 22-1021_Luo_ReadData.R) and seperate it into its components

### Info
# I am able to split this matrix into samples by cell ID because the cells are prefixed with the sample ID
# For instance, I can select columns to include that have "ATC_LJ" in the Cell ID for the ATC_LJ sample. 

### Load packages
library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)

### Load matrix as saved in 22-1031_ReadData.R
raw_counts <- readRDS(file = "data_in_use/thyroid_raw_count_matrix.rds")

# Change the name of column X to "GENE"
raw_counts <- raw_counts %>% dplyr::rename("GENE" = "X")

# Print out the first few samples within raw_counts
head(raw_counts)

### Sort out samples containing ATC_LJ. Note: also add GENE so that the GENE column is NOT removed
ATC_LJ_raw_counts <- raw_counts[, grepl("ATC_LJ|GENE", names(raw_counts))]

# Check to make sure sorted out
head(ATC_LJ_raw_counts)

# Remove ATC_LJ from the raw counts and make this as a new restricted matrix
raw_counts_restricted <- raw_counts %>% select(-contains("ATC_LJ"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is ATC_MSQ

### Sourt out samples containing ATC_MSQ
ATC_MSQ_raw_counts <- raw_counts[, grepl("ATC_MSQ|GENE", names(raw_counts))]

# Check to make sure sorted out
head(ATC_MSQ_raw_counts)

# Remove ATC_MSQ from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("ATC_MSQ"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is ATC_WYF

### Sourt out samples containing ATC_WYF
ATC_WYF_raw_counts <- raw_counts[, grepl("ATC_WYF|GENE", names(raw_counts))]

# Check to make sure sorted out
head(ATC_WYF_raw_counts)

# Remove ATC_WYF from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("ATC_WYF"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is NOM_XTZ

### Sourt out samples containing NOM_XTZ
NOM_XTZ_raw_counts <- raw_counts[, grepl("NOM_XTZ|GENE", names(raw_counts))]

# Check to make sure sorted out
head(NOM_XTZ_raw_counts)

# Remove NOM_XTZ from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("NOM_XTZ"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is PTC_WJL1

### Sort out samples containing PTC_WJL1
PTC_WJL1_raw_counts <- raw_counts[, grepl("PTC_WJL1|GENE", names(raw_counts))]

# Check to make sure sorted out
head(PTC_WJL1_raw_counts)

# Remove PTC_WJL1 from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("PTC_WJL1"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is PTC_WJL2

### Sort out samples containing PTC_WJL2
PTC_WJL2_raw_counts <- raw_counts[, grepl("PTC_WJL2|GENE", names(raw_counts))]

# Check to make sure sorted out
head(PTC_WJL2_raw_counts)

# Remove PTC_WJL2 from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("PTC_WJL2"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is PTC_XHY1

### Sort out samples containing PTC_XHY1
PTC_XHY1_raw_counts <- raw_counts[, grepl("PTC_XHY1|GENE", names(raw_counts))]

# Check to make sure sorted out
head(PTC_XHY1_raw_counts)

# Remove PTC_XHY1 from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("PTC_XHY1"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is PTC_PTC_XHY2

### Sort out samples containing PTC_XHY2
PTC_XHY2_raw_counts <- raw_counts[, grepl("PTC_XHY2|GENE", names(raw_counts))]

# Check to make sure sorted out
head(PTC_XHY2_raw_counts)

# Remove PTC_XHY2 from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("PTC_XHY2"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is PTC_XTZ

### Sort out samples containing PTC_XTZ
PTC_XTZ_raw_counts <- raw_counts[, grepl("PTC_XTZ|GENE", names(raw_counts))]

# Check to make sure sorted out
head(PTC_XTZ_raw_counts)

# Remove PTC_XTZ from the raw counts restricted matrix 
raw_counts_restricted <- raw_counts_restricted %>% select(-contains("PTC_XTZ"))

# Check to see what next name is 
head(raw_counts_restricted) # Next name is just the gene names column...all done! 

### Save individual files as RDS objects
saveRDS(ATC_LJ_raw_counts, "data_in_use/ATC_LJ_raw_counts.RDS") # This will be ATC1
saveRDS(ATC_MSQ_raw_counts, "data_in_use/ATC_MSQ_raw_counts.RDS") # This will be ATC2
saveRDS(ATC_WYF_raw_counts, "data_in_use/ATC_WYF_raw_counts.RDS") # This will be ATC3
saveRDS(PTC_WJL1_raw_counts, "data_in_use/PTC_WJL1_raw_counts.RDS") # This will be PTC1 (part 1)
saveRDS(PTC_WJL2_raw_counts, "data_in_use/PTC_WJL2_raw_counts.RDS") # This is PTC1 (part 2)
saveRDS(PTC_XHY1_raw_counts, "data_in_use/PTC_XHY1_raw_counts.RDS") # This will be PTC2 (part 1)
saveRDS(PTC_XHY2_raw_counts, "data_in_use/PTC_XHY2_raw_counts.RDS") # This is PTC2 (part 2)
saveRDS(PTC_XTZ_raw_counts, "data_in_use/PTC_XTZ_raw_counts.RDS") # This will be PTC3
saveRDS(NOM_XTZ_raw_counts, "data_in_use/NOM_XTZ_raw_counts.RDS") # This is normal thyroid

# cleaning up 
rm(list = ls())
