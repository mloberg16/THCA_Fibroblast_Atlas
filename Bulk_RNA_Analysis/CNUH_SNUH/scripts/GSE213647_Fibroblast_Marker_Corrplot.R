### Author: Matthew Aaron Loberg
### Date: December 13, 2024
### Script: GSE213647_GSVA_Fibroblast_Marker_Corrplot.R
### Source Script Name: 24-1213_GSE213647_GSVA_Fibroblast_Marker_Corrplot.R

# Goal: Corrplot of ssGSVA scores for stromal populations for single-cell RNA-seq
# Will also add in MAP score as well as marker genes

# This bulk RNA sequencing data from Lee et al. was TPM normalized and subset to protein coding genes on 24-0910 by Hua-Chang

# I then downloaded Hua-Chang's TPM file 
# It is within the 2024_Bulk_RNA_Deconvolution_Analysis "data_in_use" folder saved as follows:
# "GSE213647_thyroid_cancer.TPM.protein_coding.tsv"


##### Load packages #####
library(tidyverse)

##### Load TPM data #####
# read in file + format with genes as columns, including a GEO_ID column as well for sample identification
tpm <- read.table(file = 'data_in_use/GSE213647_thyroid_cancer.TPM.protein_coding.tsv', sep = '\t', header = TRUE)
tpm_t <- t(tpm)
colnames(tpm_t) <- tpm_t[2,]
tpm_t <- tpm_t[5:nrow(tpm_t),]
tpm_t <- as.data.frame(tpm_t) # matrix format was causing issues with adding SampleID variable so I am making it into a data frame
tpm_t$GEO_ID <- rownames(tpm_t)
tpm_t$GEO_ID <- substr(tpm_t$GEO_ID, 1, 10)

# Read in meta data
GSE213647_meta <- read.csv(file = "data_in_use/GSE213647_Meta_Data.csv")
GSE213647_GSVA <- readRDS(file = "data_in_use/deconvolution_scores/24-1213_GSE213647_Thyroid_Fibroblast_GSVA_Scores_MAP.RDS")
# Change GSE_ID column to exclude extra characters
GSE213647_GSVA$GEO_ID <- substr(GSE213647_GSVA$GEO_ID, 1, 10)
# Merge meta and GSVA
GSE213647_Merged <- merge(GSE213647_meta, GSE213647_GSVA)

# Check # of diagnoses 
table(GSE213647_Merged$Diagnosis)

# Merge meta and GSVA
GSE213647_Merged <- merge(GSE213647_Merged, tpm_t)

# Cleaning up
rm(tpm, tpm_t, GSE213647_meta)

# Check # of diagnoses 
table(GSE213647_Merged$Diagnosis)



GSE213647_Merged$COL1A1 <- as.numeric(GSE213647_Merged$COL1A1)
GSE213647_Merged$FAP <- as.numeric(GSE213647_Merged$FAP)
GSE213647_Merged$ACTA2 <- as.numeric(GSE213647_Merged$ACTA2)
GSE213647_Merged$RGS5 <- as.numeric(GSE213647_Merged$RGS5)
GSE213647_Merged$POSTN <- as.numeric(GSE213647_Merged$POSTN)
GSE213647_Merged$APOD <- as.numeric(GSE213647_Merged$APOD)



correlation_data <- GSE213647_Merged[c("myCAF", "iCAF", "APOE_CAF", "iPVCAF", "dPVCAF", "FAP", "COL1A1", "ACTA2", "RGS5", "POSTN", "APOD", "MAP")] # myCAF populations
corr_mat=cor(correlation_data,method="s") #create Spearman correlation matrix
corr_mat_p <- cor.mtest(correlation_data, method = "s", conf.level = 0.95)

outputdir <- "outputs/24-1213_GSE213647_Fibroblast_Corrplot/"
dir.create(outputdir)

png(file = paste0(outputdir, "24-1213_GSE213647_Fibroblast_Corrplot.png"), res = 300, height = 1920, width = 1920)
corrplot(corr_mat, 
         p.mat = corr_mat_p$p, 
         sig.level = c(.001, .01, .05), 
         pch.cex = 0.9,
         insig = 'label_sig',
         pch.col = "black",
         method = "circle",
         order = "hclust", 
         #addCoef.col = "black",
         #tl.pos = 'd',
         tl.col = "black",
         addrect = 2,
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()



