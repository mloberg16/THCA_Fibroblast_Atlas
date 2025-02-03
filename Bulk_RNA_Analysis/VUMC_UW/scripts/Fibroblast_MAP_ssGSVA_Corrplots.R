### Author: Matthew Aaron Loberg
### Date: June 5, 2024
### Script: Fibroblast_MAP_ssGSVA_Corrplots.R
### Source Script Name: 24-1213_CGRNA_Fibroblast_MAP_GSVA_Corrplot.R

### Goal: 
# Make spearman's correlation corrplots with corrplot package
# Do both for entire malignant cohort and for the ATC cohort
# For the ATC cohort, include MxIF staining data scored by VWL

### Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)
library(corrplot)

### Set WD (or open script as a project within defined directory)
# Set your working directory if needed (I do not as working within 2024_Bulk_RNA_Deconvolution_Analysis R Project)
# setwd("_") # Input in the "_" the working directory that you wish to assign

### Load Data
ClinicalData <- readRDS(file = "data_in_use/22-0606_CleanedMergedData_DESeq2NormalizedReads.rds")

# Read in stromal subclustering ssGSVA data
gsva_results_fibroblast <- readRDS(file = "data_in_use/deconvolution_scores/24-0903_Thyroid_Fibroblast_GSVA_Scores.RDS")
gsva_results_fibroblast$RNA.ID <- sub('.', '', gsva_results_fibroblast$RNA.ID) # Get rid of the 'X" at the start of GSVA_Results

# read in MAP score ssGSVA data (MAP score is from the Weiss lab 2023 Cell Genomics sequencing paper)
# Link to paper: https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00223-9
gsva_results_MAP <- readRDS(file = "data_in_use/deconvolution_scores/24-1210_CGRNA_MAP_GSVA_Score.RDS")
gsva_results_MAP$RNA.ID <- sub('.', '', gsva_results_MAP$RNA.ID) # Get rid of the 'X" at the start of GSVA_Results

#gsva_results_MAP <- gsva_results_MAP %>% dplyr::rename("MAP_gsva" = "MAP")


# Merge normalized counts and clinical data
cohort <- merge(x = ClinicalData, y = gsva_results_fibroblast, by.x = "RNA.ID", by.y = "RNA.ID")
rm(ClinicalData, gsva_results_fibroblast)

cohort <- cohort %>% merge(gsva_results_MAP)

cohort$BRS <- as.numeric(cohort$BRS)



##### MALIGNANT COHORT #####
malignant <- cohort %>% subset(Diagnosis != "FA" & 
                               Diagnosis != "HA" &
                               Diagnosis != "HT" &
                               Diagnosis != "MNG")

correlation_data <- malignant[c("myCAF", "iCAF", "APOE_CAF", "iPVCAF", "dPVCAF", "FAP", "COL1A1", "ACTA2", "RGS5", "POSTN", "APOD", "MAP")] # myCAF populations
corr_mat=cor(correlation_data,method="s") #create Spearman correlation matrix
corr_mat_p <- cor.mtest(correlation_data, method = "s", conf.level = 0.95)
write_csv(as.data.frame(corr_mat), file = "data_in_use/24-1213_malignant_corr_mat.csv", col_names = TRUE)
write_csv(as.data.frame(corr_mat_p$p), file = "data_in_use/24-1213_malignant_corr_mat_p.csv", col_names = TRUE)
outputdir <- "outputs/24-1213_CGRNA_Corrplot/"
dir.create(outputdir)

png(file = paste0(outputdir, "24-1213_malignant_Cohort_Corrplot.png"), res = 300, height = 1920, width = 1920)
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

### Create ATC specific cohort
ATCs <- cohort %>% subset(Diagnosis == "ATC")

####### ADD IN MxIF DATA + merge with ATC cohort
staining_data <- read.csv(file = "data_in_use/24-0812_PanCK_FAP_Staining_Analysis.csv")

staining_data$Overall_FAP_Pattern_Simplified <- staining_data$Overall.FAP.pattern..Fibroblast.1vs.tumor.2.vs.none.0.
staining_data$Tumor_PanCK <- staining_data$Tumor.Cell.PanCK.Intensity..0.3.
staining_data$Fibroblast_FAP <- staining_data$Fibroblast.FAP.Intensity..0.3.
staining_data$Tumor_FAP <- staining_data$Tumor.Cell.FAP.Intensity..0.3.
staining_data$Fibroblast_Pattern <- staining_data$Fibroblast.FAP.pattern..e.g...tumor.adjacent.1.vs.distance.2.vs.intemixed.3.vs.none.0.

staining_data <- staining_data[c(1:33),]
staining_data <- staining_data[,c("George.IP..", "Tumor_PanCK", "Fibroblast_FAP", "Tumor_FAP")]

staining_data <- staining_data %>% dplyr::rename("IP" = "George.IP..")

ATCs <- ATCs %>% merge(staining_data)

correlation_data <- ATCs[c("myCAF", "iCAF", "APOE_CAF", "iPVCAF", "dPVCAF", "MAP", "BRS", "Tumor_PanCK", "Fibroblast_FAP", "Tumor_FAP")] # 
corr_mat <- cor(correlation_data,method="s") #create Spearman correlation matrix

png(file = paste0(outputdir, "24-1210_ATC_Cohort_Corrplot_Staining_Data.png"), res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black",
         tl.col = "black",
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()

# v2
png(file = paste0(outputdir, "24-1210_ATC_Cohort_Corrplot_Staining_Data_lower.png"), res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "circle",
         order = "hclust", 
         #addCoef.col = "black",
         tl.col = "black",
         addrect = 3,
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()

# v2 with BRAF score
ATCs$BRAF_Score <- ATCs$BRS*-1
correlation_data <- ATCs[c("myCAF", "iCAF", "APOE_CAF", "iPVCAF", "dPVCAF", "MAP", "BRS", "BRAF_Score", "Tumor_PanCK", "Fibroblast_FAP", "Tumor_FAP")] # myCAF populations
corr_mat <- cor(correlation_data,method="s") #create Spearman correlation matrix
png(file = paste0(outputdir, "24-1210_ATC_Cohort_Corrplot_Staining_Data_V2_BRAF_Score.png"), res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "circle",
         order = "hclust", 
         #addCoef.col = "black",
         tl.col = "black",
         addrect = 3,
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()


# v2 with BRAF score + p-values
ATCs$BRAF_Score <- ATCs$BRS*-1
correlation_data <- ATCs[c("myCAF", "iCAF", "APOE_CAF", "iPVCAF", "dPVCAF", "MAP", "BRS", "BRAF_Score", "Tumor_PanCK", "Fibroblast_FAP", "Tumor_FAP")] # myCAF populations
corr_mat <- cor(correlation_data,method="s") #create Spearman correlation matrix
corr_mat_p <- cor.mtest(correlation_data, method = "s", conf.level = 0.95)
write_csv(as.data.frame(corr_mat), file = "data_in_use/24-1213_ATC_corr_mat.csv", col_names = TRUE)
write_csv(as.data.frame(corr_mat_p$p), file = "data_in_use/24-1213_ATC_corr_mat_p.csv", col_names = TRUE)

png(file = paste0(outputdir, "24-1210_ATC_Cohort_Corrplot_Staining_Data_V2_BRAF_Score_P_Values.png"), res = 300, height = 1920, width = 1920)
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
         addrect = 3,
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()
