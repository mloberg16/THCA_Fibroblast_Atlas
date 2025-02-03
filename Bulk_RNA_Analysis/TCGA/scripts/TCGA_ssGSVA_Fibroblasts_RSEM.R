### Author: Matthew Aaron Loberg
### Date: November 14, 2024
### Script: TCGA_ssGSEA_Fibroblasts_RSEM.R
### Source Script Name: 24-1114_TCGA_ssGSEA_Fibroblasts_RSEM.R

# Goal: Run  ssGSEA for TCGA PTCs
# In this script, I will run it for my INTERNAL CAF scores that were generated from scRNA-seq data

# Going to install the following package: 
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#1_Quick_start
# package is GSVA
# Install as below (will comment out after installed): 
# BiocManager::install("GSVA")

# Additional Info 
# This is a copy over of 24-0911_GSE213647_ssGSEA_Fibroblasts_TPM_EXTERNAL.R
# I am using that as a template for this script on TCGA data

##### 24-1018 UPDATE 
# Adding in "TCGA_PTC_fibrosis_tsv.txt" from Maxime Tarabichi

##### 24-1114 Update
# Looking at CXCL12, POSTN, FAP, ACTA2, RGS5, APOE, others by BRS status - FAP and POSTN making it into the paper Fig. 3

##### Load packages #####
library(GSVA)
library(tidyverse)

##### Load TPM data #####
# read in file 
TCGA_RSEM <- read.table(file = 'data_in_use/data_RNA_Seq_v2_expression_median.txt', sep = '\t', header = TRUE)
TCGA_RSEM_rownames <- TCGA_RSEM$Hugo_Symbol # save the gene names columns as the rownames
TCGA_RSEM_matrix <- as.matrix(TCGA_RSEM[,3:ncol(TCGA_RSEM)]) # make into matrix with just the numerical data
rownames(TCGA_RSEM_matrix) <- TCGA_RSEM_rownames # set the rownames of t_matrix to the stored gene names
rm(TCGA_RSEM, TCGA_RSEM_rownames) # remove TCGA_RSEM and TCGA_RSEM_rownames (only keeping the matrix)

# My Gene Sets
geneset_readdirs <- list("data_in_use/gene_sets/24-0903_APOE_PVL.txt",
                         "data_in_use/gene_sets/24-0903_vSMC.txt",
                         "data_in_use/gene_sets/24-0903_Pericyte.txt",
                         "data_in_use/gene_sets/24-0903_myCAF.txt",
                         "data_in_use/gene_sets/24-0903_iCAF.txt")

gene_sets <- list()
for(i in 1:length(geneset_readdirs)){
  placeholder <- data.table::fread(geneset_readdirs[[i]], header = FALSE, sep = NULL)
  gene_sets[[i]] <- placeholder$V1
}
names(gene_sets) <- c("APOE_PVL",
                      "vSMC",
                      "Pericyte",
                      "myCAF",
                      "iCAF")

# Followinng this tutorial to set up ssGSEA: 
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#1_Quick_start

# First make the param object
#gsvaPar <- gsvaParam(fpkm_matrix, gene_sets)
#gsvaPar <- GSVA::gsvaParam(fpkm_matrix, gene_sets)

# Next run GSVA using method = "ssgsea"; 24-0903 change - just running default paramaters
gsva_results  <- GSVA::gsva(expr = TCGA_RSEM_matrix,
                            gset.idx.list = gene_sets,
                            #method = "ssgsea", # NOT RUNNING ssgsea, just default GSVA this time
                            #ssgsea.norm = TRUE,# NOT RUNNING ssgsea, just default GSVA this time
                            verbose = TRUE)
gsva_results_transformed <- t(gsva_results)
TCGA_ID <- rownames(gsva_results_transformed)
gsva_results_transformed <- gsva_results_transformed %>% as_tibble()
gsva_results_transformed$TCGA_ID <- TCGA_ID
rm(TCGA_ID)

# Save the GSVA Results
saveRDS(gsva_results_transformed, file = "data_in_use/deconvolution_scores/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Scores.RDS")

# Cleaning Up 
rm(list = ls())

##### Visualization #####

### Data formatting ###
# Read in TCGA meta data 
# Load Data
TCGA_Patient_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE) # Load patient data file
TCGA_Patient_Data <- as_tibble(TCGA_Patient_Data) # Make patient data into tibble for easier data cleaning

TCGA_Sample_Data <- read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE) # Load sample data file
TCGA_Sample_Data <- as_tibble(TCGA_Sample_Data) # Make sample data into tibble for easier data cleaning

# Restrict data to the columns needed
# This probably isn't entirely necessary, but I'm doing it to familiarize myself with the variables and think about the things that would be worthwhile to look at
# For patient data, I want the following columns: Patient ID, Tumor Status, Overall survival (OS) status, OS months, disease free status (DFS), DFS months
# TCGA_Patient_Data_Restricted <- TCGA_Patient_Data[c("PATIENT_ID", "TUMOR_STATUS", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS")]
# 22-0405 update: Not sure yet what columns I want. I need to think of which columns might be predictors of aggressive disease.
# For now working with all columns/non-restricted file

# For sample data, I want the following columns:
# 22-0405 update: Not sure yet what columns I want. I need to think of which columns might be predictors of aggressive disease.
# TCGA_Sample_Data_Restricted <- TCGA_Sample_Data[c("PATIENT_ID", "SAMPLE_ID", "BRAFV600E_RAS", "EXTRATHYROIDAL_EXTENSION", "DIFFERENTIATION_SCORE", "ERK_SCORE", "DISEASE_STAGE", "BRAFV600E_RAF_SCORE", "RNASEQ_DATA", "SAMPLE_TYPE")]
# For now working with all columns/non-restricted file

# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)
# TCGA_Sample_Data_Restricted$PATIENT_ID <- substr(TCGA_Sample_Data_Restricted$PATIENT_ID, 1, nchar(TCGA_Sample_Data_Restricted$PATIENT_ID)-3) # Old line of code from when I had a restricted file
# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

#Patient_Sample_Combined <- Patient_Sample_Combined %>% subset(RNASEQ_DATA == "Yes")
# Read in GSVA results
TCGA_GSVA <- readRDS(file = "data_in_use/deconvolution_scores/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Scores.RDS")

TCGA_GSVA$TCGA_ID <- gsub("\\.", "-", TCGA_GSVA$TCGA_ID)
TCGA_GSVA <- TCGA_GSVA %>% dplyr::rename("SAMPLE_ID" = "TCGA_ID")

TCGA_merged <- merge(Patient_Sample_Combined, TCGA_GSVA)
rm(TCGA_GSVA, TCGA_Patient_Data, TCGA_Sample_Data, Patient_Sample_Combined)

##### Load TCGA RSEM data and format for merge with meta data #####
# read in file 
TCGA_RSEM <- read.table(file = 'data_in_use/data_RNA_Seq_v2_expression_median.txt', sep = '\t', header = TRUE)
Gene_Symbols <- TCGA_RSEM$Hugo_Symbol # save the gene names columns as the rownames
TCGA_RSEM_matrix <- as.matrix(TCGA_RSEM[,3:ncol(TCGA_RSEM)]) # make into matrix with just the numerical data
TCGA_Transposed <- t(TCGA_RSEM_matrix)
SAMPLE_ID <- rownames(TCGA_Transposed)
TCGA_Transposed <- as_tibble(TCGA_Transposed)
colnames(TCGA_Transposed) <- Gene_Symbols # Make Gene Symbols the column names
TCGA_Transposed$SAMPLE_ID <- SAMPLE_ID # Set sample IDs
TCGA_Transposed$SAMPLE_ID <- gsub("\\.", "-", TCGA_Transposed$SAMPLE_ID)

TCGA_merged <- merge(TCGA_merged, TCGA_Transposed)

# Add in Maxime fibrosis scores (24-1018 UPDATE)
TCGA_Fibrosis_Scores <- read.table(file = "data_in_use/TCGA_PTC_fibrosis_tsv.txt", sep = '\t', header = TRUE) 
# Need to remove the "A" at the end for compatability to merge
TCGA_Fibrosis_Scores$tcga_id <- TCGA_Fibrosis_Scores$tcga_id %>% substr(1, nchar(TCGA_Fibrosis_Scores$tcga_id)-1)
# Need to rename "tcga_id" for merge compatibility 
TCGA_Fibrosis_Scores <- TCGA_Fibrosis_Scores %>% dplyr::rename("SAMPLE_ID" = "tcga_id")

TCGA_merged_fibrosis <- merge(TCGA_merged, TCGA_Fibrosis_Scores)

TCGA_Fibrosis_csv <- TCGA_merged_fibrosis[c("SAMPLE_ID", "PATIENT_ID", "myCAF", "fibrosis_pc")]
write_csv(TCGA_Fibrosis_csv, file = "data_in_use/24-1114_TCGA_myCAF_Fibrosis_scores.csv")

TCGA_merged_BRS <- TCGA_merged %>% subset(BRAFV600E_RAS != "")
### Plot CAF Scores by Diagnosis

# myCAF - compare BRAF-like vs RAS-like
max(TCGA_merged_BRS$myCAF) # shows max of 0.7571946
min(TCGA_merged_BRS$myCAF) # shows min of -0.5933045
plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "myCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.7, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1016_BRS_Split_myCAF.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(myCAF ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value < 2.2e-16
pairwise.wilcox.test(TCGA_merged_BRS$myCAF, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # <2e-16

# iCAF - compare BRAF-like vs RAS-like
max(TCGA_merged_BRS$iCAF) # shows max of 0.6731989
min(TCGA_merged_BRS$iCAF) # shows min of -0.5514999
plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "iCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.6, .7)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1016_BRS_Split_iCAF.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(iCAF ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value < 2.2e-16
pairwise.wilcox.test(TCGA_merged_BRS$iCAF, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # <2e-16


# APOE_PVL - compare BRAF-like vs RAS-like
max(TCGA_merged_BRS$APOE_PVL) # shows max of 0.6760419
min(TCGA_merged_BRS$APOE_PVL) # shows min of -0.6944996
plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, APOE_PVL)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "APOE_PVL Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.7, .7)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1016_BRS_Split_APOE_PVL.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(APOE_PVL ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 0.1407
pairwise.wilcox.test(TCGA_merged_BRS$APOE_PVL, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # 0.14


# Pericyte - compare BRAF-like vs RAS-like
max(TCGA_merged_BRS$Pericyte) # shows max of 0.7459793
min(TCGA_merged_BRS$Pericyte) # shows min of -0.6608335
plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, Pericyte)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "Pericyte Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.7, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1016_BRS_Split_Pericyte.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(Pericyte ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 0.002105
pairwise.wilcox.test(TCGA_merged_BRS$Pericyte, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # 0.0021

# vSMC - compare BRAF-like vs RAS-like
max(TCGA_merged_BRS$vSMC) # shows max of 0.7459793
min(TCGA_merged_BRS$vSMC) # shows min of -0.6608335
plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, vSMC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "vSMC Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .6),
                     limits = c(-.7, .7)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1016_BRS_Split_vSMC.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(vSMC ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 0.3151
pairwise.wilcox.test(TCGA_merged_BRS$vSMC, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # 0.32

##### 24-1114 update is adding plots for FAP and POSTN #####
# 
TCGA_merged_BRS$FAP_log2 <- log2(TCGA_merged_BRS$FAP+1)
max(TCGA_merged_BRS$FAP_log2) # max is 11.46379
min(TCGA_merged_BRS$FAP_log2) # min is 1.182502

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, FAP_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(FAP)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12),
                     limits = c(0, 12)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2FAP.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(FAP_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value < 2.2e-16
pairwise.wilcox.test(TCGA_merged_BRS$FAP_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # < 2.2e-16

# POSTN
TCGA_merged_BRS$POSTN_log2  <- log2(TCGA_merged_BRS$POSTN+1)
max(TCGA_merged_BRS$POSTN_log2) # max is 15.79301
min(TCGA_merged_BRS$POSTN_log2) # min is 3.225182

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, POSTN_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(POSTN)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(0, 16)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2POSTN.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(POSTN_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 7.109e-09
pairwise.wilcox.test(TCGA_merged_BRS$POSTN_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # p = 7.1e-09

# ACTA2
TCGA_merged_BRS$ACTA2_log2  <- log2(TCGA_merged_BRS$ACTA2+1)
max(TCGA_merged_BRS$ACTA2_log2) # max is 14.82804
min(TCGA_merged_BRS$ACTA2_log2) # min is 8.568928

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, ACTA2_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(ACTA2)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(0, 15)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2ACTA2.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(ACTA2_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 0.5847
pairwise.wilcox.test(TCGA_merged_BRS$ACTA2_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # p = 0.59


# RGS5
TCGA_merged_BRS$RGS5_log2  <- log2(TCGA_merged_BRS$RGS5+1)
max(TCGA_merged_BRS$RGS5_log2) # max is 15.45353
min(TCGA_merged_BRS$RGS5_log2) # min is 8.342499

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, RGS5_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(RGS5)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(0, 16)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2RGS5.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(RGS5_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 1.455e-08
pairwise.wilcox.test(TCGA_merged_BRS$RGS5_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # p = 1.5e-08

# CXCL12
TCGA_merged_BRS$CXCL12_log2  <- log2(TCGA_merged_BRS$CXCL12+1)
max(TCGA_merged_BRS$CXCL12_log2) # max is 12.51361
min(TCGA_merged_BRS$CXCL12_log2) # min is 5.809101

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, CXCL12_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(CXCL12)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12),
                     limits = c(0, 13)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2CXCL12.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(CXCL12_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 8.849e-05
pairwise.wilcox.test(TCGA_merged_BRS$CXCL12_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # p = 8.9e-05

# APOE
TCGA_merged_BRS$APOE_log2  <- log2(TCGA_merged_BRS$APOE+1)
max(TCGA_merged_BRS$APOE_log2) # max is 12.51361
min(TCGA_merged_BRS$APOE_log2) # min is 5.809101

plot <- ggplot(TCGA_merged_BRS, aes(BRAFV600E_RAS, APOE_log2)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = BRAFV600E_RAS),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "indianred1")) +
  labs (x = "BRS", y = "log2(APOE)") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="BRS", limits = c("Ras-like", "Braf-like")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(0, 17)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1114_TCGA_Thyroid_Fibroblast_GSVA_Deconvolution/24-1114_BRS_Split_Log2APOE.png", 
       width = 3, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(APOE_log2 ~ BRAFV600E_RAS, 
             data = TCGA_merged_BRS) # p-value = 2.798e-15
pairwise.wilcox.test(TCGA_merged_BRS$APOE_log2, 
                     TCGA_merged_BRS$BRAFV600E_RAS, 
                     p.adjust.method = "bonferroni") # p = 2.8e-15

library(corrplot)
#TCGA_merged_fibrosis <- TCGA_merged_fibrosis %>% subset(!is.na(purity_adjusted_fibrosis))
TCGA_merged_fibrosis <- TCGA_merged_fibrosis %>% dplyr::rename("fibrosis" = "fibrosis_pc")
correlation_data <- TCGA_merged_fibrosis[c("fibrosis", "myCAF", "iCAF", "APOE_PVL", "Pericyte", "vSMC", "FAP", "COL1A1", "ACTA2", "RGS5", "POSTN", "CXCL12")] # myCAF populations
corr_mat=cor(correlation_data,method="s") #create Spearman correlation matrix

png(file = "outputs/24-1112_TCGA_CAF_ssGSEA/24-1112_TCGA_CAF_ssGSEA_fibrosis_correlation.png", res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black",
         tl.col = "black",
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()

# BRAF-RAS score added
TCGA_merged_fibrosis_BRS <- TCGA_merged_fibrosis %>% subset(!is.na(BRAFV600E_RAF_SCORE))
TCGA_merged_fibrosis_BRS <- TCGA_merged_fibrosis_BRS %>% dplyr::rename("BRS" = "BRAFV600E_RAF_SCORE")

correlation_data <- TCGA_merged_fibrosis_BRS[c("fibrosis_pc", "myCAF", "iCAF", "APOE_PVL", "Pericyte", "vSMC", "FAP", "COL1A1", "ACTA2", "RGS5", "POSTN", "CXCL12", "BRS")] # myCAF populations
corr_mat=cor(correlation_data,method="s") #create Spearman correlation matrix

png(file = "outputs/24-1018_TCGA_CAF_ssGSEA/24-1018_TCGA_CAF_ssGSEA_fibrosis_correlation_BRS.png", res = 300, height = 1920, width = 1920)
corrplot(corr_mat, method = "color",
         type = "upper", order = "hclust", 
         addCoef.col = "black",
         tl.col = "black",
         col = rev(COL2("RdBu", 200))) # takes the RdBu scale that is used in the packaged and flips it)

dev.off()


correlation <- ggplot(TCGA_merged_fibrosis, 
                           aes(fibrosis_pc, myCAF)) + 
  geom_point(alpha = 0.7) + 
  theme_classic()
ggsave(file = "outputs/24-1018_TCGA_CAF_ssGSEA/24-1018_fibrosis_pc_myCAF_Correlation.png", 
       correlation, dpi = 600, height = 4.5, width = 4.5)
correlation <- cor(TCGA_merged_fibrosis$fibrosis_pc, TCGA_merged_fibrosis$myCAF, method = "spearman")
print(correlation)
cor_test_result <- cor.test(TCGA_merged_fibrosis$fibrosis_pc, TCGA_merged_fibrosis$myCAF, method = "spearman")
print(cor_test_result)
