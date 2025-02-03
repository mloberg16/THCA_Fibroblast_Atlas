### Author: Matthew Aaron Loberg
### Date: October 15, 2024
### Script: GSE213647_ssGSVA_Fibroblast_BoxPlots.R
### Source Script Name: 24-1015_GSE213647_ssGSEA_Fibroblasts.R

# Goal: Run  ssGSEA for GSE213647 from Lee et al. 2024 Nature Communications
# In this script, I will run it for my INTERNAL CAF scores that were generated from scRNA-seq data

# This bulk RNA sequencing data from Lee et al. was TPM normalized and subset to protein coding genes on 24-0910 by Hua-Chang

# Going to install the following package: 
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#1_Quick_start
# package is GSVA
# Install as below (will comment out after installed): 
# BiocManager::install("GSVA")

# 24-1015 Update
# JUST UPDATING visualization colors (from 24-0911 script)

##### Load packages #####
library(GSVA)
library(tidyverse)


##### Load TPM data #####
# read in file 
tpm <- read.table(file = 'data_in_use/GSE213647_thyroid_cancer.TPM.protein_coding.tsv', sep = '\t', header = TRUE)
tpm_rownames <- tpm$gene_name # save the gene names columns as the rownames
tpm_matrix <- as.matrix(tpm[,5:ncol(tpm)]) # make into matrix with just the numerical data
rownames(tpm_matrix) <- tpm_rownames # set the rownames of tpm_matrix to the stored gene names
rm(tpm, tpm_rownames) # remove tpm and tpm_rownames (only keeping the matrix)
 
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
gsva_results  <- GSVA::gsva(expr = tpm_matrix,
                            gset.idx.list = gene_sets,
                            #method = "ssgsea", # NOT RUNNING ssgsea, just default GSVA this time
                            #ssgsea.norm = TRUE,# NOT RUNNING ssgsea, just default GSVA this time
                            verbose = TRUE)
gsva_results_transformed <- t(gsva_results)
GEO_ID <- rownames(gsva_results_transformed)
gsva_results_transformed <- gsva_results_transformed %>% as_tibble()
gsva_results_transformed$GEO_ID <- GEO_ID
rm(GEO_ID)

# Save the GSVA Results
saveRDS(gsva_results_transformed, file = "data_in_use/deconvolution_scores/24-0911_GSE213647_Thyroid_Fibroblast_GSVA_Scores.RDS")

# Cleaning Up 
rm(list = ls())

##### Visualization #####

### Data formatting ###
# Read in meta data 
GSE213647_meta <- read.csv(file = "data_in_use/GSE213647_Meta_Data.csv")
# Read in GSVA results
GSE213647_GSVA <- readRDS(file = "data_in_use/deconvolution_scores/24-0911_GSE213647_Thyroid_Fibroblast_GSVA_Scores.RDS")
# Change GSE_ID column to exclude extra characters
GSE213647_GSVA$GEO_ID <- substr(GSE213647_GSVA$GEO_ID, 1, 10)
# Merge meta and GSVA
GSE213647_Merged <- merge(GSE213647_meta, GSE213647_GSVA)

# Check # of diagnoses 
table(GSE213647_Merged$Diagnosis)


### Plot CAF Scores by Diagnosis

# myCAF - all diagnoses
max(GSE213647_Merged$myCAF) # shows max of 0.7977652
min(GSE213647_Merged$myCAF) # shows min of -0.6541204
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c(
                               "#0075DC",        # ATC color
                               "#808080",        # Normal color
                               "paleturquoise2", # PDTC color
                               "#FFCC99"         # PTC color
                               )) + 
  labs (x = "Diagnosis", y = "myCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.7, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_All_Diagnoses_myCAF.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(myCAF ~ Diagnosis, 
             data = GSE213647_Merged) # p-value < 2.2e-16
pairwise.wilcox.test(GSE213647_Merged$myCAF, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 8.3e-07; ATC/PTC: 1.00; PTC/normal: <2e-16; PDTC/ATC: 0.92; PDTC/PTC: 0.29

# myCAF excluded to Normal, PTC, ATC
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99")) +
  labs (x = "Diagnosis", y = "myCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.7, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_PDTC_Excluded_myCAF.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(myCAF ~ Diagnosis, 
             data = temp) # p-value < 2.2e-16
pairwise.wilcox.test(temp$myCAF, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 4.2e-07; ATC/PTC: 0.81; PTC/normal: <2e-16
rm(temp)


# iCAF - all diagnoses
max(GSE213647_Merged$iCAF) # shows max of 0.7296673
min(GSE213647_Merged$iCAF) # shows min of -0.6933978
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c(
                               "#0075DC",        # ATC color
                               "#808080",        # Normal color
                               "paleturquoise2", # PDTC color
                               "#FFCC99"         # PTC color
  )) + 
  labs (x = "Diagnosis", y = "iCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.75, .75)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_All_Diagnoses_iCAF.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(iCAF ~ Diagnosis, 
             data = GSE213647_Merged) # p-value: 1.985e-08
pairwise.wilcox.test(GSE213647_Merged$iCAF, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.0084; ATC/PTC: 0.5406; PTC/normal: 4.8e-06; PDTC/ATC: 0.1911; PDTC/PTC: 0.0396

# iCAF excluded to Normal, PTC, ATC
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99")) +
  labs (x = "Diagnosis", y = "iCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.75, .75)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_PDTC_Excluded_iCAF.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(iCAF ~ Diagnosis, 
             data = temp) # p-value: 3.094e-07
pairwise.wilcox.test(temp$iCAF, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.0042; ATC/PTC: 0.2703; PTC/normal: 2.4e-06
rm(temp)

# Pericyte - all diagnoses
max(GSE213647_Merged$Pericyte) # shows max of 0.7726496
min(GSE213647_Merged$Pericyte) # shows min of -0.7165331
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, Pericyte)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c(
                               "#0075DC",        # ATC color
                               "#808080",        # Normal color
                               "paleturquoise2", # PDTC color
                               "#FFCC99"         # PTC color
  )) + 
  labs (x = "Diagnosis", y = "Pericyte Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.75, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_All_Diagnoses_Pericyte.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(Pericyte ~ Diagnosis, 
             data = GSE213647_Merged) # p-value: 3.414e-06
pairwise.wilcox.test(GSE213647_Merged$Pericyte, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.5427; ATC/PTC: 0.0055; PTC/normal: 8.3e-05; PDTC/ATC: 1.00; PDTC/PTC: 0.5846

# Pericyte excluded to Normal, PTC, ATC
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, Pericyte)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99")) +
  labs (x = "Diagnosis", y = "Pericyte Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.75, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_PDTC_Excluded_Pericyte.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(Pericyte ~ Diagnosis, 
             data = temp) # p-value: 2.177e-06
pairwise.wilcox.test(temp$Pericyte, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.2714; ATC/PTC: 0.0027; PTC/normal: 4.2e-05
rm(temp)


# vSMC - all diagnoses
max(GSE213647_Merged$vSMC) # shows max of 0.7267613
min(GSE213647_Merged$vSMC) # shows min of -0.6509528
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, vSMC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c(
                               "#0075DC",        # ATC color
                               "#808080",        # Normal color
                               "paleturquoise2", # PDTC color
                               "#FFCC99"         # PTC color
  )) + 
  labs (x = "Diagnosis", y = "vSMC Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.7, .75)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_All_Diagnoses_vSMC.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(vSMC ~ Diagnosis, 
             data = GSE213647_Merged) # p-value: 0.0002727
pairwise.wilcox.test(GSE213647_Merged$vSMC, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.0426; ATC/PTC: 0.0027; PTC/normal: 0.4775; PDTC/ATC: 1.00; PDTC/PTC: 0.0773

# vSMC excluded to Normal, PTC, ATC
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, vSMC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99")) +
  labs (x = "Diagnosis", y = "vSMC Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.7, .75)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_PDTC_Excluded_vSMC.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(vSMC ~ Diagnosis, 
             data = temp) # p-value: 0.001254
pairwise.wilcox.test(temp$vSMC, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.0213; ATC/PTC: 0.0014; PTC/normal: 0.2388
rm(temp)


# APOE_PVL - all diagnoses
max(GSE213647_Merged$APOE_PVL) # shows max of 0.7439274
min(GSE213647_Merged$APOE_PVL) # shows min of -0.75037
plot <- ggplot(GSE213647_Merged, aes(Diagnosis, APOE_PVL)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c(
                               "#0075DC",        # ATC color
                               "#808080",        # Normal color
                               "paleturquoise2", # PDTC color
                               "#FFCC99"         # PTC color
  )) + 
  labs (x = "Diagnosis", y = "APOE_PVL Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.8, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_All_Diagnoses_APOE_PVL.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
kruskal.test(APOE_PVL ~ Diagnosis, 
             data = GSE213647_Merged) # p-value: 2.427e-06
pairwise.wilcox.test(GSE213647_Merged$APOE_PVL, 
                     GSE213647_Merged$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.51; ATC/PTC: 1.00; PTC/normal: 1.8e-06; PDTC/ATC: 1.00; PDTC/PTC: 0.92

# APOE_PVL excluded to Normal, PTC, ATC
plot <- ggplot(GSE213647_Merged %>% subset(Diagnosis != "PDTC"), aes(Diagnosis, APOE_PVL)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.3, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99")) +
  labs (x = "Diagnosis", y = "APOE_PVL Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Normal", "PTC", "ATC")) +
  scale_y_continuous(breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.8, .8)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1015_GSE213647_Thyroid_Fibroblast_GSVA_Deconvolution/24-1015_PDTC_Excluded_APOE_PVL.png", 
       width = 4, height = 3, 
       plot, dpi = 600, create.dir = TRUE)
temp <- GSE213647_Merged %>% subset(Diagnosis != "PDTC")
kruskal.test(APOE_PVL ~ Diagnosis, 
             data = temp) # p-value: 1.289e-06
pairwise.wilcox.test(temp$APOE_PVL, 
                     temp$Diagnosis, 
                     p.adjust.method = "bonferroni") # ATC/Normal: 0.26; ATC/PTC: 1.00; PTC/normal: 8.8e-07
rm(temp)

