### Author: Matthew Aaron Loberg
### Date: 24-1130
### Script: Stromal_subclustering_ssGSVA_Boxplots_by_Diagnosis.R
### Source Script Name: 24-1130_Progression_Series_Fibroblasts_GSVA_Scores.R

# Goal: Plot new GSVA scores (generated from TPM normalized data with GSVA package; see 24-0829_ssGSEA_Fibroblasts_TPM.R)
# Plotting will be by diagnosis/other metrics

# Note on script origin: 
# Script is modified from "24-0417_Progression_Series_Fibroblasts.R"
# For this modification, I also copied the object "22-0606_CleanedMergedData_DESeq2NormalizedReads.rds" into the "data_in_use folder for this project

### 24-1016 Update
# This is an update to "24-0903_Progression_Series_Fibroblast_GSVA_Scores.R"
# I am just copying that script over and updating several colors to match paper color formatting
# Also going to plot FAP and POSTN

### 24-1130 Update
# Adding red dotted line to the scores plots

# Load necessary packages:
library(tidyverse)
library(cowplot)
library(RColorBrewer)

####### Loading and Setting Up Clinical Data with DESeq2 normalized counts ###################
#GeneExpression <- read.table(file = "data_in_use/21-1214_DESeq2_Normalized_RNA_Counts.txt", sep = '\t', header = TRUE)
#GeneExpression$RNA.ID <- sub('.', '', GeneExpression$RNA.ID)
#ClinicalData <- read_csv(file = "data_in_use/VUMC.cohort.GX_4-26-22.csv")
#MergedData <- as_tibble(merge(ClinicalData, GeneExpression, by.x = "RNA.ID", by.y = "RNA.ID", all = T))
#CleanedMergedData <- MergedData %>% drop_na(FAP) # Drop all data without RNA reads
#rm(GeneExpression, ClinicalData, MergedData)
#saveRDS(CleanedMergedData, file = "data_in_use/22-0606_CleanedMergedData_DESeq2NormalizedReads.rds")
#rm(CleanedMergedData)

##### In future, can skip to here #####
CleanedMergedData <- readRDS(file = "data_in_use/22-0606_CleanedMergedData_DESeq2NormalizedReads.rds")

##### Load in GSVA Results (Internal) #####
GSVA_Results <- readRDS(file = "data_in_use/deconvolution_scores/24-0903_Thyroid_Fibroblast_GSVA_Scores.RDS")
GSVA_Results$RNA.ID <- sub('.', '', GSVA_Results$RNA.ID) # Get rid of the 'X" at the start of GSVA_Results

##### Merged together #####
cohort <- merge(x = CleanedMergedData, y = GSVA_Results, by.x = "RNA.ID", by.y = "RNA.ID")

##### Load in GSVA Results (External) #####
GSVA_Results <- readRDS(file = "data_in_use/deconvolution_scores/24-0905_Thyroid_Fibroblast_GSVA_Scores_External.RDS")
GSVA_Results$RNA.ID <- sub('.', '', GSVA_Results$RNA.ID) # Get rid of the 'X" at the start of GSVA_Results

##### Merged together #####
cohort <- merge(x = cohort, y = GSVA_Results, by.x = "RNA.ID", by.y = "RNA.ID")

#### Establish a column for labeling as Follicular, Papillary, or Transformed ####
cohort$Category <- "Benign"
# Anything NOT follicular, papillary, or transformed will have a category of benign

# Follicular lesions
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "NIFTP" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "FC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "HC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "EFVPTC"){
        cohort$Category[i] <- "Follicular"
  }
}

# Papillary lesions
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "IFVPTC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PTC"){
        cohort$Category[i] <- "Papillary"
  }
}

# Transformed lesion
for(i in 1:nrow(cohort)){
  if(cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "PDTC" |
     cohort$Diagnosis.with.iFVPTC.and.eFVPTC[i] == "ATC"){
        cohort$Category[i] <- "Transformed"
  }
}

#### Create a new variable that will be the following groups: Benign, RAS-like WDTC, BRAF-like WDTC, PDTC, ATC
cohort$Diagnosis_Simplified <- cohort$Diagnosis
cohort$BRS <- as.numeric(cohort$BRS)
for(i in 1:nrow(cohort)){
  if(cohort$Category[i] == "Papillary" | cohort$Category[i] == "Follicular"){
    if(cohort$BRS[i] >= 0){
      cohort$Diagnosis_Simplified[i] <- "RAS-Like WDTC"
    }
    else if(cohort$BRS[i] < 0){
      cohort$Diagnosis_Simplified[i] <- "BRAF-Like WDTC"
    }
  }
  if(cohort$Category[i] == "Benign"){
    cohort$Diagnosis_Simplified[i] <- "Benign"
  }
}
table(cohort$Diagnosis_Simplified)

########## Plots ###############
### myCAF GSVA score
max(cohort$myCAF) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) 
min(cohort$myCAF) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis.with.iFVPTC.and.eFVPTC, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Category),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "indianred1", "#FFCC99", "#0075DC")) +
  labs (x = "Diagnosis", y = "myCAF Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("MNG", "FA", "HA", "HT", "FC", "HC", "NIFTP", "EFVPTC", "IFVPTC", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_myCAF.png", 
       width = 7.5, height = 3, 
       plot, dpi = 600)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis_Simplified, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis_Simplified),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99", "paleturquoise2", "indianred1")) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  labs (x = "Diagnosis", y = "myCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "RAS-Like WDTC", "BRAF-Like WDTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Simplified_myCAF.png", 
       width = 5, height = 3, 
       plot, dpi = 600)

### Statistics
kruskal.test(myCAF ~ Diagnosis_Simplified, data = cohort) # shouldn't be necessary for only two groups - use if more than two groups
pairwise.wilcox.test(cohort$myCAF, cohort$Diagnosis_Simplified, 
                     p.adjust.method = "BH")

kruskal.test(iCAF ~ Diagnosis_Simplified, data = cohort) # shouldn't be necessary for only two groups - use if more than two groups
pairwise.wilcox.test(cohort$iCAF, cohort$Diagnosis_Simplified, 
                     p.adjust.method = "BH")

kruskal.test(Pericyte ~ Diagnosis_Simplified, data = cohort) # shouldn't be necessary for only two groups - use if more than two groups
pairwise.wilcox.test(cohort$Pericyte, cohort$Diagnosis_Simplified, 
                     p.adjust.method = "BH")

kruskal.test(vSMC ~ Diagnosis_Simplified, data = cohort) # shouldn't be necessary for only two groups - use if more than two groups
pairwise.wilcox.test(cohort$vSMC, cohort$Diagnosis_Simplified, 
                     p.adjust.method = "BH")

kruskal.test(APOE_PVL~ Diagnosis_Simplified, data = cohort) # shouldn't be necessary for only two groups - use if more than two groups
pairwise.wilcox.test(cohort$APOE_PVL, cohort$Diagnosis_Simplified,
                     p.adjust.method = "BH")

kruskal.test(log2(FAP+1) ~ Diagnosis_Simplified, data = cohort)
kruskal.test(log2(POSTN+1) ~ Diagnosis_Simplified, data = cohort)




### iCAF GSVA score
max(cohort$iCAF) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) 
min(cohort$iCAF) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis.with.iFVPTC.and.eFVPTC, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Category),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "indianred1", "#FFCC99", "#0075DC")) +
  labs (x = "Diagnosis", y = "iCAF Score") +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("MNG", "FA", "HA", "HT", "FC", "HC", "NIFTP", "EFVPTC", "IFVPTC", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_iCAF.png", 
       width = 7.5, height = 3, 
       plot, dpi = 600)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis_Simplified, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis_Simplified),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99", "paleturquoise2", "indianred1")) +
  labs (x = "Diagnosis", y = "iCAF Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "RAS-Like WDTC", "BRAF-Like WDTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Simplified_iCAF.png", 
       width = 5, height = 3, 
       plot, dpi = 600)





### vSMC GSVA score
max(cohort$vSMC) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) 
min(cohort$vSMC) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis.with.iFVPTC.and.eFVPTC, vSMC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Category),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "indianred1", "#FFCC99", "#0075DC")) +
  labs (x = "Diagnosis", y = "vSMC Score") +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("MNG", "FA", "HA", "HT", "FC", "HC", "NIFTP", "EFVPTC", "IFVPTC", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_vSMC.png", 
       width = 7.5, height = 3, 
       plot, dpi = 600)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis_Simplified, vSMC)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis_Simplified),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99", "paleturquoise2", "indianred1")) +
  labs (x = "Diagnosis", y = "vSMC Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "RAS-Like WDTC", "BRAF-Like WDTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.65, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Simplified_vSMC.png", 
       width = 5, height = 3, 
       plot, dpi = 600)


### Pericyte GSVA score
max(cohort$Pericyte) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) (14.59)
min(cohort$Pericyte) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  (0)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis.with.iFVPTC.and.eFVPTC, Pericyte)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Category),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "indianred1", "#FFCC99", "#0075DC")) +
  labs (x = "Diagnosis", y = "Pericyte Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("MNG", "FA", "HA", "HT", "FC", "HC", "NIFTP", "EFVPTC", "IFVPTC", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.85, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Pericyte.png", 
       width = 7.5, height = 3, 
       plot, dpi = 600)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis_Simplified, Pericyte)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis_Simplified),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99", "paleturquoise2", "indianred1")) +
  labs (x = "Diagnosis", y = "Pericyte Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "RAS-Like WDTC", "BRAF-Like WDTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8),
                     limits = c(-.85, .85)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Simplified_Pericyte.png", 
       width = 5, height = 3, 
       plot, dpi = 600)

### APOE GSVA score
max(cohort$APOE_PVL) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) (14.59)
min(cohort$APOE_PVL) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  (0)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis.with.iFVPTC.and.eFVPTC, APOE_PVL)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Category),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "indianred1", "#FFCC99", "#0075DC")) +
  labs (x = "Diagnosis", y = "APOE CAF Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("MNG", "FA", "HA", "HT", "FC", "HC", "NIFTP", "EFVPTC", "IFVPTC", "PTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-8, -.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.85, .65)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_APOE_PVL.png", 
       width = 7.5, height = 3, 
       plot, dpi = 600)

# use the above 2 values to set the breaks for the plot
plot <- ggplot(cohort, aes(Diagnosis_Simplified, APOE_PVL)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Diagnosis_Simplified),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#0075DC", "#808080", "#FFCC99", "paleturquoise2", "indianred1")) +
  labs (x = "Diagnosis", y = "APOE_PVL Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "RAS-Like WDTC", "BRAF-Like WDTC", "PDTC", "ATC")) +
  scale_y_continuous(breaks = c(-8, -.6, -.4, -.2, 0, .2, .4, .6),
                     limits = c(-.85, .65)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1130_Thyroid_Fibroblast_GSVA_Deconvolution/24-1130_All_Diagnoses_Simplified_APOE_PVL.png", 
       width = 5, height = 3, 
       plot, dpi = 600)

