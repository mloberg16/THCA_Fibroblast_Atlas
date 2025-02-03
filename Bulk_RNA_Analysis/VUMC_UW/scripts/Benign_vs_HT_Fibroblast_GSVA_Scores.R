### Author: Matthew Aaron Loberg
### Date: 24-1207
### Script: Benign_vs_HT_Fibroblast_GSVA_Scores.R
### Source Script Name: 24-1207_Benign_vs_HT_Fibroblast_GSVA_Scores.R

# Goal: Plot new GSVA scores to compare benign lesions vs benign lesions classified as HT (Hashimoto's thyroiditis)
# See if Hashimoto's has a different stromal profile than other benign lesions

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
cohort <- merge(x = cohort, y = GSVA_Results, by.x = "RNA.ID", by.y = "RNA.ID")

##### Breakdown of benign into two subcohorts: benign and (MNG, HA, FA) and HT
table(cohort$Diagnosis)
Benign_cohort <- cohort %>% subset(Diagnosis == "FA" | 
                            Diagnosis == "HA" |
                            Diagnosis == "HT" |
                            Diagnosis == "MNG")
table(Benign_cohort$Diagnosis)

Benign_cohort$Benign_Breakdown <- ifelse(Benign_cohort$Diagnosis == "HT", "HT", "Benign")
table(Benign_cohort$Benign_Breakdown)

########## Plots ###############
### myCAF GSVA score
max(Benign_cohort$myCAF) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) 
min(Benign_cohort$myCAF) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  

# use the above 2 values to set the breaks for the plot
plot <- ggplot(Benign_cohort, aes(Benign_Breakdown, myCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Benign_Breakdown),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "#FE8F42")) +
  labs (x = "Diagnosis", y = "myCAF Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "HT")) +
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-.6, .5)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1207_Benign_vs_HT_Fibroblast_GSVA/24-1207_Benign_vs_HT_myCAF.png", 
       width = 2, height = 3, create.dir = TRUE, 
       plot, dpi = 600)
pairwise.wilcox.test(Benign_cohort$myCAF, Benign_cohort$Benign_Breakdown, 
                     p.adjust.method = "bonferroni")


### iCAF GSVA score
max(Benign_cohort$iCAF) # Use this value to set the plot max value (e.g., if 8.5, would set as 9) 
min(Benign_cohort$iCAF) # Use this value to set the plot min value (e.g., if 0.5 would set as 0)  

# use the above 2 values to set the breaks for the plot
plot <- ggplot(Benign_cohort, aes(Benign_Breakdown, iCAF)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Benign_Breakdown),
               alpha = 0.7, 
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5, 
              alpha = 0.7,
              show.legend = FALSE) +
  scale_fill_manual(values = c("#808080", "#FE8F42")) +
  labs (x = "Diagnosis", y = "iCAF Score") + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Diagnosis", limits = c("Benign", "HT")) +
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-.5, .75)) 
# Note, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/24-1207_Benign_vs_HT_Fibroblast_GSVA/24-1207_Benign_vs_HT_iCAF.png", 
       width = 2, height = 3, create.dir = TRUE, 
       plot, dpi = 600)
pairwise.wilcox.test(Benign_cohort$iCAF, Benign_cohort$Benign_Breakdown, 
                     p.adjust.method = "bonferroni")
