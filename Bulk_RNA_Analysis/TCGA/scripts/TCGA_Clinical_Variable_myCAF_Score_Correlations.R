### Author: Matthew Aaron Loberg
### Date: November 30, 2024
### Script: TCGA_Clinical_Variable_myCAF_Score_Correlations.R
### Source Script Name: 24-1130_TCGA_Survival_Extrathyroidal_Extension_myCAFs.R

### Goal: 
# Look at associations of clinical variables with myCAF score.
# Clinical variables to look at: 
# extrathyroidal extension (yes/no), LN metastasis (yes/no), distant metastasis (yes/no), 2009 ATA risk score (low or int/high), histology (follicular-variant, classical papillary, diffuse sclerosing/tall cell)

### Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot) # Using cowplot due to ggsave issues with survival curves

######################### Load Data ###########################
# TCGA patient data
TCGA_Patient_Data <- as_tibble(read.table(file = "data_in_use/22-0107_TCGA_data_clinical_patient.txt", sep = '\t', header = TRUE)) # Load patient data file and make into tibble for easier cleaning

# TCGA Sample data
TCGA_Sample_Data <- as_tibble(read.table(file = "data_in_use/22-0107_TCGA_data_clinical_sample.txt", sep = '\t', header = TRUE)) # Load sample data file and make into tibble for easier data cleaning

# TCGA fibroblast deconvolution
TCGA_Fibroblast_Deconvolution <- readRDS(file = "data_in_use/deconvolution_scores/24-1016_TCGA_Thyroid_Fibroblast_GSVA_Scores.RDS")

############################## Data Cleaning for analysis ###########################
# PATIEND_ID has an "01" appended to the end of each sample. Here, I will remove it from each sample. 
TCGA_Sample_Data$PATIENT_ID <- substr(TCGA_Sample_Data$PATIENT_ID, 1, nchar(TCGA_Sample_Data$PATIENT_ID)-3)

# Combine sample data and patient data
Patient_Sample_Combined <- as_tibble(merge(TCGA_Sample_Data, TCGA_Patient_Data, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all = T))

# subset to just samples with associated RNA sequencing data
Patient_Sample_Combined <- Patient_Sample_Combined %>% subset(RNASEQ_DATA == "Yes")

# Clean up environment - only leave the data that I need (Patient_Sample_Combined)
remove(list = c("TCGA_Patient_Data", "TCGA_Sample_Data"))

# DFS_STATUS needs to be numeric, so separate column
Patient_Sample_Combined <- separate(data = Patient_Sample_Combined, col = DFS_STATUS, into = c("DFS_Status_Numeric", "DFS_Status_Char"), sep = ":")
Patient_Sample_Combined$DFS_Status_Numeric <- as.numeric(Patient_Sample_Combined$DFS_Status_Numeric) # Change to numeric

# OS_STATUS needs to be numeric, so separate column
Patient_Sample_Combined <-separate(data = Patient_Sample_Combined, col = OS_STATUS, into = c("OS_Status_Numeric", "OS_Status_Char"), sep = ":")
Patient_Sample_Combined$OS_Status_Numeric <- as.numeric(Patient_Sample_Combined$OS_Status_Numeric) # Change to numeric

# Add fibroblast scores to Patient_Sample_Combined
# First need to modify TCGA.ID of fibroblast scores to be "-" instead of "." and then rename it
TCGA_Fibroblast_Deconvolution$TCGA_ID <- gsub("\\.", "-", TCGA_Fibroblast_Deconvolution$TCGA_ID)
TCGA_Fibroblast_Deconvolution <- TCGA_Fibroblast_Deconvolution %>% dplyr::rename("SAMPLE_ID" = "TCGA_ID")

Patient_Sample_Combined <- Patient_Sample_Combined %>% merge(TCGA_Fibroblast_Deconvolution)
rm(TCGA_Fibroblast_Deconvolution)



#################### DISEASE STAGE ################################

## Disease stage in depth 
# Will start by subsetting just to lesions with disease stage data
Disease_Stage_Tibble <- Patient_Sample_Combined %>% subset(DISEASE_STAGE != "" & DISEASE_STAGE != "[Not Available]")
# Rename Stage IVA and Stage IVC as just stage IV
for(i in 1:nrow(Disease_Stage_Tibble)){
  if(Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage IVA" |
     Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage IVC"){
    Disease_Stage_Tibble$DISEASE_STAGE[i] <- "Stage IV"
  }
}
Disease_Stage_Tibble$DISEASE_STAGE_Simplified <- Disease_Stage_Tibble$DISEASE_STAGE
for(i in 1:nrow(Disease_Stage_Tibble)){
  if(Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage I" |
     Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage II"){
    Disease_Stage_Tibble$DISEASE_STAGE_Simplified[i] <- "Stage I/II"
  }
  else if(Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage III" |
          Disease_Stage_Tibble$DISEASE_STAGE[i] == "Stage IV"){
    Disease_Stage_Tibble$DISEASE_STAGE_Simplified[i] <- "Stage III/IV"
  }
}
table(Disease_Stage_Tibble$DISEASE_STAGE_Simplified)
# Shows that there are 309 stage I/II and 143 stage III/IV
# Run stats 
#res.aov <- aov(myCAF ~ DISEASE_STAGE, data = Disease_Stage_Tibble)
#summary(res.aov)
#TukeyHSD(res.aov)
# Make the plot with restricted groups
plot <- ggplot(Disease_Stage_Tibble, aes(DISEASE_STAGE_Simplified, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = DISEASE_STAGE_Simplified),
               alpha = 0.2,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c(

                                "#00BFC4", # Stage I/II Color
                                "purple")) + # Stage III/IV Color
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Disease Stage", limits = c("Stage I/II", "Stage III/IV")) + 
  labs(y = "LUM Score") + 
  ggtitle("LUM Score\nby Disease Stage") +
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-.6, .8)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(    
        panel.border = element_rect(colour = "black", size = 3, fill = NA),
        plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
        axis.title = element_text(size = 30, color = "black", face = "bold"),
        axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_Disease_Stage.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for myCAF score by Disease Stage
kruskal.test(myCAF ~ DISEASE_STAGE_Simplified, data = Disease_Stage_Tibble) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Disease_Stage_Tibble$myCAF, Disease_Stage_Tibble$DISEASE_STAGE_Simplified, # this returned 0.2475
                     p.adjust.method = "BH")
# cleaning up disease stage
rm(Disease_Stage_Tibble)


################# RISK GROUP #####################
## Risk Group light
# Risk groups with myCAF plotted
plot <- ggplot(Patient_Sample_Combined, aes(RISK_GROUP, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = RISK_GROUP),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(color = DFS_Status_Char),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) 
plot

# Risk group fancy
# Updated 22-0425
Patient_Sample_Combined_RISK <- Patient_Sample_Combined %>% subset(RISK_GROUP != "")
plot <- ggplot(Patient_Sample_Combined_RISK, aes(RISK_GROUP, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = RISK_GROUP),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c("lightgrey", # High Risk group color
                               "darkgrey", # Intermediate risk group color
                               "black")) + # Low risk group color 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Risk Group", limits = c("Low", "Intermediate", "High")) + 
  labs(y = "SFRP2 Score") + 
  ggtitle("SFRP2 Score\nby Risk Group") +
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 2)  + 
  theme(
        panel.border = element_rect(colour = "black", size = 3, fill = NA),
        plot.title = element_text(size = 30, color = "black", face = "bold", hjust = 0.5, vjust = 1.8),
        axis.title = element_text(size = 30, color = "black", face = "bold"),
        axis.text = element_text(size = 25, color = "black", face = "bold"))
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_Risk_Group.png",
       plot, 
       width = 8, height = 7, dpi = 600)

# Stats for myCAF score by Risk group
kruskal.test(myCAF ~ RISK_GROUP, data = Patient_Sample_Combined_RISK) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Patient_Sample_Combined_RISK$myCAF, Patient_Sample_Combined_RISK$RISK_GROUP, # this returned 0.2475
                     p.adjust.method = "BH")


# Risk group simplified
Patient_Sample_Combined_RISK$RISK_GROUP_SIMPLIFIED <- ifelse(Patient_Sample_Combined_RISK$RISK_GROUP == "Low", "Low", "Int/High")

# csv for JCI
csv <- Patient_Sample_Combined_RISK[c("PATIENT_ID", "SAMPLE_ID", "RISK_GROUP_SIMPLIFIED", "myCAF")]
write_csv(csv, file = "data_in_use/24-1130_TCGA_Risk_Group_JCI.csv")

# Plot Risk group simplified
plot <- ggplot(Patient_Sample_Combined_RISK, aes(RISK_GROUP_SIMPLIFIED, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = RISK_GROUP_SIMPLIFIED),
               alpha = 0.2,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c("black", # Intermediate/high Risk group color
                                          "lightgrey")) + # Low risk group color 
  scale_x_discrete(name = "Risk Group", limits = c("Low", "Int/High")) + 
  labs(y = "myCAF Score") + 
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = .75)  + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) 
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_Risk_Group_simplified.png",
       plot, 
       width = 3, height = 3, dpi = 600)

# stats
kruskal.test(myCAF ~ RISK_GROUP_SIMPLIFIED, data = Patient_Sample_Combined_RISK) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Patient_Sample_Combined_RISK$myCAF, Patient_Sample_Combined_RISK$RISK_GROUP_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "BH")

# Cleaining up from risk 
rm(Patient_Sample_Combined_RISK)

################ Extrathyroidal Extension ####################
plot <- ggplot(Patient_Sample_Combined, aes(EXTRATHYROIDAL_EXTENSION, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = EXTRATHYROIDAL_EXTENSION),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Extrathyroidal Extension", limits = c("None", "Minimal (T3)", "Moderate/Advanced (T4a)", "Very Advanced (T4b)")) + 
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red") 
plot

### Extrathyroidal Modified ### 
# I am going to simplfy the EXTRATHYROIDAL EXTENSION category 
# Specifically, I will merge "Moderate/Advanced (T4a)" with "Very Advanced (T4b)" to form "Moderate/Advanced (T4)
Patient_Sample_Combined$Extrathyroidal_Modified <- Patient_Sample_Combined$EXTRATHYROIDAL_EXTENSION
for(i in 1:nrow(Patient_Sample_Combined)){
  if(Patient_Sample_Combined$Extrathyroidal_Modified[i] == "Moderate/Advanced (T4a)"){
    Patient_Sample_Combined$Extrathyroidal_Modified[i] <- "T4"
  }
  else if(Patient_Sample_Combined$Extrathyroidal_Modified[i] == "Very Advanced (T4b)"){
    Patient_Sample_Combined$Extrathyroidal_Modified[i] <- "T4"
  }
  else if(Patient_Sample_Combined$Extrathyroidal_Modified[i] == "Minimal (T3)"){
    Patient_Sample_Combined$Extrathyroidal_Modified[i] <- "T3"
  }
}

# Extrathyroidal Modified Extension Stats
Extrathyroidal_Stats <- Patient_Sample_Combined %>% subset(Extrathyroidal_Modified == "None" |
                                                             Extrathyroidal_Modified == "T3" |
                                                             Extrathyroidal_Modified == "T4")

# Extrathyroidal Simplified
Extrathyroidal_Stats$Extrathyroidal_Simplified <- Extrathyroidal_Stats$Extrathyroidal_Modified
for(i in 1:nrow(Extrathyroidal_Stats)){
  if(Extrathyroidal_Stats$Extrathyroidal_Simplified[i] == "None"){
    Extrathyroidal_Stats$Extrathyroidal_Simplified[i] <- "NO"
  }
  else if(Extrathyroidal_Stats$Extrathyroidal_Simplified[i] == "T3" |
          Extrathyroidal_Stats$Extrathyroidal_Simplified[i] == "T4"){
    Extrathyroidal_Stats$Extrathyroidal_Simplified[i] <- "YES"
  }
}
table(Extrathyroidal_Stats$Extrathyroidal_Simplified)
# NO = 309; YES = 132

# Now I will plot the new extrathyroidal modified category 

# Save csv for JCI
csv <- Extrathyroidal_Stats[c("PATIENT_ID", "SAMPLE_ID", "Extrathyroidal_Simplified", "myCAF")]
write_csv(csv, file = "data_in_use/24-1130_TCGA_Extrathyroidal_JCI.csv")

plot <- ggplot(Extrathyroidal_Stats, aes(Extrathyroidal_Simplified, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = Extrathyroidal_Simplified),
               alpha = 0.2,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c("lightgrey", "black")) + 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "Extrathy. Ext.", limits = c("NO", "YES")) + 
  labs(y = "myCAF Score") + 
  theme_classic() + 
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75)  + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) 
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_Extrathyroidal_simplified.png",
       plot, 
       width = 3, height = 3, dpi = 600)

# Stats for myCAF score by Extrathyroidal Extension
kruskal.test(myCAF ~ Extrathyroidal_Modified, data = Extrathyroidal_Stats) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Extrathyroidal_Stats$myCAF, Extrathyroidal_Stats$Extrathyroidal_Modified, # this returned 0.2475
                     p.adjust.method = "BH")

# Simplified
kruskal.test(myCAF ~ Extrathyroidal_Simplified, data = Extrathyroidal_Stats) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Extrathyroidal_Stats$myCAF, Extrathyroidal_Stats$Extrathyroidal_Simplified, # this returned 0.2475
                     p.adjust.method = "BH")
# Cleaning up from extrathyroidal extension
rm(Extrathyroidal_Stats)

############################### PTC HISTOLOGY SUBTYPES ###################################

### Altering histological type 
# Here I will look at histological type and how it compares to both the myCAF score and how it compares to BRS
Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX <- Patient_Sample_Combined$HISTOLOGICAL_TYPE
for(i in 1:nrow(Patient_Sample_Combined)){
  if(Patient_Sample_Combined$HISTOLOGICAL_TYPE[i] == ""){
    Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Unknown\nNot Specified"
  }
  else if(Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse Sclerosing Variant" |
          Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "Diffuse sclerosing variant" |
          Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "Papillary carcinoma, diffuse sclerosing variant"){
    Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Diffuse\nSclerosing"
  }
  else if(Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "cribriform morular" |
          Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "cribiform morular"){
    Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Cribriform\nMorular"
  }
  else if(Patient_Sample_Combined$HISTOLOGICAL_TYPE_OTHER[i] == "encapsulated follicular"){
    Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Follicular" 
  }
  else if(Patient_Sample_Combined$HISTOLOGICAL_TYPE[i] == "Tall Cell"){
    Patient_Sample_Combined$HISTOLOGICAL_TYPE_COMPLEX[i] <- "Tall\nCell"
  }
}

# Subset the data to just the following groups: 
# Follicular 
# Diffuse\nSclerosing
# Tall Cell
# Classical
# Unknown\nNot Specified
Histology_Plot_Cohort <- Patient_Sample_Combined %>% subset(HISTOLOGICAL_TYPE_COMPLEX == "Classical" |
                                                            HISTOLOGICAL_TYPE_COMPLEX == "Diffuse\nSclerosing" |
                                                            HISTOLOGICAL_TYPE_COMPLEX == "Follicular" | 
                                                            HISTOLOGICAL_TYPE_COMPLEX == "Tall\nCell")
                                                            #HISTOLOGICAL_TYPE_COMPLEX == "Unknown\nNot Specified")
# Histological Type
plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_COMPLEX, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = HISTOLOGICAL_TYPE_COMPLEX),
               alpha = 0.4,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = c(0.8),
              show.legend = FALSE) +
  scale_colour_manual(values = c("black", "red", "black")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "Diffuse\nSclerosing", "Tall\nCell")) + 
  labs(y = "myCAF Score") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) 
plot

### Make a new plot combining diffuse sclerosing and tall cell
Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED <- Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX
for(i in 1:nrow(Histology_Plot_Cohort)){
  if(Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Diffuse\nSclerosing" |
     Histology_Plot_Cohort$HISTOLOGICAL_TYPE_COMPLEX[i] == "Tall\nCell"){
    Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED[i] <- "DST"
  }
}

# save csv for JCI
csv <- Histology_Plot_Cohort[c("PATIENT_ID", "SAMPLE_ID", "HISTOLOGICAL_TYPE_SIMPLIFIED", "myCAF")]
write_csv(csv, file = "data_in_use/24-1130_TCGA_Histology_JCI.csv")

plot <- ggplot(Histology_Plot_Cohort, aes(HISTOLOGICAL_TYPE_SIMPLIFIED, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = HISTOLOGICAL_TYPE_SIMPLIFIED),
               alpha = 0.3,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = c(0.8),
              show.legend = FALSE) +
  scale_fill_manual(values = c("#FFCC99", "black", "red")) + 
  scale_x_discrete(name = "Histology", limits = c("Follicular", "Classical", "DST")) + 
  labs(y = "myCAF Score") + 
  #ggtitle("SFRP2 Score\nby Histology") +
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  theme_classic() + 
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) 
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_TCGA_Histology_Breakdown_by_myCAF.png",
       plot, 
       width = 4, height = 3, dpi = 600)

# Stats for myCAF score by HISTOLOGICAL_TYPE_SIMPLIFIED
kruskal.test(myCAF ~ HISTOLOGICAL_TYPE_SIMPLIFIED, data = Histology_Plot_Cohort) # shouldn't be necessary for only two groups
pairwise.wilcox.test(Histology_Plot_Cohort$myCAF, Histology_Plot_Cohort$HISTOLOGICAL_TYPE_SIMPLIFIED, # this returned 0.2475
                     p.adjust.method = "BH")
rm(Histology_Plot_Cohort, plot)

######################### LYMPH NODE ASSOCIATION #############################
table(Patient_Sample_Combined$PATH_N_STAGE)
LN_Plot_Cohort <- Patient_Sample_Combined %>% subset(PATH_N_STAGE == "N0" |
                                                     PATH_N_STAGE == "N1" |
                                                     PATH_N_STAGE == "N1a" |
                                                     PATH_N_STAGE == "N1b")

LN_Plot_Cohort$LN_Binary <- ifelse(LN_Plot_Cohort$PATH_N_STAGE == "N0", "NO", "YES") # set LN_Binary as either NO (N0) or YES (N1, N1a, N1b)

# save csv for JCI
csv <- LN_Plot_Cohort[c("PATIENT_ID", "SAMPLE_ID", "LN_Binary", "myCAF")]
write_csv(csv, file = "data_in_use/24-1130_TCGA_LN_Mets_JCI.csv")

# Plot binary category and perform stats
plot <- ggplot(LN_Plot_Cohort, aes(LN_Binary, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = LN_Binary),
               alpha = 0.2,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c("lightgrey", "blue")) + 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "LN", limits = c("NO", "YES")) + 
  labs(y = "myCAF Score") + 
  theme_classic() + 
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75)  + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) 
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_LN_Spread.png",
       plot, 
       width = 3, height = 3, dpi = 600)

# LN stats
kruskal.test(myCAF ~ LN_Binary, data = LN_Plot_Cohort) # shouldn't be necessary for only two groups
pairwise.wilcox.test(LN_Plot_Cohort$myCAF, LN_Plot_Cohort$LN_Binary, #
                     p.adjust.method = "BH")
# Cleaning up from extrathyroidal extension
rm(LN_Plot_Cohort)


######################### DISTANT MET ASSOCIATION #############################
table(Patient_Sample_Combined$PATH_M_STAGE)
DM_Plot_Cohort <- Patient_Sample_Combined %>% subset(PATH_M_STAGE == "M0" |
                                                     PATH_M_STAGE == "M1")

# save csv for JCI
csv <- DM_Plot_Cohort[c("PATIENT_ID", "SAMPLE_ID", "PATH_M_STAGE", "myCAF")]
write_csv(csv, file = "data_in_use/24-1130_TCGA_Distant_Mets_JCI.csv")

# Plot binary category and perform stats
plot <- ggplot(DM_Plot_Cohort, aes(PATH_M_STAGE, myCAF)) + 
  geom_boxplot(outlier.size = -1,
               aes(fill = PATH_M_STAGE),
               alpha = 0.2,
               show.legend = FALSE) + 
  geom_jitter(aes(),
              position = position_jitter(width = 0.1, height = 0),
              size = 1.5,
              alpha = 0.7,
              show.legend = FALSE) + 
  scale_fill_manual(values = c("lightgrey", "red")) + 
  # scale_colour_manual(values = c("black", "red")) + 
  scale_x_discrete(name = "LN", limits = c("M0", "M1")) + 
  labs(y = "myCAF Score") + 
  theme_classic() + 
  scale_y_continuous(breaks = c(-.5, 0, .5),
                     limits = c(-0.6, .8)) +
  geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75)  + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 16)) 
plot
ggsave("outputs/24-1130_myCAF_TCGA_Survival_Extrathyroidal/24-1130_myCAF_by_Distant_Met_Spread.png",
       plot, 
       width = 3, height = 3, dpi = 600)

# LN stats
kruskal.test(myCAF ~ PATH_M_STAGE, data = DM_Plot_Cohort) # shouldn't be necessary for only two groups
pairwise.wilcox.test(DM_Plot_Cohort$myCAF, DM_Plot_Cohort$PATH_M_STAGE, #
                     p.adjust.method = "BH") # 0.43
# Cleaning up from extrathyroidal extension
rm(DM_Plot_Cohort)





