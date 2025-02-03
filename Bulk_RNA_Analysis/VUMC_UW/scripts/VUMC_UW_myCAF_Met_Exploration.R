### Author: Matthew Aaron Loberg
### Date: November 30, 2024
### Script: VUMC_UW_myCAF_Met_Exploration.R
### Source Script Name: 24-1130_CGRNA_myCAF_Met_Exploration.R

### Goal: 
# Look at associations of lymph node met/distant met, etc. with myCAF score in the VUMC/UW cohort

### Load required packages
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)

### Set WD (or open script as a project within defined directory)
# Set your working directory if needed (I do not as working within 2024_Bulk_RNA_Deconvolution_Analysis R Project)
# setwd("_") # Input in the "_" the working directory that you wish to assign

### Load Data
# Note, while doing this, I decided to switch to the most recent version of the clinical file
# My concern was that previous versions of the clinical file may be outdated with regard to some of the met categories
# I switched from "VUMC.cohort.GX_9-8-22_v2.csv" to the most recent in George's folder, which is
# I have now updated to my own most recent: 5-31-24
ClinicalData <- as_tibble(read_csv("data_in_use/clinical_cohort_data/VUMC.cohort.MAL_5-31-24.csv")) # All Data

# Read in myCAF score data
gsva_results <- readRDS(file = "data_in_use/deconvolution_scores/24-0903_Thyroid_Fibroblast_GSVA_Scores.RDS")
gsva_results$RNA.ID <- sub('.', '', gsva_results$RNA.ID) # Get rid of the 'X" at the start of GSVA_Results


# Merge normalized counts and clinical data
cohort <- merge(x = ClinicalData, y = gsva_results, by.x = "RNA.ID", by.y = "RNA.ID")
rm(ClinicalData, gsva_results)

# Change "Distant.met.patient" to "YES" or "NO"
for(i in 1:nrow(cohort)){
  if(!is.na(cohort$Distant.met.patient[i])){
    if(cohort$Distant.met.patient[i] == 0){
      cohort$Distant.met.patient[i] <- "NO"
    }
    else if(cohort$Distant.met.patient[i] == 1){
      cohort$Distant.met.patient[i] <- "YES"
    }
  }
  if(!is.na(cohort$Local.met.patient[i])){
    if(cohort$Local.met.patient[i] == "0"){
      cohort$Local.met.patient[i] <- "NO"
    }
    else if(cohort$Local.met.patient[i] == "1"){
      cohort$Local.met.patient[i] <- "YES"
    }
  }
  
  if(!is.na(cohort$Extrathyroidal.extension[i])){
    if(cohort$Extrathyroidal.extension[i] == "E" | 
       cohort$Extrathyroidal.extension[i] == "M" |
       cohort$Extrathyroidal.extension[i] == "P"){
          cohort$Extrathyroidal.extension[i] <- "YES"
    }
    else if(cohort$Extrathyroidal.extension[i] == "X"){
          cohort$Extrathyroidal.extension[i] <- "NO"
    }
    else if(cohort$Extrathyroidal.extension[i] == "I" | 
            cohort$Extrathyroidal.extension[i] == "NR"){
          cohort$Extrathyroidal.extension[i] <- "Unknown"
    }
  }
  
}

####### 24-1024 UPDATE #######
# Removing local disease samples that are NOT actually primary
cohort <- cohort %>% subset(Location.type.detailed != "Localdisease" | PFS_Sample != "NO")


# Make a function for all of the plots that I want to make
MetPlots <- function(cohort, output_folder){
    
    # Distant met yes/no for PRIMARIES only
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" & !is.na(Distant.met.patient)), aes(Distant.met.patient, myCAF)) +
      geom_boxplot(outlier.size = -1, 
                   aes(fill = Distant.met.patient),
                   alpha = 0.4, 
                   show.legend = FALSE) + 
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, 
                  alpha = 0.7,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("lightgrey", "red")) +
      geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
      labs (x ="Met patient", y = "MyCAF Score") + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) +
      scale_x_discrete(name ="Distant Met Patient", limits = c("NO", "YES")) +
      scale_y_continuous(breaks = c(-.5, 0, .5),
                         limits = c(-0.6, .85))
    ggsave(filename = file.path(output_folder, "Distant_Met_Patient_LocalDisease_Only.png"),
           plot, height = 3, width = 3, dpi = 600)
    
    test <- cohort %>% subset(Location.type == "Primary" & !is.na(Distant.met.patient))
    kruskal <- kruskal.test(myCAF ~ Distant.met.patient, data = test)
    pairwise <- pairwise.wilcox.test(test$myCAF, test$Distant.met.patient, 
                                     p.adjust.method = "bonferroni")
    rm(test)
    
    # Local met yes/no for PRIMARIES only
    plot <- ggplot(cohort %>% subset(Location.type == "Primary"), aes(Local.met.patient, myCAF)) +
      geom_boxplot(outlier.size = -1, 
                   aes(fill = Local.met.patient),
                   alpha = 0.4, 
                   show.legend = FALSE) + 
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, 
                  alpha = 0.7,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("lightgrey", "blue")) +
      geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
      labs (x ="Met patient", y = "myCAF") + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) +
      scale_x_discrete(name ="Local Met Patient", limits = c("NO", "YES")) +
      scale_y_continuous(breaks = c(-.5, 0, .5),
                         limits = c(-0.6, .85))
    ggsave(filename = file.path(output_folder, "Local_Met_Patient_LocalDisease_Only.png"),
           plot, height = 3, width = 3, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary")
    kruskal <- kruskal %>% append(kruskal.test(myCAF ~ Local.met.patient, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(test$myCAF, test$Local.met.patient, 
                                                         p.adjust.method = "bonferroni"))
    rm(test)
    
    # Samples by Location.type (for primary and LN met only)
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet"), aes(Location.type, myCAF)) +
                     geom_boxplot(outlier.size = -1, 
                                  aes(fill = Location.type),
                                  alpha = 0.4, 
                                  show.legend = FALSE) + 
                     geom_jitter(aes(),
                                 position = position_jitter(width = 0.1, height = 0),
                                 size = 1.5, 
                                 alpha = 0.7,
                                 show.legend = FALSE) +
                     scale_fill_manual(values = c("lightgrey", "blue")) + # can modify here to define colors
                     geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
                     labs (x ="Location Simplified", y = "myCAF Score") + 
                     theme_classic() + 
                     theme(
                       plot.title = element_text(face = "bold", hjust = 0.5), 
                       axis.title.x = element_text(face = "bold", size = 14.5), 
                       axis.title.y = element_text(face = "bold", size = 12),
                       axis.text.x = element_text(face = "bold", size = 10.5),
                       axis.text.y = element_text(face = "bold", size = 16)) +
                     scale_x_discrete(name ="Location Simplified ", limits = c("Primary", "Localmet")) +
      scale_y_continuous(breaks = c(-.5, 0, .5),
                         limits = c(-0.6, .85))
    ggsave(filename = file.path(output_folder, "Location_Simplified.png"),
           plot, height = 3, width = 3, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet")
    kruskal <- kruskal %>% append(kruskal.test(myCAF ~ Location.type, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(test$myCAF, test$Location.type, 
                                                         p.adjust.method = "bonferroni"))
    
    # Samples by Location.type (for primary, LN, and distant met only)
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet" | Location.type == "Distantmet"), aes(Location.type, myCAF)) +
      geom_boxplot(outlier.size = -1, 
                   aes(fill = Location.type),
                   alpha = 0.4, 
                   show.legend = FALSE) + 
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, 
                  alpha = 0.7,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("red", "blue", "lightgrey")) + # can modify here to define colors
      geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
      labs (x ="Location Simplified", y = "myCAF Score") + 
      theme_classic() + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) +
      scale_x_discrete(name ="Location Simplified ", limits = c("Primary", "Localmet", "Distantmet")) +
      scale_y_continuous(breaks = c(-.5, 0, .5),
                         limits = c(-0.6, .85))
    ggsave(filename = file.path(output_folder, "Location_Prim_LN_Dist_Simplified.png"),
           plot, height = 3, width = 4, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary" | Location.type == "Localmet" | Location.type == "Distantmet")
    kruskal <- kruskal %>% append(kruskal.test(myCAF ~ Location.type, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(test$myCAF, test$Location.type, 
                                                         p.adjust.method = "bonferroni"))
    
    # Plot extrathyroidal extension
    plot <- ggplot(cohort %>% subset(Location.type == "Primary" & 
                                     !is.na(Extrathyroidal.extension) & 
                                     Extrathyroidal.extension != "Unknown"), 
                   aes(Extrathyroidal.extension, myCAF)) +
      geom_boxplot (outlier.size = -1,
                    aes(fill = Extrathyroidal.extension),
                    alpha = 0.4,
                    show.legend = FALSE) +
      geom_jitter(aes(),
                  position = position_jitter(width = 0.1, height = 0),
                  size = 1.5, alpha = 0.7, show.legend = FALSE) +
      scale_fill_manual(values = c("lightgrey", "black")) + # can modify here to define colors
      geom_hline(yintercept = 0.0, linetype = 2, colour = "red", size = 0.75) + 
      labs(x = "Extrathy. Ext.", y = "myCAF Score") +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 14.5), 
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 16)) + 
      scale_x_discrete(name = "Extrathy. Ext.", limits = c("NO", "YES")) +
      scale_y_continuous(breaks = c(-.5, 0, .5),
                         limits = c(-0.6, .85))
    ggsave(filename = file.path(output_folder, "Extrathyroidal_Extension.png"),
           plot, height = 3, width = 3, dpi = 600)
    test <- cohort %>% subset(Location.type == "Primary" & !is.na(Extrathyroidal.extension) & Extrathyroidal.extension != "Unknown")
    kruskal <- kruskal %>% append(kruskal.test(myCAF ~ Extrathyroidal.extension, data = test)) 
    pairwise <- pairwise %>% append(pairwise.wilcox.test(test$myCAF, test$Extrathyroidal.extension, 
                                                         p.adjust.method = "bonferroni"))
    
    # Return(NULL)
    return(list(kruskal, pairwise))
}

# Create an output directory: 
dir.create("outputs/24-1130_CGRNA_myCAF_Met_Plots")

# MetPlots malignant cohort (NIFTP Excluded)
output_folder <- "outputs/24-1130_CGRNA_myCAF_Met_Plots/Malignant_NIFTP_Excluded/"
dir.create(output_folder)
stats <- MetPlots(cohort = cohort %>% subset(Diagnosis != "MNG" & Diagnosis != "FA" &
                                             Diagnosis != "HA" & Diagnosis != "HT" & Diagnosis != "NIFTP"),
                            output_folder = output_folder)
saveRDS(stats, file.path(output_folder, "24-1130_Malignant_Samples_NIFTP_Excluded_Stats_List.RDS"))


# MetPlots WDTC Cohort (NIFTP Excluded)
output_folder <- "outputs/24-1130_CGRNA_myCAF_Met_Plots/WDTC_NIFTP_Excluded/"
dir.create(output_folder)
stats <- MetPlots(cohort = cohort %>% subset(Diagnosis != "MNG" & Diagnosis != "FA" &
                                             Diagnosis != "HA" & Diagnosis != "HT" & 
                                             Diagnosis != "PDTC" & Diagnosis != "ATC" & Diagnosis != "NIFTP"),
                  output_folder = output_folder)
saveRDS(stats, file.path(output_folder, "24-1130_WDTC_NIFTP_Excluded_Stats_List.RDS"))

##### Cleaning UP ######
rm(list = ls())
