### Author: Matthew Aaron Loberg
### Date: November 24, 2024
### Script: MERINGUE_PTC_myCAF_Pericyte_Spatial_Correlation_Analysis.R
### Source Script Name: 24-1124_MERINGUE_PTC_myCAF_Pericyte_Spatial_Correlation_Analysis.R

### Adapted form the following script
# Script: 24-1123_MERINGUE_PTC_pEMT_myCAF_Pericyte_Spatial_Correlation_Analysis.R

### Goal: Spatial correlation with MERINGUE on PTCs (Thy5, Thy7, Thy15, Thy16, Thy17)
# The goal here is to look at the localization of pericytes in PTC
# I want to know if the pericyte score is localizing within tumors (e.g., within fibrovascular cores of tumors) or on the edge of tumor (e.g., stroma)
# Also ALL 8 pediatric PTCs (Peds01 - Peds08)
# Also including PTC/ATC mixed specimens (Thy1, Thy5, Thy6, Thy10, Thy11)
# For stats purposes to look at the most pure sample, I will subset down to adult PTCs (Thy7, Thy15, Thy16, Thy17) and Peds PTCs (Peds01-Peds08)

# In this script, I am using the 24-1001 RCTD distinction with NO pEMT designation (Just PTC or ATC)

# See package info here for MERINGUE spatial cross-correlation: 
# https://github.com/JEFworks-Lab/MERINGUE

### Load packages
library(corrplot) # for correlation plots
library(MERINGUE) # for spatial correlation statistic
library(Seurat) # Load Seurat objects
library(tidyverse) # ggsave, others

### read directory list
# The read directories will specifically be PTC Visium samples
readdirs <- c("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy1.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy5.RDS",
              "Data_in_USe/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy6.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy7.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy10.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy11.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy15.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy16_V2.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy17.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds01.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds02.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds03.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds04.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds05.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds06.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds07.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Peds08.RDS")

save_vector <- c("1", "5", "6", "7", "10", "11", "15", "16", "17", "01", "02", "03", "04", "05", "06", "07", "08")

spatial_SO_list <- list()
for(i in 1:length(readdirs)){
  spatial_SO_list[[i]] <- readRDS(file = readdirs[i])
}

filterDist <- c(rep(7.5, times = 17))
outputdir = c(rep("outputs/MERINGUE/Thy", times = 9), rep("outputs/MERINGUE/Peds", times = 8))
sample <- c(rep("Thy", times = 9), rep("Peds", times = 8))

correlations <- list()
for(i in 1:length(spatial_SO_list)){
  # Pull out data in a new data.frame (correlation_data_frame) to perform a correlation test on 
  correlation_data_frame <- data.frame(
      PTC_RCTD = spatial_SO_list[[i]]$PTC_RCTD,
      myCAF_RCTD = spatial_SO_list[[i]]$myCAF_RCTD,
      Pericyte_RCTD = spatial_SO_list[[i]]$Pericyte_RCTD
  )
  
  # make correlation_data_frame into a matrix + transpose it for optimal formatting for running MERINGUE
  mat <- t(as.matrix(correlation_data_frame))
  
  # pull out position data
  pos <- GetTissueCoordinates(spatial_SO_list[[i]])
  
  # Get spatial neighbors from position data
  w <- getSpatialNeighbors(pos, filterDist = filterDist[i])
  
  # Save a network of spatial neighbors (done prior, not repeating)
  #png(file = paste0(outputdir[i], save_vector[i], "/plotNetwork.png"), 
  #                  res = 300, height = 1920, width = 1920)
  #print(plotNetwork(pos, w))
  #dev.off()
  
  # Get spatial patterns
  I <- getSpatialPatterns(mat, w)
  
  results.filter <- filterSpatialPatterns(mat = mat, # Commenting out because I want to view ALL of the correlations
                                           I = I,
                                           w = w,
                                           adjustPv = TRUE,
                                           alpha = 0.05,
                                           minPercentCells = 0.05,
                                           verbose = TRUE)
  
  scc <- spatialCrossCorMatrix(mat = mat, 
                               weight = w)
  
  png(file = paste0(outputdir[i], save_vector[i], "/24-1124_PTC_myCAF_Pericyte_corrplot.png"), 
                    res = 300, height = 1920, width = 1920)
  print(corrplot(scc, 
                 method = "color", 
                 order = "hclust",
                 col = rev(COL2("RdBu", 200)), # takes the RdBu scale that is used in the packaged and flips it
                 tl.col = "black",
                 outline = TRUE) )
  dev.off()
  
  # save a data.frame of correlations for each sample -> will use for downstream stats
  correlations[[i]] <- data.frame(
       Sample = c(paste0(sample[i], save_vector[i])),
       PTC_Pericyte = c(scc[1,3]),
       myCAF_Pericyte = c(scc[2,3]),
       PTC_myCAF = c(scc[1,2])
  )

}

correlations_overview <- rbind(correlations[[1]],
                               correlations[[2]],
                               correlations[[3]],
                               correlations[[4]],
                               correlations[[5]],
                               correlations[[6]],
                               correlations[[7]],
                               correlations[[8]],
                               correlations[[9]],
                               correlations[[10]],
                               correlations[[11]],
                               correlations[[12]],
                               correlations[[13]],
                               correlations[[14]],
                               correlations[[15]],
                               correlations[[16]],
                               correlations[[17]])

# Need to pivot longer the correlations_overview
Correlations_overview_plotting <- correlations_overview %>% 
  pivot_longer(
    cols = starts_with("PTC_Pericyte") | starts_with("myCAF_Pericyte") | starts_with("PTC_myCAF"),
    names_to = c("Correlation_Group"),
    values_to = "Spatial_Correlation"
  )

# Plot Pericyte correlations
Pericyte_Correlations <- Correlations_overview_plotting %>% subset(Correlation_Group == "PTC_Pericyte" | Correlation_Group == "myCAF_Pericyte")


# ATCs excluded (see goals section)
Pericyte_Correlations_ATC_Excluded <- Pericyte_Correlations %>% subset(Sample != "Thy1" &
                                                                   Sample != "Thy5" &
                                                                   Sample != "Thy6" &
                                                                   Sample != "Thy10" &
                                                                   Sample != "Thy11")

max(Pericyte_Correlations_ATC_Excluded$Spatial_Correlation)
min(Pericyte_Correlations_ATC_Excluded$Spatial_Correlation)
# Go with max .4, min .4
plot <- ggplot(Pericyte_Correlations_ATC_Excluded, aes(Correlation_Group, Spatial_Correlation)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Correlation_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             #position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  geom_line(aes(group = Sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("#783FC1", "#FFCC99")) +
  labs (x = "Tumor", y = "Pericyte Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  #scale_x_discrete(name ="Diagnosis", limits = c("Not")) +
  scale_y_continuous(breaks = c(-.4, -.2, 0, .2, .4),
                     limits = c(-.4, .40)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/Pericyte_Correlations/24-1124_Pericyte_PTC_myCAF_Correlation_Boxplots_ATC_Exclude.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

# stats
Pericyte_ATC_Exclude_wide <- pivot_wider(Pericyte_Correlations_ATC_Excluded, names_from = Correlation_Group, values_from = Spatial_Correlation)
t.test(Pericyte_ATC_Exclude_wide$PTC_Pericyte, Pericyte_ATC_Exclude_wide$myCAF_Pericyte, paired = TRUE) 

## Could exclude additional based on histology: 
# Peds01, Peds02, etc. that do NOT have large PTC areas (or only pEMT in vessels); will be conservative and leave in for now
