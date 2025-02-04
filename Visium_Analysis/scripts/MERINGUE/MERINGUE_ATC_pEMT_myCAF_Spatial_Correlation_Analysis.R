### Author: Matthew Aaron Loberg
### Date: December 25, 2024
### Script: MERINGUE_ATC_pEMT_myCAF_Spatial_Correlation_Analysis.R
### Source Script Name: 24-1225_MERINGUE_ATC_pEMT_myCAF_Spatial_Correlation_Analysis.R

# Goal: Spatial correlation with MERINGUE on ATCs; looking at pEMT/myCAF interaction

# See package info here: 
# https://github.com/JEFworks-Lab/MERINGUE

# 24-1225 Update
# Repeating for ATC samples looking at pEMT-ATC and ATC spatial correlation with myCAF
# This script was adapted from 24-1123_MERINGUE_PTC_pEMT_myCAF_iPVCAF_Spatial_Correlation_Analysis

# Load packages
library(corrplot) # for correlation plots
library(MERINGUE) # for spatial correlation statistic
library(Seurat) # Load Seurat objects
library(tidyverse) # ggsave, others

# read directory list
# The read directories will specifically be PTC Visium samples
readdirs <- c("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy1.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy4.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy5.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy6.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy8.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy10.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy11.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1130_Broad_Lu_ATC09_pEMT/Thy12.RDS")

save_vector <- c("1", "4", "5", "6", "8", "10", "11", "12")

spatial_SO_list <- list()
for(i in 1:length(readdirs)){
  spatial_SO_list[[i]] <- readRDS(file = readdirs[i])
}

filterDist <- c(rep(7.5, times = 8))
outputdir = c(rep("outputs/MERINGUE/Thy", times = 8))
sample <- c(rep("Thy", times = 8))

correlations <- list()
for(i in 1:length(spatial_SO_list)){
  # Pull out data in a new data.frame (correlation_data_frame) to perform a correlation test on 
  correlation_data_frame <- data.frame(
      ATC_RCTD = spatial_SO_list[[i]]$ATC_RCTD,
      pEMT_RCTD = spatial_SO_list[[i]]$pEMT_ATC_RCTD,
      myCAF_RCTD = spatial_SO_list[[i]]$myCAF_RCTD
  )
  
  # make correlation_data_frame into a matrix + transpose it for optimal formatting for running MERINGUE
  mat <- t(as.matrix(correlation_data_frame))
  
  # pull out position data
  pos <- GetTissueCoordinates(spatial_SO_list[[i]])
  
  # Get spatial neighbors from position data
  w <- getSpatialNeighbors(pos, filterDist = filterDist[i])
  
  # Save a network of spatial neighbors
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
  
  png(file = paste0(outputdir[i], save_vector[i], "/24-1225_myCAF_ATC_pEMT_ATC_corrplot.png"), 
                    res = 300, height = 1920, width = 1920)
  print(corrplot(scc, 
                 method = "color", 
                 order = "hclust",
                 col = rev(COL2("RdBu", 200)), # takes the RdBu scale that is used in the packaged and flips it
                 tl.col = "black",
                 outline = TRUE) )
  dev.off()
  
  # make a data.frame of 
  correlations[[i]] <- data.frame(
       Sample = c(paste0(sample[i], save_vector[i])),
       ATC_pEMT = c(scc[1,2]),
       ATC_myCAF = c(scc[1,3]),
       pEMT_myCAF = c(scc[2,3])
  )

}

correlations_overview <- rbind(correlations[[1]],
                               correlations[[2]],
                               correlations[[3]],
                               correlations[[4]],
                               correlations[[5]],
                               correlations[[6]],
                               correlations[[7]],
                               correlations[[8]])

# Need to pivot longer the correlations_overview
Correlations_overview_plotting <- correlations_overview %>% 
  pivot_longer(
    cols = starts_with("ATC_myCAF") | starts_with("pEMT_myCAF"),
    names_to = c("Correlation_Group"),
    values_to = "Spatial_Correlation"
  )

# Save as .csv for JCI
write_csv(Correlations_overview_plotting, file = "Data_in_Use/24-1225_MERINGUE_ATC_pEMTATC_myCAF_Correlations.csv")

# Plot correlations
max(Correlations_overview_plotting$Spatial_Correlation)
min(Correlations_overview_plotting$Spatial_Correlation)
# Go with max .4, min .4
plot <- ggplot(Correlations_overview_plotting %>% subset(Sample != "Thy6"), aes(Correlation_Group, Spatial_Correlation)) +
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
  scale_fill_manual(values = c("#0075DC", "darkgrey")) +
  labs (x = "Tumor", y = "myCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  #scale_x_discrete(name ="Diagnosis", limits = c("Not")) +
  scale_y_continuous(breaks = c(-.2, 0, .2, .4),
                     limits = c(-.3, .5)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/ATC_pEMT_myCAF_Correlations/24-1225_ATC_pEMT_myCAF_Correlation_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)



# stats
wide_stats <- pivot_wider(Correlations_overview_plotting %>% subset(Sample != "Thy6"), names_from = Correlation_Group, values_from = Spatial_Correlation)
t.test(wide_stats$ATC_myCAF, wide_stats$pEMT_myCAF, paired = TRUE) # two-sided t-test: 0.03706
