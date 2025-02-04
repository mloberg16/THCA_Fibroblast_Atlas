### Author: Matthew Aaron Loberg
### Date: November 23, 2024
### Script: MERINGUE_PTC_pEMT_myCAF_Spatial_Correlation_Analysis.R
### Source Script Name: 24-1123_MERINGUE_PTC_pEMT_myCAF_iPVCAF_Spatial_Correlation_Analysis.R

# Goal: Spatial correlation with MERINGUE on PTCs (Thy7, Thy15, Thy16, Thy17)
# Also ALL 8 pediatric PTCs
# Goal is to look at how the pEMT-PTC deconvolution aligns with myCAF relative to the PTC deconvolution w/ no pEMT

# See package info here: 
# https://github.com/JEFworks-Lab/MERINGUE

# Load packages
library(corrplot) # for correlation plots
library(MERINGUE) # for spatial correlation statistic
library(Seurat) # Load Seurat objects
library(tidyverse) # ggsave, others

# read directory list
# The read directories will specifically be PTC Visium samples (note: Thy1, Thy5, Thy10, Thy11 contain ATC as well and were later excluded)
readdirs <- c("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy1.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy5.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy7.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy10.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy11.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy15.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy16_V2.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy17.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds01.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds02.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds03.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds04.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds05.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds06.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds07.RDS",
              "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds08.RDS")

save_vector <- c("1", "5", "7", "10", "11", "15", "16", "17", "01", "02", "03", "04", "05", "06", "07", "08")

spatial_SO_list <- list()
for(i in 1:length(readdirs)){
  spatial_SO_list[[i]] <- readRDS(file = readdirs[i])
}

filterDist <- c(rep(7.5, times = 16))
outputdir = c(rep("outputs/MERINGUE/Thy", times = 8), rep("outputs/MERINGUE/Peds", times = 8))
sample <- c(rep("Thy", times = 8), rep("Peds", times = 8))

correlations <- list()
for(i in 1:length(spatial_SO_list)){
  # Pull out data in a new data.frame (correlation_data_frame) to perform a correlation test on 
  correlation_data_frame <- data.frame(
      PTC_RCTD = spatial_SO_list[[i]]$PTC_RCTD,
      pEMT_RCTD = spatial_SO_list[[i]]$pEMT_PTC_RCTD,
      myCAF_RCTD = spatial_SO_list[[i]]$myCAF_RCTD
  )
  
  # make correlation_data_frame into a matrix + transpose it for optimal formatting for running MERINGUE
  mat <- t(as.matrix(correlation_data_frame))
  
  # pull out position data
  pos <- GetTissueCoordinates(spatial_SO_list[[i]])
  
  # Get spatial neighbors from position data
  w <- getSpatialNeighbors(pos, filterDist = filterDist[i])
  
  # Save a network of spatial neighbors
  png(file = paste0(outputdir[i], save_vector[i], "/plotNetwork.png"), 
                    res = 300, height = 1920, width = 1920)
  print(plotNetwork(pos, w))
  dev.off()
  
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
  
  png(file = paste0(outputdir[i], save_vector[i], "/corrplot.png"), 
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
       PTC_myCAF = c(scc[1,3]),
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
                               correlations[[8]],
                               correlations[[9]],
                               correlations[[10]],
                               correlations[[11]],
                               correlations[[12]],
                               correlations[[13]],
                               correlations[[14]],
                               correlations[[15]],
                               correlations[[16]])

# Need to pivot longer the correlations_overview
Correlations_overview_plotting <- correlations_overview %>% 
  pivot_longer(
    cols = starts_with(starts_with("pEMT_myCAF") | starts_with("PTC_myCAF"),
    names_to = c("Correlation_Group"),
    values_to = "Spatial_Correlation"
  )

# Plot myCAF correlations
myCAF_Correlations <- Correlations_overview_plotting %>% subset(Correlation_Group == "PTC_myCAF" | Correlation_Group == "pEMT_myCAF")

max(myCAF_Correlations$Spatial_Correlation)
min(myCAF_Correlations$Spatial_Correlation)
# Go with max -.7, .5
plot <- ggplot(myCAF_Correlations, aes(Correlation_Group, Spatial_Correlation)) +
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
  scale_fill_manual(values = c("grey", "#990000")) +
  labs (x = "Tumor", y = "myCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Tumor", limits = c("PTC_myCAF", "pEMT_myCAF")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4),
                     limits = c(-.7, .50)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/myCAF_Correlations/24-1123_myCAF_PTC_pEMT_Correlation_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)


# No connected line
plot <- ggplot(myCAF_Correlations, aes(Correlation_Group, Spatial_Correlation)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Correlation_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  #geom_line(aes(group = Sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("grey", "#990000")) +
  labs (x = "Tumor", y = "myCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Tumor", limits = c("PTC_myCAF", "pEMT_myCAF")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4),
                     limits = c(-.7, .50)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/myCAF_Correlations/24-1123_myCAF_PTC_pEMT_Correlation_Boxplots_no_connections.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)


# ATCs excluded
myCAF_Correlations_ATC_Excluded <- myCAF_Correlations %>% subset(Sample != "Thy1" &
                                                                 Sample != "Thy5" &
                                                                 Sample != "Thy10" &
                                                                 Sample != "Thy11")


# save as .csv for JCI
write_csv(myCAF_Correlations_ATC_Excluded, file = "Data_in_Use/24-1123_MERINGUE_PTC_pEMT_myCAF_Correlations.csv")

max(myCAF_Correlations_ATC_Excluded$Spatial_Correlation)
min(myCAF_Correlations_ATC_Excluded$Spatial_Correlation)
plot <- ggplot(myCAF_Correlations_ATC_Excluded, aes(Correlation_Group, Spatial_Correlation)) +
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
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Tumor", limits = c("PTC_myCAF", "pEMT_myCAF")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4),
                     limits = c(-.7, .50)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/myCAF_Correlations/24-1123_myCAF_PTC_pEMT_Correlation_Boxplots_ATC_Exclude.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)


# No connected line
plot <- ggplot(myCAF_Correlations_ATC_Excluded, aes(Correlation_Group, Spatial_Correlation)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Correlation_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  #geom_line(aes(group = Sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Tumor", limits = c("PTC_myCAF", "pEMT_myCAF")) +
  scale_y_continuous(breaks = c(-.6, -.4, -.2, 0, .2, .4),
                     limits = c(-.7, .50)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/MERINGUE/myCAF_Correlations/24-1123_myCAF_PTC_pEMT_Correlation_Boxplots_ATC_Exclude_no_connections.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

# stats
myCAF_ATC_Exclude_wide <- pivot_wider(myCAF_Correlations_ATC_Excluded, names_from = Correlation_Group, values_from = Spatial_Correlation)
t.test(myCAF_ATC_Exclude_wide$pEMT_myCAF, myCAF_ATC_Exclude_wide$PTC_myCAF, paired = TRUE) # two-sided t-test: 0.006445
