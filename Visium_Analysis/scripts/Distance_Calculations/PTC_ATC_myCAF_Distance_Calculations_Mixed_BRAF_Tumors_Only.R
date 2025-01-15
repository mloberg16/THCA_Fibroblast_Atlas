### Author: Matthew Aaron Loberg
### Date: December 02, 2024
### Script: PTC_ATC_myCAF_Distance_Calculations_Mixed_BRAF_Tumors_Only.R
### Source Script Name: 24-1202_PTC_ATC_myCAF_Distance_Calculations_Mixed_BRAF_Tumors_Only.R

# Goal:
# Calculate distance between ATC regions and myCAFs and PTC regions and myCAFs in tumors with BRAFV600E and mixed PTC/ATC morphologies

# Load packages
library(Seurat)
library(tidyverse)
source("scripts/Function_Scripts/24-1124_distance_comparison_function.R")

SO_readdirs <- c("Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy1.RDS",     # BRAFV600E PTC/ATC mix
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy5.RDS",     # BRAFV600E PTC/ATC mix
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy6.RDS",     # BRAFV600E PTC/ATC mix
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy10.RDS",    # BRAFV600E PTC/ATC mix
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1001_iCAF2_Excluded/Thy11.RDS")    # BRAFV600E PTC/ATC mix

samples <- c("Thy1", "Thy5", "Thy6", "Thy10", "Thy11")

mean_median_min_distances <- list()
for(i in 1:length(SO_readdirs)){
  
  # Load in spatial seurat object (SO) from readdirs
  spatial_SO <- readRDS(SO_readdirs[i])
  
  # Load in distance matrix for sample
  distance_matrix <- readRDS(file = paste0("Data_in_Use/Distance_Matrices/", samples[i], "_Distance_Matrix.RDS"))
  
  proportions <- data.frame(
      "ATC" = c(spatial_SO$ATC_RCTD),
      "PTC" = c(spatial_SO$PTC_RCTD),
      "myCAF" = c(spatial_SO$myCAF_RCTD)
  )
  
  proportions$Tumor <- proportions$PTC + proportions$ATC
  
  PTC_myCAF_distance <- distance_comparison(distances = distance_matrix,
                                             proportions = proportions,
                                             topic1 = "PTC",
                                             topic2 = "myCAF",
                                             t1_thresh = 0.1,
                                             t2_thresh = 0.1)
  
  ATC_myCAF_distance <- distance_comparison(distances = distance_matrix,
                                            proportions = proportions,
                                            topic1 = "ATC",
                                            topic2 = "myCAF",
                                            t1_thresh = 0.1,
                                            t2_thresh = 0.1)
  
  combined_distance <- rbind(PTC_myCAF_distance, ATC_myCAF_distance)
  
  distance_violin <- ggplot(combined_distance, aes(x = barcodes, y = min, fill = barcodes)) +
    geom_violin(scale = "width") +
    theme_classic() +
    scale_fill_manual(values = c("#0075DC", "#FFCC99")) +
    scale_x_discrete(name ="Deconvolution", limits = c("PTC vs. myCAF", "ATC vs. myCAF"))
    
    
  ggsave(file = paste0("outputs/Distance_Plots/", samples[i], "/PTC_ATC_myCAF_Min_Distance_Violin.png"),
         distance_violin, height = 5, width = 6, dpi = 600)
  
  mean_median_min_distances[[i]] <- data.frame(
      "sample" = c(samples[i]),
      "PTC_myCAF_mean" = c(mean(PTC_myCAF_distance$min)),
      "ATC_myCAF_mean" = c(mean(ATC_myCAF_distance$min)),
      "PTC_myCAF_median" = c(median(PTC_myCAF_distance$min)),
      "ATC_myCAF_median" = c(median(ATC_myCAF_distance$min))
  )
}

distances_overview <- rbind(mean_median_min_distances[[1]],
                            mean_median_min_distances[[2]],
                            mean_median_min_distances[[3]],
                            mean_median_min_distances[[4]],
                            mean_median_min_distances[[5]])

distances_overview_plotting <- distances_overview %>% 
  pivot_longer(
    cols = starts_with("PTC_myCAF_mean") | starts_with("PTC_myCAF_median") | starts_with("ATC_myCAF_mean") | starts_with("ATC_myCAF_median"),
    names_to = c("Distance_Group"),
    values_to = "Distance"
  )

median_distances_plotting <- distances_overview_plotting %>% subset(Distance_Group == "PTC_myCAF_median" | Distance_Group == "ATC_myCAF_median")

plot <- ggplot(median_distances_plotting, aes(Distance_Group, Distance)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Distance_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             #position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  geom_line(aes(group = sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("#0075DC", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Distance") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_median", "ATC_myCAF_median")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 105)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_ATC_myCAF/24-1202_PTC_ATC_myCAF_Median_Min_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

######### MEAN plots ################
mean_distances_plotting <- distances_overview_plotting %>% subset(Distance_Group == "PTC_myCAF_mean" | Distance_Group == "ATC_myCAF_mean")

plot <- ggplot(mean_distances_plotting, aes(Distance_Group, Distance)) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Distance_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             #position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  geom_line(aes(group = sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("#0075DC", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Distance") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_mean", "ATC_myCAF_mean")) +
  scale_y_continuous(breaks = c(0, 30, 60, 90, 120),
                     limits = c(0, 130)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_ATC_myCAF/24-1202_PTC_ATC_myCAF_Mean_Min_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

# Log2 of above
plot <- ggplot(mean_distances_plotting, aes(Distance_Group, log2(Distance+1))) +
  geom_boxplot(outlier.size = -1, 
               aes(fill = Distance_Group),
               alpha = 0.9, 
               show.legend = FALSE) + 
  geom_point(aes(),
             #position = position_jitter(width = 0.1, height = 0),
             size = 2, 
             alpha = 0.7,
             show.legend = FALSE) +
  geom_line(aes(group = sample), color = "black", linewidth = 0.5, alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("#0075DC", "#FFCC99")) +
  labs (x = "Tumor", y = "myCAF Distance") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_mean", "ATC_myCAF_mean")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6),
                     limits = c(0, 7.1)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_ATC_myCAF/24-1202_PTC_ATC_myCAF_Log2_Mean_Min_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

distance_wide <- pivot_wider(mean_distances_plotting, names_from = Distance_Group, values_from = Distance)
t.test(log2(distance_wide$PTC_myCAF_mean+1), log2(distance_wide$ATC_myCAF_mean+1), paired = TRUE) # two-sided t-test: 0.006445
