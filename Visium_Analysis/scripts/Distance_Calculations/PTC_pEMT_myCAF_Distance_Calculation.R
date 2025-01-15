### Author: Matthew Aaron Loberg
### Date: November 24, 2024
### Script: PTC_pEMT_myCAF_Distance_Calculations.R
### Source Script Name: 24-1124_PTC_pEMT_myCAF_Distance_Calculations.R

# Goal:
# PTC Samples pEMT PTC vs myCAF distance/PTC vs myCAF distance
# Looking to see if there is a difference in the distance between the hypothesized leading-edge pEMT PTC population and myCAFs compared to the broader PTC population and myCAFs
# Hypothesis is that pEMT PTC is closer to myCAFs than PTC
# Excluding samples with mixed ATC in this analysis--only samples with pure PTC component

# Load packages
library(Seurat)
library(tidyverse)

# Load source script for w/ function for distance comparison: distance_comparison()
source("scripts/Function_Scripts/24-1124_distance_comparison_function.R")

SO_readdirs <- c(#"Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy1.RDS",    # PTC/ATC mix - excluded
                 #"Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy5.RDS",    # PTC/ATC mix - excluded
                 #"Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy6.RDS",    # PTC/ATC mix - excluded
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy7.RDS",     # pure PTC
                 #"Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy10.RDS",   # PTC/ATC mix - excluded
                 #"Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy11.RDS",   # PTC/ATC mix - excluded
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy15.RDS",    # Pure PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy16_V2.RDS", # Pure PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Thy17.RDS",    # Pure PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds01.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds02.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds03.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds04.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds05.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds06.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds07.RDS",   # Peds PTC
                 "Data_in_Use/Pre_Processed_with_RCTD_Annotations/24-1104_iCAF2_Excluded/Peds08.RDS")   # Peds PTC

samples <- c("Thy7", "Thy15", "Thy16", "Thy17", 
             "Peds01", "Peds02", "Peds03", "Peds04", 
             "Peds05", "Peds06", "Peds07", "Peds08")

mean_median_min_distances <- list()
for(i in 1:length(SO_readdirs)){
  
  # Load in spatial seurat object (SO) from readdirs
  spatial_SO <- readRDS(SO_readdirs[i])
  
  # Load in distance matrix for sample
  distance_matrix <- readRDS(file = paste0("Data_in_Use/Distance_Matrices/", samples[i], "_Distance_Matrix.RDS"))

  # Create proportions data frame for thresholding barcodes for distance comparisons
  proportions <- data.frame(
      "pEMT" = c(spatial_SO$pEMT_PTC_RCTD),
      "PTC" = c(spatial_SO$PTC_RCTD),
      "myCAF" = c(spatial_SO$myCAF_RCTD)
  )

  # Calculate min. distance of spots with at least 10% pEMT PTC phenotype from spots with at least 10% myCAF phenotype
  pEMT_myCAF_distance <- distance_comparison(distances = distance_matrix,
                                             proportions = proportions,
                                             topic1 = "pEMT",
                                             topic2 = "myCAF",
                                             t1_thresh = 0.1,
                                             t2_thresh = 0.1)

  # Calculate min. distance of spots with at least 10% PTC phenotype from spots with at least 10% myCAF phenotype
  PTC_myCAF_distance <- distance_comparison(distances = distance_matrix,
                                            proportions = proportions,
                                            topic1 = "PTC",
                                            topic2 = "myCAF",
                                            t1_thresh = 0.1,
                                            t2_thresh = 0.1)
  
  combined_distance <- rbind(pEMT_myCAF_distance, PTC_myCAF_distance)
  
  # Violin plot to compare distances 
  distance_violin <- ggplot(combined_distance, aes(x = barcodes, y = min, fill = barcodes)) +
    geom_violin(scale = "width") +
    theme_classic() +
    scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
    scale_x_discrete(name ="Deconvolution", limits = c("PTC vs. myCAF", "pEMT vs. myCAF"))
    
    
  # Save violin plot
  ggsave(file = paste0("outputs/Distance_Plots/", samples[i], "/PTC_pEMT_myCAF_Min_Distance_Violin.png"),
         distance_violin, height = 5, width = 6, dpi = 600)
  
  # For each sample calculate the mean and the median of the minimum distance from myCAF for pEMT PTC and PTC
  mean_median_min_distances[[i]] <- data.frame(
      "sample" = c(samples[i]),
      "pEMT_myCAF_mean" = c(mean(pEMT_myCAF_distance$min)),
      "PTC_myCAF_mean" = c(mean(PTC_myCAF_distance$min)),
      "pEMT_myCAF_median" = c(median(pEMT_myCAF_distance$min)),
      "PTC_myCAF_median" = c(median(PTC_myCAF_distance$min))
  )
}

# Combine the distances from all of the samples 
distances_overview <- rbind(mean_median_min_distances[[1]],
                            mean_median_min_distances[[2]],
                            mean_median_min_distances[[3]],
                            mean_median_min_distances[[4]],
                            mean_median_min_distances[[5]],
                            mean_median_min_distances[[6]],
                            mean_median_min_distances[[7]],
                            mean_median_min_distances[[8]],
                            mean_median_min_distances[[9]],
                            mean_median_min_distances[[10]],
                            mean_median_min_distances[[11]],
                            mean_median_min_distances[[12]])

# Expand data frame for plotting
distances_overview_plotting <- distances_overview %>% 
  pivot_longer(
    cols = starts_with("pEMT_myCAF_mean") | starts_with("pEMT_myCAF_median") | starts_with("PTC_myCAF_mean") | starts_with("PTC_myCAF_median"),
    names_to = c("Distance_Group"),
    values_to = "Distance"
  )

# plots based on median
median_distances_plotting <- distances_overview_plotting %>% subset(Distance_Group == "pEMT_myCAF_median" | Distance_Group == "PTC_myCAF_median")

# make plot 
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
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "iPVCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_median", "pEMT_myCAF_median")) +
  scale_y_continuous(breaks = c(0, 300, 600, 900, 1200, 1500),
                     limits = c(0, 1550)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_pEMT_myCAF/24-1124_PTC_pEMT_myCAF_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

# Plot with distance as log2 with a pseudocount to avoid log2(0)
plot <- ggplot(median_distances_plotting, aes(Distance_Group, log2(Distance+1))) +
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
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "iPVCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_median", "pEMT_myCAF_median")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 11)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_pEMT_myCAF/24-1124_PTC_pEMT_myCAF_Log2_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)


######### MEAN plots ################
# Repeat of median plots above but using mean instead of median
# mean does a better job of capturing the variability 
mean_distances_plotting <- distances_overview_plotting %>% subset(Distance_Group == "pEMT_myCAF_mean" | Distance_Group == "PTC_myCAF_mean")

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
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "iPVCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_mean", "pEMT_myCAF_mean")) +
  scale_y_continuous(breaks = c(0, 400, 800, 1200, 1600),
                     limits = c(0, 1700)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_pEMT_myCAF/24-1124_PTC_pEMT_myCAF_Mean_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)

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
  scale_fill_manual(values = c("darkgrey", "#FFCC99")) +
  labs (x = "Tumor", y = "iPVCAF Correlation") + 
  theme_classic() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 14.5), 
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 16)) +
  scale_x_discrete(name ="Deconvolution", limits = c("PTC_myCAF_mean", "pEMT_myCAF_mean")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 11)) 
# Not, max and min values need to be restricted for each new gene that is plotted
ggsave("outputs/Distance_Plots/PTC_pEMT_myCAF/24-1124_PTC_pEMT_myCAF_Log2_Mean_Distance_Boxplots.png", 
       width = 2.75, height = 3, 
       plot, dpi = 600, create.dir = TRUE)



# stats
distance_wide <- pivot_wider(mean_distances_plotting, names_from = Distance_Group, values_from = Distance)
t.test(distance_wide$pEMT_myCAF_mean, distance_wide$PTC_myCAF_mean, paired = TRUE) 
