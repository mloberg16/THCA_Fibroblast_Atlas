# Author: Matthew Aaron Loberg
# Date: November 24, 2024
# Script: "24-1124_distance_comparison_function.R"

# Goal: 
# Write a function for comparing distances between deconvolved populations (e.g., RCTD, STDeconvolve, other)

# More info: 
# adapted from Andres_Distance_Calcs.Rmd

distance_comparison <- function(distances, proportions, topic1, topic2, t1_thresh, t2_thresh){
  
  # Filter the desired topics for barcode proportions greater than threshold and only keep barcode and topic column
  filtered_props_1 <- proportions[proportions[,topic1] > t1_thresh, ]
  rownames_holder <- rownames(filtered_props_1)
  filtered_props_1 <- data.frame(filtered_props_1[, c(topic1)])
  rownames(filtered_props_1) <- rownames_holder
  rm(rownames_holder)
  colnames(filtered_props_1) <- c(topic1)
  print(filtered_props_1)
  
  # Filter the desired topics for barcode proportions greater than threshold and only keep barcode and topic column
  filtered_props_2 <- proportions[proportions[,topic2] > t2_thresh, ]
  rownames_holder <- rownames(filtered_props_2)
  filtered_props_2 <- data.frame(filtered_props_2[, c(topic2)])
  rownames(filtered_props_2) <- rownames_holder
  rm(rownames_holder)
  colnames(filtered_props_2) <- c(topic2)
  print(filtered_props_2)
  
  # Create a new matrix that contains the first topic as rows and the second topic as columns
  dist_filtered <- distances
  dist_filtered <- dist_filtered[rownames(dist_filtered) %in% rownames(filtered_props_1), ]
  
  # Filter columns
  dist_filtered <- dist_filtered[, colnames(dist_filtered) %in% rownames(filtered_props_2)]
  
  print(dist_filtered)
  
  # Apply min function to each row (1=rows, 2=columns)
  min_df <- apply(dist_filtered, 1, function(x) min(x)) 
  
  # Convert to dataframe
  min_df <- as.data.frame(min_df)
  colnames(min_df) <- c("min")
  
  comp_str <- paste0(topic1, " vs. ", topic2)
  min_df$barcodes <- comp_str
  print(min_df)
  
  return(min_df)
}
