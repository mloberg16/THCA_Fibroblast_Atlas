# Author: Matthew Aaron Loberg
# Date: November 23, 2024
# Script: Visium_Distance_Matrix_Calculations.R
# Source Script Name: 24-1123_Visium_Distance_Matrix_Calculations.R

# Credit to Andres Ocampo for help in developing this distance matrix code 

### Goal: 
# For each visium object, I want a distance matrix (how far each barcoded area is from every other barcoded area)
# This will be calculated as the Euclidean distance using pythagorean theorem 
# This distance matrix can be used to subset the distances between spots of interest

### More iinfo
# For reference, see "Andres_Distance_Calcs.RMD"
# This is an idea that we had during Andres's lab rotation

# Load Packages
library(Seurat)

# readdirs
readdirs <- list(
  "Data_in_Use/2021_JHU_Data/Thy1_Processed/22-0913_Thy1_Raw_PreProcessed.rds",
  "Data_in_Use/2021_JHU_Data/Thy2_Processed/22-1129_Thy2_Raw_PreProcessed.rds",
  "Data_in_Use/2021_JHU_Data/Thy3_Processed/22-1129_Thy3_Raw_PreProcessed.rds",
  "Data_in_Use/2021_JHU_Data/Thy4_Processed/22-1129_Thy4_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy5_Processed/22-0825_Thy5_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy6_Processed/22-0915_Thy6_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy7_Processed/22-1122_Thy7_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy8_Processed/22-1122_Thy8_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy9_ManualAlign_Processed/22-1003_Thy9_ManualAlign_Raw_PreProcessed.rds", # This is manual align from Lana Olson
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy10_Processed/22-0915_Thy10_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy11_Processed/22-0905_Thy11_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy12_Processed/22-0915_Thy12_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy13_Processed/22-0915_Thy13_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy14_Processed/22-1129_Thy14_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy15_Processed/22-1129_Thy15_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy16_Processed/24-1121_Thy16_Raw_PreProcessed.rds", # This is manual align from Lana Olson
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy17_Processed/22-1129_Thy17_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy18_Processed/22-1129_Thy18_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy19_Processed/22-1129_Thy19_Raw_PreProcessed.rds",
  "Data_in_Use/August_2022_VANTAGE_Visium_Run/Processed_Outputs/Thy20_Processed/22-1129_Thy20_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds01v2/24-1104_Peds01_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds02v2/24-1104_Peds02_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds03v2/24-1104_Peds03_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds04v2/24-1104_Peds04_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds05v2/24-1104_Peds05_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds06v2/24-1104_Peds06_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds07v2/24-1104_Peds07_Raw_PreProcessed.rds",
  "Data_in_Use/Belcher_Peds_Visium/Processed/Peds08v2/24-1104_Peds08_Raw_PreProcessed.rds"
)

samplenames <- c("Thy1", "Thy2", "Thy3", "Thy4", "Thy5", "Thy6",
                 "Thy7", "Thy8", "Thy9", "Thy10", "Thy11", "Thy12",
                 "Thy13", "Thy14", "Thy15", "Thy16", "Thy17", "Thy18", "Thy19", "Thy20", "Peds01", "Peds02",
                 "Peds03", "Peds04", "Peds05", "Peds06", "Peds07", "Peds08")

spot_diameter_fullres <- c(35.67556561084767,   # Thy1 spot_diameter_fullres
                           35.67948420101959,   # Thy2 spot_diameter_fullres
                           35.684497508174154,  # Thy3 spot_diameter_fullres
                           35.689101703753444,  # Thy4 spot_diameter_fullres
                           236.90587837361966,  # Thy5 spot_diameter_fullres
                           235.44108564074304,  # Thy6 spot_diameter_fullres
                           239.13277820640744,  # Thy7 spot_diameter_fullres
                           238.70588835658745,  # Thy8 spot_diameter_fullres
                           236.64687,           # Thy9 manual align spot_diameter_fullres
                           235.47597500610541,  # Thy10 spot_diameter_fullres
                           235.4675536841221,   # Thy11 spot_diameter_fullres
                           241.8773846930639,   # Thy12 spot_diameter_fullres
                           235.4538460374921,   # Thy13 spot_diameter_fullres
                           235.50181073012482,  # Thy14 spot_diameter_fullres
                           235.46194252361263,  # Thy15 spot_diameter_fullres
                           235.21776,           # Thy16 manual align spot_diameter_fullres
                           235.47463457815473,  # Thy17 spot_diameter_fullres
                           235.4781329036873,   # Thy18 spot_diameter_fullres
                           235.4551843553312,   # Thy19 spot_diameter_fullres
                           235.48634643598638,  # Thy20 spot_diameter_fullres
                           235.46887353654128,  # Peds01 v2 spot_diameter_fullres
                           235.6411295425445,   # Peds02 v2 spot_diameter_fullres
                           235.53693237192695,  # Peds03 v3 spot_diameter_fullres
                           235.6634359013668,   # Peds04 v4 spot_diameter_fullres
                           235.57711517211794,  # Peds05 v5 spot_diameter_fullres
                           235.6559465769829,   # Peds06 v6 spot_diameter_fullres
                           235.53617784020585,  # Peds07 v7 spot_diameter_fullres
                           235.6667708600977)   # Peds08 v8 spot_diameter_fullres     


for(i in 1:length(readdirs)){
  
  # Read in spatial SO
  spatial_SO <- readRDS(file = readdirs[[i]])
  
  # Pull out position data
  # scale = NULL REQUIRED to get original pixels from full-resolution image
  # Need pixels in order to convert to microns; allows for exact distance measurement
  pos <- GetTissueCoordinates(spatial_SO, scale = NULL) # returns positions in pixels
  pos$imagerow <- pos$imagerow*(65/spot_diameter_fullres[i]) # converts pixels to microns
  pos$imagecol <- pos$imagecol*(65/spot_diameter_fullres[i]) # converts pixels to microns
  
  # Initialize distance matrix
  dist_mat <- matrix(NA, nrow = nrow(pos), ncol = nrow(pos))
  
  # Nested Loop through rows and columns of position (pos) df
  for(n in 1:nrow(pos)){
    for(j in 1:nrow(pos)){ 
      if(n == j){
        dist_mat[n,j] <- 0 
      } else {
        point1 <- as.list(pos[n,])
        point2 <- as.list(pos[j,])
        dist <- sqrt((point1$imagerow - point2$imagerow)^2 + (point1$imagecol - point2$imagecol)^2)
        dist_mat[n,j] <- dist
        dist_mat[j,n] <- dist
      }
    }
  }
  colnames(dist_mat) <- rownames(pos)
  rownames(dist_mat) <- rownames(pos)
  saveRDS(dist_mat, file = paste0("Data_in_Use/Distance_Matrices/", samplenames[i], "_Distance_Matrix.RDS"))
}

##### Cleaning up
rm(list = ls())

##### Session Info
sessionInfo()
