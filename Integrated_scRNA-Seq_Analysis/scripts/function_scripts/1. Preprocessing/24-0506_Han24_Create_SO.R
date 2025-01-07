### Author: Matthew Aaron Loberg
### Date: 24-0506
### Script: 24-0506_Han24_Create_SO.R

# Goal: Read in Han et al. count data and create SO

# Update info:
# 24-0506 Update of 23-1127_Pu_Create_SO.R for Han et al. 2024 dataset

# Function definition
Han_Create_SO <- function(readdir, Cell_IDs){

  # Read in data and format for creating a Seurat object
  counts <- Seurat::Read10X(readdir)
  colnames(counts) <- paste(paste0("Han24_", substr(readdir, 46, nchar(readdir))), colnames(counts), sep = '_')
  SO <- CreateSeuratObject(counts = counts,
                           project = "Han_etal_2024",
                           min.cells = 5,
                           min.features = 200,
                           orig.ident = paste("Han24", substr(readdir, nchar(readdir)-6, nchar(readdir)), sep = '_'))

  # Return seurat object
  return(SO)
  
}
