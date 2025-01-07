# 23-1127_Pu_Create_SO.R
# Author: Matthew Aaron Loberg
# Take a readdir from the Pu et al. data set and return a Seurat Object

# 23-1127 Update
# Changing nchar(readdirs from -6 to 61)

# Define the Pu_Create_SO function
Pu_Create_SO <- function(readdir, Cell_IDs){
  # Read in data and format for creating a Seurat object
  counts <- Seurat::Read10X(readdir)
  colnames(counts) <- paste(paste0("Pu_", substr(readdir, 61, nchar(readdir))), colnames(counts), sep = '_')
  SO <- CreateSeuratObject(counts = counts,
                           project = "Pu_etal_2021",
                           min.cells = 5,
                           min.features = 200,
                           orig.ident = paste("Pu", substr(readdir, nchar(readdir)-6, nchar(readdir)), sep = '_'))
  return(SO)
}
