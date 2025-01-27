### Author: Matthew Aaron Loberg
### Date: 24-0201
### Script: 24-0201_Wang_Create_SO.R

# Goal: Read in Wang et al. count data and create SO

# 24-0201 Wang Create Seurat Object function
Wang_Create_SO <- function(readdir){

  # Read in data and format for creating a Seurat object
  counts <- Seurat::Read10X_h5(readdir)
  colnames(counts) <- paste(paste0("Wang_", substr(readdir, 100, nchar(readdir)-3)), colnames(counts), sep = '_')
  SO <- CreateSeuratObject(counts = counts,
                           project = "Wang_etal_2022",
                           min.cells = 5,
                           min.features = 200,
                           orig.ident = paste("Wang", substr(readdir, 100, nchar(readdir)-3), sep = '_'))

  # Return Seruat Object (SO)
  return(SO)

}
