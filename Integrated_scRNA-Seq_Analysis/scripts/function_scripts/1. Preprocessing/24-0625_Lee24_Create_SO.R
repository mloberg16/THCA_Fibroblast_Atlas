### Author: Matthew Aaron Loberg
### Date: 24-0625
### Script: 24-0625_Lee24_Create_SO.R

# Goal: Read in Lee et al. count data and create SO

# Update info:
# 24-0625 Update of 24-0506_Han_Create_SO.R for Lee et al. 2024 dataset

Lee_Create_SO <- function(readdir, Cell_IDs){

  # Read in data and format for creating a Seurat object
  counts <- as.data.frame(readr::read_table(readdir))
  rownames(counts) <- counts$gene_name
  counts <- counts[,2:ncol(counts)]
  SO <- CreateSeuratObject(counts = counts,
                           project = "Lee_etal_2024",
                           min.cells = 5,
                           min.features = 200,
                           orig.ident = paste("Lee24", substr(readdir, 52, nchar(readdir)-13), sep = '_'))

  # Set paper meta data for each sample to "Lee24"
  SO$Paper <- "Lee24"

  # Label Histology meta data as either "PTC" or "ATC" based on start of orig.ident label
  if(grepl("PT", SO$orig.ident[1], fixed = TRUE)){
    SO$Histology <- "PTC"
  }
  else if(grepl("AT", SO$orig.ident[1], fixed = TRUE)){
    SO$Histology <- "ATC"
  }

  # Return seurat object (SO)
  return(SO)

}
