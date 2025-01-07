# Author: Matthew Aron Loberg
# Date: August 16, 2024
# Script: 24-0816_Seurat_Basic_QC.R

# 24-0625 Update
# Note: this is an update of the "23-1009_Seurat_Basic_QC.R" script; the only difference is the use of dir.create(outputdir) at the start of the script

# 24-0814 Update
# This is an update of 24-0625_Seurat_Basic_QC.R
# I am changing the nFeature_RNA subset to be > 500

# 24-0816 Update
# This is an update of 24-0814_Seurat_Basic_QC.R
# I am cahnging the nFeature_RNA subset to be >= 500 instead of > 500

Seurat_Basic_QC <- function(SO, outputdir){
  # Create the outputdir
  dir.create(outputdir)

  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  # Add QC Stats for percent mitochondrial reads
  SO[["percent.mt"]] <- Seurat::PercentageFeatureSet(SO, pattern = "^MT-")

  # Visualize QC metrics as a violin plot
  QC_Violin <- Seurat::VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(outputdir, "QC_Violin.png"),
         QC_Violin, width = 5, height = 5, dpi = 600)

  # Visiualize QC feature correlations
  plot1 <- Seurat::FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- Seurat::FeatureScatter(SO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(file.path(outputdir, "QC_Correlations.png"),
         plot1 + plot2, width = 10, height = 5, dpi = 600)

  # Subset down to cells that meet QC cutoffs
  # Changing here to >= 500
  SO <- subset(SO, subset = nFeature_RNA > 200 & nCount_RNA >= 500 & percent.mt < 10)
  # Note that nFeature_RNA < 2500 is also often done to reduce doublets
  # This is maybe less accepted as good practice overall
  # Alternate options exist (e.g., doublet finder)
  # Here I will just do the cutoff of 200 features and 10% mitochondrial; I have later added the nCount_RNA >= 500 to get rid of some low quality cells.

  # Getting rid of cells with high mitochondrial percents generates features that are present in less than 5 of the remaining cells
  # This is problematic for running SCTransform and having mismatching features counts between SCTransform and RNA assays
  # I need to subset down to features that are in 5 cells
  # I cannot find a way to do this other than creating a new Seurat object ...
  SO_Cleaned <- CreateSeuratObject(counts = SO@assays$RNA@counts,
                                   min.cells = 5,
                                   min.features = 200,
                                   meta.data = SO@meta.data)

  ### Visualize QC again after cleaning
  # Visualize QC metrics as a violin plot
  QC_Violin <- Seurat::VlnPlot(SO_Cleaned, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(outputdir, "QC_Violin_PostCleaning.png"),
         QC_Violin, width = 5, height = 5, dpi = 600)

  # Visiualize QC feature correlations
  plot1 <- Seurat::FeatureScatter(SO_Cleaned, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- Seurat::FeatureScatter(SO_Cleaned, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(file.path(outputdir, "QC_Correlations_PostCleaning.png"),
         plot1 + plot2, width = 10, height = 5, dpi = 600)

  return(SO_Cleaned)
}
