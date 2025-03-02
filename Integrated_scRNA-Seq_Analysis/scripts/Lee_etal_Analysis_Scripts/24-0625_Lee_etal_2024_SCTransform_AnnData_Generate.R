# Author: Matthew Aaron Loberg
# Date: June 25, 2024
# Script: 24-0625_Lee_etal_2024_SCTransform_AnnData_Generate

#### INFORMATION on SCRIPT ####
# Creating SCTransformed Seurat Objects and AnnData objects for Lee et al. 2024
# Paper link: https://www.nature.com/articles/s41467-024-45366-0
# GEO link for data download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232237
# I will use these to merge with other samples + run FastMNN (eventually after further exploration)

#### HISTORY OF THIS SCRIPT STARTING WITH Pu et al. -> Han et al. -> Lee et al. ####
# 23-1108 Update
# NOT RUNNING SoupX
# Seurat version 4 is REQUIRED for use of SeuratDisk

# 23-1121 Update
# Running ALL Pu Samples this time ... including lymph nodes, mets, etc.
# Note, script finished + run on 23-1122, so output files are from 23-1122

# 23-1127 Update
# I am re-running this with a new Pu Create SO function that will name Cell IDs without overlapping identities

# 24-0506 Update
# I am adapting the 23-1127 Pu et al. script to the new Han et al. 2024 ATC (BRAF mutant) samples

# 24-0624 Update
# I am adapting the 24-0506 version for Han et al. 2024 to use for Lee et al. 2024

### Load packages
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SoupX)
library(celldex) # Provides access to several reference datasets ("Pokedex for Cell Types) - celldex source: https://bioconductor.org/packages/3.17/data/experiment/vignettes/celldex/inst/doc/userguide.html
library(SingleR) # Load in SingleR - make sure to install the Bioconductor version as there is an old version on Github
# Packages required for scDblFinder
library(SingleCellExperiment)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
# scDblFinder
library(scDblFinder)

# Run source function scripts
source("function_scripts/24-0625_Lee24_Create_SO.R") # New Create SO function that fixes the error of overlapping names
source("function_scripts/24-0625_Seurat_Basic_QC.R")
source("function_scripts/24-0625_DoubletFinder_Function.R")
source("function_scripts/23-1026_SingleR_Prediction_Function.R")
#source("function_scripts/23-1031_AmbientRNA_Function.R") # No longer running SoupX, so will NOT run this (for now)

# Read Directories in as a list
# Read Directories
Lee_readdirs <- list("data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324266_PT10.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324267_PT12.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324268_PT3.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324269_PT5.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324270_PT7.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324271_PT8.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324272_PT9.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324273_AT9.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324274_AT13.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324275_AT16.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324276_AT17.count.tsv.gz",
                     "data_in_use/Lee_etal_2024_scRNA/RAW_DATA/GSM7324277_AT20.count.tsv.gz")

# Create Seurat Object in a list for each directory
Lee_SOs <- list()
for(i in 1:length(Lee_readdirs)){
  Lee_SOs[[i]] <- Lee_Create_SO(readdir = Lee_readdirs[[i]])
}
# Will get the following warning: "Feature names cannot have underscores ('_'), replacing with dashes ('-')
# This warning is NOT a problem
# Note: I did NOT get this warning this time

# Cleaning up
rm(Lee_readdirs)

# Explore Lee_SOs to make sure it is what I want - confirmed, commenting out for now
#Lee_SOs[[1]]
#Lee_SOs[[2]]$orig.ident
#Lee_SOs[[3]]$orig.ident
#Lee_SOs[[4]]$orig.ident
#Lee_SOs[[4]]$Paper
#Lee_SOs[[4]]$Histology
# Seems that orig.ident, paper, and histology were successfully defined

#### Doublet Detection ####
dir.create("outputs/Lee_etal_2024_Analysis_Outputs")
outputdir <- "outputs/Lee_etal_2024_Analysis_Outputs/24-0625_scDblFinder/"
dir.create(outputdir)
# Run doublet detection
for(i in 1:length(Lee_SOs)){
  Lee_SOs[[i]] <- Doublet_Detection(Lee_SOs[[i]], outputdir = paste0(outputdir, Lee_SOs[[i]]$orig.ident[1]))
}
rm(outputdir, i)

# Statistics for doublets + subsetting
Lee_Doublet_Info <- list()
for(i in 1:length(Lee_SOs)){
  Lee_Doublet_Info[[i]] <- table(Lee_SOs[[i]]$scDblFinder.class)
  Lee_SOs[[i]] <- Lee_SOs[[i]] %>% subset(scDblFinder.class == "singlet")
  Lee_SOs[[i]]$scDblFinder.class <- NULL
}
Lee_Doublet_Info
saveRDS(Lee_Doublet_Info, file = "data_in_use/Lee_etal_2024_scRNA/24-0625_Lee24_Doublet_Info.rds")
rm(Lee_Doublet_Info, i)

#### COMMENTING OUT - NOT RUNNING SoupX ####
#### SoupX Ambient RNA Detection ####
#outputdir <- "outputs/Pu_etal_Analysis_Outputs/23-1031_SoupX/"
#for(i in 1:length(Pu_PTCs)){
#  Pu_PTCs[[i]] <- AmbientRNA_Processing(Pu_PTCs[[i]], outputdir = paste0(outputdir, Pu_PTCs[[i]]$orig.ident[1]))
#}
#rm(outputdir, i)
# A note on the Rho values:
# PTC01_T: 0.04;
# PTC_02_T: 0.20;
# PTC_03_T: 0.09;
# PTC_05_T: 0.10;
# PTC_08_T: 0.01;
# PTC_09_T: 0.01;
# PTC-10_T: 0.02;
# The following warning: In sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w], :  'giveCsparse' is deprecated; setting repr="T" for you
# Should NOT need to worry about the warning

# Can check and see corrected counts
# Pu_PTCs[[1]]@assays$RNA@counts # Commenting out for now - do not want this to run every time

#### Basic QC ####
# QC plots + mitochondrial subsetting
QC_OutputDir <- "outputs/Lee_etal_2024_Analysis_Outputs/24-0625_All_Sample_QC/"
dir.create(QC_OutputDir)
for(i in 1:length(Lee_SOs)){
  Lee_SOs[[i]] <- Seurat_Basic_QC(Lee_SOs[[i]], outputdir = paste0(QC_OutputDir,Lee_SOs[[i]]$orig.ident[1]))
}
rm(QC_OutputDir, i)

#### Normalization with SCTransform ####
# Run SCTransform for all Lee SOs using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Lee_SOs <- lapply(X = Lee_SOs,
                  FUN = SCTransform,
                  variable.features.n = 7000,
                  return.only.var.genes = FALSE,
                  vst.flavor = "v2")

### Run SingleR on the Lee SOs list
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Lee_etal_2024_scRNA/24-0625_SingleR_Predictions/" # Set savedir for the tables with prediction info
dir.create(savedir)
outputdir <- "outputs/Lee_etal_2024_Analysis_Outputs/24-0625_SingleR_UMAPs/" # set the outputdir for sample umaps
dir.create(outputdir)
Lee_SOs <- SingleR_Predictions(SO_List = Lee_SOs, outputdir = outputdir, savedir = savedir)

### Add in MetaData about samples from the paper
# There is virtually NO meta data on these samples outside of histology (e.g., PTC or ATC)

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Lee_SOs)){
  print(length(rownames(GetAssayData(Lee_SOs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Lee_SOs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Lee_SOs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Lee_SOs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Lee_SOs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad
# I am going to save to the PC to save space in my OneDrive. I will format this the same way that the 2023_Integrated_scRNA-Seq_Aanalysis project is formatted

# Save SCTransformed data for each of the Han 2024 ATCs
savedir <- "data_in_use/Lee_etal_2024_scRNA/Processed_Data/Individual_Samples/"
dir.create(savedir)
for(i in 1:length(Lee_SOs)){
  savedir_temp <- paste0(savedir, Lee_SOs[[i]]$orig.ident[1])
  dir.create(savedir_temp)
  SaveH5Seurat(Lee_SOs[[i]], filename = paste0(savedir_temp, "/24-0625_SCTransformed.h5Seurat"), overwrite = TRUE)
  #Convert(paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE) # To save space, I will NOT be saving the AnnData object here ... can always come back and run in future
}
rm(list = ls())


# Session info
sessionInfo()
# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
#
# Matrix products: default
#
#
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8
#
# time zone: America/Chicago
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] scDblFinder_1.14.0          SingleCellExperiment_1.22.0 SingleR_2.2.0               celldex_1.10.1
# [5] SummarizedExperiment_1.30.2 Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4
# [9] IRanges_2.34.1              S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3
# [13] matrixStats_1.1.0           SoupX_1.6.2                 SeuratDisk_0.0.0.9020       SeuratObject_4.1.4
# [17] Seurat_4.4.0                lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1
# [21] dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0
# [25] tibble_3.2.1                ggplot2_3.5.0               tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21              splines_4.3.1                 later_1.3.1                   BiocIO_1.10.0
# [5] bitops_1.0-7                  filelock_1.0.2                polyclip_1.10-6               XML_3.99-0.15
# [9] lifecycle_1.0.4               edgeR_3.42.4                  globals_0.16.2                lattice_0.21-8
# [13] hdf5r_1.3.8                   MASS_7.3-60                   magrittr_2.0.3                limma_3.56.2
# [17] plotly_4.10.3                 yaml_2.3.7                    metapod_1.8.0                 httpuv_1.6.12
# [21] sctransform_0.4.1             sp_2.1-1                      spatstat.sparse_3.0-3         reticulate_1.34.0
# [25] cowplot_1.1.1                 pbapply_1.7-2                 DBI_1.1.3                     RColorBrewer_1.1-3
# [29] abind_1.4-5                   zlibbioc_1.46.0               Rtsne_0.16                    RCurl_1.98-1.13
# [33] rappdirs_0.3.3                GenomeInfoDbData_1.2.10       ggrepel_0.9.4                 irlba_2.3.5.1
# [37] listenv_0.9.0                 spatstat.utils_3.0-4          goftest_1.2-3                 dqrng_0.3.1
# [41] spatstat.random_3.2-1         fitdistrplus_1.1-11           parallelly_1.36.0             DelayedMatrixStats_1.22.6
# [45] leiden_0.4.3.1                codetools_0.2-19              DelayedArray_0.26.7           scuttle_1.10.3
# [49] tidyselect_1.2.0              viridis_0.6.4                 ScaledMatrix_1.8.1            BiocFileCache_2.8.0
# [53] spatstat.explore_3.2-5        GenomicAlignments_1.36.0      jsonlite_1.8.7                BiocNeighbors_1.18.0
# [57] ellipsis_0.3.2                progressr_0.14.0              scater_1.28.0                 ggridges_0.5.4
# [61] survival_3.5-5                tools_4.3.1                   ica_1.0-3                     Rcpp_1.0.11
# [65] glue_1.6.2                    gridExtra_2.3                 withr_2.5.2                   BiocManager_1.30.22
# [69] fastmap_1.1.1                 bluster_1.10.0                fansi_1.0.5                   rsvd_1.0.5
# [73] digest_0.6.33                 timechange_0.2.0              R6_2.5.1                      mime_0.12
# [77] colorspace_2.1-0              scattermore_1.2               tensor_1.5                    spatstat.data_3.0-3
# [81] RSQLite_2.3.3                 utf8_1.2.4                    generics_0.1.3                data.table_1.14.8
# [85] rtracklayer_1.60.1            httr_1.4.7                    htmlwidgets_1.6.2             S4Arrays_1.0.6
# [89] uwot_0.1.16                   pkgconfig_2.0.3               gtable_0.3.4                  blob_1.2.4
# [93] lmtest_0.9-40                 XVector_0.40.0                htmltools_0.5.7               scales_1.3.0
# [97] png_0.1-8                     scran_1.28.2                  rstudioapi_0.15.0             rjson_0.2.21
# [101] tzdb_0.4.0                    reshape2_1.4.4                nlme_3.1-162                  curl_5.1.0
# [105] zoo_1.8-12                    cachem_1.0.8                  BiocVersion_3.17.1            KernSmooth_2.23-21
# [109] vipor_0.4.5                   parallel_4.3.1                miniUI_0.1.1.1                AnnotationDbi_1.62.2
# [113] restfulr_0.0.15               pillar_1.9.0                  grid_4.3.1                    vctrs_0.6.4
# [117] RANN_2.6.1                    promises_1.2.1                BiocSingular_1.16.0           beachmat_2.16.0
# [121] dbplyr_2.3.4                  xtable_1.8-4                  cluster_2.1.4                 beeswarm_0.4.0
# [125] locfit_1.5-9.8                Rsamtools_2.16.0              cli_3.6.1                     compiler_4.3.1
# [129] rlang_1.1.2                   crayon_1.5.2                  future.apply_1.11.0           ggbeeswarm_0.7.2
# [133] plyr_1.8.9                    stringi_1.8.1                 BiocParallel_1.34.2           viridisLite_0.4.2
# [137] deldir_1.0-9                  munsell_0.5.0                 Biostrings_2.68.1             lazyeval_0.2.2
# [141] spatstat.geom_3.2-7           Matrix_1.6-1                  ExperimentHub_2.8.1           hms_1.1.3
# [145] patchwork_1.2.0               sparseMatrixStats_1.12.2      bit64_4.0.5                   future_1.33.0
# [149] statmod_1.5.0                 KEGGREST_1.40.1               shiny_1.8.0                   interactiveDisplayBase_1.38.0
# [153] AnnotationHub_3.8.0           ROCR_1.0-11                   igraph_2.0.3                  memoise_2.0.1
# [157] xgboost_1.7.5.1               bit_4.0.5
