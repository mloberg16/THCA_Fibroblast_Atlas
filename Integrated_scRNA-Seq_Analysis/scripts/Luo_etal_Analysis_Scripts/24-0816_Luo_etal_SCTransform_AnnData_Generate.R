### Author: Matthew Aaron Loberg
### Date: August 16, 2024
### Script: 24-0816_Luo__etal_SCTransform_AnnData_Generate.R

# Goal:
# Creating SCTransformed Seurat Objects and AnnData objects for Luo et al.

# 23-1031 Update
# Incorporate Doublet Finder
# Incorporate SoupX
# SingleR of Myeloid, CD4, CD8, and HPCA (Fine + Main)

# 23-1108 Update
# NOT RUNNING SoupX

# 24-0814 Update
# Running for just the following: PTC_WJL; ATC_LJ; PTC_XHY
# These are the samples that have cells with < 500 nCount_RNA
# I want to re-run this with the new 24-0814_Seurat_Basic_QC.R and remove these cells from the analysis
# Since it is ALREADY 5/9 of the samples, I am just going to re-run it for ALL of them with new QC
# Also using 24-0625 Doublet Finder function (updated from 23-1025) to AVOID outputir problems
# ALso using 24-0814 SingleR Prediction Function to AVOID outputdir problems

# 24-0816 Update
# Running with new Basic QC function -> switch to nCount_RNA >= 500
# For simplicity/consistency, I just went ahead and ran this script for all of Luo et al. samples
# They are now all pre-processed identically through this script

# Load packages
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
source("function_scripts/24-0816_Seurat_Basic_QC.R")
source("function_scripts/24-0625_DoubletFinder_Function.R")
source("function_scripts/24-0814_SingleR_Prediction_Function.R")
#source("function_scripts/23-1031_AmbientRNA_Function.R")

# Read Directories in as a list
# Read Directories
# Note, change: I am going to read WJL1 and WJL2 SEPARATE (same for XHY1 and XHY2). However, I am going to set them as the same "orig.ident" and give them different "Identifier" variables for the purpose of downstream integration
# These read directories are to previously generated Seurat objects
# See the following pre-processing scripts for how Seurat objects were generated: 
# 22-1031_Luo_ReadData.R -> 22-1031_Luo_Seperate_Matrix.R -> 22-1101_Luo_Create_SOs.R
Luo_SO_readdirs <- list("data_in_use/Luo_etal_2021_ATC_scRNA/ATC_WYF_ATC1_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/ATC_MSQ_ATC2_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/ATC_LJ_ATC3_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/PTC_WJL1_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/PTC_WJL2_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/PTC_XHY1_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/PTC_XHY2_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/PTC_XTZ_SO.RDS",
                        "data_in_use/Luo_etal_2021_ATC_scRNA/NOM_XTZ_SO.RDS")

# Read in Seurat Object in a list for each directory
Luo_SOs <- list()
for(i in 1:length(Luo_SO_readdirs)){
  Luo_SOs[[i]] <- readRDS(Luo_SO_readdirs[[i]])
}

### Metadata Updates ###
# Look at each of these - note, they already have sample prepended to cell IDs, so that is ideal!
# Will make sure orig.ident is what I want
# Will also add a "Histology" column
# Will also add a "Paper" column
# Histology and paper columns are for future subsetting when integrated with all data later on
# The paper also has mutation, sex, and Hashimotos metadata for each sample; I will add that here as well. 

# Luo ATC WYF meta data
Luo_SOs[[1]]
Luo_SOs[[1]]$orig.ident
Luo_SOs[[1]]$orig.ident <- "Luo_ATC_WYF"
Luo_SOs[[1]]$Identifier <- "Luo_ATC_WYF"
Luo_SOs[[1]]$Histology <- "ATC"
Luo_SOs[[1]]$Paper <- "Luo"
Luo_SOs[[1]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[1]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[1]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[1]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[1]]$Age <- 78
Luo_SOs[[1]]$Sex <- "F"
Luo_SOs[[1]]$Concurrent_HT <- "YES"

# Luo ATC MSQ meta data
Luo_SOs[[2]]
Luo_SOs[[2]]$orig.ident
Luo_SOs[[2]]$orig.ident <- "Luo_ATC_MSQ"
Luo_SOs[[2]]$Identifier <- "Luo_ATC_MSQ"
Luo_SOs[[2]]$Histology <- "ATC"
Luo_SOs[[2]]$Paper <- "Luo"
Luo_SOs[[2]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[2]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[2]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[2]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[2]]$Age <- 78
Luo_SOs[[2]]$Sex <- "F"
Luo_SOs[[2]]$Concurrent_HT <- "YES"

# Luo ATC LJ meta data
Luo_SOs[[3]]
Luo_SOs[[3]]$orig.ident
Luo_SOs[[3]]$orig.ident <- "Luo_ATC_LJ"
Luo_SOs[[3]]$Identifier <- "Luo_ATC_LJ"
Luo_SOs[[3]]$Histology <- "ATC"
Luo_SOs[[3]]$Paper <- "Luo"
Luo_SOs[[3]]$BRAFV600E <- "YES_Luo"
Luo_SOs[[3]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[3]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[3]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[3]]$Age <- 69
Luo_SOs[[3]]$Sex <- "M"
Luo_SOs[[3]]$Concurrent_HT <- "NO"

# Luo PTC WJL sample #1 meta data
Luo_SOs[[4]]
Luo_SOs[[4]]$orig.ident
Luo_SOs[[4]]$orig.ident <- "Luo_PTC_WJL"
Luo_SOs[[4]]$Identifier <- "Luo_PTC_WJL1"
Luo_SOs[[4]]$Histology <- "PTC"
Luo_SOs[[4]]$Paper <- "Luo"
Luo_SOs[[4]]$BRAFV600E <- "YES_Luo"
Luo_SOs[[4]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[4]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[4]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[4]]$Age <- 56
Luo_SOs[[4]]$Sex <- "M"
Luo_SOs[[4]]$Concurrent_HT <- "NO"

# Luo PTC WJL sample #2 meta data
Luo_SOs[[5]]
Luo_SOs[[5]]$orig.ident
Luo_SOs[[5]]$orig.ident <- "Luo_PTC_WJL"
Luo_SOs[[5]]$Identifier <- "Luo_PTC_WJL2"
Luo_SOs[[5]]$Histology <- "PTC"
Luo_SOs[[5]]$Paper <- "Luo"
Luo_SOs[[5]]$BRAFV600E <- "YES_Luo"
Luo_SOs[[5]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[5]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[5]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[5]]$Age <- 56
Luo_SOs[[5]]$Sex <- "M"
Luo_SOs[[5]]$Concurrent_HT <- "NO"

# Luo PTC XHY sample #1 meta data
Luo_SOs[[6]]
Luo_SOs[[6]]$orig.ident
Luo_SOs[[6]]$orig.ident <- "Luo_PTC_XHY"
Luo_SOs[[6]]$Identifier <- "Luo_PTC_XHY1"
Luo_SOs[[6]]$Histology <- "PTC"
Luo_SOs[[6]]$Paper <- "Luo"
Luo_SOs[[6]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[6]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[6]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[6]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[6]]$Age <- 64
Luo_SOs[[6]]$Sex <- "F"
Luo_SOs[[6]]$Concurrent_HT <- "NO"

# Luo PTC XHY sample #2 meta data
Luo_SOs[[7]]
Luo_SOs[[7]]$orig.ident
Luo_SOs[[7]]$orig.ident <- "Luo_PTC_XHY"
Luo_SOs[[7]]$Identifier <- "Luo_PTC_XHY2"
Luo_SOs[[7]]$Histology <- "PTC"
Luo_SOs[[7]]$Paper <- "Luo"
Luo_SOs[[7]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[7]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[7]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[7]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[7]]$Age <- 64
Luo_SOs[[7]]$Sex <- "F"
Luo_SOs[[7]]$Concurrent_HT <- "NO"

# Luo PTC XTZ meta data
Luo_SOs[[8]]
Luo_SOs[[8]]$orig.ident
Luo_SOs[[8]]$orig.ident <- "Luo_PTC_XTZ"
Luo_SOs[[8]]$Identifier <- "Luo_PTC_XTZ"
Luo_SOs[[8]]$Histology <- "PTC"
Luo_SOs[[8]]$Paper <- "Luo"
Luo_SOs[[8]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[8]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[8]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[8]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[8]]$Age <- 66
Luo_SOs[[8]]$Sex <- "F"
Luo_SOs[[8]]$Concurrent_HT <- "YES"

# Luo Nom XTZ meta data
Luo_SOs[[9]]
Luo_SOs[[9]]$orig.ident
Luo_SOs[[9]]$orig.ident <- "Luo_NOM_XTZ"
Luo_SOs[[9]]$Identifier <- "Luo_NOM_XTZ"
Luo_SOs[[9]]$Histology <- "Normal"
Luo_SOs[[9]]$Paper <- "Luo"
Luo_SOs[[9]]$BRAFV600E <- "NO_Luo"
Luo_SOs[[9]]$TERT_MUT <- "NO_Luo"
Luo_SOs[[9]]$TP53_MUT <- "NO_Luo"
Luo_SOs[[9]]$RAS_MUT <- "NO_Luo"
Luo_SOs[[9]]$Age <- 66
Luo_SOs[[9]]$Sex <- "F"
Luo_SOs[[9]]$Concurrent_HT <- "YES"

#### Doublet Detection ####
outputdir <- "outputs/Luo_etal_Individual_Analysis/24-0816_scDblFinder/"
dir.create(outputdir)
# Run doublet detection
for(i in 1:length(Luo_SOs)){
  Luo_SOs[[i]] <- Doublet_Detection(Luo_SOs[[i]], outputdir = paste0(outputdir, Luo_SOs[[i]]$Identifier[1]))
}
rm(outputdir, i)

# Statistics for doublets + subsetting
Luo_Doublet_Info <- list()
for(i in 1:length(Luo_SOs)){
  Luo_Doublet_Info[[i]] <- table(Luo_SOs[[i]]$scDblFinder.class)
  Luo_SOs[[i]] <- Luo_SOs[[i]] %>% subset(scDblFinder.class == "singlet")
  Luo_SOs[[i]]$scDblFinder.class <- NULL
}
saveRDS(Luo_Doublet_Info, file = "data_in_use/Luo_etal_2021_ATC_scRNA/24-0816_Lu_Doublet_Info.rds")
rm(Luo_Doublet_Info, i)

#### Commented out - NOT RUNNING SoupX ####
#### SoupX Ambient RNA Detection ####
# outputdir <- "outputs/Luo_etal_Individual_Analysis/23-1031_SoupX/"
# for(i in 1:length(Luo_SOs)){
#   Luo_SOs[[i]] <- AmbientRNA_Processing(Luo_SOs[[i]], outputdir = paste0(outputdir, Luo_SOs[[i]]$Identifier[1]))
# }
# rm(outputdir, i)
# A note on the Rho values:
# [INSERT INFO]
# The following warning: In sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w], :  'giveCsparse' is deprecated; setting repr="T" for you
# Should NOT need to worry about the warning

# Can check and see corrected counts
# Luo_ATCs[[1]]@assays$RNA@counts # Commenting out for now - do not want this to run every time


# QC plots + mitochondrial subsetting
QC_OutputDir <- "outputs/Luo_etal_Individual_Analysis/24-0816_All_Sample_QC/"
dir.create(QC_OutputDir)
for(i in 1:length(Luo_SOs)){
  Luo_SOs[[i]] <- Seurat_Basic_QC(Luo_SOs[[i]], outputdir = paste0(QC_OutputDir,Luo_SOs[[i]]$Identifier[1]))
}

# Run SCTransform for all Luo_Sos using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Luo_SOs <- lapply(X = Luo_SOs,
                  FUN = SCTransform,
                  variable.features.n = 7000,
                  return.only.var.genes = FALSE,
                  vst.flavor = "v2")

#### SingleR Annotations ####
### Run SingleR on the Lu_SOs list
### But First ... I would like to RUN SingleR to do some deconvolution to add to this data!
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Luo_etal_2021_ATC_scRNA/24-0816_SingleR_Predictions/" # Set savedir for the tables with prediction info
dir.create(savedir)
outputdir <- "outputs/Luo_etal_Individual_Analysis/24-0816_SingleR_UMAPs/" # set the outputdir for sample umaps
dir.create(outputdir)
Luo_SOs <- SingleR_Predictions(SO_List = Luo_SOs, outputdir = outputdir, savedir = savedir)

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Luo_SOs)){
  print(length(rownames(GetAssayData(Luo_SOs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Luo_SOs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Luo_SOs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Luo_SOs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Luo_SOs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad

# Save SCTransformed data for each of the Luo_SOs
savedir <- "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Luo_etal_2021_ATC_scRNA/Processed_Data/Individual_Samples/"
for(i in 1:length(Luo_SOs)){
  savedir_temp <- paste0(savedir, Luo_SOs[[i]]$Identifier[1])
  SaveH5Seurat(Luo_SOs[[i]], filename = paste0(savedir_temp, "/24-0816_SCTransformed.h5Seurat"), overwrite = TRUE)
  #Convert(paste0(savedir_temp, "/24-0816_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE)
}

##### Session Info #####
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
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C
# [5] LC_TIME=English_United States.utf8
#
# time zone: America/Chicago
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] scDblFinder_1.14.0          SingleCellExperiment_1.22.0 SingleR_2.2.0
# [4] celldex_1.10.1              SummarizedExperiment_1.30.2 Biobase_2.60.0
# [7] GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1
# [10] S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3
# [13] matrixStats_1.1.0           SoupX_1.6.2                 SeuratDisk_0.0.0.9020
# [16] SeuratObject_4.1.4          Seurat_4.4.0                lubridate_1.9.3
# [19] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4
# [22] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0
# [25] tibble_3.2.1                ggplot2_3.5.0               tidyverse_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21              splines_4.3.1                 later_1.3.1
# [4] BiocIO_1.10.0                 bitops_1.0-7                  filelock_1.0.2
# [7] polyclip_1.10-6               XML_3.99-0.15                 lifecycle_1.0.4
# [10] edgeR_3.42.4                  globals_0.16.2                lattice_0.21-8
# [13] hdf5r_1.3.8                   MASS_7.3-60                   magrittr_2.0.3
# [16] limma_3.56.2                  plotly_4.10.3                 yaml_2.3.7
# [19] metapod_1.8.0                 httpuv_1.6.12                 sctransform_0.4.1
# [22] sp_2.1-1                      spatstat.sparse_3.0-3         reticulate_1.34.0
# [25] cowplot_1.1.1                 pbapply_1.7-2                 DBI_1.1.3
# [28] RColorBrewer_1.1-3            abind_1.4-5                   zlibbioc_1.46.0
# [31] Rtsne_0.16                    RCurl_1.98-1.13               rappdirs_0.3.3
# [34] GenomeInfoDbData_1.2.10       ggrepel_0.9.4                 irlba_2.3.5.1
# [37] listenv_0.9.0                 spatstat.utils_3.0-4          goftest_1.2-3
# [40] dqrng_0.3.1                   spatstat.random_3.2-1         fitdistrplus_1.1-11
# [43] parallelly_1.36.0             DelayedMatrixStats_1.22.6     leiden_0.4.3.1
# [46] codetools_0.2-19              DelayedArray_0.26.7           scuttle_1.10.3
# [49] tidyselect_1.2.0              viridis_0.6.4                 ScaledMatrix_1.8.1
# [52] BiocFileCache_2.8.0           spatstat.explore_3.2-5        GenomicAlignments_1.36.0
# [55] jsonlite_1.8.7                BiocNeighbors_1.18.0          ellipsis_0.3.2
# [58] progressr_0.14.0              scater_1.28.0                 ggridges_0.5.4
# [61] survival_3.5-5                tools_4.3.1                   ica_1.0-3
# [64] Rcpp_1.0.11                   glue_1.6.2                    gridExtra_2.3
# [67] withr_2.5.2                   BiocManager_1.30.22           fastmap_1.1.1
# [70] bluster_1.10.0                fansi_1.0.5                   rsvd_1.0.5
# [73] digest_0.6.33                 timechange_0.2.0              R6_2.5.1
# [76] mime_0.12                     colorspace_2.1-0              scattermore_1.2
# [79] tensor_1.5                    spatstat.data_3.0-3           RSQLite_2.3.3
# [82] utf8_1.2.4                    generics_0.1.3                data.table_1.14.8
# [85] rtracklayer_1.60.1            httr_1.4.7                    htmlwidgets_1.6.2
# [88] S4Arrays_1.0.6                uwot_0.1.16                   pkgconfig_2.0.3
# [91] gtable_0.3.4                  blob_1.2.4                    lmtest_0.9-40
# [94] XVector_0.40.0                htmltools_0.5.7               scales_1.3.0
# [97] png_0.1-8                     scran_1.28.2                  rstudioapi_0.15.0
# [100] rjson_0.2.21                  tzdb_0.4.0                    reshape2_1.4.4
# [103] nlme_3.1-162                  curl_5.1.0                    zoo_1.8-12
# [106] cachem_1.0.8                  BiocVersion_3.17.1            KernSmooth_2.23-21
# [109] vipor_0.4.5                   parallel_4.3.1                miniUI_0.1.1.1
# [112] AnnotationDbi_1.62.2          restfulr_0.0.15               pillar_1.9.0
# [115] grid_4.3.1                    vctrs_0.6.4                   RANN_2.6.1
# [118] promises_1.2.1                BiocSingular_1.16.0           beachmat_2.16.0
# [121] dbplyr_2.3.4                  xtable_1.8-4                  cluster_2.1.4
# [124] beeswarm_0.4.0                locfit_1.5-9.8                Rsamtools_2.16.0
# [127] cli_3.6.1                     compiler_4.3.1                rlang_1.1.2
# [130] crayon_1.5.2                  future.apply_1.11.0           ggbeeswarm_0.7.2
# [133] plyr_1.8.9                    stringi_1.8.1                 BiocParallel_1.34.2
# [136] viridisLite_0.4.2             deldir_1.0-9                  munsell_0.5.0
# [139] Biostrings_2.68.1             lazyeval_0.2.2                spatstat.geom_3.2-7
# [142] Matrix_1.6-1                  ExperimentHub_2.8.1           hms_1.1.3
# [145] patchwork_1.2.0               sparseMatrixStats_1.12.2      bit64_4.0.5
# [148] future_1.33.0                 statmod_1.5.0                 KEGGREST_1.40.1
# [151] shiny_1.8.0                   interactiveDisplayBase_1.38.0 AnnotationHub_3.8.0
# [154] ROCR_1.0-11                   igraph_2.0.3                  memoise_2.0.1
# [157] xgboost_1.7.5.1               bit_4.0.5
