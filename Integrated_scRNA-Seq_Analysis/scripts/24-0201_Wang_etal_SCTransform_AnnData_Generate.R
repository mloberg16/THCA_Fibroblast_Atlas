# Author: Matthew Aaron Loberg
# Date: February 1, 2024
# Script: 24-0201_Wang_etal_SCTransform_AnnData_Generate.R

# Goal: Creating SCTransformed Seurat h5ad objects for all of the Wang et al. Thyroid scRNA samples
# I will use these samples to merge with all other thyroid scRNA samples + run FastMNN to build my atlas

##### Old script info from Pu script #####
# Creating SCTransformed Seurat Objects and AnnData objects for Pu et al.
# I will use these to merge with other samples + run FastMNN

# 23-1108 Update
# NOT RUNNING SoupX
# Seurat version 4 is REQUIRED for use of SeuratDisk

# 23-1121 Update
# Running ALL Pu Samples this time ... including lymph nodes, mets, etc.
# Note, script finished + run on 23-1122, so output files are from 23-1122

# 23-1127 Update
# I am re-running this with a new Pu Create SO function that will name Cell IDs without overlapping identities

##### New script info running on Wang et al. #####
# 24-0201 Update
# Adapting 23-1127_Pu_etal_SCTransform_AnnData_Generate.R to run on Wang et al. Samples

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
source("function_scripts/24-0201_Wang_Create_SO.R") # New Create SO function that fixes the error of overlapping names
source("function_scripts/23-1009_Seurat_Basic_QC.R")
source("function_scripts/23-1025_DoubletFinder_Function.R")
source("function_scripts/23-1026_SingleR_Prediction_Function.R")
#source("function_scripts/23-1031_AmbientRNA_Function.R") # No longer running SoupX, so will NOT run this (for now)

# Read Directories in as a list
# Read Directories
Wang_readdirs <- list("~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743021_T1L.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743022_T1R.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743023_T2L.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743024_T2R.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743025_T3L.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743026_T3R.h5",
                      "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/GSE191288_RAW/GSM5743027_NT.h5")

# Create Seurat Object in a list for each directory
Wang_PTCs <- list()
for(i in 1:length(Wang_readdirs)){
  Wang_PTCs[[i]] <- Wang_Create_SO(readdir = Wang_readdirs[[i]])
  # Note sure why orig.ident is not working in the function but I will run it myself here
  Wang_PTCs[[i]]$orig.ident <- paste("Wang", substr(Wang_readdirs[[i]], 100, nchar(Wang_readdirs[[i]])-3), sep = '_')
  Wang_PTCs[[i]]$Paper <- "Wang"
}
# Cleaning up - remove readdirs
rm(Wang_readdirs)

# Explore Wang_PTCs to make sure it is what I want - confirmed, commenting out for now
# Wang_PTCs[[1]]
# head(Wang_PTCs[[1]]$orig.ident)
# head(Wang_PTCs[[2]]$orig.ident)
# head(Wang_PTCs[[3]]$orig.ident)
# head(Wang_PTCs[[4]]$orig.ident)
# head(Wang_PTCs[[5]]$orig.ident)
# head(Wang_PTCs[[6]]$orig.ident)
# head(Wang_PTCs[[7]]$orig.ident)

#### Doublet Detection ####
outputdir <- "outputs/Wang_etal_Analysis_Outputs/24-0201_scDblFinder/"
# Run doublet detection
for(i in 1:length(Wang_PTCs)){
  Wang_PTCs[[i]] <- Doublet_Detection(Wang_PTCs[[i]], outputdir = paste0(outputdir, Wang_PTCs[[i]]$orig.ident[1]))
}
rm(outputdir, i) # Cleaning up

# Statistics for doublets + subsetting
Wang_Doublet_Info <- list()
for(i in 1:length(Wang_PTCs)){
  Wang_Doublet_Info[[i]] <- table(Wang_PTCs[[i]]$scDblFinder.class)
  Wang_PTCs[[i]] <- Wang_PTCs[[i]] %>% subset(scDblFinder.class == "singlet")
  Wang_PTCs[[i]]$scDblFinder.class <- NULL # Remove scDbleFinder.class from meta-data
}
saveRDS(Wang_Doublet_Info, file = "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/24-0201_Wang_Doublet_Info.rds")
rm(Wang_Doublet_Info, i)

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
QC_OutputDir <- "outputs/Wang_etal_Analysis_Outputs/24-0201_All_Sample_QC/"
for(i in 1:length(Wang_PTCs)){
  Wang_PTCs[[i]] <- Seurat_Basic_QC(Wang_PTCs[[i]], outputdir = paste0(QC_OutputDir,Wang_PTCs[[i]]$orig.ident[1]))
}
rm(QC_OutputDir, i)

#### Normalization with SCTransform ####
# Run SCTransform for all Pu PTCs using vst.flavor = v2
# Returning more variable features than default (7000 vs 3000)
# Returning ALL genes (not just variable genes) -> doing for purpose of integration later/have SCT values for all genes
Wang_PTCs <- lapply(X = Wang_PTCs,
                    FUN = SCTransform,
                    variable.features.n = 7000,
                    return.only.var.genes = FALSE,
                    vst.flavor = "v2")

### Run SingleR on the Wang_PTCs list do add cell type annotations to the data
# Source: SingleR book - http://bioconductor.org/books/release/SingleRBook/
# https://github.com/dviraran/SingleR/issues/115 - see here for adding SingleR labels back to Seurat object metadata
savedir <- "data_in_use/Wang_etal_2022_PTC_scRNA/24-0201_SingleR_Predictions/" # Set savedir for the tables with prediction info
outputdir <- "outputs/Wang_etal_Analysis_Outputs/24-0201_SingleR_UMAPs/" # set the outputdir for sample umaps
Wang_PTCs <- SingleR_Predictions(SO_List = Wang_PTCs, outputdir = outputdir, savedir = savedir)

### Add in MetaData about samples from the paper
# Driver Mutation
Wang_PTCs[[1]]$Driver_Mut <- "RET/FARP1" # This is Wang_T1L
Wang_PTCs[[2]]$Driver_Mut <- "BRAFV600E" # This is Wang_T1R
Wang_PTCs[[3]]$Driver_Mut <- "BRAFV600E" # This is Wang_T2L
Wang_PTCs[[4]]$Driver_Mut <- "BRAFV600E" # This is Wang_T2R
Wang_PTCs[[5]]$Driver_Mut <- "BRAFV600E" # This is Wang_T3L
Wang_PTCs[[6]]$Driver_Mut <- "BRAFV600E"  # This is Wang_T3R
Wang_PTCs[[7]]$Driver_Mut <- "Normal" # This is Wang_NT

# BRAF Mutation Status
Wang_PTCs[[1]]$BRAFV600E <- "NO" # This is Wang_T1L
Wang_PTCs[[2]]$BRAFV600E <- "YES" # This is Wang_T1R
Wang_PTCs[[3]]$BRAFV600E <- "YES" # This is Wang_T2L
Wang_PTCs[[4]]$BRAFV600E <- "YES" # This is Wang_T2R
Wang_PTCs[[5]]$BRAFV600E <- "YES" # This is Wang_T3L
Wang_PTCs[[6]]$BRAFV600E <- "YES"  # This is Wang_T3R
Wang_PTCs[[7]]$BRAFV600E <- "NO" # This is Wang_NT

#RET Fusion Status
Wang_PTCs[[1]]$RET <- "YES" # This is Wang_T1L
Wang_PTCs[[2]]$RET <- "NO" # This is Wang_T1R
Wang_PTCs[[3]]$RET <- "NO" # This is Wang_T2L
Wang_PTCs[[4]]$RET <- "NO" # This is Wang_T2R
Wang_PTCs[[5]]$RET <- "NO" # This is Wang_T3L
Wang_PTCs[[6]]$RET <- "NO"  # This is Wang_T3R
Wang_PTCs[[7]]$RET <- "NO" # This is Wang_NT

# TERT Mutation Status (note: TERT mutations in patients 4 (subcutaneous met) and 6 (lymph node))
Wang_PTCs[[1]]$TERT_MUT <- "Unknown" # This is Wang_T1L
Wang_PTCs[[2]]$TERT_MUT <- "Unknown" # This is Wang_T1R
Wang_PTCs[[3]]$TERT_MUT <- "Unknown" # This is Wang_T2L
Wang_PTCs[[4]]$TERT_MUT <- "Unknown" # This is Wang_T2R
Wang_PTCs[[5]]$TERT_MUT <- "Unknown" # This is Wang_T3L
Wang_PTCs[[6]]$TERT_MUT <- "Unknown"  # This is Wang_T3R
Wang_PTCs[[7]]$TERT_MUT <- "NO" # This is Wang_NT

# RAS Mutation Status (note: Sequenced NRAS, KRAS, HRAS - no mutations in any of them)
Wang_PTCs[[1]]$RAS_MUT <- "NO" # This is Wang_T1L
Wang_PTCs[[2]]$RAS_MUT <- "NO" # This is Wang_T1R
Wang_PTCs[[3]]$RAS_MUT <- "NO" # This is Wang_T2L
Wang_PTCs[[4]]$RAS_MUT <- "NO" # This is Wang_T2R
Wang_PTCs[[5]]$RAS_MUT <- "NO" # This is Wang_T3L
Wang_PTCs[[6]]$RAS_MUT <- "NO"  # This is Wang_T3R
Wang_PTCs[[7]]$RAS_MUT <- "NO" # This is Wang_NT

# Location
Wang_PTCs[[1]]$Location <- "Thyroid" # This is Wang_T1L
Wang_PTCs[[2]]$Location <- "Thyroid" # This is Wang_T1R
Wang_PTCs[[3]]$Location <- "Thyroid" # This is Wang_T2L
Wang_PTCs[[4]]$Location <- "Thyroid" # This is Wang_T2R
Wang_PTCs[[5]]$Location <- "Thyroid" # This is Wang_T3L
Wang_PTCs[[6]]$Location <- "Thyroid"  # This is Wang_T3R
Wang_PTCs[[7]]$Location <- "Thyroid" # This is Wang_NT

# Histology
Wang_PTCs[[1]]$Histology <- "PTC" # This is Wang_T1L
Wang_PTCs[[2]]$Histology <- "PTC" # This is Wang_T1R
Wang_PTCs[[3]]$Histology <- "PTC" # This is Wang_T2L
Wang_PTCs[[4]]$Histology <- "PTC" # This is Wang_T2R
Wang_PTCs[[5]]$Histology <- "PTC" # This is Wang_T3L
Wang_PTCs[[6]]$Histology <- "PTC"  # This is Wang_T3R
Wang_PTCs[[7]]$Histology <- "Normal" # This is Wang_NT

# Histology_Subtype
Wang_PTCs[[1]]$Histology_Subtype <- "Unknown" # This is Wang_T1L
Wang_PTCs[[2]]$Histology_Subtype <- "Unknown" # This is Wang_T1R
Wang_PTCs[[3]]$Histology_Subtype <- "Unknown" # This is Wang_T2L
Wang_PTCs[[4]]$Histology_Subtype <- "Unknown" # This is Wang_T2R
Wang_PTCs[[5]]$Histology_Subtype <- "Unknown" # This is Wang_T3L
Wang_PTCs[[6]]$Histology_Subtype <- "Unknown"  # This is Wang_T3R
Wang_PTCs[[7]]$Histology_Subtype <- "Normal" # This is Wang_NT

# TNM Stage
Wang_PTCs[[1]]$TNM_Stage <- "N0" # This is Wang_T1L
Wang_PTCs[[2]]$TNM_Stage <- "N0" # This is Wang_T1R
Wang_PTCs[[3]]$TNM_Stage <- "N1b" # This is Wang_T2L
Wang_PTCs[[4]]$TNM_Stage <- "N1b" # This is Wang_T2R
Wang_PTCs[[5]]$TNM_Stage <- "N1a" # This is Wang_T3L
Wang_PTCs[[6]]$TNM_Stage <- "N1a"  # This is Wang_T3R
Wang_PTCs[[7]]$TNM_Stage <- "Normal" # This is Wang_NT

# Concomitant Hashimoto's
Wang_PTCs[[1]]$Hashimotos <- "CLT" # This is Wang_T1L
Wang_PTCs[[2]]$Hashimotos <- "Unknown" # This is Wang_T1R
Wang_PTCs[[3]]$Hashimotos <- "CLT" # This is Wang_T2L
Wang_PTCs[[4]]$Hashimotos <- "Unknown" # This is Wang_T2R
Wang_PTCs[[5]]$Hashimotos <- "Unknown" # This is Wang_T3L
Wang_PTCs[[6]]$Hashimotos <- "Unknown"  # This is Wang_T3R
Wang_PTCs[[7]]$Hashimotos <- "Unknown" # This is Wang_NT

# Treatment History? (Binary: YES or NO)
Wang_PTCs[[1]]$Tx_Hx_Binary <- "Unknown" # This is Wang_T1L
Wang_PTCs[[2]]$Tx_Hx_Binary <- "Unknown" # This is Wang_T1R
Wang_PTCs[[3]]$Tx_Hx_Binary <- "Unknown" # This is Wang_T2L
Wang_PTCs[[4]]$Tx_Hx_Binary <- "Unknown" # This is Wang_T2R
Wang_PTCs[[5]]$Tx_Hx_Binary <- "Unknown" # This is Wang_T3L
Wang_PTCs[[6]]$Tx_Hx_Binary <- "Unknown"  # This is Wang_T3R
Wang_PTCs[[7]]$Tx_Hx_Binary <- "Unknown" # This is Wang_NT

# Treatment History? (Descriptive)
Wang_PTCs[[1]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_T1L
Wang_PTCs[[2]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_T1R
Wang_PTCs[[3]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_T2L
Wang_PTCs[[4]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_T2R
Wang_PTCs[[5]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_T3L
Wang_PTCs[[6]]$Tx_Hx_Descriptive <- "Unknown"  # This is Wang_T3R
Wang_PTCs[[7]]$Tx_Hx_Descriptive <- "Unknown" # This is Wang_NT

# The convert function for h5ad can be affected by the gene length of assays
# Check to make sure that the gene length is the same between assays
for(i in 1:length(Wang_PTCs)){
  print(length(rownames(GetAssayData(Wang_PTCs[[i]], assay = "RNA", slot = "counts"))))
  print(length(rownames(GetAssayData(Wang_PTCs[[i]], assay = "RNA", slot = "data"))))
  print(length(rownames(GetAssayData(Wang_PTCs[[i]], assay = "SCT", slot = "counts"))))
  print(length(rownames(GetAssayData(Wang_PTCs[[i]], assay = "SCT", slot = "data"))))
  print(length(rownames(GetAssayData(Wang_PTCs[[i]], assay = "SCT", slot = "scale.data"))))
}

# Now that we have confirmed that the above is all true, we can proceed with saving/converting the data to h5ad
# I am going to save to the PC to save space in my OneDrive. I will format this the same way that the 2023_Integrated_scRNA-Seq_Aanalysis project is formatted

# Save SCTransformed data for each of the Pu PTCs
savedir <- "~/2023_Integrated_scRNA-Seq_Analysis/data_in_use/Wang_etal_2022_PTC_scRNA/Processed_Data/Individual_Samples/"
for(i in 1:length(Wang_PTCs)){
  savedir_temp <- paste0(savedir, Wang_PTCs[[i]]$orig.ident[1])
  SaveH5Seurat(Wang_PTCs[[i]], filename = paste0(savedir_temp, "/24-0201_SCTransformed.h5Seurat"), overwrite = TRUE) # Only need h5Seurat for now
  #Convert(paste0(savedir_temp, "/23-1108_SCTransformed.h5Seurat"), dest = "h5ad", overwrite = TRUE) # To save space, I will NOT be saving the AnnData object here ... can always come back and run in future
}
rm(list = ls())

##### SESSION INFO #####
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
#   [1] scDblFinder_1.14.0          SingleCellExperiment_1.22.0 SingleR_2.2.0               celldex_1.10.1              SummarizedExperiment_1.30.2
# [6] Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1              S4Vectors_0.38.2
# [11] BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.1.0           SoupX_1.6.2                 SeuratDisk_0.0.0.9020
# [16] SeuratObject_4.1.4          Seurat_4.4.0                lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1
# [21] dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1
# [26] ggplot2_3.5.0               tidyverse_2.0.0
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
