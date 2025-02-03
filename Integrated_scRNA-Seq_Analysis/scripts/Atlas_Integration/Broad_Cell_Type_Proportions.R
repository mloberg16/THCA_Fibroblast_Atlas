### Author: Matthew Aaron Loberg
### Date: Novermber 28, 2024
### Script: Broad_Cell_Type_Proportions.R
### Source Script Name: 24-1128_Broad_Cell_Type_Proportions.R

### Goal:
# Make proportion plots for broad cell types for Figure 1/supplement 1
# Will do for each individual orig.ident and for a more simplified by paper/histology plot

### Original Script date: October 9, 2024

### 24-1128 Update
# Repleacing Lee24 normal as Hong et al.
# The Lee24 normal are from the following paper:
# Hong et al. 2024 Endocrinology: "Single Cell Analysis of Human Thyroid Reveals the Transcriptional Signatures of Aging"

### Load packages:
library(Seurat)
library(tidyverse)

### Load Data
Merged_SO_FastMNN <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_AFTER_FastMNN_3000.RDS")

###### INDIVIDUAL SAMPLES ######
# Set the colors
cols <- c("#F0A0FF", "#0075DC", "#993F00", "#9900CC", "#191919", "#005C31", "#2BCE48", "#808080", "#FFCC99", "#94FFB5")

# Pull out preportion data
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final
cells_number <- as.data.frame(table(Merged_SO_FastMNN@active.ident, Merged_SO_FastMNN@meta.data$orig.ident))
plot <- ggplot(cells_number,
               aes(x = Var2, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill",
           width = .8,
           colour = "black") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/BarPlots"
dir.create(savedir)
ggsave(file.path(savedir, "RNA_snn_res.0.6_Final_Labels_by_orig.ident.png"),
       plot,
       height = 5, width = 15, dpi = 600)

### Make the above plot with the x-axis scaled specifically
plot <- ggplot(cells_number,
               aes(x = Var2, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill",
           width = .8,
           colour = "black") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(name ="Sample",
                   limits = c("N3-GEX", "Thy01", "Thy04", "Thy05", "Thy06", "Thy10", "Thy15", # Lee Normal (Hong et al.)
                              "BLANK SPACE",
                              "PT3", "PT5", "PT7", "PT8", "PT9", "PT10", "PT12", # Lee PTC
                              "BLANK SPACE",
                              "AT9", "AT13", "AT16", "AT17", "AT20", # Lee ATC

                              "BLANK SPACE",
                              "Han24_ATC34", "Han24_ATC35", "Han24_ATC36", "Han24_ATC37", # Han ATCs
                              "BLANK SPACE",

                              "NORM03_Lu", "NORM07_Lu", "NORM18_Lu", "NORM19_Lu", "NORM20_Lu", # Lu normal/paratumor
                              "BLANK SPACE",
                              "PTC01_Lu", "PTC02_Lu", "PTC03_Lu", "PTC04_Lu", "PTC05_Lu", "PTC06_Lu", "PTC07_Lu", # Lu PTC
                              "BLANK SPACE",
                              "ATC08_Lu", "ATC09_Lu", "ATC10_Lu", "ATC11_Lu", "ATC12_Lu", "ATC13_Lu", "ATC14_Lu", "ATC15_Lu", "ATC17_Lu", "ATC18_Lu", # Lu ATC
                              "BLANK SPACE",

                              "Luo_NOM_XTZ", # Luo Normal (paratumor)
                              "BLANK SPACE",
                              "Luo_PTC_WJL", "Luo_PTC_XHY", "Luo_PTC_XTZ", # Luo PTCs
                              "BLANK SPACE",
                              "Luo_ATC_LJ", "Luo_ATC_WYF", "Luo_ATC_MSQ", # Luo ATCs
                              "BLANK SPACE",

                              "Pu_PTC01_P", "Pu_PTC02_P", "Pu_PTC03_P", "Pu_PTC05_P", "Pu_PTC08_P", "Pu_PTC09_P", # Pu normal/paratumor
                              "BLANK SPACE",
                              "Pu_PTC01_T", "Pu_PTC02_T", "Pu_PTC02_LeftLN", "Pu_PTC03_T", "Pu_PTC03_LeftLN", "Pu_PTC03_RightLN", "Pu_PTC04_SC", "Pu_PTC05_T", "Pu_PTC05_RightLN", "Pu_PTC06_RightLN", "Pu_PTC07_RightLN", "Pu_PTC08_T", "Pu_PTC09_T", "Pu_PTC10_RightLN", "Pu_PTC11_RightLN", "Pu_PTC11_SC", # Pu PTC primary tumors/LNs/distant mets
                              "BLANK SPACE",

                              "Wang_NT", "Wang_T1L", "Wang_T1R", "Wang_T2L", "Wang_T2R", "Wang_T3L", "Wang_T3R" # Wang et al. samples
                              ))

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/BarPlots"
ggsave(file.path(savedir, "RNA_snn_res.0.6_Final_Labels_by_orig.ident_discrete_x.png"),
       plot,
       height = 5, width = 15, dpi = 600)

# Reordered
# Input desired order:
desired_order <- c("NK/T", "Plasma", "B_Cell", "pDC", "Myeloid", "Endothelial", "Fibroblast", "ATC", "PTC", "Thyrocyte")

# Reorder Var1 in cells_number
cells_number$Var1 <- factor(cells_number$Var1, levels = desired_order)

# reorder cols
cols <- c("#F0A0FF", # NK/T color
          "#005C31", # Plasma color
          "#993F00", # B_Cell color
          "#94FFB5", # pDC color (note, this is a new color, switch with PTC from prior)
          "#191919", # Myeloid color
          "#2BCE48", # Endothelial color
          "#9900CC", # Fibroblast color
          "#0075DC", # ATC color
          "#FFCC99", # PTC color (note, this is a new color, switched with pDC from prior, but I think it pop better)
          "#808080") # Normal/paratumor color

# Plot again
plot <- ggplot(cells_number,
               aes(x = Var2, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill",
           width = .8,
           colour = "black") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(name ="Sample",
                   limits = c("N3-GEX", "Thy01", "Thy04", "Thy05", "Thy06", "Thy10", "Thy15", # Hong Normal
                              "BLANK SPACE",
                              "PT3", "PT5", "PT7", "PT8", "PT9", "PT10", "PT12", # Lee PTCs
                              "BLANK SPACE",
                              "AT9", "AT13", "AT16", "AT17", "AT20", # Lee ATC

                              "BLANK SPACE",
                              "Han24_ATC34", "Han24_ATC35", "Han24_ATC36", "Han24_ATC37", # Han ATCs
                              "BLANK SPACE",

                              "NORM03_Lu", "NORM07_Lu", "NORM18_Lu", "NORM19_Lu", "NORM20_Lu", # Lu normal/paratumor
                              "BLANK SPACE",
                              "PTC01_Lu", "PTC02_Lu", "PTC03_Lu", "PTC04_Lu", "PTC05_Lu", "PTC06_Lu", "PTC07_Lu", # Lu PTC
                              "BLANK SPACE",
                              "ATC08_Lu", "ATC09_Lu", "ATC10_Lu", "ATC11_Lu", "ATC12_Lu", "ATC13_Lu", "ATC14_Lu", "ATC15_Lu", "ATC17_Lu", "ATC18_Lu", # Lu ATC
                              "BLANK SPACE",

                              "Luo_NOM_XTZ", # Luo Normal (paratumor)
                              "BLANK SPACE",
                              "Luo_PTC_WJL", "Luo_PTC_XHY", "Luo_PTC_XTZ", # Luo PTCs
                              "BLANK SPACE",
                              "Luo_ATC_LJ", "Luo_ATC_WYF", "Luo_ATC_MSQ", # Luo ATCs
                              "BLANK SPACE",

                              "Pu_PTC01_P", "Pu_PTC02_P", "Pu_PTC03_P", "Pu_PTC05_P", "Pu_PTC08_P", "Pu_PTC09_P", # Pu normal/paratumor
                              "BLANK SPACE",
                              "Pu_PTC01_T", "Pu_PTC02_T", "Pu_PTC02_LeftLN", "Pu_PTC03_T", "Pu_PTC03_LeftLN", "Pu_PTC03_RightLN", "Pu_PTC04_SC", "Pu_PTC05_T", "Pu_PTC05_RightLN", "Pu_PTC06_RightLN", "Pu_PTC07_RightLN", "Pu_PTC08_T", "Pu_PTC09_T", "Pu_PTC10_RightLN", "Pu_PTC11_RightLN", "Pu_PTC11_SC", # Pu PTC primary tumors/LNs/distant mets
                              "BLANK SPACE",

                              "Wang_NT",
                              "BLANK SPACE",
                              "Wang_T1L", "Wang_T1R", "Wang_T2L", "Wang_T2R", "Wang_T3L", "Wang_T3R" # Wang samples
                   ))

savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/BarPlots"
ggsave(file.path(savedir, "RNA_snn_res.0.6_Final_Labels_by_orig.ident_discrete_x_ordered_y.png"),
       plot,
       height = 5, width = 15, dpi = 600)

###### BY PAPER + HISTOLOGY COMBINED ######
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.6_ClusterLabels_Final
cells_number <- as.data.frame(table(Merged_SO_FastMNN@active.ident, Merged_SO_FastMNN@meta.data$Paper_Histology_Combined))

# Input desired order:
desired_order <- c("NK/T", "Plasma", "B_Cell", "pDC", "Myeloid", "Endothelial", "Fibroblast", "ATC", "PTC", "Thyrocyte")

# Reorder Var1 in cells_number
cells_number$Var1 <- factor(cells_number$Var1, levels = desired_order)

# reorder cols
cols <- c("#F0A0FF", # NK/T color
                   "#005C31", # Plasma color
                   "#993F00", # B_Cell color
                   "#94FFB5", # pDC color (note, this is a new color, switch with PTC from prior)
                   "#191919", # Myeloid color
                   "#2BCE48", # Endothelial color
                   "#9900CC", # Fibroblast color
                   "#0075DC", # ATC color
                   "#FFCC99", # PTC color (note, this is a new color, switched with pDC from prior, but I think it pop better)
                   "#808080") # Normal/paratumor color


plot <- ggplot(cells_number,
               aes(x = Var2, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill",
           width = .8,
           colour = "black") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = c("Lee_Normal",
                              "EMPTY SPACE",
                              "Lee_PTCs",
                              "Lee_ATCs",
                              "EMPTY SPACE",
                              "Han_ATCs",
                              "EMPTY SPACE",
                              "Lu_Paratumor",
                              "Lu_PTCs",
                              "Lu_ATCs",
                              "EMPTY SPACE",
                              "Luo_Paratumor",
                              "Luo_PTCs",
                              "Luo_ATCs",
                              "EMPTY SPACE",
                              "Pu_Paratumor",
                              "Pu_PTCs",
                              "EMPTY SPACE",
                              "Wang_Paratumor",
                              "Wang_PTCs"))
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/BarPlots"
dir.create(savedir)
ggsave(file.path(savedir, "RNA_snn_res.0.6_Final_Labels_by_Paper_Histo_V3.png"),
       plot,
       height = 5, width = 10.5, dpi = 600)


### Session Info
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.4        tidyr_1.3.0
# [8] tibble_3.2.1       ggplot2_3.5.0      tidyverse_2.0.0    SeuratObject_4.1.4 Seurat_4.4.0
#
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.7         magrittr_2.0.3         spatstat.utils_3.1-2
# [6] ggbeeswarm_0.7.2       vctrs_0.6.4            ROCR_1.0-11            spatstat.explore_3.2-5 htmltools_0.5.7
# [11] sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-21     htmlwidgets_1.6.2      ica_1.0-3
# [16] plyr_1.8.9             plotly_4.10.3          zoo_1.8-12             igraph_2.0.3           mime_0.12
# [21] lifecycle_1.0.4        pkgconfig_2.0.3        Matrix_1.6-1           R6_2.5.1               fastmap_1.1.1
# [26] fitdistrplus_1.1-11    future_1.33.0          shiny_1.8.0            digest_0.6.33          colorspace_2.1-0
# [31] patchwork_1.2.0        tensor_1.5             irlba_2.3.5.1          progressr_0.14.0       timechange_0.2.0
# [36] fansi_1.0.5            spatstat.sparse_3.0-3  httr_1.4.7             polyclip_1.10-6        abind_1.4-5
# [41] compiler_4.3.1         withr_2.5.2            MASS_7.3-60            tools_4.3.1            vipor_0.4.5
# [46] lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.12          future.apply_1.11.0    goftest_1.2-3
# [51] glue_1.6.2             nlme_3.1-162           promises_1.2.1         grid_4.3.1             Rtsne_0.16
# [56] cluster_2.1.4          reshape2_1.4.4         generics_0.1.3         gtable_0.3.4           spatstat.data_3.0-3
# [61] tzdb_0.4.0             hms_1.1.3              data.table_1.14.8      sp_2.1-1               utf8_1.2.4
# [66] spatstat.geom_3.2-7    RcppAnnoy_0.0.21       ggrepel_0.9.4          RANN_2.6.1             pillar_1.9.0
# [71] later_1.3.1            splines_4.3.1          lattice_0.21-8         survival_3.5-5         deldir_1.0-9
# [76] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3          scattermore_1.2
# [81] matrixStats_1.1.0      stringi_1.8.1          lazyeval_0.2.2         codetools_0.2-19       cli_3.6.1
# [86] uwot_0.1.16            xtable_1.8-4           reticulate_1.34.0      munsell_0.5.0          Rcpp_1.0.11
# [91] globals_0.16.2         spatstat.random_3.2-1  png_0.1-8              ggrastr_1.0.2          parallel_4.3.1
# [96] ellipsis_0.3.2         listenv_0.9.0          viridisLite_0.4.2      scales_1.3.0           ggridges_0.5.4
# [101] leiden_0.4.3.1         rlang_1.1.2            cowplot_1.1.1
