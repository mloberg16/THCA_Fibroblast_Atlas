### Author: Matthew Aaron Loberg
### Date: November 28, 2024
### Script: Fibroblast_Subclustering_FastMNN_3000_Cell_Type_Proportions.R
### Source Script Name: 24-1128_Fibroblast_Subclustering_FastMNN_3000_Cell_Type_Proportions.R

##### Goal: #####
# Make the following plots for the fibroblast subclustering object:
# 1. By paper proportions
# 2. By paper + histology combined proportions
# 3. By individual sample proportions

### 24-1128 Update
# This is an update of 24-1014_FastMNN_Fibroblast_Subclustering_Cell_Type_Proportions.R
# In this update, I am adding a designation for HONG et al. paper (previously had as Lee et al. normal)

##### Load Packages #####
library(Seurat)
library(tidyverse)

##### Load Data #####
Merged_SO_FastMNN <- readRDS(file = "~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

##### Proportions by "Paper_Simplified" #####
# Use Cohen et al. code as a reference
# Pull out paper cluster proportions
# Also used this stack overflow post as reference: https://bioinformatics.stackexchange.com/questions/11150/percentage-of-each-cluster-in-seurat

# First, need to make Paper_Simplified from Paper
# The 24-1128 update of Paper_Simplified incorporates the addition of Hong et al. in place of Lee normals
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$Paper
levels(Idents(Merged_SO_FastMNN))
Merged_SO_FastMNN$Paper_Simplified <- Idents(Merged_SO_FastMNN)
levels(Merged_SO_FastMNN$Paper_Simplified)
levels(Merged_SO_FastMNN$Paper_Simplified) <- c("Lee", "Han", "Lu", "Luo", "Pu", "Hong", "Wang")
levels(Merged_SO_FastMNN$Paper_Simplified)

Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.2
cells_number <- as.data.frame(table(Merged_SO_FastMNN@active.ident, Merged_SO_FastMNN@meta.data$Paper_Simplified))

desired_order <- c("3", "0", "4", "2", "1", "5")
cols <- c("#990000", "#0075DC", "#426600", "#783FC1", "#FE8F42", "#808080")
cells_number$Var1 <- factor(cells_number$Var1, levels = desired_order)
levels(cells_number$Var1)

# If percentage is needed specifically
#cells_number$Percentage <- with(cells_number, Freq / ave(Freq, Var2, FUN = sum) * 100)
# Save as well
#write.csv(cells_number, file = "data_in_use/Integrated_Data/24-1128_Stromal_Cell_Type_Proportions_by_Paper.csv")


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
  scale_x_discrete(name = "Sample", limits = c("Hong", "Lee", "Han", "Lu", "Luo", "Pu", "Wang"))
savedir = "outputs/Fibroblast_subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/BarPlots"
dir.create(savedir)
ggsave(file.path(savedir, "Fibroblast_RNA_snn_res.0.2_Proportion_by_Paper_V2.png"),
       plot,
       height = 5, width = 4, dpi = 600)

##### Proportions by "Paper_Histology_Combined" #####

# First, need to make "Paper_Histology_Combined" variable
# This is a more in depth proportion by paper simplified that splits by histology

# Lee ATCs
Lee_ATCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Lee" & Histology_Simplified == "ATC")
Lee_ATCs$Paper_Histology_Combined <- "Lee_ATCs"
Lee_ATCs_Annotation <- as.data.frame(Lee_ATCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Lee_ATCs$Paper_Histology_Combined")
rm(Lee_ATCs)

# Lee PTCs
Lee_PTCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Lee" & Histology_Simplified == "PTC")
Lee_PTCs$Paper_Histology_Combined <- "Lee_PTCs"
Lee_PTCs_Annotation <- as.data.frame(Lee_PTCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Lee_PTCs$Paper_Histology_Combined")
rm(Lee_PTCs)

# Lee Normal - replacing with HONG NORMAL (24-1128)
Hong_Normal <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Hong" & Histology_Simplified == "Paratumor/Normal")
Hong_Normal$Paper_Histology_Combined <- "Hong_Normal"
Hong_Normal_Annotation <- as.data.frame(Hong_Normal$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Hong_Normal$Paper_Histology_Combined")
rm(Hong_Normal)

# Lu ATCs
Lu_ATCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Lu" & Histology_Simplified == "ATC")
Lu_ATCs$Paper_Histology_Combined <- "Lu_ATCs"
Lu_ATCs_Annotation <- as.data.frame(Lu_ATCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Lu_ATCs$Paper_Histology_Combined")
rm(Lu_ATCs)

# Lu PTCs
Lu_PTCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Lu" & Histology_Simplified == "PTC")
Lu_PTCs$Paper_Histology_Combined <- "Lu_PTCs"
Lu_PTCs_Annotation <- as.data.frame(Lu_PTCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Lu_PTCs$Paper_Histology_Combined")
rm(Lu_PTCs)

# Lu Normal
Lu_Normal <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Lu" & Histology_Simplified == "Paratumor/Normal")
Lu_Normal$Paper_Histology_Combined <- "Lu_Paratumor"
Lu_Normal_Annotation <- as.data.frame(Lu_Normal$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Lu_Normal$Paper_Histology_Combined")
rm(Lu_Normal)

# Luo ATCs
Luo_ATCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Luo" & Histology_Simplified == "ATC")
Luo_ATCs$Paper_Histology_Combined <- "Luo_ATCs"
Luo_ATCs_Annotation <- as.data.frame(Luo_ATCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Luo_ATCs$Paper_Histology_Combined")
rm(Luo_ATCs)

# Luo PTCs
Luo_PTCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Luo" & Histology_Simplified == "PTC")
Luo_PTCs$Paper_Histology_Combined <- "Luo_PTCs"
Luo_PTCs_Annotation <- as.data.frame(Luo_PTCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Luo_PTCs$Paper_Histology_Combined")
rm(Luo_PTCs)

# Luo Normal
Luo_Normal <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Luo" & Histology_Simplified == "Paratumor/Normal")
Luo_Normal$Paper_Histology_Combined <- "Luo_Paratumor"
Luo_Normal_Annotation <- as.data.frame(Luo_Normal$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Luo_Normal$Paper_Histology_Combined")
rm(Luo_Normal)

# Han ATCs
Han_ATCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Han" & Histology_Simplified == "ATC")
Han_ATCs$Paper_Histology_Combined <- "Han_ATCs"
Han_ATCs_Annotation <- as.data.frame(Han_ATCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Han_ATCs$Paper_Histology_Combined")
rm(Han_ATCs)

# Pu PTCs
Pu_PTCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Pu" & Histology_Simplified == "PTC")
Pu_PTCs$Paper_Histology_Combined <- "Pu_PTCs"
Pu_PTCs_Annotation <- as.data.frame(Pu_PTCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Pu_PTCs$Paper_Histology_Combined")
rm(Pu_PTCs)

# Pu Normal
Pu_Normal <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Pu" & Histology_Simplified == "Paratumor/Normal")
Pu_Normal$Paper_Histology_Combined <- "Pu_Paratumor"
Pu_Normal_Annotation <- as.data.frame(Pu_Normal$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Pu_Normal$Paper_Histology_Combined")
rm(Pu_Normal)

# Wang PTCs
Wang_PTCs <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Wang" & Histology_Simplified == "PTC")
Wang_PTCs$Paper_Histology_Combined <- "Wang_PTCs"
Wang_PTCs_Annotation <- as.data.frame(Wang_PTCs$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Wang_PTCs$Paper_Histology_Combined")
rm(Wang_PTCs)

# Wang Normal
Wang_Normal <- Merged_SO_FastMNN %>% subset(Paper_Simplified == "Wang" & Histology_Simplified == "Paratumor/Normal")
Wang_Normal$Paper_Histology_Combined <- "Wang_Paratumor"
Wang_Normal_Annotation <- as.data.frame(Wang_Normal$Paper_Histology_Combined) %>% dplyr::rename("Paper_Histology_Combined" = "Wang_Normal$Paper_Histology_Combined")
rm(Wang_Normal)

# rbind all
Paper_Histology_Combined <- rbind(Lee_ATCs_Annotation,
                                  Lee_PTCs_Annotation,
                                  Hong_Normal_Annotation,
                                  Lu_ATCs_Annotation,
                                  Lu_PTCs_Annotation,
                                  Lu_Normal_Annotation,
                                  Luo_ATCs_Annotation,
                                  Luo_PTCs_Annotation,
                                  Luo_Normal_Annotation,
                                  Han_ATCs_Annotation,
                                  Pu_PTCs_Annotation,
                                  Pu_Normal_Annotation,
                                  Wang_PTCs_Annotation,
                                  Wang_Normal_Annotation)

# Add meta data
Merged_SO_FastMNN$Paper_Histology_Combined <- Paper_Histology_Combined

# Cleaning up
rm(Paper_Histology_Combined, Lee_ATCs_Annotation, Lee_PTCs_Annotation, Hong_Normal_Annotation,
   Lu_ATCs_Annotation, Lu_PTCs_Annotation, Lu_Normal_Annotation, Luo_ATCs_Annotation, Luo_PTCs_Annotation,
   Luo_Normal_Annotation, Han_ATCs_Annotation, Pu_PTCs_Annotation, Pu_Normal_Annotation, Wang_PTCs_Annotation, Wang_Normal_Annotation)

### Now make Paper_Histology_Combined plot
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.2
cells_number <- as.data.frame(table(Merged_SO_FastMNN@active.ident, Merged_SO_FastMNN@meta.data$Paper_Histology_Combined))
cells_number$Var1 <- as.character(cells_number$Var1)

# If percentage is needed specifically
#cells_number$Percentage <- with(cells_number, Freq / ave(Freq, Var2, FUN = sum) * 100)
# Save as well
#write.csv(cells_number, file = "data_in_use/Integrated_Data/24-1128_Stromal_Cell_Type_Proportions_by_Paper_Histology_Combined.csv")

desired_order <- c("3", "0", "4", "2", "1", "5")
cols <- c("#990000", "#0075DC", "#426600", "#783FC1", "#FE8F42", "#808080")
cells_number$Var1 <- factor(cells_number$Var1, levels = desired_order)
levels(cells_number$Var1)

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
  scale_x_discrete(limits = c("Hong_Normal",
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

savedir = "outputs/Fibroblast_subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/BarPlots"
dir.create(savedir)
ggsave(file.path(savedir, "Fibroblast_RNA_snn_res.0.2_Prop_by_Paper_Histology_Combined_V2.png"),
       plot,
       height = 5, width = 10, dpi = 600)


###### INDIVIDUAL SAMPLES ######

# Pull out preportion data
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$RNA_snn_res.0.2
cells_number <- as.data.frame(table(Merged_SO_FastMNN@active.ident, Merged_SO_FastMNN@meta.data$orig.ident))
desired_order <- c("3", "0", "4", "2", "1", "5")
cols <- c("#990000", "#0075DC", "#426600", "#783FC1", "#FE8F42", "#808080")
cells_number$Var1 <- factor(cells_number$Var1, levels = desired_order)
levels(cells_number$Var1)

# If percentage is needed specifically
#cells_number$Percentage <- with(cells_number, Freq / ave(Freq, Var2, FUN = sum) * 100)
# Save as well
#write.csv(cells_number, file = "data_in_use/Integrated_Data/24-1128_Stromal_Cell_Type_Proportions_by_Sample.csv")

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
                 limits = c("N3-GEX", "Thy01", "Thy04", "Thy05", "Thy06", "Thy10", "Thy15", # Lee Normal
                            "BLANK SPACE",
                            "PT3", "PT5", "PT7", "PT8", "PT9", "PT10", "PT12",
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
                            "Luo_PTC_WJL", "Luo_PTC_XHY", "Luo_PTC_XTZ",
                            "BLANK SPACE",
                            "Luo_ATC_LJ", "Luo_ATC_WYF", "Luo_ATC_MSQ", # Luo ATCs
                            "BLANK SPACE",

                            "Pu_PTC01_P", "Pu_PTC02_P", "Pu_PTC03_P", "Pu_PTC05_P", "Pu_PTC08_P", "Pu_PTC09_P", # Pu normal/paratumor
                            "BLANK SPACE",
                            "Pu_PTC01_T", "Pu_PTC02_T", "Pu_PTC02_LeftLN", "Pu_PTC03_T", "Pu_PTC03_LeftLN", "Pu_PTC03_RightLN", "Pu_PTC04_SC", "Pu_PTC05_T", "Pu_PTC05_RightLN", "Pu_PTC06_RightLN", "Pu_PTC07_RightLN", "Pu_PTC08_T", "Pu_PTC09_T", "Pu_PTC10_RightLN", "Pu_PTC11_RightLN", "Pu_PTC11_SC", # Pu PTC primary tumors/LNs/distant mets
                            "BLANK SPACE",

                            "Wang_NT", "Wang_T1L", "Wang_T1R", "Wang_T2L", "Wang_T2R", "Wang_T3L", "Wang_T3R"
                 ))

savedir = "outputs/Fibroblast_subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/BarPlots"
dir.create(savedir)
ggsave(file.path(savedir, "Fibroblast_RNA_snn_res.0.2_Proportion_by_orig.ident_V2.png"),
       plot,
       height = 5, width = 15, dpi = 600)


##### OUTPUT TOTAL FIBROBLAST NUMBERS #####
Idents(Merged_SO_FastMNN) <- Merged_SO_FastMNN$orig.ident
IDs <- levels(Merged_SO_FastMNN)

total_fibroblast <- list()
for(i in 1:length(IDs)){
  total_fibroblast[[i]] <- 0
  for(n in 1:nrow(cells_number)){
    if(as.character(cells_number$Var2[n]) == IDs[i]){
      total_fibroblast[[i]] <- total_fibroblast[[i]] + cells_number$Freq[n]
      #print(IDs[i])
    }
  }
  print(paste0(IDs[i], ": ", total_fibroblast[[i]], " Fibroblasts"))
}

##### OUTPUT IDs WITH myCAFs #####
cells_number_2 <- cells_number %>% subset(Var1 == 2)
cells_number_2 <- cells_number_2 %>% subset(Freq > 0)
print(cells_number_2$Var2)


### SESSION INFO ###
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
