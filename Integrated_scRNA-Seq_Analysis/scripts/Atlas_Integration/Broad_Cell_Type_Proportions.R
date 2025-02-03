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
