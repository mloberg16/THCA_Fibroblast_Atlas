### Author: Matthew Aaron Loberg
### Date: December 13, 2024
### Script: FastMNN_Fibroblast_Subclustering_3000_PROGENy_Scores.R
###Source Script Name: 24-1213_FastMNN_Fibroblast_Subclustering_3000_PROGENy_Scores.R

### Goal ###
# Load in the Fibroblast subclustering object from 24-0821
# Calculate PROGENy scores to look at key activated signaling pathways

### Additional INFO ###
# See this paper: https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00319-2
# They used the PROGENy method
# Their GitHub of this is here: https://github.com/aliceygao/pan-Fibroblast/blob/main/Figure3.R
# I will follow their Figure 3 PROGENy GitHub
# PROGENy is a bioconductor package
# Here is a link to the PDF for PROGENy manual: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://bioconductor.org/packages/release/bioc/manuals/progeny/man/progeny.pdf

### PROGENy Installation ###
# Only do ONCE
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("progeny")
# Note: hit update none for the package updates

##### Load packages #####
library(Seurat)
library(SeuratWrappers)
library(tidyverse) # for ggsave
library(RColorBrewer) # For plot customization
library(progeny)

##### READ RDS Seurat Object of stromal cell subclustering #####
Merged_SO_FastMNN <- readRDS(file = "~/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS")

##### RUNNING PROGENy according to the paper's GitHub #####
# See link here: https://github.com/aliceygao/pan-Fibroblast/blob/main/Figure3.R
metaData <- Merged_SO_FastMNN@meta.data

### Make a CellsClusters df for annotation
CellsClusters <- data.frame(
  Cell = colnames(Merged_SO_FastMNN),
  CellType = as.character(Merged_SO_FastMNN@meta.data$CAF_Labels),
  stringsAsFactors = FALSE
)

### Run progeny
progeny_results <- progeny(Merged_SO_FastMNN, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
# scale the results
progeny_results <- Seurat::ScaleData(progeny_results, assay = "progeny") # scale the results
# Extract the results to a data frame
progeny_scores_df <- as.data.frame(t(GetAssayData(progeny_results, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

### Join the progeny scores df and the CellsClusters df
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

### Summarize progeny scores
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
sub_progeny_scores <- summarized_progeny_scores #%>% filter(Pathway %in% c("TGFb", "WNT", "PI3K", "Hypoxia", "TNFa", "NFkB"))

### Subset to JUST myCAF, iCAF, APOE+_CAF, dPVCAF, iPVCAF (e.g., remove iCAF2)
cluster.lineages <- c("myCAF", "iCAF", "APOE+_CAF", "dPVCAF", "iPVCAF")
sub_progeny_scores <- sub_progeny_scores %>% filter(CellType %in% cluster.lineages)

### Setting variables for plot
# Make cell type a factor
sub_progeny_scores$CellType <- factor(sub_progeny_scores$CellType, levels = rev(cluster.lineages))

# Set up "group" variable similar to their paper
sub_progeny_scores$group <- "PVL"
sub_progeny_scores$group[sub_progeny_scores$CellType %in% c("myCAF", "iCAF")] <- "CAF"
sub_progeny_scores$group <- factor(sub_progeny_scores$group, levels = c("CAF", "PVL"))

# groupMean variable calculation
groupMean <- sub_progeny_scores %>%
  dplyr::group_by(Pathway) %>%
  dplyr::mutate(MeanV = mean(avg)) %>%
  dplyr::distinct(Pathway, .keep_all = TRUE) %>%
  dplyr::select(Pathway, MeanV)

# make plot
g1 <- ggplot(sub_progeny_scores, aes(x = CellType, y = avg, color = CellType)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("myCAF" = "#783FC1",
                                "iCAF" = "#FE8F42",
                                "dPVCAF" = "#990000",
                                "iPVCAF" = "#0075DC",
                                "APOE+_CAF" = "#426600")) +
  geom_hline(data = groupMean, aes(yintercept = MeanV), linetype = "dashed") +
  coord_flip() +
  theme_gray(base_size = 12) +
  theme(axis.text.x = element_text(size = 7)) +
  facet_grid(group ~ Pathway, scales = "free")

# This is the ggsave from the paper
# ggsave(sprintf("%s/Figure3D.pdf", figdir), g1, width = 6, height = 4, useDingbats = F)
# If want to save, make own ggsave with directory specific command

ggsave("outputs/Fibroblast_Subclustering/24-0821_Fibroblast_Subclustering/FastMNN_3000/PROGENy_Plots/24-1213_PROGENy_Plot.png",
       g1, width = 16, height = 4, dpi = 600, create.dir = TRUE)


### Cleaning up 
rm(list = ls()
