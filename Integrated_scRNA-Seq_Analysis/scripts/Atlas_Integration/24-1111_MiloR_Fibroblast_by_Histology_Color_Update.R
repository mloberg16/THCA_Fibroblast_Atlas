### Matthew Loberg
### November 11, 2024
### Script: 24-1111_MiloR_Fibroblast_by_Histology.R

##### GOAL: #####
# Run differential abundance testing using Milo
# Will compare ATC vs PTC; ATC vs paratumor/normal; PTC vs paratumor/normal
# Interested in how epithelial populations change (thyrocyte, PTC, ATC) and how stromal populations change

##### Info: #####
# This is adapted from 24-0321 scripts: "24-0321_MiloR_Attempt_1.R"
# That was my first attempt at MiloR
# Here, I will apply it to the 24-0819 Integrated atlas


# Wrote this code using the following tutorial:
# https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#5_Finding_markers_of_DA_populations

## 24-1018 Update
# Updating colors to reflect colors used throughout the paper

## 24-111
# Finalizing color update after VW agreed with what I had done


###### Load packages #####
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

##### Load data #####
Merged_SO_FastMNN <- readRDS(file = "~/24-0819_Merged_SOs_scRNA_AFTER_FastMNN_3000.RDS")

##### Data Pre-Processing #####

# Rename histology simplified - "Paratumor/Normal" was NOT syntactically valid
levels(Merged_SO_FastMNN$Histology_Simplified) # Checking levels of Histology Simplified
levels(Merged_SO_FastMNN$Histology_Simplified) <- c("ATC", "PTC", "Para") # Replace levels
levels(Merged_SO_FastMNN$Histology_Simplified) # Check again that levels are replaced

# MiloR requires objects to be in single cell experiment format instead of seurat format
# Create single cell object from seurat object
sce_object <- as.SingleCellExperiment(Merged_SO_FastMNN)
sce_object
plotReducedDim(sce_object, colour_by="Histology_Simplified", dimred = "UMAP") # just making sure it looks as expected

# Going to go ahead and remove Merged_SO_FastMNN for space purposes
rm(Merged_SO_FastMNN)

##### Generate Milo Object #####
# Create Milo object
milo_object <- Milo(sce_object)
rm(sce_object) # remove sce object to save space
milo_object
milo_object <- buildGraph(milo_object, k = 30, d = 30, reduced.dim = "MNN")
milo_object <- makeNhoods(milo_object, prop = 0.1, k = 30, d=30, refined = TRUE, refinement_scheme="graph", reduced_dims = "MNN")
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
dir.create(savedir)
ggsave(file.path(savedir, "24-1111_Milo_NHoodSize.png"),
       plotNhoodSizeHist(milo_object),
       height = 5, width = 8, dpi = 600)
milo_object <- countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")

# Design for Differential Abundance testing
design <- data.frame(colData(milo_object))[,c("orig.ident", "Histology_Simplified", "Paper")]
design <- distinct(design)
rownames(design) <- design$orig.ident
design <- design[colnames(nhoodCounts(milo_object)), , drop=FALSE]
table(design$Histology_Simplified)

############### FIRST DA COMPARISON: ATC samples vs PTC samples ############################

# define first contrast
contrast.1 <- c("Histology_SimplifiedATC - Histology_SimplifiedPTC")
# See GitHub comment here: NO NEED TO DO calcNhoodDistance:
# https://github.com/MarioniLab/miloR/issues/293
# As long as refinement_scheme="graph" is included in "makeNhoods" command
# Example code with for miloR_graph_refinement: https://gist.github.com/emdann/05b4438415e5f5f1ec7fa6bc92d01742
#milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = "MNN")

# generate da_results using testNhoods function
da_results <- testNhoods(milo_object,
                         design = ~ 0 + Histology_Simplified,
                         design.df = design,
                         model.contrasts = contrast.1,
                         fdr.weighting="graph-overlap",
                         reduced.dim = "MNN")

# table of da_results
table(da_results$SpatialFDR < 0.1)

# bulkdNhoodGraph
milo_object <- buildNhoodGraph(milo_object)

# plotUMAP(milo_object, colour_by="Histology_Simplified")
# plotNhoodGraphDA(milo_object, da_results, alpha=0.1) # Careful - this plot takes forever to render

da_results <- annotateNhoods(milo_object, da_results, coldata_col = "Broad_Labels_CAFs_Included")
head(da_results)

ggplot(da_results, aes(Broad_Labels_CAFs_Included_fraction)) + geom_histogram(bins=50)
da_results$Broad_Labels_CAFs_Included <- ifelse(da_results$Broad_Labels_CAFs_Included_fraction < 0.7, "Mixed", da_results$Broad_Labels_CAFs_Included)

# Test plot to see what it looks like WITHOUT modification
plotDAbeeswarm(da_results, group.by = "Broad_Labels_CAFs_Included")

# Make da_results_modified to be just the categories of interest
# Will start with just Thyrocyte/PTC/ATC (da_results_modified_thyroid)
da_results_modified_thyroid <- da_results %>% subset(Broad_Labels_CAFs_Included == "Thyrocyte" |
                                                     Broad_Labels_CAFs_Included == "PTC" |
                                                     Broad_Labels_CAFs_Included == "ATC")

# Check levels for graphing order
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_thyroid$Broad_Labels_CAFs_Included <- factor(da_results_modified_thyroid$Broad_Labels_CAFs_Included,
                                                                 levels = c("ATC",
                                                                            "PTC",
                                                                            "Thyrocyte"))
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # check levels

# Testing to make sure this looks how I want it to
plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included")

# Save the plot
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
# trace(plotDAbeeswarm, edit = TRUE) # edit colors here in scale_color_gradient2 (low and high colors)
# I have the above line of code commented out; can run it if plot color change is desired
ggsave(file.path(savedir, "24-1111_Milo_Thyroid_DAbeeswarm_ATCvsPTC.png"),
       plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 5, width = 10, dpi = 600)

# Option to saveRDS here; 24-1015 update NOT doing this
#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results_modified.RDS")

##### MODIFY FOR FIBROBLASTS #####

# Subset down to just fibroblast groups
da_results_modified_fibroblasts <- da_results %>% subset(Broad_Labels_CAFs_Included == "myCAF" |
                                                           Broad_Labels_CAFs_Included == "iCAF" |
                                                           #Broad_Labels_CAFs_Included == "iCAF2" | # Not including iCAF2
                                                           Broad_Labels_CAFs_Included == "iPVCAF" |
                                                           Broad_Labels_CAFs_Included == "dPVCAF" |
                                                           Broad_Labels_CAFs_Included == "APOE+_CAF")

# Check levels for graphing order
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_fibroblasts$Broad_Labels_CAFs_Included <- factor(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included,
                                                                     levels = c(#"iCAF2", # Not including iCAF2
                                                                                "iCAF",
                                                                                "myCAF",
                                                                                "APOE+_CAF",
                                                                                "iPVCAF",
                                                                                "dPVCAF"))
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # check levels

# Now make the plot
savedir = savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
ggsave(file.path(savedir, "24-1111_Milo_Fibroblast_DAbeeswarm_ATCvsPTC.png"),
       plotDAbeeswarm(da_results_modified_fibroblasts, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 7, width = 10, dpi = 600)

# option to saveRDS - not doing that here
#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results_modified.RDS")


############## Repeating above for ATC samples vs Normal/Paratumor samples #################

contrast.1 <- c("Histology_SimplifiedATC - Histology_SimplifiedPara")
# See GitHub comment here: NO NEED TO DO calcNhoodDistance:
# https://github.com/MarioniLab/miloR/issues/293
# As long as refinement_scheme="graph" is included in "makeNhoods" command
# Example code with for miloR_graph_refinement: https://gist.github.com/emdann/05b4438415e5f5f1ec7fa6bc92d01742
#milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = "MNN")
da_results <- testNhoods(milo_object,
                         design = ~ 0 + Histology_Simplified,
                         design.df = design,
                         model.contrasts = contrast.1,
                         fdr.weighting="graph-overlap",
                         reduced.dim = "MNN")
table(da_results$SpatialFDR < 0.1)
milo_object <- buildNhoodGraph(milo_object)

# plotUMAP(milo_object, colour_by="Histology_Simplified")
# plotNhoodGraphDA(milo_object, da_results, alpha=0.1) # Careful - this plot takes forever to render

da_results <- annotateNhoods(milo_object, da_results, coldata_col = "Broad_Labels_CAFs_Included")
head(da_results)

ggplot(da_results, aes(Broad_Labels_CAFs_Included_fraction)) + geom_histogram(bins=50)
da_results$Broad_Labels_CAFs_Included <- ifelse(da_results$Broad_Labels_CAFs_Included_fraction < 0.7, "Mixed", da_results$Broad_Labels_CAFs_Included)

# Test plot to see what it looks like WITHOUT modification
# trace(plotDAbeeswarm, edit = TRUE) # adjust colors as desired within scale_gradient2
# I have the above line of code commented out; can run it if plot color change is desired
plotDAbeeswarm(da_results, group.by = "Broad_Labels_CAFs_Included")

##### MODIFY FOR Thyrocyte/Tumor populations #####

# modify to just groups of interest
da_results_modified_thyroid <- da_results %>% subset(Broad_Labels_CAFs_Included == "Thyrocyte" |
                                                     Broad_Labels_CAFs_Included == "PTC" |
                                                     Broad_Labels_CAFs_Included == "ATC")

# Check levels for graphing order
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_thyroid$Broad_Labels_CAFs_Included <- factor(da_results_modified_thyroid$Broad_Labels_CAFs_Included,
                                                                 levels = c("ATC",
                                                                            "PTC",
                                                                            "Thyrocyte"))
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # check levels

# Testing to make sure this looks how I want it to
plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included")

# Save the plot
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
ggsave(file.path(savedir, "24-111_Milo_Thyroid_DAbeeswarm_ATCvsNormal.png"),
       plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 5, width = 10, dpi = 600)

# Option to saveRDS here; 24-1015 NOT doing this
#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results_modified.RDS")


##### MODIFY FOR FIBROBLASTS #####

# Subset down to just fibroblast groups
da_results_modified_fibroblasts <- da_results %>% subset(Broad_Labels_CAFs_Included == "myCAF" |
                                                         Broad_Labels_CAFs_Included == "iCAF" |
                                                         #Broad_Labels_CAFs_Included == "iCAF2" |
                                                         Broad_Labels_CAFs_Included == "iPVCAF" |
                                                         Broad_Labels_CAFs_Included == "dPVCAF" |
                                                         Broad_Labels_CAFs_Included == "APOE+_CAF")

# Check levels for graphing order
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_fibroblasts$Broad_Labels_CAFs_Included <- factor(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included,
                                                                     levels = c(#"iCAF2",
                                                                                "iCAF",
                                                                                "myCAF",
                                                                                "APOE+_CAF",
                                                                                "iPVCAF",
                                                                                "dPVCAF"))
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # check levels

# Now make the plot
savedir = savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
ggsave(file.path(savedir, "24-1111_Milo_Fibroblast_DAbeeswarm_ATCvsNormal.png"),
       plotDAbeeswarm(da_results_modified_fibroblasts, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 7, width = 10, dpi = 600)

#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results_modified.RDS")


############## Repeating above for PTC samples vs Normal/paratumor samples #################

contrast.1 <- c("Histology_SimplifiedPTC - Histology_SimplifiedPara")
# See GitHub comment here: NO NEED TO DO calcNhoodDistance:
# https://github.com/MarioniLab/miloR/issues/293
# As long as refinement_scheme="graph" is included in "makeNhoods" command
# Example code with for miloR_graph_refinement: https://gist.github.com/emdann/05b4438415e5f5f1ec7fa6bc92d01742
#milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = "MNN")
da_results <- testNhoods(milo_object,
                         design = ~ 0 + Histology_Simplified,
                         design.df = design,
                         model.contrasts = contrast.1,
                         fdr.weighting="graph-overlap",
                         reduced.dim = "MNN")
table(da_results$SpatialFDR < 0.1)
milo_object <- buildNhoodGraph(milo_object)

# plotUMAP(milo_object, colour_by="Histology_Simplified")
# plotNhoodGraphDA(milo_object, da_results, alpha=0.1) # Careful - this plot takes forever to render

da_results <- annotateNhoods(milo_object, da_results, coldata_col = "Broad_Labels_CAFs_Included")
head(da_results)

ggplot(da_results, aes(Broad_Labels_CAFs_Included_fraction)) + geom_histogram(bins=50)
da_results$Broad_Labels_CAFs_Included <- ifelse(da_results$Broad_Labels_CAFs_Included_fraction < 0.7, "Mixed", da_results$Broad_Labels_CAFs_Included)

# Test plot to see what it looks like WITHOUT modification
# trace(plotDAbeeswarm, edit = TRUE)
# I have the above line of code commented out; can run it if plot color change is desired
plotDAbeeswarm(da_results, group.by = "Broad_Labels_CAFs_Included")

##### MODIFY FOR Thyrocyte/Tumor populations #####

da_results_modified_thyroid <- da_results %>% subset(Broad_Labels_CAFs_Included == "Thyrocyte" |
                                                       Broad_Labels_CAFs_Included == "PTC" |
                                                       Broad_Labels_CAFs_Included == "ATC")

# Check levels for graphing order
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_thyroid$Broad_Labels_CAFs_Included <- factor(da_results_modified_thyroid$Broad_Labels_CAFs_Included,
                                                                 levels = c("ATC",
                                                                            "PTC",
                                                                            "Thyrocyte"))
levels(da_results_modified_thyroid$Broad_Labels_CAFs_Included) # check levels

# Testing to make sure this looks how I want it to
plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included")

# Save the plot
savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
ggsave(file.path(savedir, "24-1111_Milo_Thyroid_DAbeeswarm_PTCvsNormal.png"),
       plotDAbeeswarm(da_results_modified_thyroid, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 5, width = 10, dpi = 600)

# Option to saveRDS here; 24-1015 NOT doing this
#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_PTC_da_results_modified.RDS")


##### MODIFY FOR FIBROBLASTS #####

# Subset down to just fibroblast groups
da_results_modified_fibroblasts <- da_results %>% subset(Broad_Labels_CAFs_Included == "myCAF" |
                                                           Broad_Labels_CAFs_Included == "iCAF" |
                                                           #Broad_Labels_CAFs_Included == "iCAF2" |
                                                           Broad_Labels_CAFs_Included == "iPVCAF" |
                                                           Broad_Labels_CAFs_Included == "dPVCAF" |
                                                           Broad_Labels_CAFs_Included == "APOE+_CAF")

# Check levels for graphing order
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # No levels - must create levels
# make into a factor w/ levels
da_results_modified_fibroblasts$Broad_Labels_CAFs_Included <- factor(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included,
                                                                     levels = c(#"iCAF2",
                                                                                "iCAF",
                                                                                "myCAF",
                                                                                "APOE+_CAF",
                                                                                "iPVCAF",
                                                                                "dPVCAF"))
levels(da_results_modified_fibroblasts$Broad_Labels_CAFs_Included) # check levels

# Now make the plot
savedir = savedir = "outputs/24-0819_Atlas_Integration/FastMNN_3000/MiloPlots"
ggsave(file.path(savedir, "24-1111_Milo_Fibroblast_DAbeeswarm_PTCvsNormal.png"),
       plotDAbeeswarm(da_results_modified_fibroblasts, group.by = "Broad_Labels_CAFs_Included") +
         theme(
           axis.title.x = element_text(face = "bold"),
           axis.title.y = element_blank(),
           axis.text.y = element_text(size = 15)),
       height = 7, width = 10, dpi = 600)

#saveRDS(da_results, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results.RDS")
#saveRDS(da_results_modified, file = "data_in_use/Integrated_Data/24-0321_ATC_Normal_da_results_modified.RDS")

########## CLEANING UP #############
rm(list = ls())

########## SESSION INFO ############
sessionInfo()

