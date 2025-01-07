# Author: Matthew Loberg
# Date: August 21, 2024
# Script: 24-0821_AddFibroblastModuleScores.R

### Function Info
# In this function, I will add module scores from fibroblast gene sets to a given seurat object
# The fibroblast gene sets are from the following papers:
# Cohen et al. 2024 Nature Communications (CAF-S4 gene set) (pericyte-like fibroblasts)
# Kieffer et al. 2020 Cancer Discovery (Normal fibroblast, CAF-S1, and CAF-S1 subset gene sets)

### 24-0821 Update
# This is updated from the 24-0812 version of this script
# I added Wu et al. scores to this version of the script

# Function
AddFibroblastScores <- function(SO){

  ### Broad Subsets

  # CAF S4 Broad Subset:
  CAFS4 <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Cohen_etal_CAFS4.txt", header = FALSE, sep = NULL)
  CAFS4_list <- list(CAFS4$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = CAFS4_list,
                                      name = "CAFS4_Mod_Score")
  rm(CAFS4, CAFS4_list)

  # CAF S1 Broad Subset:
  CAFS1 <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Kieffer_etal_CAFS1.txt", header = FALSE, sep = NULL)
  CAFS1_list <- list(CAFS1$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = CAFS1_list,
                                      name = "CAFS1_Mod_Score")
  rm(CAFS1, CAFS1_list)

  # Normal Fibroblast Broad Subset:
  Normal <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Kieffer_etal_Normal_Fibroblast.txt", header = FALSE, sep = NULL)
  Normal_list <- list(Normal$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = Normal_list,
                                      name = "Normal_Fibroblast_Mod_Score")
  rm(Normal, Normal_list)

  ### CAF S1 Breast Cancer Subsets

  # Kieffer Detox iCAF (CAF-S1 subset)
  Detox_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_Detox_iCAF.txt", header = FALSE, sep = NULL)
  Detox_iCAF_list <- list(Detox_iCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = Detox_iCAF_list,
                                      name = "Detox_iCAF_Mod_Score")
  rm(Detox_iCAF, Detox_iCAF_list)

  # Kieffer ecm myCAF (CAF-S1 subset)
  ecm_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_ecm_myCAF.txt", header = FALSE, sep = NULL)
  ecm_myCAF_list <- list(ecm_myCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = ecm_myCAF_list,
                                      name = "ecm_myCAF_Mod_Score")
  rm(ecm_myCAF, ecm_myCAF_list)

  # Kieffer IFNg iCAF (CAF-S1 subset)
  IFNg_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_IFNg_iCAF.txt", header = FALSE, sep = NULL)
  IFNg_iCAF_list <- list(IFNg_iCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = IFNg_iCAF_list,
                                      name = "IFNg_iCAF_Mod_Score")
  rm(IFNg_iCAF, IFNg_iCAF_list)

  # Kieffer IL iCAF (CAF-S1 subset)
  IL_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_IL_iCAF.txt", header = FALSE, sep = NULL)
  IL_iCAF_list <- list(IL_iCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = IL_iCAF_list,
                                      name = "IL_iCAF_Mod_Score")
  rm(IL_iCAF, IL_iCAF_list)

  # Kieffer TGFB myCAF (CAF-S1 subset)
  TGFB_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_TGFB_myCAF.txt", header = FALSE, sep = NULL)
  TGFB_myCAF_list <- list(TGFB_myCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = TGFB_myCAF_list,
                                      name = "TGFB_myCAF_Mod_Score")
  rm(TGFB_myCAF, TGFB_myCAF_list)

  # Kieffer TGFB myCAF (CAF-S1 subset)
  Wound_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_Wound_myCAF.txt", header = FALSE, sep = NULL)
  Wound_myCAF_list <- list(Wound_myCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = Wound_myCAF_list,
                                      name = "Wound_myCAF_Mod_Score")
  rm(Wound_myCAF, Wound_myCAF_list)

  ### CAF-S1 pancreatic cancer subsets (from Elyada et al.)

  # Elyada Human iCAF (iCAF subset)
  iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_iCAF.txt", header = FALSE, sep = NULL)
  iCAF_list <- list(iCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = iCAF_list,
                                      name = "Elyada_iCAF_Mod_Score")
  rm(iCAF, iCAF_list)

  # Elyada Human myCAF (myCAF subset)
  myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_myCAF.txt", header = FALSE, sep = NULL)
  myCAF_list <- list(myCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = myCAF_list,
                                      name = "Elyada_myCAF_Mod_Score")
  rm(myCAF, myCAF_list)

  # Elyada Mouse apCAF (apCAF subset)
  apCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_apCAF_Mouse.txt", header = FALSE, sep = NULL)
  apCAF_list <- list(apCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = apCAF_list,
                                      name = "Elyada_apCAF_Mouse_Mod_Score")
  rm(apCAF, apCAF_list)

  ### Wu et al. CAF subsets (TN Breast Cancer)

  # Wu human myCAF (myCAF subset)
  myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_myCAF.txt", header = FALSE, sep = NULL)
  myCAF_list <- list(myCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = myCAF_list,
                                      name = "Wu_TN_myCAF_Mod_Score")
  rm(myCAF, myCAF_list)

  # Wu human iCAF (iCAF subset)
  iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_iCAF.txt", header = FALSE, sep = NULL)
  iCAF_list <- list(iCAF$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = iCAF_list,
                                      name = "Wu_TN_iCAF_Mod_Score")
  rm(iCAF, iCAF_list)

  # Wu human dPVL (CAF-S4 subset)
  dPVL <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_dPVL.txt", header = FALSE, sep = NULL)
  dPVL_list <- list(dPVL$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = dPVL_list,
                                      name = "Wu_TN_dPVL_Mod_Score")
  rm(dPVL, dPVL_list)

  # Wu human imPVL (CAF-S4 subset)
  imPVL <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_imPVL.txt", header = FALSE, sep = NULL)
  imPVL_list <- list(imPVL$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = imPVL_list,
                                      name = "Wu_TN_imPVL_Mod_Score")
  rm(imPVL, imPVL_list)

  return(SO)
}
