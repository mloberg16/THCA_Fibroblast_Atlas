# Author: Matthew Loberg
# Date: November 25, 2024
# Script: 24-1125_AddFibroblastModuleScores.R

### Function Info
# In this function, I will add module scores from fibroblast gene sets to a given seurat object
# The fibroblast gene sets are from the following papers:
# Cohen et al. 2024 Nature Communications (CAF-S4 gene set) (pericyte-like fibroblasts)
# Kieffer et al. 2020 Cancer Discovery (Normal fibroblast, CAF-S1, and CAF-S1 subset gene sets)

### 24-0821 Update
# This is updated from the 24-0812 version of this script
# I added Wu et al. scores to this version of the script

### 24-1125 Update
# This is updated from the 24-0821 version of this script
# I added Cords et al. and Hornburg et al. data sets to this list

# Function
AddFibroblastScores <- function(SO){

  ### Broad Subsets

  # CAF S4 Broad Subset:
  CAFS4 <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Cohen_etal_CAFS4.txt", header = FALSE, sep = NULL)
  CAFS4_list <- list(CAFS4$V1)
  SO$CAFS4_Mod_Score1 <- NULL # removing in case I've already done this for this object it won't make two versions
  SO <- SO %>% Seurat::AddModuleScore(features = CAFS4_list,
                                      name = "CAFS4_Mod_Score")
  rm(CAFS4, CAFS4_list)

  # CAF S1 Broad Subset:
  CAFS1 <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Kieffer_etal_CAFS1.txt", header = FALSE, sep = NULL)
  CAFS1_list <- list(CAFS1$V1)
  SO$CAFS1_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = CAFS1_list,
                                      name = "CAFS1_Mod_Score")
  rm(CAFS1, CAFS1_list)

  # Normal Fibroblast Broad Subset:
  Normal <- data.table::fread("data_in_use/Fibroblast_Gene_Sets/24-0812_Kieffer_etal_Normal_Fibroblast.txt", header = FALSE, sep = NULL)
  Normal_list <- list(Normal$V1)
  SO$Normal_Fibroblast_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Normal_list,
                                      name = "Normal_Fibroblast_Mod_Score")
  rm(Normal, Normal_list)

  ### CAF S1 Breast Cancer Subsets

  # Kieffer Detox iCAF (CAF-S1 subset)
  Detox_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_Detox_iCAF.txt", header = FALSE, sep = NULL)
  Detox_iCAF_list <- list(Detox_iCAF$V1)
  SO$Detox_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Detox_iCAF_list,
                                      name = "Detox_iCAF_Mod_Score")
  rm(Detox_iCAF, Detox_iCAF_list)

  # Kieffer ecm myCAF (CAF-S1 subset)
  ecm_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_ecm_myCAF.txt", header = FALSE, sep = NULL)
  ecm_myCAF_list <- list(ecm_myCAF$V1)
  SO$ecm_myCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = ecm_myCAF_list,
                                      name = "ecm_myCAF_Mod_Score")
  rm(ecm_myCAF, ecm_myCAF_list)

  # Kieffer IFNg iCAF (CAF-S1 subset)
  IFNg_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_IFNg_iCAF.txt", header = FALSE, sep = NULL)
  IFNg_iCAF_list <- list(IFNg_iCAF$V1)
  SO$IFNg_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = IFNg_iCAF_list,
                                      name = "IFNg_iCAF_Mod_Score")
  rm(IFNg_iCAF, IFNg_iCAF_list)

  # Kieffer IL iCAF (CAF-S1 subset)
  IL_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_IL_iCAF.txt", header = FALSE, sep = NULL)
  IL_iCAF_list <- list(IL_iCAF$V1)
  SO$IL_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = IL_iCAF_list,
                                      name = "IL_iCAF_Mod_Score")
  rm(IL_iCAF, IL_iCAF_list)

  # Kieffer TGFB myCAF (CAF-S1 subset)
  TGFB_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_TGFB_myCAF.txt", header = FALSE, sep = NULL)
  TGFB_myCAF_list <- list(TGFB_myCAF$V1)
  SO$TGFB_myCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = TGFB_myCAF_list,
                                      name = "TGFB_myCAF_Mod_Score")
  rm(TGFB_myCAF, TGFB_myCAF_list)

  # Kieffer TGFB myCAF (CAF-S1 subset)
  Wound_myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0812_Kieffer_etal_Wound_myCAF.txt", header = FALSE, sep = NULL)
  Wound_myCAF_list <- list(Wound_myCAF$V1)
  SO$Wound_myCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Wound_myCAF_list,
                                      name = "Wound_myCAF_Mod_Score")
  rm(Wound_myCAF, Wound_myCAF_list)

  ### CAF-S1 pancreatic cancer subsets

  # Elyada Human iCAF (iCAF subset)
  iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_iCAF.txt", header = FALSE, sep = NULL)
  iCAF_list <- list(iCAF$V1)
  SO$Elyada_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = iCAF_list,
                                      name = "Elyada_iCAF_Mod_Score")
  rm(iCAF, iCAF_list)

  # Elyada Human myCAF (myCAF subset)
  myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_myCAF.txt", header = FALSE, sep = NULL)
  myCAF_list <- list(myCAF$V1)
  SO$Elyada_myCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = myCAF_list,
                                      name = "Elyada_myCAF_Mod_Score")
  rm(myCAF, myCAF_list)

  # Elyada Mouse apCAF (apCAF subset)
  apCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0813_Elyada_etal_apCAF_Mouse.txt", header = FALSE, sep = NULL)
  apCAF_list <- list(apCAF$V1)
  SO$Elyada_apCAF_Mouse_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = apCAF_list,
                                      name = "Elyada_apCAF_Mouse_Mod_Score")
  rm(apCAF, apCAF_list)

  ### Wu et al. CAF subsets (TN Breast Cancer)

  # Wu human myCAF (myCAF subset)
  myCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_myCAF.txt", header = FALSE, sep = NULL)
  myCAF_list <- list(myCAF$V1)
  SO$Wu_TN_myCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = myCAF_list,
                                      name = "Wu_TN_myCAF_Mod_Score")
  rm(myCAF, myCAF_list)

  # Wu human iCAF (iCAF subset)
  iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_iCAF.txt", header = FALSE, sep = NULL)
  iCAF_list <- list(iCAF$V1)
  SO$Wu_TN_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = iCAF_list,
                                      name = "Wu_TN_iCAF_Mod_Score")
  rm(iCAF, iCAF_list)

  # Wu human dPVL (CAF-S4 subset)
  dPVL <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_dPVL.txt", header = FALSE, sep = NULL)
  dPVL_list <- list(dPVL$V1)
  SO$Wu_TN_dPVL_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = dPVL_list,
                                      name = "Wu_TN_dPVL_Mod_Score")
  rm(dPVL, dPVL_list)

  # Wu human imPVL (CAF-S4 subset)
  imPVL <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0816_Wu_etal_imPVL.txt", header = FALSE, sep = NULL)
  imPVL_list <- list(imPVL$V1)
  SO$Wu_TN_imPVL_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = imPVL_list,
                                      name = "Wu_TN_imPVL_Mod_Score")
  rm(imPVL, imPVL_list)

  ### Cords et al. CAF subsets (Breast Cancer? - yes, but PANCANCER)

  # Cords apCAF subset
  Cords_apCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_apCAF.txt", header = FALSE, sep = NULL)
  Cords_apCAF_list <- list(Cords_apCAF$V1)
  SO$Cords_apCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_apCAF_list,
                                      name = "Cords_apCAF_Mod_Score")
  rm(Cords_apCAF, Cords_apCAF_list)

  # Cords dCAF subset
  Cords_dCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_dCAF.txt", header = FALSE, sep = NULL)
  Cords_dCAF_list <- list(Cords_dCAF$V1)
  SO$Cords_dCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_dCAF_list,
                                      name = "Cords_dCAF_Mod_Score")
  rm(Cords_dCAF, Cords_dCAF_list)

  # Cords hsp_tCAF subset
  Cords_hsp_tCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_hsp_tCAF.txt", header = FALSE, sep = NULL)
  Cords_hsp_tCAF_list <- list(Cords_hsp_tCAF$V1)
  SO$Cords_hsp_tCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_hsp_tCAF_list,
                                      name = "Cords_hsp_tCAF_Mod_Score")
  rm(Cords_hsp_tCAF, Cords_hsp_tCAF_list)

  # Cords iCAF subset
  Cords_iCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_iCAF.txt", header = FALSE, sep = NULL)
  Cords_iCAF_list <- list(Cords_iCAF$V1)
  SO$Cords_iCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_iCAF_list,
                                      name = "Cords_iCAF_Mod_Score")
  rm(Cords_iCAF, Cords_iCAF_list)

  # Cords ifnCAF subset
  Cords_ifnCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_ifnCAF.txt", header = FALSE, sep = NULL)
  Cords_ifnCAF_list <- list(Cords_ifnCAF$V1)
  SO$Cords_ifnCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_ifnCAF_list,
                                      name = "Cords_ifnCAF_Mod_Score")
  rm(Cords_ifnCAF, Cords_ifnCAF_list)

  # Cords mCAF subset
  Cords_mCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_mCAF.txt", header = FALSE, sep = NULL)
  Cords_mCAF_list <- list(Cords_mCAF$V1)
  SO$Cords_mCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_mCAF_list,
                                      name = "Cords_mCAF_Mod_Score")
  rm(Cords_mCAF, Cords_mCAF_list)

  # Cords Pericyte subset
  Cords_Pericyte <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_Pericyte.txt", header = FALSE, sep = NULL)
  Cords_Pericyte_list <- list(Cords_Pericyte$V1)
  SO$Cords_Pericyte_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_Pericyte_list,
                                      name = "Cords_Pericyte_Mod_Score")
  rm(Cords_Pericyte, Cords_Pericyte_list)

  # Cords rCAF subset
  Cords_rCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_rCAF.txt", header = FALSE, sep = NULL)
  Cords_rCAF_list <- list(Cords_rCAF$V1)
  SO$Cords_rCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_rCAF_list,
                                      name = "Cords_rCAF_Mod_Score")
  rm(Cords_rCAF, Cords_rCAF_list)

  # Cords tCAF subset
  Cords_tCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_tCAF.txt", header = FALSE, sep = NULL)
  Cords_tCAF_list <- list(Cords_tCAF$V1)
  SO$Cords_tCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_tCAF_list,
                                      name = "Cords_tCAF_Mod_Score")
  rm(Cords_tCAF, Cords_tCAF_list)

  # Cords vCAF subset
  Cords_vCAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-0919_Cords_etal_vCAF.txt", header = FALSE, sep = NULL)
  Cords_vCAF_list <- list(Cords_vCAF$V1)
  SO$Cords_vCAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Cords_vCAF_list,
                                      name = "Cords_vCAF_Mod_Score")
  rm(Cords_vCAF, Cords_vCAF_list)

  ### Hornburg et al. CAF subsets (Ovarian Cancer)

  # Hornburg IL1_CAF subset
  Hornburg_IL1_CAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-1125_Hornburg_etal_IL1_CAF.txt", header = FALSE, sep = NULL)
  Hornburg_IL1_CAF_list <- list(Hornburg_IL1_CAF$V1)
  SO$Hornburg_IL1_CAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Hornburg_IL1_CAF_list,
                                      name = "Hornburg_IL1_CAF_Mod_Score")
  rm(Hornburg_IL1_CAF, Hornburg_IL1_CAF_list)

  # Hornburg TGFB_CAF subset
  Hornburg_TGFB_CAF <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-1125_Hornburg_etal_TGFB_CAF.txt", header = FALSE, sep = NULL)
  Hornburg_TGFB_CAF_list <- list(Hornburg_TGFB_CAF$V1)
  SO$Hornburg_TGFB_CAF_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Hornburg_TGFB_CAF_list,
                                      name = "Hornburg_TGFB_CAF_Mod_Score")
  rm(Hornburg_TGFB_CAF, Hornburg_TGFB_CAF_list)

  # Hornburg Proliferative subset
  Hornburg_Proliferative <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-1125_Hornburg_etal_Proliferative.txt", header = FALSE, sep = NULL)
  Hornburg_Proliferative_list <- list(Hornburg_Proliferative$V1)
  SO$Hornburg_Proliferative_Mod_Score1 <- NULL
  SO <- SO %>% Seurat::AddModuleScore(features = Hornburg_Proliferative_list,
                                      name = "Hornburg_Proliferative_Mod_Score")
  rm(Hornburg_Proliferative, Hornburg_Proliferative_list)

  ### Iesato et al. CAF subsets (Fibroblast and Pericyte - Thyroid paper from Carmelo Nucera)
  # https://academic-oup-com.proxy.library.vanderbilt.edu/jcem/article/106/12/3569/6327813

  # Iesato Pericyte subset
  Iesato_Pericyte <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-1125_Iesato_etal_Pericyte.txt", header = FALSE, sep = NULL)
  Iesato_Pericyte_list <- list(Iesato_Pericyte$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = Iesato_Pericyte_list,
                                      name = "Iesato_Pericyte_Mod_Score")
  rm(Iesato_Pericyte, Iesato_Pericyte_list)

  # Iesato Fibroblast subset
  Iesato_Fibroblast <- data.table::fread("data_in_use/Fibroblast_gene_sets/24-1125_Iesato_etal_Fibroblast.txt", header = FALSE, sep = NULL)
  Iesato_Fibroblast_list <- list(Iesato_Fibroblast$V1)
  SO <- SO %>% Seurat::AddModuleScore(features = Iesato_Fibroblast_list,
                                      name = "Iesato_Fibroblast_Mod_Score")
  rm(Iesato_Fibroblast, Iesato_Fibroblast_list)

  return(SO)
}
