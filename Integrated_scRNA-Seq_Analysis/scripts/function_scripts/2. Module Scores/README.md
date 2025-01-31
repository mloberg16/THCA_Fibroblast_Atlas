The AddFibroblastModuleScores function takes a Seurat Object as input and returns a Seurat Object with fibroblast module scores added. The Fibroblast module scores are from fibroblast gene sets from published papers (see below). 

This function has three iterations as additional gene sets were added:  
- 24-0821_AddFibroblastModuleScores.R (original function)  
- 24-1125_AddFibroblastModuleScores.R (edited to add gene sets from Cords et al. and Hornburg et al.)  
- 24-1226_AddFibroblasTModuleScores.R (edited to add gene sets from Belle et al., Dinh et al., and Ma et al.)

**Fibroblast papers for public gene sets used to generate module scores:** 

Cords et al.:  
https://www.nature.com/articles/s41467-023-39762-1

Elyada et al.: Pancreatic cancer  
https://aacrjournals.org/cancerdiscovery/article/9/8/1102/42174/Cross-Species-Single-Cell-Analysis-of-Pancreatic

Wu et al.: Triple negative breast cancer  
https://www.embopress.org/doi/full/10.15252/embj.2019104063

Kieffer et al.: Breast cancer  
https://aacrjournals.org/cancerdiscovery/article/10/9/1330/2752/Single-Cell-Analysis-Reveals-Fibroblast-Clusters

Cohen et al.:  
https://www.nature.com/articles/s41467-024-44886-z

Hornburg et al.: Ovarian Cancer  
https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00212-9

Belle et al.: Pancreatic cancer  
https://aacrjournals.org/cancerdiscovery/article/14/7/1324/746125/Senescence-Defines-a-Distinct-Subset-of

Dinh et al.: Esophageal squamous cell carcinoma  
https://www.nature.com/articles/s41467-021-27599-5

Ma et al.: Pan-cancer  
https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-023-01876-x

Iesato et al.:  
https://academic.oup.com/jcem/article/106/12/3569/6327813?login=true
