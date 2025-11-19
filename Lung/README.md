# STExplorer GWR/SA Analysis

GWR: Geographically Weighted Regression
SA: Spatial Autocorrelation

A repository for the analysis of the Lung dataset (from https://doi.org/10.1038/s41588-024-01819-2) with the GWR/SA methods from the STExplorer package. The scripts in this subfolder generate the plots used in main figure 5 and supplementary figure .

# Folder structure
```
Lung
├── scripts
│	├── analysis_1_setUp.R 			 (set up workspace and load required packages, metadata, and data)
│   ├── analysis_2_preprocessing.R 	 (preprocessing, QC, and normalisation)
│   ├── analysis_3_ModuleScores.R 	 (calculate SenMayom, ECM, EMT, TGFB module scores per spot)
│   ├── analysis_4_MainAnalysis.R 	 (main analysis)
│   ├── analysis_5_additionalPlots.R (additonal annotation and other plots)
│   └── analysis_6_reloadImages.R 	 (script to reload images into the SFE objects after saving)
├── data 
│	├── hs_misc
│	│	├── GeneSet_ECM.xlsx			(gene set used for calculating the fibrotic score)
│	│	├── GeneSet_EMT.xlsx			(gene set used for calculating the epith-to-mesench score)
│	│	├── GeneSet_SenMayo_human.xlsx  (gene set used for calculating the senescent score)
│	│	└── GeneSet_TGFB.txt			(gene set used for calculating the responce to TGFB score)
│	└── data_source_publication.pdf
└── README.md
```
