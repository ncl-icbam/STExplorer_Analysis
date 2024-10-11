library(STExplorer)
library(readr)
library(ggplot2)
library(svglite)
library(magick)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(AnnotationHub)

#library(SpatialFeatureExperiment)

setwd("P:/STExplorer/Scripts/Code_test")

## Data import and preparation of the MSFE object
msfe <- MetaSpatialFeatureExperiment()

sampleDir <- c("P:/STExplorer/Prostate/test_data/Visium_Human_Prostate/Patient1/H1_4",
               "P:/STExplorer/Prostate/test_data/Visium_Human_Prostate/Patient1/H2_5")
               
sampleNames <- c("H1_4","H2_5")
names(sampleDir) <- sampleNames
  
for (i in seq_along(sampleNames)) {
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i]))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures")))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures/GWPCA")))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures/GSEA")))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures/FGC")))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures/FGC/Automated")))
  dir.create(paste0("P:/STExplorer/Scripts/Code_test/",sampleNames[i],paste0("/Figures/Spatial_expression")))
}
               
for (i in seq_along(sampleNames)) {
  
  msfe <- addSFE(msfe,
                 read10xVisiumSFE(samples = sampleDir[i], 
                                  sample_id = sampleNames[i], 
                                  type = "HDF5", 
                                  data = "filtered", 
                                  images = "lowres", 
                                  style = "W", 
                                  zero.policy = TRUE))
}

ground_truth <- read_table("P:/STExplorer/Prostate/test_data/Visium_Human_Prostate/Patient1/spotZonationGroup.txt")
gTruth_list <- list(H1_4 = ground_truth[ground_truth$sample_id == "H1_4",],
                    H2_5 = ground_truth[ground_truth$sample_id == "H2_5",])
                    
## Spot-level Quality Control

## Mark a subset of mitochondrial genes
is_mito <- getSubset(msfe,
                      sample_id = TRUE,
                      subset = "(^MT-)|(^mt-)",
                      set = "rowData",
                      col_name = "symbol")

for (i in seq_along(sampleNames)) {
    message("Working on sample: ", sampleNames[i])
    ## Add location-related statistics
    msfe <- addPerLocQC(msfe,
                        sample_id = sampleNames[i],
                        gTruth = gTruth_list[[i]],
                        assay = "counts",
                        MARGIN = 2,
                        subsets = list(mito = is_mito[[i]]))
    message("\tAdded location-related statistics")
    
    ## Add geometries
    msfe <- addGeometries(msfe,
                          samples = sampleDir[i],
                          sample_id = sampleNames[i],
                          res = "fullres",
                          flipped=FALSE)
    message("\tAdded geometries")
    
    ## Add gene/feature-related statistics
    msfe <- addPerGeneQC(msfe,
                         sample_id = sampleNames[i],
                         assay = "counts",
                         version = NULL,
                         mirror = NULL,
                         add=c("zeroexpr", "exprstats"))
    message("\tAdded gene/feature-related statistics")
  }

  ## Set QC thresholds

for (i in seq_along(sampleNames)) {
  msfe@sfe_data[[i]] <- setQCthresh_LibSize(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.15)), max_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.99)))
  msfe@sfe_data[[i]] <- setQCthresh_Mito(msfe@sfe_data[[i]], sample_id = TRUE, min_t = NA, max_t = quantile(msfe@sfe_data[[i]]@colData$subsets_mito_percent, probs = c(.99), na.rm = TRUE))
  msfe@sfe_data[[i]] <- setQCthresh_GenesExpr(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.15)), max_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.99)))
  
  ## Set the combined filtering threshold using the QC metrics
  msfe@sfe_data[[i]] <- setQCtoDiscard_loc(msfe@sfe_data[[i]], sample_id = TRUE, filters = TRUE)}
  
## Remove combined set of low-quality spots
msfe <- applyQCthresh_loc(msfe, sample_id = TRUE)

## Normalisation of counts

### Compute library size factors
msfe <- computeLibSizeFactors(msfe)

### Log-tranformation of counts
msfe <- normaliseCounts(msfe)

## Gene-level Quality Control

## Calculate the mean of log counts over the number of locations a gene is present
msfe <- perGeneLogMean(msfe,sample_id = TRUE)

## Zero expression genes
msfe <- setQCthresh_ZeroExpr(msfe,sample_id = TRUE)

## Lowly expressed (noise?!) genes
msfe <- setQCthresh_LowLogMean(msfe,threshold = 1,sample_id = TRUE)

## Remove mitochondrial and other genes

for (i in seq_along(sampleNames)) {
msfe@sfe_data[[i]] <- setQCthresh_custom(msfe@sfe_data[[i]], MARGIN = 1, qcMetric = is_mito[[i]])}

## QC discard Features
## Set the combined filtering threshold using the QC metrics
msfe <- setQCtoDiscard_feat(msfe, filters = TRUE,sample_id = TRUE)

## FEATURE SELECTION
## Apply gene-level QC threshold
msfe <- applyQCthresh_feat(msfe,sample_id = TRUE)

## Fit mean-variance relationship
dec <- modelGeneVariance(msfe, sample_id = TRUE, method = "Var")
## Select top HVGs

top_10_hvgs <- getTopHighVarGenes(dec,
                                  var.field = "bio",
                                  prop = 0.1,
                                  var.threshold = 0,
                                  fdr.threshold = 0.1)

### Spatially variable genes (SVGs)

## Add a neighbour graph using a weighted distance matrix
msfe <- addSpatialNeighGraphs(msfe, sample_id = TRUE, type = "knearneigh", style = "W", distMod = "raw", k = 6)

## Calculate a simple distance matrix
msfe <- addDistMat(msfe, p = 2)

## Geographically Weighted Principal Components Analysis (GWPCA)

## Set the number of components to be retained
k = 20
## Set the kernel to be used
kernel = "gaussian"
## Set the Minkowski distance power: p = 2 --> Euclidean
p = 2
## Is the bandwidth adaptive?: No because spots are fixed
adaptive = FALSE
## Cross-Validate GWPCA?
cv = TRUE
## Calculate PCA scores?
scores = FALSE
## Run a robust GWPCA?
robust = FALSE

  #Setting the GWPCA parameters
 
  bw_6 <- 6*msfe@sfe_data[["H1_4"]]@metadata[["spotDiameter"]][["H1_4"]][["spot_diameter_fullres"]]
  vars_10_hvgs <- top_10_hvgs[["H1_4"]]
 
  #Perform GWPCA for the top 10 HVGs
  my.cl <- makeClusterGWPCA(spec=7,type = "PSOCK") #(for Windows with parallel computing)
  pcagw_top10_output <- gwpcaSTE(sfe = msfe@sfe_data[["H1_4"]],
                                 assay = "logcounts",
                                 vars = vars_10_hvgs,
                                 p = p,
                                 k = k,
                                 bw = bw_6,
                                 kernel = kernel,
                                 adaptive = adaptive,
                                 scores = scores,
                                 robust = robust,
                                 cv = cv,
                                 future = FALSE,
                                 strategy = "cluster",
                                 workers = my.cl,
                                 verbose = FALSE)
  
  #Extract the GWPCA leading genes for the top 10 HVGs
  pcagw_top10_LeadingGene_single <- gwpca_LeadingGene(gwpca = pcagw_top10_output,
                                                             m_sfe = msfe@sfe_data[["H1_4"]],
                                                             pc_nos = 1:4,
                                                             sample_id="H1_4",
                                                             type = "single",
                                                             names = "gene_names")
  plotGWPCA_leadingG(gwpca = pcagw_top10_LeadingGene_single,
                     comps = 1,
                     type = "single",
                     arrange = FALSE)
  
  ggsave("P:/STExplorer/Scripts/Code_test/H1_4/Figures/GWPCA/H1_4_top_10_HVG_GWPCA_leading_single.svg",
         device = "svg",
         width = 15,
         height = 8,
         units = "in",
         dpi = 300)
  
  ## Calculate the GWPCA PTV for multiple Components
  pcagw_top10_output <- gwpca_PropVar(gwpca = pcagw_top10_output, n_comp = 2:10, m_sfe = msfe@sfe_data[["H1_4"]])

  ## Map GWPCA PTV 
  plotGWPCA_ptv(gwpca = pcagw_top10_output,
                comps = 1:6,
                type = "map")
  
  ggsave("P:/STExplorer/Scripts/Code_test/H1_4/Figures/GWPCA/H1_4_top_10_HVG_GWPCA_PTV_map.svg",
         device = "svg",
         width = 15,
         height = 8,
         units = "in",
         dpi = 300)

### GSEA Functional BP clustering 
msigdb <- getMSigDBData("Homo sapiens")
t2g_BP <- getTerm2Gene(msig_data = msigdb, cat = "C5", subcat = "GO:BP")

gsea_map_BP <- gwpca_FunctionalClustering(gwpca = pcagw_top10_output,
                                                   pc = 1,
                                                   genes_no = 2,
                                                   NES = 1.5,
                                                   minGSSize = 5,
                                                   pvalueCutoff = 0.25,
                                                   TERM2GENE = t2g_BP,
                                                   pAdjustMethod = "fdr",
                                                   scoreType = "std",
                                                   nPermSimple = 10000,
                                                   mc.cores = 10,type = "PSOCK")
  
gsea_plot_BP <- plotGWPCA_FuncCLust(gsea_map_BP, count = 8, legend = "right",legend.title="GO_BP_GSEA")

ggsave("P:/STExplorer/Scripts/Code_test/H1_4/Figures/GSEA/H1_4_top_10_HVG_GWPCA_gsea_map_BP.svg", gsea_plot_BP,
       device = "svg",
       width = 15,
       height = 8,
       units = "in",
       dpi = 300)

  ## Fuzzy Geographically Weighted Clustering (FGWC)
  
  ## Find optimum number of Factors
best_k_nmf <- fgwc_nmfFactorNumber(m_sfe = msfe,
                                 sample_id = "H2_5",
                                 assay = "logcounts",
                                 top_hvgs = top_10_hvgs[["H2_5"]],
                                 k_range = seq(2, 5, 1),
                                 n_cores = 1,
                                 do_plot = FALSE,
                                 seed = 1,
                                 loss = "mse",
                                 max.iter = 250)

  ## Run NMF
sfe_nmf <- fgwc_nmf(m_sfe = msfe,
                                sample_id = "H2_5",
                                top_hvgs = top_10_hvgs[["H2_5"]],
                                ncomponents = best_k_nmf[["k"]])

  ## Set number of FGWC clusters
 
  fgwc_param <- fgwc_params(algorithm = "classic",
                                      ncluster = 4)

# Run FGWC ----

  fgwc_list <- fgwcSTE(m_sfe = msfe,
                            sample_id = "H2_5",
                            data = sfe_nmf,
                            dMetric = "euclidean",
                            algorithm = "classic",
                            parameters = fgwc_param)

# Plot FGWC clusters ----

  ## Plot single FGWC clusters
  p1 <- plotFGWC_singleMap(fgwc = fgwc_list,
                           m_sfe = msfe,
                           sample_id = "H2_5")
  
  p2 <- plotQC_spotsAnnotation(msfe,
                               sample_id = "H2_5",
                               type = "hex")
  
  print(p1 + p2)
  
  ggplot2::ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_single_NMF.svg",
         device = "svg",
         width = 15,
         height = 8,
         units = "in",
         dpi = 300)
  
# Plot multi FGWC clusters ----

  ## Plot multi clusters
  p2 <- plotFGWC_multiMap(fgwc = fgwc_list,
                          m_sfe = msfe,
                          sample_id = "H2_5")
  
  print(p2)
  ggplot2::ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_multi_NMF.svg",
                  device = "svg",
                  width = 15,
                  height = 8,
                  units = "in",
                  dpi = 300)

# Plot FGWC pie-doughnuts ----

  
  ## Matching FGWC clusters and annotation of locations
  plotFGWC_pie(fgwc = fgwc_list,
               m_sfe = msfe,
               sample_id = "H2_5",
               mapping = aes(pie = cluster, donut = annotation),donutLabelSize = 10,pieLabelSize = 9)
  ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_pieClustAnn_NMF.svg",
                  device = "svg",
                  width = 15,
                  height = 8,
                  units = "in",
                  dpi = 300)
  
  ## Matching FGWC clusters and NMF factors
  plotFGWC_pie(fgwc = fgwc_list,
               m_sfe = msfe,
               sample_id = "H2_5",
               mapping = aes(pie = cluster, donut = factors),donutLabelSize = 10,pieLabelSize = 9)
  ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_pieClustFact_NMF.svg",
                  device = "svg",
                  width = 15,
                  height = 8,
                  units = "in",
                  dpi = 300)

  ## Plot a factor heatmap alongside spot annotations and/or FGWC clusters
  
  p3 <- plotFGWC_nmfFactorsHeatmap(fgwc = fgwc_list,
                             loc_annot = "annotation",
                             order_rows = "annotation")
  ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_factorsHeatMap_NMF_annotation.svg",p3,
         device = "svg",
         width = 15,
         height = 8,
         units = "in",
         dpi = 300)
  
  p4 <- plotFGWC_nmfFactorsHeatmap(fgwc = fgwc_list,
                             loc_annot = "cluster",
                             order_rows = "cluster")
  ggsave("P:/STExplorer/Scripts/Code_test/H2_5/Figures/FGC/Automated/H2_5_factorsHeatMap_NMF_cluster.svg",p4,
         device = "svg",
         width = 15,
         height = 8,
         units = "in",
         dpi = 300)
  
  # Plot the expression of selected genes on tissue slices
  
  library(terra)
  
  plotGeneExpression_list <- list()
  
  selected = as.list(c("ENSG00000125144","ENSG00000116133","ENSG00000118523","ENSG00000183036","ENSG00000164687","ENSG00000130529","ENSG00000148346","ENSG00000081041"))
  names(selected) = c("ENSG00000125144","ENSG00000116133","ENSG00000118523","ENSG00000183036","ENSG00000164687","ENSG00000130529","ENSG00000148346","ENSG00000081041")
  
  for (k in selected){
    plotGeneExpression_list[["H2_5"]][[k]] <- plotGeneExpression(m_sfe = msfe@sfe_data[["H2_5"]],
                                                                 genes = selected[[k]],
                                                                 sample_id = "H2_5",
                                                                 assay = "logcounts",
                                                                 minmax = c(2, 10),
                                                                 type = "hex",
                                                                 res = "lowres",
                                                                 fill_args = list(option = "viridis",na.value = "grey"
                                                                 ),alpha = 0.3)
    
    ggsave(paste0("P:/STExplorer/Scripts/Code_test/H2_5/Figures/Spatial_expression/",k,paste0("_plotGeneExpression.svg")),plotGeneExpression_list[["H2_5"]][[k]],
           device = "svg",
           width = 15,
           height = 8,
           units = "in",
           dpi = 300)}
  
  save.image(file='STExplorer_code.RData')
  