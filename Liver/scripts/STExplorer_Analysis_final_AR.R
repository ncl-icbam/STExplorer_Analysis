knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", eval = FALSE
)

## To load the packages use:
library(STExplorer)
library(readr)
library(SpatialFeatureExperiment)
library(tidyverse)
library(scran)
library(scater)
library(sf)
library(spdep)
library(GWmodel)
library(tidyterra)
library(ggplot2)
library(igraph)
library(pheatmap)
library(ggExtra)
library(future)
library(doFuture)
library(foreach)
library(progressr)
library(parallel)
library(cols4all)
library(pheatmap)
library(RColorBrewer)
library(readxl)


## Set working directory to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Create the MSFE object
msfe <- MetaSpatialFeatureExperiment()

## Prepare vectors with the paths to the data folders
sampleDir <- c("./data/Human_Liver_Healthy_JBO018_Results","./data/Human_Liver_Steatotic_JBO019_Results")

sampleNames <- c("JBO018", "JBO019")

names(sampleDir) <- sampleNames

## Load sfe objects inside the msfe
for (i in seq_along(sampleNames)) {
  message("Adding sample: ", sampleNames[i])
  msfe <- addSFE(msfe,
                 read10xVisiumSFE(samples = sampleDir[i],
                                  sample_id = sampleNames[i],
                                  type = "HDF5",
                                  data = "filtered",
                                  images = "lowres",
                                  style = "W",
                                  zero.policy = TRUE))
}


ground_truth <- read.table("./data/metadata/annot_humanVisium.csv", sep =",", header = 1)

ground_truth$sample_id <- gsub("JBO", "JBO0", ground_truth$sample)
ground_truth$Barcode = gsub("\\_.*","",ground_truth$spot )
ground_truth$annotation = ground_truth$zonationGroup

gTruth_list <- list(JBO018 = ground_truth[ground_truth$sample_id == "JBO018",],
                    JBO019 = ground_truth[ground_truth$sample_id == "JBO019",])

str(gTruth_list)

## Dataset dimensions before the filtering
for (id in sampleNames) {
  message("Sample: ", id)
  print(dim(msfe@sfe_data[[id]]))
}

## Mark a subset of mitochondrial genes
is_mito <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "(^MT-)|(^mt-)",
                     set = "rowData",
                     col_name = "symbol")

for (id in sampleNames) {
  message("Working on sample: ", id)
  ## Add location-related statistics
  msfe <- addPerLocQC(msfe,
                      sample_id = id,
                      gTruth = gTruth_list[[id]],
                      assay = "counts",
                      MARGIN = 2,
                      subsets = list(mito = is_mito[[id]]))
  message("\tAdded location-related statistics")

  ## Add geometries
  msfe <- addGeometries(msfe,
                        samples = sampleDir[id],
                        sample_id = id,
                        res = "fullres",
                        flipped = TRUE)
  message("\tAdded geometries")

  ## Add gene/feature-related statistics
  msfe <- addPerGeneQC(msfe,
                       sample_id = id,
                       assay = "counts",
                       version = NULL,
                       mirror = NULL,
                       add = c("zeroexpr", "exprstats"))
  message("\tAdded gene/feature-related statistics")
}

#### FILTER JBO018 ####
sfe <- getSFE(msfe, "JBO018")

## Plot spatial coordinates with annotations
plotQC_spotsAnnotation(sfe, type = "spot", sample_id = NULL)
plotQC_spotsAnnotation(sfe, type = "hex", sample_id = "JBO018")

plotQC_tissueImg(sfe, res = "lowres", type = "spot", sample_id = NULL, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "hex", sample_id = "JBO018", annotate = TRUE, alpha = 0.3)

## Density and histogram of library sizes
plotQC_hist(sfe, metric = "libsize")

## Map the library sizes
plotQC_map(sfe, metric = "libsize")

## Select threshold
sfe <- setQCthresh_LibSize(sfe, sample_id = TRUE,
                           min_t = quantile(sfe@colData$sum, probs = c(.05)),
                           max_t = quantile(sfe@colData$sum, probs = c(1)))

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "libsize", sample_id = TRUE)

## Density and histogram of expressed genes
plotQC_hist(sfe, metric = "detected")

## Map the library sizes
plotQC_map(sfe, metric = "detected")

## Select threshold
sfe <- setQCthresh_GenesExpr(sfe, sample_id = TRUE,
                             min_t = quantile(sfe@colData$detected, probs = c(.05)),
                             max_t = quantile(sfe@colData$detected, probs = c(1)))

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "detected", sample_id = TRUE)

## Density and histogram of percentage of mitochondrial expression
plotQC_hist(sfe, metric = "mito")
plotQC_hist(sfe, metric = "mito", limits = c(NA, 15))

## Map the library sizes
plotQC_map(sfe, metric = "mito")

## Select threshold
sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, min_t = NA, max_t = 15)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "mito", sample_id = TRUE)

## Select locations without annotation
sfe <- setQCthresh_NAs(sfe, sample_id = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "NAs")

## Set the combined filtering threshold using the QC metrics
sfe <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "discard", sample_id = TRUE)

## Remove combined set of low-quality spots
sfe <- applyQCthresh_loc(sfe, sample_id = TRUE)

# Show annotation of spots after filtering
plotQC_spotsAnnotation(sfe, type = "hex")

dim(sfe)

msfe <- addSFE(msfe, sfe)

#### FILTER JBO019 ####
sfe <- getSFE(msfe, "JBO019")

## Plot spatial coordinates with annotations
plotQC_spotsAnnotation(sfe, type = "spot", sample_id = NULL)
plotQC_spotsAnnotation(sfe, type = "hex", sample_id = "JBO019")

plotQC_tissueImg(sfe, res = "lowres", type = "spot", sample_id = NULL, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "hex", sample_id = "JBO019", annotate = TRUE, alpha = 0.3)

## Density and histogram of library sizes
plotQC_hist(sfe, metric = "libsize")

## Map the library sizes
plotQC_map(sfe, metric = "libsize")

## Select threshold
sfe <- setQCthresh_LibSize(sfe, sample_id = TRUE,
                           min_t = quantile(sfe@colData$sum, probs = c(.05)),
                           max_t = quantile(sfe@colData$sum, probs = c(1)))

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "libsize", sample_id = TRUE)

## Density and histogram of expressed genes
plotQC_hist(sfe, metric = "detected")

## Map the library sizes
plotQC_map(sfe, metric = "detected")

## Select threshold
sfe <- setQCthresh_GenesExpr(sfe, sample_id = TRUE,
                             min_t = quantile(sfe@colData$detected, probs = c(.05)),
                             max_t = quantile(sfe@colData$detected, probs = c(1)))

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "detected", sample_id = TRUE)

## Density and histogram of percentage of mitochondrial expression
plotQC_hist(sfe, metric = "mito")
plotQC_hist(sfe, metric = "mito", limits = c(NA, 15))

## Map the library sizes
plotQC_map(sfe, metric = "mito")

## Select threshold
sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, min_t = NA, max_t = 15)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "mito", sample_id = TRUE)

## Select locations without annotation
sfe <- setQCthresh_NAs(sfe, sample_id = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "NAs")

## Set the combined filtering threshold using the QC metrics
sfe <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "discard", sample_id = TRUE)

## Remove combined set of low-quality spots
sfe <- applyQCthresh_loc(sfe, sample_id = TRUE)

# Show annotation of spots after filtering
plotQC_spotsAnnotation(sfe, type = "hex")

dim(sfe)

msfe <- addSFE(msfe, sfe)

#### WORK WITH FILTERED MSFE OBJECT ####
## Calculate library size factors
msfe <- computeLibSizeFactors(msfe)

## Calculate logcounts using library size factors
msfe <- normaliseCounts(msfe)

## Calculate the mean of log counts over the number of locations a gene is present
msfe <- perGeneLogMean(msfe)

# ## Zero expression genes
msfe <- setQCthresh_ZeroExpr(msfe)

## Fit mean-variance relationship
dec <- modelGeneVariance(msfe, sample_id = TRUE, method = "Var")

## Select top HVGs
top_hvgs <- getTopHighVarGenes(dec,
                               var.field = "bio",
                               prop = 0.5,
                               var.threshold = 0,
                               fdr.threshold = 0.05)

## Examine overlap in HVGs
library(ggVennDiagram)
x <- list(JBO019=top_hvgs$JBO019,
          JBO018=top_hvgs$JBO018)

library(ggplot2)
ggVennDiagram(x,label = "count") + scale_fill_gradient(low="grey90",high = "red")

## Visualize mean-variance relationship
plotGeneVariance(dec = dec, hvgs = top_hvgs)

## Add a neighbour graph using a weighted distance matrix
msfe <- addSpatialNeighGraphs(msfe, sample_id = TRUE, type = "knearneigh", style = "W", distMod = "raw", k = 6)

## Calculate a simple distance matrix
msfe <- addDistMat(msfe, p = 2)

## Plot the neighbours graph
plotNeighbourGraph(msfe, sample_id = TRUE,
                   res = "lowres", plotImage = TRUE)

# Prepare lists ----
samples <- names(msfe@sfe_data)

best_k_nmf <- list()
sfe_nmf_list <- list()
best_k_fgwc <- list()
fgwc_param_list <- list()
fgwc_list <- list()

marker_heatmap_list <- list()
subHeatmap_list <- list()

# Find optimum factor number ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Find optimum number of Factors
  result <- fgwc_nmfFactorNumber(m_sfe = msfe,
                                 sample_id = s,
                                 assay = "logcounts",
                                 top_hvgs = top_hvgs[[s]],
                                 k_range = seq(2, 10, 1),
                                 n_cores = 1,
                                 do_plot = FALSE,
                                 seed = 1,
                                 loss = "mse",
                                 max.iter = 250)

  best_k_nmf[[s]] <- result

  print(plotFGWC_factorSelection(result))

  ## Housekeeping
  rm(result)
}

# Calculate NMF ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Run NMF
  sfe_nmf_list[[s]] <- fgwc_nmf(m_sfe = msfe,
                                sample_id = s,
                                top_hvgs = top_hvgs[[s]],
                                ncomponents = 3)
}

for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Find best number of clusters
  fgwc_param_list[[s]] <- fgwc_params(algorithm = "classic", ncluster = 5)

  best_k_fgwc[[s]] <- fgwc_findOptimumK(fgwc_in = sfe_nmf_list[[s]],
                                        k_range = 2:10,
                                        index_type = "FPC",
                                        elbow_method = "knee",
                                        m_sfe = msfe,
                                        sample_id = s,
                                        algorithm = "classic",
                                        parameters = fgwc_param_list[[s]])

  ## update the parameters input
  fgwc_param_list[[s]] <- fgwc_params(algorithm = "classic",
                                      ncluster = best_k_fgwc[[s]])
}


# Run FGWC ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Run FGWC
  fgwc_list[[s]] <- fgwcSTE(m_sfe = msfe,
                            sample_id = s,
                            data = sfe_nmf_list[[s]],
                            dMetric = "euclidean",
                            algorithm = "classic",
                            parameters = fgwc_param_list[[s]])
}

# Reorder clusters for visualisation
JBO018_names = factor(fgwc_list$JBO018$cluster)
JBO018_names = factor(JBO018_names,
                      levels = c(3,4,5,1,2),
                      labels = c(1,2,3,4,5))
#fgwc_list$JBO018$cluster = JBO018_names
fgwc_list$JBO018$finaldata$cluster = fgwc_list$JBO018$cluster


JBO019_names = factor(fgwc_list$JBO019$cluster)
JBO019_names = factor(JBO019_names,
                      levels = c(2,3,5,4,1),
                      labels = c(1,2,3,4,5))
#fgwc_list$JBO019$cluster = JBO019_names
fgwc_list$JBO019$finaldata$cluster = fgwc_list$JBO019$cluster

# Plot single clusters ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Plot single clusters
  p1 <- plotFGWC_singleMap(fgwc = fgwc_list[[s]],
                           m_sfe = msfe,
                           sample_id = s)

  p2 <- plotQC_spotsAnnotation(msfe,
                               sample_id = s,
                               type = "hex")

  print(p1 + p2)
}

# Plot clusters and annotation
p1 = plotQC_spotsAnnotation(msfe@sfe_data$JBO018, type = "hex",colours = c("#6fafd6","#FFFFBF", "#FDAE61", "#D7191C"))
p2 = plotFGWC_singleMap(fgwc = fgwc_list[[1]],m_sfe = msfe, sample_id = sampleNames[1], colours = c( "#6fafd6","#ABDDA4", "#FFFFBF","#FDAE61","#D7191C"))
p1 + p2

p1 = plotQC_spotsAnnotation(msfe@sfe_data$JBO019, type = "hex",colours = c("#6fafd6","#FFFFBF", "#FDAE61", "#D7191C"))
p2 = plotFGWC_singleMap(fgwc = fgwc_list[[2]],m_sfe = msfe, sample_id = sampleNames[2], colours = c( "#6fafd6","#ABDDA4", "#FFFFBF","#FDAE61","#D7191C"))
p1 + p2

# Plot multi clusters ----
for (s in samples) {
  message("# ---------------------- #\n",
          "Working on sample: ", s)
  ## Plot multi clusters
  p1 <- plotFGWC_multiMap(fgwc = fgwc_list[[s]],
                          m_sfe = msfe,
                          sample_id = s)

  print(p1)
}


# Plot memberships in violins ----
for (s in samples) {
  ## Plot memberships in violins
  p1 <- plotFGWC_multiViolin(fgwc = fgwc_list[[s]],
                             m_sfe = msfe,
                             sample_id = s)

  print(p1)
}


# Plot FGWC pie-doughnuts ----
for (s in samples) {
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m

  ## Matching FGWC clusters and annotation of locations
  plotFGWC_pie(fgwc = fgwc_list[[s]],
               m_sfe = msfe,
               sample_id = s,
               mapping = aes(pie = cluster, donut = annotation),
               pieColours =  c( "#6fafd6","#ABDDA4", "#FFFFBF","#FDAE61","#D7191C"))

  ## Matching FGWC clusters and NMF factors
  plotFGWC_pie(fgwc = fgwc_list[[s]],
               m_sfe = msfe,
               sample_id = s,
               mapping = aes(pie = cluster, donut = factors),
               pieColours =  c( "#6fafd6","#ABDDA4", "#FFFFBF","#FDAE61","#D7191C"))
}


# Plot map of factor scores ----
for (s in samples) {
  ## Plot a map of the NMF factor scores for each location
  p1 <- plotFGWC_nmfFactorsMap(nmf = sfe_nmf_list[[s]],
                               m_sfe = msfe,
                               sample_id = s)

  print(p1)

  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m

}

# Plot a factor heatmap alongside spot annotations and/or FGWC clusters
for (s in samples) {
  plotFGWC_nmfFactorsHeatmap(fgwc = fgwc_list[[s]],
                             loc_annot = "both",
                             order_rows = "cluster")
  ncomp = ncol(fgwc_list[[s]]$membership)
  ncluster = fgwc_list[[s]]$call$ncluster
  a = fgwc_list[[s]]$call$a
  m = fgwc_list[[s]]$call$m
}

# GWPCA analysis ----
gsea_map_list <- list(JBO018 = list(),
                      JBO019 = list())
## Select the sample you would like to perform a GWPCA analysis
sfe <- getSFE(msfe, "JBO018")
## Get the gene names that are going to be evaluated
vars = top_hvgs[["JBO018"]]
## Set a fixed bandwidth
## bw is an important parameter as it defines the neighbourhood for which the
##  PCA will be calculated. The distance is measured in ultra-high resolution
##  image pixels. The default is 3x the diameter of the Visium spot. Make sure
##  to adjust it if it is too large or too small for your setting.
bw = 3*sfe@metadata[["spotDiameter"]][["JBO018"]][["spot_diameter_fullres"]]
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
## Make a cluster for parallel computing (otherwise GWPCA is slow!)
my.cl <- makeClusterGWPCA(type = "FORK")

## Run GWPCA
pcagw <- gwpcaSTE(sfe = sfe,
                  assay = "logcounts",
                  vars = vars,
                  p = p,
                  k = k,
                  bw = bw,
                  kernel = kernel,
                  adaptive = adaptive,
                  scores = scores,
                  robust = robust,
                  cv = cv,
                  future = FALSE,
                  strategy = "cluster",
                  workers = my.cl,
                  verbose = FALSE)


plotGWPCA_global(gwpca = pcagw,
                 comps = 1:10,
                 type = "scree",
                 point_args = list(size = 3, colour = "red"),
                 line_args = list(linewidth = 1, colour = "dodgerblue"))

## Extract leading genes
pcagw <- gwpca_LeadingGene(gwpca = pcagw,
                           m_sfe = sfe,
                           pc_nos = 1:4,
                           type = "single",
                           names = "gene_names")

pcagw <- gwpca_LeadingGene(gwpca = pcagw,
                           m_sfe = sfe,
                           pc_nos = 1:4,
                           genes_n = 4,
                           type = "multi",
                           method = "membership",
                           names = "gene_names")

## Plot leading genes
plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1:2,
                   type = "single",
                   arrange = FALSE)

plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1,
                   type = "multi",
                   arrange = FALSE)

### Plot multi type (extra parameters)
plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1,
                   type = "multi",
                   arrange = FALSE,
                   legend.position = "bottom",
                   cutoff = 12,
                   size = 8)

### Functional clustering
## Steatosis - no enrichment in this sample
t2g = read.csv("./data/t2g_files/steatosis_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO018"]][["steatosis"]] <- gsea_map

## Fibrosis
t2g = read.csv("./data/t2g_files/fibrosis_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO018"]][["fibrosis"]] <- gsea_map

##  NAS > 4 - not corrected for F stage
t2g = read.csv("./data/t2g_files/NAS_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO018"]][["nas"]] <- gsea_map

## Metabolism of lipids
msigdb <- getMSigDBData("Homo sapiens")
reactome_t2g <- getTerm2Gene(msig_data = msigdb, cat = "C2", subcat = "CP:REACTOME")
t2g = reactome_t2g[reactome_t2g$term == "REACTOME_METABOLISM_OF_LIPIDS",]

gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO018"]][["lipid_metab"]] <- gsea_map

## Select the sample you would like to perform a GWPCA analysis
sfe <- getSFE(msfe, "JBO019")
## Get the gene names that are going to be evaluated
vars = top_hvgs[["JBO019"]]
## Set a fixed bandwidth
## bw is an important parameter as it defines the neighbourhood for which the
##  PCA will be calculated. The distance is measured in ultra-high resolution
##  image pixels. The default is 3x the diameter of the Visium spot. Make sure
##  to adjust it if it is too large or too small for your setting.
bw = 3*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
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
scores = TRUE
## Run a robust GWPCA?
robust = FALSE
## Make a cluster for parallel computing (otherwise GWPCA is slow!)
my.cl <- makeClusterGWPCA(type = "FORK")

# Run GWPCA
pcagw <- gwpcaSTE(sfe = sfe,
                  assay = "logcounts",
                  vars = vars,
                  p = p,
                  k = k,
                  bw = bw,
                  kernel = kernel,
                  adaptive = adaptive,
                  scores = scores,
                  robust = robust,
                  cv = cv,
                  future = FALSE,
                  strategy = "cluster",
                  workers = my.cl,
                  verbose = FALSE)


plotGWPCA_global(gwpca = pcagw,
                 comps = 1:10,
                 type = "scree",
                 point_args = list(size = 3, colour = "red"),
                 line_args = list(linewidth = 1, colour = "dodgerblue"))

## Extract leading genes
pcagw <- gwpca_LeadingGene(gwpca = pcagw,
                           m_sfe = sfe,
                           pc_nos = 1:4,
                           type = "single",
                           names = "gene_names")

pcagw <- gwpca_LeadingGene(gwpca = pcagw,
                           m_sfe = sfe,
                           pc_nos = 1:4,
                           genes_n = 4,
                           type = "multi",
                           method = "membership",
                           names = "gene_names")

## Plot leading genes
plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1:2,
                   type = "single",
                   arrange = FALSE)

plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1,
                   type = "multi",
                   arrange = FALSE)

### Plot multi type (extra parameters)
plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1,
                   type = "multi",
                   arrange = FALSE,
                   legend.position = "bottom",
                   cutoff = 12,
                   size = 8)

### Functional clustering
## Steatosis
t2g = read.csv("./data/t2g_files/steatosis_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 2,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO019"]][["steatosis"]] <- gsea_map

## Fibrosis
t2g = read.csv("./data/t2g_files/fibrosis_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO019"]][["fibrosis"]] <- gsea_map

##  NAS > 4 - not corrected for F stage
t2g = read.csv("./data/t2g_files/NAS_RNA-Seq_t2g.csv")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO019"]][["nas"]] <- gsea_map

## Metabolism of lipids
msigdb <- getMSigDBData("Homo sapiens")
reactome_t2g <- getTerm2Gene(msig_data = msigdb, cat = "C2", subcat = "CP:REACTOME")
t2g = reactome_t2g[reactome_t2g$term == "REACTOME_METABOLISM_OF_LIPIDS",]

gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways")
gsea_map_list[["JBO019"]][["lipid_metab"]] <- gsea_map

## Expression of marker genes

p1 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO018,"CYP3A4",assay = "logcounts") # hepatocytes
p2 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO018,"RGS5",assay = "logcounts") #HSC - Hepatic stellate cells (immune)
p3 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO018,"MFAP4",assay = "logcounts") #FB - fibroblast
p4 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO018,"MYH11",assay = "logcounts") #VSMC - Vascular Smooth Muscle Cell (blood vessel)
p1 + p2 + p3 + p4


p1 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO019,"CYP3A4",assay = "logcounts") # hepatocytes
p2 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO019,"RGS5",assay = "logcounts") #HSC - Hepatic stellate cells (immune)
p3 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO019,"MFAP4",assay = "logcounts") #FB - fibroblast
p4 = STExplorer::plotGeneExpression(msfe@sfe_data$JBO019,"MYH11",assay = "logcounts") #VSMC - Vascular Smooth Muscle Cell (blood vessel)
p1 + p2 + p3 + p4


#save.image("/analysis_objs/011124_AR.Rdata")
#load("/analysis_objs/011124_AR.Rdata")

# Additional figures plots ----

## NMF Metagene Signatures
for (s in samples) {
  p <- plotFGWC_nmfMetagenesHeatmap(fgwc_list[[s]])

  print(p)
}

## Pericentral cluster leading genes
leaders_list_single <- list()
leaders_list_multi <- list()
clust2_leaders_single <- list()
clust2_leaders_multi <- list()

for (s in samples) {
  # clust2 <- ifelse(s == "JBO019", 3, 4)
  pcagw <- get(ifelse(s == "JBO019", "pcagw_19", "pcagw_18"))
  clusters <- fgwc_list[[s]][["finaldata"]] %>%
    dplyr::filter(cluster == 2)
  barcodes <- rownames(clusters)

  dt_single <- pcagw$leadingGeneSingle
  rownames(dt_single) <- rownames(pcagw[["loadings"]])

  dt_multi <- pcagw$leadingGeneMulti
  rownames(dt_multi) <- rownames(pcagw[["loadings"]])

  leaders_list_single[[s]] <- dt_single # store exports
  leaders_list_multi[[s]] <- dt_multi # store exports

  clust2_leaders_single[[s]] <- dt_single[rownames(dt_single) %in% barcodes,]
  clust2_leaders_multi[[s]] <- dt_multi[rownames(dt_multi) %in% barcodes,]

  for (pc in paste0("PC", 1:2)) {
    cols_single <- getColours(length(unique(clust2_leaders_single[[s]][,pc])))
    # cols_multi <- getColours(length(unique(clust2_leaders_multi[[s]][,pc])))

    p1 <- ggplot(clust2_leaders_single[[s]],
                 aes(fill = .data[[pc]])) +
      geom_sf(aes(geometry = .data[["geometry"]][["geometry"]])) +
      scale_fill_manual(values = cols_single) +
      labs(title = s,
           fill = paste0(pc,"\n", "Leading\nGene")) +
      theme_void()

    p2 <- ggplot(clust2_leaders_multi[[s]],
                 aes(fill = .data[[pc]])) +
      geom_sf(aes(geometry = .data[["geometry"]][["geometry"]])) +
      # scale_fill_manual(values = cols_multi) +
      labs(title = s,
           fill = paste0(pc,"\n", "Leading\nGenes")) +
      theme_void()

    print(p1)
    print(p2)
  }
}

## Functional clustering plots
for (s in samples) {
  for (c in c("nas", "fibrosis")) {
    gsea_map <- gsea_map_list[[s]][[c]] # %>%
      # filter(p.adjust < 0.1)
    print(plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right", legend.title = "Pathways"))
  }
}

## Multi leaders clustering based on original or absolute leading scores
## JBO019
data_list <- vector(mode = "list", length = 2)
names(data_list) <- c("abs", "original")
for (sort in c("abs", "original")) {
  pcagw_19 <- gwpca_LeadingGene(gwpca = pcagw_19,
                                m_sfe = sfe,
                                pc_nos = 1:4,
                                genes_n = 10,
                                type = "multi",
                                method = "membership",
                                sort = sort,
                                names = "gene_names")

  # create dataset
  data_list[[sort]] <- pcagw_19$leadingGeneMulti
}

## character string clustering
# Load required library
# install.packages("stringdist")
library(stringdist)
library(dbscan)
# library(igraph)
# library(qgraph)

epsilon_abs <- list("1" = c(0.05,0.10,0.15),
                    "2" = c(0.1,0.15,0.20),
                    "3" = c(0.15,0.17,0.2),
                    "4" = c(0.17,0.20,0.23))

epsilon_org <- list("1" = c(0.05,0.08,0.12),
                    "2" = c(0.15,0.2,0.23),
                    "3" = c(0.15,0.2,0.25),
                    "4" = c(0.2,0.23,0.26))

dist_m_list <- list(abs = list(), original = list())
for (sort in c("abs", "original")) {
  data <- data_list[[sort]]
  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    dist_m_list[[sort]][[pc]] <- stringdistmatrix(data[[paste0("PC", pc)]], method = "jaccard")
  }
}

for (sort in c("abs", "original")) {
  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    distance_matrix <- dist_m_list[[sort]][[pc]]

    # DBSCAN clustering
    if (sort == "abs") {
      epsilon <- epsilon_abs
    } else {
      epsilon <- epsilon_org
    }

    dbscan::kNNdistplot(distance_matrix, k =  5)
    abline(h = epsilon[[pc]], lty = 2)
  }
}

for (sort in c("abs", "original")) {
  # DBSCAN clustering
  if (sort == "abs") {
    epsilon <- epsilon_abs
  } else {
    epsilon <- epsilon_org
  }

  data <- data_list[[sort]]

  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    distance_matrix <- dist_m_list[[sort]][[pc]]

    dbscan_19 <- vector(mode = "list", length = length(epsilon))
    names(dbscan_19) <- as.character(epsilon)

    for (i in epsilon[[pc]]) {
      i_c <- as.character(i)
      dbscan_19[[i_c]] <- dbscan(distance_matrix, eps = i, minPts = 5)
      data[[paste0("cluster_PC", pc, "_", i)]] <- dbscan_19[[i_c]]$cluster

      dt <- data %>%
        select(matches(paste0("geom|PC", pc, "_", i))) %>%
        rename("cluster" = paste0("cluster_PC", pc, "_", i))

      p <- ggplot(dt,
                  aes(fill = as.factor(cluster))) +
        geom_sf(aes(geometry = geometry$geometry)) +
        scale_fill_manual(values = getColours(max(dt$cluster) + 1)) +
        labs(fill = paste0("PC", pc, "\nCluster\n", "(eps: ", i, ")\n(", if_else(sort == "abs", "abs", "org"), ")")) +
        theme_void()

      print(p)

      ggsave(paste0("JBO019_GWPCA_", if_else(sort == "abs", "abs", "org"), "_DBSCAN_PC", pc, "_eps", i, ".pdf"),
             dpi = 300)
    }
  }
  data_list[[sort]] <- data
}
data_list_19 <- data_list


## JBO018
sfe <- getSFE(msfe, sample_id = "JBO018")
data_list <- vector(mode = "list", length = 2)
names(data_list) <- c("abs", "original")
for (sort in c("abs", "original")) {
  pcagw_18 <- gwpca_LeadingGene(gwpca = pcagw_18,
                                m_sfe = sfe,
                                pc_nos = 1:4,
                                genes_n = 10,
                                type = "multi",
                                method = "membership",
                                sort = sort,
                                names = "gene_names")

  # create dataset
  data_list[[sort]] <- pcagw_18$leadingGeneMulti
}

## character string clustering
# Load required library
# install.packages("stringdist")
library(stringdist)
library(dbscan)
# library(igraph)
# library(qgraph)

epsilon_abs <- list("1" = c(0.04,0.07,0.1),
                    "2" = c(0.07,0.1,0.12),
                    "3" = c(0.12,0.15,0.18),
                    "4" = c(0.15,0.17,0.2))

epsilon_org <- list("1" = c(0.04,0.07,0.1),
                    "2" = c(0.13,0.15,0.18),
                    "3" = c(0.15,0.18,0.22),
                    "4" = c(0.18,0.20,0.23))

dist_m_list <- list(abs = list(), original = list())
for (sort in c("abs", "original")) {
  data <- data_list[[sort]]
  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    dist_m_list[[sort]][[pc]] <- stringdistmatrix(data[[paste0("PC", pc)]], method = "jaccard")
  }
}

for (sort in c("abs", "original")) {
  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    distance_matrix <- dist_m_list[[sort]][[pc]]

    # DBSCAN clustering
    if (sort == "abs") {
      epsilon <- epsilon_abs
    } else {
      epsilon <- epsilon_org
    }

    dbscan::kNNdistplot(distance_matrix, k =  5)
    abline(h = epsilon[[pc]], lty = 2)
  }
}

for (sort in c("abs", "original")) {
  # DBSCAN clustering
  if (sort == "abs") {
    epsilon <- epsilon_abs
  } else {
    epsilon <- epsilon_org
  }

  data <- data_list[[sort]]

  for (pc in 1:4) {
    message("Working on PC: ", pc)
    # Compute Levenshtein distance matrix
    distance_matrix <- dist_m_list[[sort]][[pc]]

    dbscan_18 <- vector(mode = "list", length = length(epsilon))
    names(dbscan_18) <- as.character(epsilon)

    for (i in epsilon[[pc]]) {
      i_c <- as.character(i)
      dbscan_18[[i_c]] <- dbscan(distance_matrix, eps = i, minPts = 5)
      data[[paste0("cluster_PC", pc, "_", i)]] <- dbscan_18[[i_c]]$cluster

      dt <- data %>%
        select(matches(paste0("geom|PC", pc, "_", i))) %>%
        rename("cluster" = paste0("cluster_PC", pc, "_", i))

      p <- ggplot(dt,
                  aes(fill = as.factor(cluster))) +
        geom_sf(aes(geometry = geometry$geometry)) +
        scale_fill_manual(values = getColours(max(dt$cluster) + 1)) +
        labs(fill = paste0("PC", pc, "\nCluster\n", "(eps: ", i, ")\n(", if_else(sort == "abs", "abs", "org"), ")")) +
        theme_void()

      print(p)

      ggsave(paste0("JBO018_GWPCA_", if_else(sort == "abs", "abs", "org"), "_DBSCAN_PC", pc, "_eps", i, ".pdf"),
             dpi = 300)
    }
  }
  data_list[[sort]] <- data
}
data_list_18 <- data_list

## GWPCA gene signatures per FGWC cluster
# BiocManager::install("org.Hs.eg.db")   # Human gene annotation database
# BiocManager::install("ReactomePA")    # Reactome enrichment
# BiocManager::install("ComplexHeatmap") # heatmap visualisations
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(cols4all)
library(parallel)

# Function: Split strings, convert to Entrez IDs, and run enrichment
process_row <- function(gene_string) {
  # Split the string into a vector of gene symbols
  gene_list <- unlist(strsplit(gene_string, ";"))
  # Check for known genes with different names:
  if ("C10orf10" %in% gene_list) {
    gene_list[gene_list %in% "C10orf10"] <- "DEPP1"
  }
  if ("IGJ" %in% gene_list) {
    gene_list[gene_list %in% "IGJ"] <- "JCHAIN"
  }
  if ("SEPT5" %in% gene_list) {
    gene_list[gene_list %in% "SEPT5"] <- "SEPTIN5"
  }
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Check if conversion is successful
  if (nrow(entrez_ids) == 0) return(list(KEGG = NULL, Reactome = NULL))

  # Perform KEGG enrichment
  kegg_enrichment <- clusterProfiler::enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

  # Perform Reactome enrichment
  reactome_enrichment <- ReactomePA::enrichPathway(gene = entrez_ids$ENTREZID, organism = "human", pvalueCutoff = 0.05)

  # Return results as a list
  list(KEGG = kegg_enrichment, Reactome = reactome_enrichment)
}

# Function: Process the list to aggregate KEGG and Reactome results
process_enrichment_results <- function(dt) {
  # Extract KEGG results
  kegg_results <- lapply(names(dt), function(barcode) {
    if (!is.null(dt[[barcode]]$KEGG)) {
      df <- dt[[barcode]]$KEGG@result
      df$location <- barcode
      df$source <- "KEGG"
      return(df)
    }
    NULL
  }) %>% dplyr::bind_rows()

  # Extract Reactome results
  reactome_results <- lapply(names(dt), function(barcode) {
    if (!is.null(dt[[barcode]]$Reactome)) {
      df <- dt[[barcode]]$Reactome@result
      df$location <- barcode
      df$source <- "Reactome"
      return(df)
    }
    NULL
  }) %>% dplyr::bind_rows()

  list(KEGG = kegg_results, Reactome = reactome_results)
}

calculate_spots_per_cluster <- function(fgwc_list, sample_id) {
  # Extract the cluster assignments from the nested list
  cluster_assignments <- fgwc_list[[sample_id]][["cluster"]]

  # Create a data frame with cluster assignments
  cluster_data <- data.frame(cluster = cluster_assignments)

  # Count the number of spots in each cluster
  cluster_counts <- cluster_data %>%
    group_by(cluster) %>%
    summarise(n_spots = n(), .groups = "drop")

  # Return the data frame with cluster counts
  return(cluster_counts)
}

# Function: Summarise KEGG and Reactome results across locations
summarise_results <- function(results, clust_counts) {
  # Summarise results and normalize Frequency
  results %>%
    group_by(cluster, Description) %>%
    summarise(
      Frequency = n(),  # Count how many times the pathway is enriched
      Mean_pvalue = mean(p.adjust, na.rm = TRUE),
      Median_pvalue = median(p.adjust, na.rm = TRUE),
      Locations = paste(unique(location), collapse = ";"),
      .groups = "drop"
    ) %>%
    # Join with the number of spots per cluster
    left_join(clust_counts, by = "cluster") %>%
    # Normalise Frequency by the number of spots in the cluster
    mutate(Normalised_Frequency = Frequency / n_spots) %>%
    # Arrange by Mean_pvalue
    arrange(Mean_pvalue)
}

# Function: Heatmap of pathway frequencies
create_heatmap <- function(results,
                           color = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100),
                           order = NULL,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_colnames = FALSE,
                           ...) {
  # Create matrix of pathway frequencies across locations
  pathway_matrix <- table(results$Description, results$location)
  if (!is.null(order)) {
    pathway_matrix <- pathway_matrix[,order]
  }

  # Plot heatmap
  pheatmap::pheatmap(
    pathway_matrix,
    color = color,
    scale = "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_colnames = show_colnames,
    ...
  )
}

# Function: Bubble plot of top pathways by frequency and p-value
create_bubble_plot <- function(summary_data, title, n = 10, minmax = c("min", "max")) {
  # Subset for the top pathways (per cluster)
  if (minmax == "max") {
    top_pathways_per_cluster <- summary_data %>%
      group_by(cluster) %>%
      slice_max(order_by = Normalised_Frequency, n = n) %>%
      ungroup() %>%
      mutate(Description = factor(Description, levels = unique(Description)))
  } else {
    top_pathways_per_cluster <- summary_data %>%
      group_by(cluster) %>%
      slice_min(order_by = Mean_pvalue, n = n) %>%
      ungroup() %>%
      mutate(Description = factor(Description, levels = unique(Description)))
  }

  # Create the bubble plot
  ggplot(top_pathways_per_cluster, aes(x = cluster, y = Description, size = Normalised_Frequency, fill = -log10(Mean_pvalue))) +
    geom_point(alpha = 0.7, shape = 21, stroke = 0.5) +
    scale_size(range = c(2, 20), name = "Norm.\nFreq.") +  # Adjust bubble size
    scale_fill_viridis_c(name = "-log10\nMean\npvalue", option = "viridis") +  # Use viridis colour palette
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 12)
    ) +
    labs(
      title = title,
      x = "Cluster",
      y = "Pathway"
    )
}

# Function to create a spatial bubble plot
create_spatial_bubble_plot <- function(data, title, pathway, geoms) {
  # Filter data for the relevant pathways
  filtered_data <- data %>%
    filter(Description %in% pathway)

  p <- ggplot(geoms) +
    geom_sf(aes(geometry = geometry), fill = "white")

  # Create the plot using ggplot2 and geom_sf
  p + geom_sf(data = filtered_data,
              aes(fill = -log10(p.adjust), geometry = geometry),
              alpha = 1, show.legend = TRUE) +
    scale_fill_viridis_c(name = "-log10\nAdjusted\np-value") +  # Use viridis colour palette
    facet_wrap(~cluster) +
    theme_void() +
    labs(title = title)
}

## JBO018 - PC1 enrichment ----
sort <- "abs"
pc <- 1
df <- data_list_18[[sort]]["PC1"]

# Apply function to each row in the dataframe
## Detect number of cores
num_cores <- detectCores() - 1

## Create a cluster
cl <- makeCluster(num_cores)

## Export necessary variables and packages to the cluster
clusterExport(cl, varlist = c("df", "process_row"))
clusterEvalQ(cl, library(clusterProfiler))
clusterEvalQ(cl, library(ReactomePA))
clusterEvalQ(cl, library(org.Hs.eg.db))

## Apply process_row in parallel
df$enrichment_results <- parLapply(cl, df$PC1, process_row)

# Stop the cluster
stopCluster(cl)

# df$enrichment_results <- apply(df["PC1"], 1, process_row)
df_18 <- df

# Filter results for specific cluster
for (i in 0:5) {
  message("Working on cluster: ", i)
  # selection_v <- data_list_18[[sort]][["cluster_PC1_0.1"]] == 0
  dt <- df$enrichment_results
  names(dt) <- rownames(data_list_18[[sort]])
  if (i != 0) {
    selection_v <- fgwc_list[["JBO018"]][["cluster"]] == i
    dt <- dt[selection_v]
  } else {
    NULL
  }

  # Apply to the data
  enrichment_data <- process_enrichment_results(dt)
  spatial_coordinates <- data_list_18[[sort]]$geometry$geometry %>%
    data.frame() %>%
    mutate(barcode = rownames(data_list_18[[sort]]))
  clusters_18 <- fgwc_list[["JBO018"]][["finaldata"]]["cluster"] %>%
    rownames_to_column(var = "barcode") %>%
    cbind(., fgwc_list[["JBO018"]][["membership"]]) %>%
    rename("4" = "1", "5" = "2", "1" = "3", "2" = "4", "3" = "5")
  kegg_results <- enrichment_data$KEGG %>%
    filter(p.adjust < 0.05 & Count >= 2) %>%
    inner_join(spatial_coordinates, by = c("location" = "barcode")) %>%
    inner_join(clusters_18, by = c("location" = "barcode"))
  reactome_results <- enrichment_data$Reactome %>%
    filter(p.adjust < 0.05 & Count >= 2) %>%
    inner_join(spatial_coordinates, by = c("location" = "barcode")) %>%
    inner_join(clusters_18, by = c("location" = "barcode"))

  # Summarise KEGG and Reactome results across locations
  clust_n_18 <- calculate_spots_per_cluster(fgwc_list, "JBO018")
  kegg_summary <- summarise_results(kegg_results, clust_n_18)
  reactome_summary <- summarise_results(reactome_results, clust_n_18)

  if (i == 0) {
    # Export summaries
    write.csv(kegg_summary, file = "kegg_summary_PC1.csv", row.names = FALSE)
    write.csv(reactome_summary, file = "reactome_summary_PC1.csv", row.names = FALSE)
  }

  # Prepare cluster annotations
  kegg_spots_left <- fgwc_list[["JBO018"]]$cluster[unique(kegg_results$location)]
  if (i != 0) {
    # order_by_cl_kegg <- order(kegg_spots_left)
    # Use sym() to refer to column by name dynamically
    column_name <- sym(as.character(i))

    order_by_cl_kegg <- kegg_results %>%
      arrange(desc(!!column_name)) %>%
      select(location) %>%
      unique() %>%
      .[["location"]]
  } else {
    order_by_cl_kegg <- order(kegg_spots_left)
  }
  annot_col_kegg <- data.frame(cluster = kegg_results$cluster,
                               barcode = kegg_results$location,
                               `% in 1` = kegg_results[["1"]],
                               `% in 2` = kegg_results[["2"]],
                               `% in 3` = kegg_results[["3"]],
                               `% in 4` = kegg_results[["4"]],
                               `% in 5` = kegg_results[["5"]],
                               check.names = FALSE) %>%
    unique() %>%
    remove_rownames %>%
    column_to_rownames(var = "barcode")

  react_spots_left <- fgwc_list[["JBO018"]]$cluster[unique(reactome_results$location)]
  if (i != 0) {
    # order_by_cl_react <- order(react_spots_left)
    order_by_cl_react <- reactome_results %>%
      arrange(desc(!!column_name)) %>%
      select(location) %>%
      unique() %>%
      .[["location"]]
  } else {
    order_by_cl_react <- order(react_spots_left)
  }
  annot_col_react <- data.frame(cluster = reactome_results$cluster,
                                barcode = reactome_results$location,
                                `% in 1` = reactome_results[["1"]],
                                `% in 2` = reactome_results[["2"]],
                                `% in 3` = reactome_results[["3"]],
                                `% in 4` = reactome_results[["4"]],
                                `% in 5` = reactome_results[["5"]],
                                check.names = FALSE) %>%
    unique() %>%
    remove_rownames %>%
    column_to_rownames(var = "barcode")

  ann_colours_kegg <- list(
    cluster = c(`1` = "#6fafd6", `2` = "#ABDDA4", `3` = "#FFFFBF", `4` = "#FDAE61", `5` = "#D7191C"),
    `% in 1` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 1`), max(annot_col_kegg$`% in 1`))),
    `% in 2` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 2`), max(annot_col_kegg$`% in 2`))),
    `% in 3` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 3`), max(annot_col_kegg$`% in 3`))),
    `% in 4` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 4`), max(annot_col_kegg$`% in 4`))),
    `% in 5` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 5`), max(annot_col_kegg$`% in 5`)))
  )

  ann_colours_react <- list(
    cluster = c(`1` = "#6fafd6", `2` = "#ABDDA4", `3` = "#FFFFBF", `4` = "#FDAE61", `5` = "#D7191C"),
    `% in 1` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 1`), max(annot_col_react$`% in 1`))),
    `% in 2` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 2`), max(annot_col_react$`% in 2`))),
    `% in 3` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 3`), max(annot_col_react$`% in 3`))),
    `% in 4` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 4`), max(annot_col_react$`% in 4`))),
    `% in 5` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 5`), max(annot_col_react$`% in 5`)))
  )

  # Plot KEGG heatmap
  create_heatmap(kegg_results,
                 color = colorRampPalette(c("white", "firebrick3"))(2),
                 order = order_by_cl_kegg,
                 annotation_col = annot_col_kegg,
                 annotation_colors = ann_colours_kegg,
                 cluster_cols = TRUE,
                 filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_cl", i, ".pdf"))

  # Plot Reactome heatmap
  create_heatmap(reactome_results,
                 color = colorRampPalette(c("white", "firebrick3"))(2),
                 order = order_by_cl_react,
                 annotation_col = annot_col_react,
                 annotation_colors = ann_colours_react,
                 cluster_cols = TRUE,
                 filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_cl", i, ".pdf"))

  if (i >= 0) {
    # Run once with cluster_cols = FALSE
    create_heatmap(kegg_results,
                   color = colorRampPalette(c("white", "firebrick3"))(2),
                   order = order_by_cl_kegg,
                   annotation_col = annot_col_kegg,
                   annotation_colors = ann_colours_kegg,
                   cluster_cols = FALSE,
                   filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_cl", i, "_clFALSE.pdf"))

    create_heatmap(reactome_results,
                   color = colorRampPalette(c("white", "firebrick3"))(2),
                   order = order_by_cl_react,
                   annotation_col = annot_col_react,
                   annotation_colors = ann_colours_react,
                   cluster_cols = FALSE,
                   filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_cl", i, "_clFALSE.pdf"))
  }

  if (i == 0) {
    for (minmax in c("min", "max")) {
      # KEGG Bubble Plot
      create_bubble_plot(kegg_summary, "KEGG Pathway Enrichment in Cluster", n = 20, minmax)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_sum_cl", i, "_", minmax, "_.pdf"),
             dpi = 300, height = 11, width = 10)

      # Reactome Bubble Plot
      create_bubble_plot(reactome_summary, "Reactome Pathway Enrichment in Cluster", n = 20, minmax)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_sum_cl", i, "_", minmax, "_.pdf"),
             dpi = 300, height = 14, width = 17)
    }

    pathways_k <- c("Linoleic acid metabolism",
                    "Retinol metabolism",
                    "Alcoholic liver disease",
                    "Glycolysis - Gluconeogenesis",
                    "Vascular smooth muscle contraction")
    pathways_r <- c("Biosynthesis of DHA-derived SPMs",
                    "Xenobiotics",
                    "RA biosynthesis pathway",
                    "Defective C1GALT1C1 causes TNPS",
                    "O2-CO2 exchange in erythrocytes")

    geometries_18 <- data_list_18[[sort]]$geometry$geometry %>%
      data.frame() %>%
      mutate(barcode = rownames(data_list_18[[sort]]))

    kegg_results$Description <- gsub("/", "-", kegg_results$Description)
    reactome_results$Description <- gsub("/", "-", reactome_results$Description)

    for (path in pathways_k) {
      message(path)
      # Create spatial bubble plots for KEGG and Reactome results
      create_spatial_bubble_plot(kegg_results,
                                 path,
                                 path,
                                 geoms = geometries_18)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_spatial_", path, ".pdf"),
             dpi = 300, height = 7, width = 10)
    }

    for (path in pathways_r) {
      message(path)
      create_spatial_bubble_plot(reactome_results,
                                 path,
                                 path,
                                 geoms = geometries_18)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_spatial_", path, ".pdf"),
             dpi = 300, height = 7, width = 10)
    }
  }

}

## JBO018 - PC1-3 enrichment ----
# Run again for leading genes from a combo of PCs
sort <- "abs"
pc <- "1.2.3"
df <- data_list_18[[sort]][c("PC1", "PC2", "PC3")] %>%
  mutate(PCs = paste(PC1, PC2, PC3, sep = ";"))

# Apply function to each row in the dataframe
## Detect number of cores
num_cores <- detectCores() - 1

## Create a cluster
cl <- makeCluster(num_cores)

## Export necessary variables and packages to the cluster
clusterExport(cl, varlist = c("df", "process_row"))
clusterEvalQ(cl, library(clusterProfiler))
clusterEvalQ(cl, library(ReactomePA))
clusterEvalQ(cl, library(org.Hs.eg.db))

## Apply process_row in parallel
df$enrichment_results <- parLapply(cl, df$PCs, process_row)

# Stop the cluster
stopCluster(cl)

# df$enrichment_results <- apply(df["PC1"], 1, process_row)
df_18_1.2.3 <- df

# Filter results for specific cluster
for (i in 0:5) {
  message("Working on cluster: ", i)
  # selection_v <- data_list_18[[sort]][["cluster_PC1_0.1"]] == 0
  dt <- df$enrichment_results
  names(dt) <- rownames(data_list_18[[sort]])
  if (i != 0) {
    selection_v <- fgwc_list[["JBO018"]][["cluster"]] == i
    dt <- dt[selection_v]
  } else {
    NULL
  }

  # Apply to the data
  enrichment_data <- process_enrichment_results(dt)
  spatial_coordinates <- data_list_18[[sort]]$geometry$geometry %>%
    data.frame() %>%
    mutate(barcode = rownames(data_list_18[[sort]]))
  clusters_18 <- fgwc_list[["JBO018"]][["finaldata"]]["cluster"] %>%
    rownames_to_column(var = "barcode") %>%
    cbind(., fgwc_list[["JBO018"]][["membership"]]) %>%
    rename("4" = "1", "5" = "2", "1" = "3", "2" = "4", "3" = "5")
  kegg_results <- enrichment_data$KEGG %>%
    filter(p.adjust < 0.05 & Count >= 2) %>%
    inner_join(spatial_coordinates, by = c("location" = "barcode")) %>%
    inner_join(clusters_18, by = c("location" = "barcode"))
  reactome_results <- enrichment_data$Reactome %>%
    filter(p.adjust < 0.05 & Count >= 2) %>%
    inner_join(spatial_coordinates, by = c("location" = "barcode")) %>%
    inner_join(clusters_18, by = c("location" = "barcode"))

  # Summarise KEGG and Reactome results across locations
  clust_n_18 <- calculate_spots_per_cluster(fgwc_list, "JBO018")
  kegg_summary <- summarise_results(kegg_results, clust_n_18)
  reactome_summary <- summarise_results(reactome_results, clust_n_18)

  if (i == 0) {
    # Export summaries
    write.csv(kegg_summary, file = "kegg_summary_PC1.2.3.csv", row.names = FALSE)
    write.csv(reactome_summary, file = "reactome_summary_PC1.2.3.csv", row.names = FALSE)
  }

  # Prepare cluster annotations
  kegg_spots_left <- fgwc_list[["JBO018"]]$cluster[unique(kegg_results$location)]
  if (i != 0) {
    # order_by_cl_kegg <- order(kegg_spots_left)
    # Use sym() to refer to column by name dynamically
    column_name <- sym(as.character(i))

    order_by_cl_kegg <- kegg_results %>%
      arrange(desc(!!column_name)) %>%
      select(location) %>%
      unique() %>%
      .[["location"]]
  } else {
    order_by_cl_kegg <- order(kegg_spots_left)
  }
  annot_col_kegg <- data.frame(cluster = kegg_results$cluster,
                               barcode = kegg_results$location,
                               `% in 1` = kegg_results[["1"]],
                               `% in 2` = kegg_results[["2"]],
                               `% in 3` = kegg_results[["3"]],
                               `% in 4` = kegg_results[["4"]],
                               `% in 5` = kegg_results[["5"]],
                               check.names = FALSE) %>%
    unique() %>%
    remove_rownames %>%
    column_to_rownames(var = "barcode")

  react_spots_left <- fgwc_list[["JBO018"]]$cluster[unique(reactome_results$location)]
  if (i != 0) {
    # order_by_cl_react <- order(react_spots_left)
    order_by_cl_react <- reactome_results %>%
      arrange(desc(!!column_name)) %>%
      select(location) %>%
      unique() %>%
      .[["location"]]
  } else {
    order_by_cl_react <- order(react_spots_left)
  }
  annot_col_react <- data.frame(cluster = reactome_results$cluster,
                                barcode = reactome_results$location,
                                `% in 1` = reactome_results[["1"]],
                                `% in 2` = reactome_results[["2"]],
                                `% in 3` = reactome_results[["3"]],
                                `% in 4` = reactome_results[["4"]],
                                `% in 5` = reactome_results[["5"]],
                                check.names = FALSE) %>%
    unique() %>%
    remove_rownames %>%
    column_to_rownames(var = "barcode")

  ann_colours_kegg <- list(
    cluster = c(`1` = "#6fafd6", `2` = "#ABDDA4", `3` = "#FFFFBF", `4` = "#FDAE61", `5` = "#D7191C"),
    `% in 1` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 1`), max(annot_col_kegg$`% in 1`))),
    `% in 2` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 2`), max(annot_col_kegg$`% in 2`))),
    `% in 3` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 3`), max(annot_col_kegg$`% in 3`))),
    `% in 4` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 4`), max(annot_col_kegg$`% in 4`))),
    `% in 5` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_kegg$`% in 5`), max(annot_col_kegg$`% in 5`)))
  )

  ann_colours_react <- list(
    cluster = c(`1` = "#6fafd6", `2` = "#ABDDA4", `3` = "#FFFFBF", `4` = "#FDAE61", `5` = "#D7191C"),
    `% in 1` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 1`), max(annot_col_react$`% in 1`))),
    `% in 2` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 2`), max(annot_col_react$`% in 2`))),
    `% in 3` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 3`), max(annot_col_react$`% in 3`))),
    `% in 4` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 4`), max(annot_col_react$`% in 4`))),
    `% in 5` = c4a(palette = "matplotlib.viridis", range = c(min(annot_col_react$`% in 5`), max(annot_col_react$`% in 5`)))
  )

  # Plot KEGG heatmap
  create_heatmap(kegg_results,
                 color = colorRampPalette(c("white", "firebrick3"))(2),
                 order = order_by_cl_kegg,
                 annotation_col = annot_col_kegg,
                 annotation_colors = ann_colours_kegg,
                 cluster_cols = TRUE,
                 filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_cl", i, ".pdf"))

  # Plot Reactome heatmap
  create_heatmap(reactome_results,
                 color = colorRampPalette(c("white", "firebrick3"))(2),
                 order = order_by_cl_react,
                 annotation_col = annot_col_react,
                 annotation_colors = ann_colours_react,
                 cluster_cols = TRUE,
                 filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_cl", i, ".pdf"),
                 width = 20)

  if (i >= 0) {
    # Run once with cluster_cols = FALSE
    create_heatmap(kegg_results,
                   color = colorRampPalette(c("white", "firebrick3"))(2),
                   order = order_by_cl_kegg,
                   annotation_col = annot_col_kegg,
                   annotation_colors = ann_colours_kegg,
                   cluster_cols = FALSE,
                   filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_cl", i, "_clFALSE.pdf"))

    create_heatmap(reactome_results,
                   color = colorRampPalette(c("white", "firebrick3"))(2),
                   order = order_by_cl_react,
                   annotation_col = annot_col_react,
                   annotation_colors = ann_colours_react,
                   cluster_cols = FALSE,
                   filename = paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_cl", i, "_clFALSE.pdf"),
                   width = 20)
  }

  if (i == 0) {
    for (minmax in c("min", "max")) {
      # KEGG Bubble Plot
      create_bubble_plot(kegg_summary, "KEGG Pathway Enrichment in Cluster", n = 20, minmax)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_sum_cl", i, "_", minmax, "_.pdf"),
             dpi = 300, height = 11, width = 10)

      # Reactome Bubble Plot
      create_bubble_plot(reactome_summary, "Reactome Pathway Enrichment in Cluster", n = 40, minmax)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_sum_cl", i, "_", minmax, "_.pdf"),
             dpi = 300, height = 14, width = 17)
    }

    pathways_k <- c("Linoleic acid metabolism",
                  "Retinol metabolism",
                  "Alcoholic liver disease",
                  "Glycolysis - Gluconeogenesis",
                  "Vascular smooth muscle contraction")
    pathways_r <- c("Biosynthesis of DHA-derived SPMs",
                    "Xenobiotics",
                    "RA biosynthesis pathway",
                    "Defective C1GALT1C1 causes TNPS",
                    "O2-CO2 exchange in erythrocytes")

    geometries_18 <- data_list_18[[sort]]$geometry$geometry %>%
      data.frame() %>%
      mutate(barcode = rownames(data_list_18[[sort]]))

    kegg_results$Description <- gsub("/", "-", kegg_results$Description)
    reactome_results$Description <- gsub("/", "-", reactome_results$Description)

    for (path in pathways_k) {
      message(path)
      # Create spatial bubble plots for KEGG and Reactome results
      create_spatial_bubble_plot(kegg_results,
                                 path,
                                 path,
                                 geoms = geometries_18)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadKEGG_PC", pc, "_spatial_", path, ".pdf"),
             dpi = 300, height = 7, width = 10)
    }

    for (path in pathways_r) {
      message(path)
      create_spatial_bubble_plot(reactome_results,
                                 path,
                                 path,
                                 geoms = geometries_18)
      ggsave(paste0("JBO018_GWPCA_", dplyr::if_else(sort == "abs", "abs", "org"), "_top10LeadREACT_PC", pc, "_spatial_", path, ".pdf"),
             dpi = 300, height = 7, width = 10)
    }
  }
}

sfe_5 <- STExplorer::getSFE(msfe, "JBO018")[,fgwc_list$JBO018$cluster == 5]
STExplorer::plotGeneExpression(sfe_5, "CCL19", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "CCL21", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "FBLN1", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "MFAP4", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "ACTA2", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "ACTG2", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "LMOD1", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "MYH11", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "TPM2", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "HBA1", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "HBA2", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_5, "HBB", sample_id = "JBO018", assay = "logcounts")

sfe_4 <- STExplorer::getSFE(msfe, "JBO018")[,fgwc_list$JBO018$cluster == 4]
STExplorer::plotGeneExpression(sfe_4, "CYP2E1", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_4, "CYP3A4", sample_id = "JBO018", assay = "logcounts")
STExplorer::plotGeneExpression(sfe_4, "CYP2A6", sample_id = "JBO018", assay = "logcounts")

for (i in 1:5) {
  sfe_i <- STExplorer::getSFE(msfe, "JBO018")[,fgwc_list$JBO018$cluster == i]
  print(STExplorer::plotGeneExpression(sfe_i, "CYP2E1", sample_id = "JBO018", assay = "logcounts"))
  print(STExplorer::plotGeneExpression(sfe_i, "CYP3A4", sample_id = "JBO018", assay = "logcounts"))
  print(STExplorer::plotGeneExpression(sfe_i, "CYP2A6", sample_id = "JBO018", assay = "logcounts"))
}

