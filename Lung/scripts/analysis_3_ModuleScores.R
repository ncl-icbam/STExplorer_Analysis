## Analysis of human lung pulmonary fibrosis data  using STExplorer
# The data used here are from:
#
# Franzén, L., Olsson Lindvall, M., Hühn, M. et al.
# Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis.
# Nat Genet 56, 1725–1736 (2024).
# https://doi.org/10.1038/s41588-024-01819-2
#
#
# ---------------------------------------------------------------------------- #
#
#
# This script is for generating the SenMayo, CM, EMT, and TGFB scores.
#
#
# ---------------------------------------------------------------------------- #
## ----Generate Module Scores --------------------------------------------------
### ----Senescence -------------------------------------------------------------
# optimalBin_senM <- optimalBinSize(seu, nbin = 30)
optimalBin_senM <- 24
seu <- AddModuleScore(seu,
                      features = list(senMayo = senMayoGenes_ENSG),
                      assay = "SCT",
                      name = "senMayo",
                      nbin = optimalBin_senM)

### ----EMT --------------------------------------------------------------------
# optimalBin_emt <- optimalBinSize(seu)
optimalBin_emt <- 24
seu <- AddModuleScore(seu,
                      features = list(emt = emtGenes_ENSG),
                      assay = "SCT",
                      name = "emt",
                      nbin = optimalBin_emt)

### ----ECM --------------------------------------------------------------------
# optimalBin_ecm <- optimalBinSize(seu)
optimalBin_ecm <- 24
seu <- AddModuleScore(seu,
                      features = list(ecm = ecmGenes_ENSG),
                      assay = "SCT",
                      name = "ecm",
                      nbin = optimalBin_ecm)

### ----ECM --------------------------------------------------------------------
# optimalBin_tgfb <- optimalBinSize(seu)
optimalBin_tgfb <- 24
seu <- AddModuleScore(seu,
                      features = list(tgfb = tgfbGenes_ENSG),
                      assay = "SCT",
                      name = "tgfb",
                      nbin = optimalBin_tgfb)

## Some genes from the above feature sets are not present in the data set.
## These genes have been filtered out. However, it is good to know which are
## these genes:
cat("Genes from SenMayo gene set not present in the dataset:", sort(names(senMayoGenes_ENSG[!senMayoGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from EMT gene set not present in the dataset:", sort(names(emtGenes_ENSG[!emtGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from ECM gene set not present in the dataset:", sort(names(ecmGenes_ENSG[!ecmGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from TGFB gene set not present in the dataset:", sort(names(tgfbGenes_ENSG[!tgfbGenes_ENSG %in% rownames(seu)])), "\n")

# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
saveRDS(seu, file = "../STExplorer_Analysis/Lung/analysis_objs/seu.rds")

## ----Convert back to SFE -----------------------------------------------------
## Move the SCT counts and the senMayo/emt/ecm/tgfb scores into the sfe_multi object.
## Then split the sfe_multi object into multiple SFE objects inside a MetaSFE

## 1. First remove the genes that have been filtered out after SCTransform
to_keep <- rownames(seu)
sfe_multi <- sfe_multi[to_keep,]

## 2. Add the SCTransformed counts
assay(sfe_multi, "sct") <- seu[["SCT"]]$data

## 3. Add the senMayo/emt/ecm/tgfb
colData(sfe_multi)$senMayo <- as.vector(seu$senMayo1)
colData(sfe_multi)$emt <- as.vector(seu$emt1)
colData(sfe_multi)$ecm <- as.vector(seu$ecm1)
colData(sfe_multi)$tgfb <- as.vector(seu$tgfb1)

## 4. Remove spots without annotation (71 in total)
for (i in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[i]
  ## If the multi-sample SFE object was saved and loaded again into the
  ## environment, the images are going to be of class "PackedRasterImage". If
  ## this is the case, we need to "unwrap" the packaged images.
  is_Packed <- inherits(imgData(sfe_multi)$data[[i]], "PackedRasterImage")
  if (is_Packed) {
    img <- SpatRasterImage(unwrap(imgData(sfe_multi)$data[[i]]))
    sfe_multi@int_metadata$imgData$data[[i]] <- img
  }
}
sfe_multi <- sfe_multi[,!is.na(colData(sfe_multi)["annotation"])]

## 5. Split into single SFEs inside an MSFE object
msfe <- addMultiSFE(msfe, sfe_multi)


# ---------------------------------------------------------------------------- #
## ----Add Distance Matrices ---------------------------------------------------
## Calculate a simple distance matrix
msfe <- addDistMat(msfe, p = 2)


# ---------------------------------------------------------------------------- #
## ----Import cell abundances --------------------------------------------------
## Import the cell abundances into the colData
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Extract abundances for specific cell types -------------------------------#
  message("\tExtracting abundances...")
  ### Extract barcodes from colData because SFEs are QC'ed but cell abundances not
  spotBarcodes <- colData(msfe@sfe_data[[id]])[["Barcode"]]
  ### Extract abundances
  abundance_df <- cellAbundance[[id]] %>%
    select(c("Barcode", paste0("q05cell_abundance_w_sf_", cellTypes))) %>%
    filter(Barcode %in% spotBarcodes)
  ### Order abundances to be in the same order as the spots in the colData
  order <- match(abundance_df$Barcode, spotBarcodes)
  abundance_df <- abundance_df[order,]
  ## Import abundances into colData -------------------------------------------#
  message("\tImporting to colData...")
  ### Merge with existing colData because otherwise it gets overwritten
  to_colData <- colData(msfe@sfe_data[[id]]) %>%
    as.data.frame() %>%
    left_join(abundance_df, by = "Barcode") %>%
    DataFrame()
  ### Import
  colData(msfe@sfe_data[[id]]) <- to_colData
  colnames(msfe@sfe_data[[id]]) <- spotBarcodes
}
## Housekeeping
rm(spotBarcodes, abundance_df, order, to_colData)

# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
saveRDS(sfe_multi, "../STExplorer_Analysis/Lung/analysis_objs/sfe_multi.rds")
saveRDS(msfe, "../STExplorer_Analysis/Lung/analysis_objs/msfe.rds")

