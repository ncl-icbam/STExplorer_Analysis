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
# This script is for preprocessing, QC-ing, and normalising the data.
#
#
# ---------------------------------------------------------------------------- #
## ----QC ----------------------------------------------------------------------
### ----Calculate QC Metrics ---------------------------------------------------
## Mark a subset of mitochondrial genes
is_mito <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "(^MT-)|(^mt-)|(^MTRNR)",
                     set = "rowData",
                     col_name = "symbol")
## Mark a subset of haemoglobin genes
is_hemo <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "^HB",
                     set = "rowData",
                     col_name = "symbol")
## Mark a subset of Ribosomal genes
is_ribo <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "^RPS|^RPL",
                     set = "rowData",
                     col_name = "symbol")

for (id in sampleNames_selected) {
  message("Working on sample: ", id)
  ## Add location-related statistics
  msfe <- addPerLocQC(msfe,
                      sample_id = id,
                      gTruth = gTruth_list[[id]],
                      assay = "counts",
                      MARGIN = 2,
                      subsets = list(mito = is_mito[[id]],
                                     hemo = is_hemo[[id]],
                                     ribo = is_ribo[[id]]))
  message("\tAdded location-related statistics")

  ## Add geometries
  msfe <- addGeometries(msfe,
                        samples = sampleDir[id],
                        sample_id = id,
                        res = "fullres",
                        flipped = FALSE,
                        geoms = "both")
  message("\tAdded geometries")

  ## Add gene/feature-related statistics
  msfe <- addPerGeneQC(msfe,
                       sample_id = id,
                       assay = "counts",
                       version = NULL,
                       mirror = NULL,
                       add = "none")
  message("\tAdded gene/feature-related statistics")
}

## Check that annotations were successfully added
for (i in seq_along(sampleNames_selected)) {
  message(sampleNames_selected[i], ": ",
          "annotation" %in% colnames(colData(msfe@sfe_data[[i]])))
}


# ---------------------------------------------------------------------------- #
### ----Quality Control: Spots--------------------------------------------------
#### ----Library size ----------------------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_library_size/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Select filter
  min <- 350 # at least 350 unique molecules per spot

  ## Density and histogram of library sizes
  p1 <- plotQC_hist(sfe, metric = "libsize",
                    limits = c(min, NA),
                    hist_args = list(bins = 100)) +
    theme(axis.text.x = element_text(size = 12, angle = 45))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "libsize", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_LibSize(sfe, sample_id = TRUE,
                             min_t = min,
                             max_t = NA)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "libsize", sample_id = TRUE, type = "hex", size = 0.5)

  ## Plot tissue image
  p4 <- plotQC_tissueImg(sfe, res = "lowres", sample_id = id, type = "none")

  ## Patch, print and save
  print((p1 + p2) / (p3 + p4))

  prfx <- id
  main <- "_QC_LibSize"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3, p4)
}


# ---------------------------------------------------------------------------- #
#### ----Number of expressed genes ---------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_gene_number/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Select filter
  min <- 100 # at least 100 genes per spot

  ## Density and histogram of expressed genes
  p1 <- plotQC_hist(sfe, metric = "detected", c(min, NA))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "detected", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_GenesExpr(sfe, sample_id = TRUE,
                               min_t = min,
                               max_t = NA)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "detected", sample_id = TRUE, type = "hex", size = 0.5)

  ## Patch, print and save
  print(p1 + p2 + p3)

  prfx <- id
  main <- "_QC_GeneNumber"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3)
}


# ---------------------------------------------------------------------------- #
#### ----Percent of mito expression --------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_mito_percent/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Density and histogram of percentage of mitochondrial expression
  p1 <- plotQC_hist(sfe, metric = "mito", limits = c(NA, 15))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "mito", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, min_t = NA, max_t = 15)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "mito", sample_id = TRUE, type = "hex", size = 0.5)

  ## Patch, print and save
  print(p1 + p2 + p3)

  prfx <- id
  main <- "_QC_MitoPercent"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3)
}


# ---------------------------------------------------------------------------- #
#### ----Percent of haemoglobin expression -------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_hemo_percent/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Get the vector
  hemo_out <- list(id = sfe@colData$subsets_hemo_percent > 30)
  names(hemo_out) <- id

  ## Map the library sizes
  p1 <- plotQC_map(sfe, metric = "custom", type = "hex", size = 0.5, metric_name = "subsets_hemo_percent")

  ## Select threshold
  sfe <- setQCthresh_custom(sfe, sample_id = TRUE, MARGIN = 2, qcMetric = hemo_out)

  ## Check putative spatial patterns of removed spots
  p2 <- plotQC_filtered(sfe,
                        metric = "custom",
                        sample_id = TRUE,
                        type = "hex",
                        size = 0.5,
                        metric_name = "qc_hemo_out",
                        metric_lab = "Filtered for % Hemoglobin Expression")

  ## Patch, print and save
  print(p1 + p2)

  prfx <- id
  main <- "_QC_HemoPercent"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, hemo_out)
}


# ---------------------------------------------------------------------------- #
#### ----Remove low-quality spots ----------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_filter_spots/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Set the combined filtering threshold using the QC metrics
  sfe <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)

  ## Check putative spatial patterns of removed spots
  plotQC_filtered(sfe, metric = "discard", sample_id = TRUE, type = "hex")

  prfx <- id
  main <- "_QC_SpotsFiltered"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Remove combined set of low-quality spots
  sfe <- applyQCthresh_loc(sfe, sample_id = TRUE)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe)
}


# ---------------------------------------------------------------------------- #
### ----Quality Control: Features ----------------------------------------------

# ---------------------------------------------------------------------------- #
#### ----Mark specific gene sets -----------------------------------------------
## Genes in the XY chromosomes
xychr_genes <- unique(subset(annotation_gene, annotation_gene$chromosome %in% c("X", "Y"))$gene_name)
is_xy <- lapply(sampleNames_selected, function(x){rowData(msfe@sfe_data[[x]])[["gene_name"]] %in% xychr_genes})
names(is_xy) <- sampleNames_selected

## Keep only protein coding, IG and TR genes
keep_biotype <- c("protein_coding", grep("_gene", unique(annotation_gene$biotype), value = TRUE))
coding_genes <- annotation_gene[annotation_gene$biotype %in% keep_biotype, "gene_name"]
is_ncGene <- lapply(sampleNames_selected, function(x){!rowData(msfe@sfe_data[[x]])[["gene_name"]] %in% coding_genes$gene_name})
names(is_ncGene) <- sampleNames_selected

## Remove mitochondrial, hemoglobin, ribosomal, XY, and anything NOT protein coding or IG or TR.
for (id in sampleNames_selected) {
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_mito[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_hemo[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_ribo[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_xy[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_ncGene[[id]])
}

# ---------------------------------------------------------------------------- #
#### ----Remove low-quality genes ----------------------------------------------
## QC discard Features
## Set the combined filtering threshold using the QC metrics
msfe <- setQCtoDiscard_feat(msfe, filters = TRUE)

## Apply gene-level QC threshold
msfe <- applyQCthresh_feat(msfe)


# ---------------------------------------------------------------------------- #
## ----Normalisation of Counts -------------------------------------------------
### ----Get MultiSFE object ----------------------------------------------------
sfe_multi <- getMultiSFE(msfe, sampleNames_selected)

## Add the subject alias in the colData to regress it out downstream
sfe_multi@colData$subject_alias <- left_join(as.data.frame(sfe_multi@colData["sample_id"]),
                                             metadata,
                                             by = "sample_id")[["subject_alias"]]

## Remove genes with zero expression in all samples
sfe_multi <- setQCthresh_ZeroExpr(sfe_multi)
sfe_multi <- setQCtoDiscard_feat(sfe_multi)
sfe_multi <- applyQCthresh_feat(sfe_multi)

# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
saveRDS(sfe_multi, "../STExplorer_Analysis/Lung/analysis_objs/sfe_multi.rds")

### ----Convert to Seurat object -----------------------------------------------
seu <- as.Seurat(sfe_multi, data = NULL)


### ----SCTransform counts -----------------------------------------------------
## Normalise and regress out sampleID and donor, to remove major effects from
## technical and interindividual differences.
vars_to_reg <- c("sample_id", "subject_alias")

options(future.globals.maxSize = 8000 * 1024^2)

seu <- SCTransform(seu, vars.to.regress = vars_to_reg, assay = "originalexp", conserve.memory = TRUE)

