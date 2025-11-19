# Load R packages
library(STExplorer)
library(readr)
library(ggplot2)
library(svglite)
# library(magick)
library(grid)
library(gridExtra)
library(RColorBrewer)
# library(AnnotationHub)

## Data import and preparation of the MSFE object
msfe <- MetaSpatialFeatureExperiment()

workingDir <- c("/mnt/c/Users/Lefteris/Downloads/STMiner_test/")
sampleDir <- c("/mnt/c/Users/Lefteris/Downloads/STMiner_test/data/H1_5/")

sampleNames <- c("H1_5")
names(sampleDir) <- sampleNames

for (i in seq_along(sampleNames)) {
    msfe <- addSFE(
        msfe,
        read10xVisiumSFE(
            samples = sampleDir[i],
            sample_id = sampleNames[i],
            type = "HDF5",
            data = "filtered",
            images = "lowres",
            style = "W",
            zero.policy = TRUE
        )
    )
}

ground_truth <- read_csv(file.path(sampleDir, "H1_5_Final_Consensus_Annotations.csv"))
gTruth_list <- list(H1_5 = ground_truth)

## Spot-level Quality Control

## Mark a subset of mitochondrial genes
is_mito <- getSubset(msfe,
    sample_id = TRUE,
    subset = "(^MT-)|(^mt-)",
    set = "rowData",
    col_name = "symbol"
)

for (i in seq_along(sampleNames)) {
    message("Working on sample: ", sampleNames[i])
    ## Add location-related statistics
    msfe <- addPerLocQC(msfe,
        sample_id = sampleNames[i],
        gTruth = gTruth_list[[i]],
        assay = "counts",
        MARGIN = 2,
        subsets = list(mito = is_mito[[i]])
    )
    message("\tAdded location-related statistics")

    ## Add geometries
    msfe <- addGeometries(msfe,
        samples = sampleDir[i],
        sample_id = sampleNames[i],
        res = "fullres",
        flipped = FALSE
    )
    message("\tAdded geometries")

    ## Add gene/feature-related statistics
    msfe <- addPerGeneQC(msfe,
        sample_id = sampleNames[i],
        assay = "counts",
        version = NULL,
        mirror = NULL,
        add = c("zeroexpr", "exprstats")
    )
    message("\tAdded gene/feature-related statistics")
}

## Set QC thresholds

for (i in seq_along(sampleNames)) {
    msfe@sfe_data[[i]] <- setQCthresh_LibSize(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.15)), max_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.99)))
    msfe@sfe_data[[i]] <- setQCthresh_Mito(msfe@sfe_data[[i]], sample_id = TRUE, min_t = NA, max_t = quantile(msfe@sfe_data[[i]]@colData$subsets_mito_percent, probs = c(.99), na.rm = TRUE))
    msfe@sfe_data[[i]] <- setQCthresh_GenesExpr(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.15)), max_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.99)))

    ## Set the combined filtering threshold using the QC metrics
    msfe@sfe_data[[i]] <- setQCtoDiscard_loc(msfe@sfe_data[[i]], sample_id = TRUE, filters = TRUE)
}

## Remove combined set of low-quality spots
msfe <- applyQCthresh_loc(msfe, sample_id = TRUE)

## Normalisation of counts

### Compute library size factors
msfe <- computeLibSizeFactors(msfe)

### Log-tranformation of counts
msfe <- normaliseCounts(msfe)

## Gene-level Quality Control

## Calculate the mean of log counts over the number of locations a gene is present
msfe <- perGeneLogMean(msfe, sample_id = TRUE)

## Zero expression genes
msfe <- setQCthresh_ZeroExpr(msfe, sample_id = TRUE)

## Lowly expressed (noise?!) genes
msfe <- setQCthresh_LowLogMean(msfe, threshold = 1, sample_id = TRUE)

## Remove mitochondrial and other genes

for (i in seq_along(sampleNames)) {
    msfe@sfe_data[[i]] <- setQCthresh_custom(msfe@sfe_data[[i]], MARGIN = 1, qcMetric = is_mito[[i]])
}

## QC discard Features
## Set the combined filtering threshold using the QC metrics
msfe <- setQCtoDiscard_feat(msfe, filters = TRUE, sample_id = TRUE)

## FEATURE SELECTION
## Apply gene-level QC threshold
msfe <- applyQCthresh_feat(msfe, sample_id = TRUE)

## Fit mean-variance relationship
dec <- modelGeneVariance(msfe, sample_id = TRUE, method = "Var")

## Select top HVGs

top_10_hvgs <- getTopHighVarGenes(dec,
    var.field = "bio",
    prop = 0.1,
    var.threshold = 0,
    fdr.threshold = 0.1
)

## Add a neighbour graph using a weighted distance matrix
msfe <- addSpatialNeighGraphs(msfe, sample_id = TRUE, type = "knearneigh", style = "W", distMod = "raw", k = 6)

## Calculate a simple distance matrix
msfe <- addDistMat(msfe, p = 2)

## Fuzzy Geographically Weighted Clustering (FGWC)

## Find optimum number of Factors
best_k_nmf <- fgwc_nmfFactorNumber(
    m_sfe = msfe,
    sample_id = "H1_5",
    assay = "logcounts",
    top_hvgs = top_10_hvgs[["H1_5"]],
    k_range = seq(2, 5, 1),
    n_cores = 1,
    do_plot = FALSE,
    seed = 1,
    loss = "mse",
    max.iter = 250
)

best_k_nmf$k

## Run NMF
sfe_nmf <- fgwc_nmf(
    m_sfe = msfe,
    sample_id = "H1_5",
    top_hvgs = top_10_hvgs[["H1_5"]],
    ncomponents = best_k_nmf[["k"]]
)

## Set number of FGWC clusters

fgwc_param <- fgwc_params(
    algorithm = "classic",
    ncluster = 5
)

# Run FGWC ----
fgwc_list <- fgwcSTE(
    m_sfe = msfe,
    sample_id = "H1_5",
    data = sfe_nmf,
    dMetric = "euclidean",
    algorithm = "classic",
    parameters = fgwc_param
)


# ================================
# FGWC metrics + exports (append)
# ================================
# Packages used below (I need to install/load some of them - check)
pkgs <- c("Matrix", "data.table", "jsonlite", "pROC", "FNN", "spdep", "fclust", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(Matrix)
library(data.table)
library(jsonlite)
library(pROC)
library(FNN)
library(spdep)
library(fclust)

res_fld <- file.path(workingDir, "results")
exp_fld <- file.path(workingDir, "export")
dir.create(res_fld, showWarnings = FALSE)
dir.create(exp_fld, showWarnings = FALSE)

sid <- "H1_5"
sfe <- msfe@sfe_data[[sid]]

# ---- 1) Pull FGWC membership matrix ----
fg <- fgwc_list
U <- fg$membership
stopifnot(!is.null(U))

rownames(U) <- colnames(sfe)
colnames(U) <- paste0("cluster_", seq_len(ncol(U)))

# ---- 2) Coordinates + kNN graph (k=6 as in the graph step above) ----
# Prefer SpatialExperiment coords if present; otherwise fallback to geometry centroids
coords <- SpatialExperiment::spatialCoords(sfe)
coords <- coords[rownames(U), , drop = FALSE]

k_nn <- 6
knn <- FNN::get.knn(coords, k = k_nn)
# Build neighbour list for Moran's I and boundary detection
nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k_nn))
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

# ---- 3) Entropy per spot + visuals ----
entropy <- -rowSums(U * log(pmax(U, 1e-12)))
df_entropy <- data.frame(
    spot = rownames(U),
    x = coords[, 1],
    y = coords[, 2],
    entropy = entropy,
    row.names = NULL
)
data.table::fwrite(df_entropy, file.path(res_fld, "fgwc_entropy_per_spot.csv"))
p_ent <- ggplot(df_entropy, aes(x, y, colour = entropy)) +
    geom_point(size = 1) +
    coord_fixed() +
    scale_colour_viridis_c(name = "Membership\nentropy") +
    theme_void() +
    ggtitle(sprintf("FGWC uncertainty — %s", sid))
ggsave(file.path(res_fld, "fgwc_entropy.png"), p_ent, width = 5, height = 5, dpi = 300)

# ---- 4) Partition indices (PC, PE) ----
PC <- mean(rowSums(U^2))
PE <- -mean(rowSums(U * log(pmax(U, 1e-12))))

# ---- 5) Fuzzy silhouette (via fclust) ----
# Build a distance matrix among spots (Euclidean in the NMF embedding)
# Reuse the NMF embedding 'sfe_nmf' (spots x components)
if (is.list(sfe_nmf) && !is.null(sfe_nmf[[sid]])) {
    Xnmf <- as.matrix(sfe_nmf[[sid]])
} else {
    Xnmf <- as.matrix(sfe_nmf) # if already a matrix for H1_5
}
Xnmf <- Xnmf[rownames(U), , drop = FALSE]
D <- as.matrix(dist(Xnmf, method = "euclidean"))
# fclust::SIL.F expects U (rows=items, cols=clusters) and a dissimilarity matrix
fsil <- fclust::SIL.F(D, U, distance = TRUE)

# ---- 6) Ground truth labels + soft AUC (one-vs-rest) ----
gt <- gTruth_list[[sid]]

bar_col <- "Barcode"
annot_col <- "annotation"
gt_sub <- gt[, c(bar_col, annot_col)]
colnames(gt_sub) <- c("spot", "label")
gt_sub$spot <- as.character(gt_sub$spot)
gt_sub <- gt_sub[gt_sub$spot %in% rownames(U), , drop = FALSE]
gt_sub$label <- factor(gt_sub$label)
U_aligned <- U[gt_sub$spot, , drop = FALSE]

# Ensure alignment: rows of U_aligned correspond to gt_sub$spot
stopifnot(identical(rownames(U_aligned), gt_sub$spot))

# Build modelling frame and drop any rows with NA predictors or NA label
dfX <- as.data.frame(U_aligned)
dfX[] <- lapply(dfX, function(v) {
    v[!is.finite(v)] <- NA_real_
    v
}) # guard against infs
keep_all <- stats::complete.cases(dfX, gt_sub$label)

gt_use <- droplevels(gt_sub$label[keep_all])
X_use <- dfX[keep_all, , drop = FALSE]

# Train simple multinomial (one-vs-rest AUCs)
auc_tab <- data.frame(label = character(), AUC = numeric(), stringsAsFactors = FALSE)

for (lev in levels(gt_sub$label)[2:6]) {
    y <- as.integer(gt_use == lev)
    # Probabilistic score: best linear comb of memberships via logistic regression
    df_fit <- data.frame(y = y, X_use)
    fit <- glm(y ~ ., data = df_fit, family = binomial())
    score <- stats::predict(fit, type = "response")
    rocobj <- pROC::roc(y, score, quiet = TRUE)
    auc_tab <- rbind(auc_tab, data.frame(label = lev, AUC = as.numeric(pROC::auc(rocobj))))
}
data.table::fwrite(auc_tab, file.path(res_fld, "fgwc_soft_auc_by_label.csv"))

# ---- 7) Boundary set + localisation score + contrast index ----
# Define boundary spots: any spot with at least one neighbour of different label
lab_vec <- setNames(as.character(gt_sub$label), gt_sub$spot)
lab_vec <- lab_vec[rownames(U)]
# adjacency from nb
adj_list <- nb
is_boundary <- vapply(seq_along(adj_list), function(i) {
    me <- names(lab_vec)[i]
    ng <- adj_list[[i]]
    if (length(ng) == 0) {
        return(FALSE)
    }
    any(lab_vec[ng] != lab_vec[i], na.rm = TRUE)
}, logical(1))
names(is_boundary) <- rownames(U)

# Boundary localisation score (share of entropy mass near boundary)
H <- entropy
BLS <- sum(H[is_boundary], na.rm = TRUE) / sum(H, na.rm = TRUE)

# Boundary contrast: here we use the tumour-vs-rest score (GG4 Cribriform-vs-rest)
main_lev <- if ("GG4 Cribriform" %in% tolower(levels(gt_sub$label))) {
    levels(gt_sub$label)[which(tolower(levels(gt_sub$label)) == "GG4 Cribriform")]
} else {
    names(sort(table(gt_sub$label), decreasing = TRUE))[1]
}
# Fit logistic on memberships to get p(main_lev)
y_main <- as.integer(gt_sub$label == main_lev)
fit_main <- glm(y_main ~ ., data = data.frame(y_main = y_main, U_aligned), family = binomial())
p_main <- setNames(predict(fit_main, newdata = data.frame(U), type = "response"), rownames(U))

# Mean absolute difference of p across edges.
# We report overall and restricted to boundary edges
edge_pairs <- do.call(rbind, lapply(seq_along(nb), function(i) {
    if (length(nb[[i]]) == 0) {
        return(NULL)
    }
    cbind(i = i, j = nb[[i]])
}))
edge_pairs <- unique(t(apply(edge_pairs, 1, function(z) sort(z)))) # undirected
p_diff <- abs(p_main[edge_pairs[, 1]] - p_main[edge_pairs[, 2]])
edge_boundary <- is_boundary[edge_pairs[, 1]] | is_boundary[edge_pairs[, 2]]
BCI_all <- mean(p_diff, na.rm = TRUE)
BCI_boundary <- mean(p_diff[edge_boundary], na.rm = TRUE)

# ---- 8) Spatial autocorrelation of uncertainty (Moran's I) ----
moran_res <- spdep::moran.test(H, lw, zero.policy = TRUE)
moran_I <- if (!is.null(moran_res)) unname(moran_res$estimate[["Moran I statistic"]]) else NA_real_
moran_p <- if (!is.null(moran_res)) moran_res$p.value else NA_real_

# ---- 9) Save metrics summary ----
metrics <- data.frame(
    sample = sid,
    K = ncol(U),
    PC = PC,
    PE = PE,
    fuzzy_silhouette = fsil,
    BLS = BLS,
    BCI_all = BCI_all,
    BCI_boundary = BCI_boundary,
    moran_I_entropy = moran_I,
    moran_p_entropy = moran_p
)
data.table::fwrite(metrics, file.path(res_fld, "fgwc_metrics_summary.tsv"), sep = "\t")

# ==========================
# Export SFE for Python use
# ==========================
# We’ll write BOTH raw counts and log-normalised matrices,
# plus features, barcodes, coords and labels.

# matrices
m_counts <- SummarizedExperiment::assay(sfe, "counts") %>% as("dgCMatrix")
m_logcounts <- SummarizedExperiment::assay(sfe, "logcounts") %>% as("dgCMatrix")

Matrix::writeMM(m_counts, file.path(exp_fld, "counts.mtx"))
Matrix::writeMM(m_logcounts, file.path(exp_fld, "logcounts.mtx"))
R.utils::gzip(file.path(exp_fld, "counts.mtx"), overwrite = TRUE)
R.utils::gzip(file.path(exp_fld, "logcounts.mtx"), overwrite = TRUE)

# features (gene names)
feat <- data.frame(
    gene_id = rownames(sfe),
    gene_name = SummarizedExperiment::rowData(sfe)$gene_name,
    stringsAsFactors = FALSE
)
rd <- as.data.frame(SummarizedExperiment::rowData(sfe))
data.table::fwrite(feat, file.path(exp_fld, "features.tsv"), sep = "\t", col.names = FALSE)
R.utils::gzip(file.path(exp_fld, "features.tsv"), overwrite = TRUE)

# barcodes
bar <- data.frame(barcode = colnames(sfe))
data.table::fwrite(bar, file.path(exp_fld, "barcodes.tsv"), sep = "\t", col.names = FALSE)
R.utils::gzip(file.path(exp_fld, "barcodes.tsv"), overwrite = TRUE)

# coordinates
coords_df <- data.frame(spot = rownames(coords), x = coords[, 1], y = coords[, 2])
data.table::fwrite(coords_df, file.path(exp_fld, "coords.csv"))

# labels from ground truth
lab_out <- data.frame(spot = colnames(sfe), label = lab_vec)
data.table::fwrite(lab_out, file.path(exp_fld, "labels.csv"))

# record what we exported for convenience
writeLines(
    c(
        "counts.mtx.gz", "logcounts.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz",
        "coords.csv", "labels.csv"
    ),
    file.path(exp_fld, "manifest.txt")
)

message("FGWC metrics saved in ./results; export files for Python written to ./export")

