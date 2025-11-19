# Load R packages
library(STExplorer)
library(readr)
library(ggplot2)
library(svglite)
library(dplyr)
library(sf)
# library(magick)
library(grid)
library(gridExtra)
library(RColorBrewer)
# library(AnnotationHub)

## Data import and preparation of the MSFE object
msfe <- MetaSpatialFeatureExperiment()

workingDir <- c("/mnt/c/Users/Lefteris/Downloads/SpatialPCA_test/")
sampleDir <- c("/mnt/c/Users/Lefteris/Downloads/SpatialPCA_test/data/s151673/")

sampleNames <- c("s151673")
names(sampleDir) <- sampleNames

for (i in seq_along(sampleNames)) {
    msfe <- addSFE(
        msfe,
        read10xVisiumSFE(
            samples = sampleDir[i],
            sample_id = sampleNames[i],
            type = "sparse",
            data = "filtered",
            images = "lowres",
            style = "W",
            zero.policy = TRUE
        )
    )
}

ground_truth <- read_csv(file.path(sampleDir, "dlpfc_s151673_annotations.csv"))
gTruth_list <- list(s151673 = ground_truth)

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
    msfe@sfe_data[[i]] <- setQCthresh_LibSize(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.01)), max_t = quantile(msfe@sfe_data[[i]]@colData$sum, probs = c(.99)))
    msfe@sfe_data[[i]] <- setQCthresh_Mito(msfe@sfe_data[[i]], sample_id = TRUE, min_t = NA, max_t = quantile(msfe@sfe_data[[i]]@colData$subsets_mito_percent, probs = c(.99), na.rm = TRUE))
    msfe@sfe_data[[i]] <- setQCthresh_GenesExpr(msfe@sfe_data[[i]], sample_id = TRUE, min_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.01)), max_t = quantile(msfe@sfe_data[[i]]@colData$detected, probs = c(.99)))

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

## Calculate the mean of log counts over number of locations a gene is present
msfe <- perGeneLogMean(msfe, sample_id = TRUE)

## Zero expression genes
msfe <- setQCthresh_ZeroExpr(msfe, sample_id = TRUE)

## Lowly expressed (noise?!) genes
msfe <- setQCthresh_LowLogMean(msfe, threshold = 1, sample_id = TRUE)

## Remove mitochondrial and other genes

for (i in seq_along(sampleNames)) {
    msfe@sfe_data[[i]] <- setQCthresh_custom(msfe@sfe_data[[i]],
        MARGIN = 1,
        qcMetric = is_mito[[i]]
    )
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

top_hvgs <- getTopHighVarGenes(dec,
    var.field = "bio",
    prop = 0.1,
    var.threshold = 0,
    fdr.threshold = 0.1
)

## Add a neighbour graph using a weighted distance matrix
msfe <- addSpatialNeighGraphs(msfe,
    sample_id = TRUE,
    type = "knearneigh",
    style = "W",
    distMod = "raw",
    k = 6
)

## Calculate a simple distance matrix
msfe <- addDistMat(msfe, p = 2)


## Geographically Weighted PCA (GWPCA)

### Select parameters
sfe <- getSFE(msfe, "s151673")

## Get the gene names that are going to be evaluated
vars <- top_hvgs[["s151673"]]

## Set a fixed bandwidth
## bw is an important parameter as it defines the neighbourhood for which the
##  PCA will be calculated. The distance is measured in ultra-high resolution
##  image pixels. The default is 3x the diameter of the Visium spot. Make sure
##  to adjust it if it is too large or too small for your setting.
## Here the bw is selected as a multiple of the spot diameter. This is due to
##  using Visium as an example. You can just input a number.
bw <- 3 * sfe@metadata[["spotDiameter"]][["s151673"]][["spot_diameter_fullres"]]
## Set the number of components to be retained
k <- 20
## Set the kernel to be used
kernel <- "gaussian"
## Set the Minkowski distance power: p = 2 --> Euclidean
p <- 2
## The bandwidth is not adaptive because spots are in fixed positions
adaptive <- FALSE
## Cross-Validate GWPCA?
cv <- TRUE
## Calculate PCA scores?
scores <- TRUE
## Run a robust GWPCA?
robust <- FALSE
## Make a cluster for parallel computing (otherwise GWPCA is slow!)
my.cl <- makeClusterGWPCA(spec = 10, type = "PSOCK")

### Run GWPCA
pcagw <- gwpcaSTE(
    sfe = sfe,
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
    verbose = TRUE
)

## Extract leading genes
pcagw <- gwpca_LeadingGene(
    gwpca = pcagw,
    m_sfe = sfe,
    pc_nos = 1:4,
    type = "single",
    names = "gene_names"
)

pcagw <- gwpca_LeadingGene(
    gwpca = pcagw,
    m_sfe = sfe,
    pc_nos = 1:4,
    genes_n = 4,
    type = "multi",
    method = "membership",
    names = "gene_names"
)

#----------------------------------#
# Checkpoint
#----------------------------------#
saveRDS(pcagw, file = "./results/pcagw.rds")
str(pcagw)


#------------------------------------------------------------------------------#
# GWPCA vs SpatialPCA
#------------------------------------------------------------------------------#
# Fix analysis constants and validate parameters -------------------------------
## Dims from pcagw object ----
dl <- dim(pcagw$loadings) # [n_spots, n_genes, k]
n <- dl[1]
p <- dl[2]
k0 <- dl[3]

## Number of GWPCA components (local PCs kept per neighbourhood) ----
# We prefer the value actually used by GWPCA
k <- as.integer(pcagw$GW.arguments$k)

## Size of the union dictionary (should be >= k: 20-30) ----
K <- 150

message(sprintf("Using k = %d (locals) and K = %d (dictionary).", k, K))

## Gene set order and spot IDs
genes <- dimnames(pcagw$loadings)[[2]]
spot_ids <- dimnames(pcagw$loadings)[[1]]

## Structural checks on eigenvalues and scores matrices ----
# local principal variances (eigenvalues): should be n x k
if (!is.matrix(pcagw$local.PV) || !all(dim(pcagw$local.PV) == c(n, k))) {
    stop("pcagw$local.PV must be an n×k matrix of local principal variances.")
}

## Export spatial coordinates ----
coords <- pcagw$SDF@coords
# Ensure the row order matches spot_ids (should be... but double-check)
if (!identical(rownames(coords), spot_ids)) {
    warning("Row names of coords don’t match spot_ids; will align by row order.")
}
## Get global PCA info ----
# (not used as reference, but handy for QC later on)
pca_global <- list(
    loadings = pcagw$pca$loadings, # genes × PCs
    sdev     = pcagw$pca$sdev,
    scores   = pcagw$pca$scores, # spots × PCs
    center   = pcagw$pca$center,
    scale    = pcagw$pca$scale
)

## Set a small config list ----
# We will use it to get parameters for later steps
gw_cfg <- list(
    n = n,
    p = p,
    k = k,
    K = K,
    genes = genes,
    spot_ids = spot_ids,
    coords = coords
)

# Step 1: Gather local bases (V_i) and eigenvalues (Lambda_i) ------------------
## Export parameters from gw_cfg ----
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k

## Extract local loading matrices Vi (genes × k for each spot) ----
# This is already stored in pcagw object as a 3-D array: [spot, gene, pc]
V_arr <- pcagw$loadings[, , 1:k, drop = FALSE] # dims: n × p × k

# Add dimnames (it should have them already, but let's be sure)
dimnames(V_arr) <- list(
    spot  = gw_cfg$spot_ids,
    gene  = gw_cfg$genes,
    pc    = paste0("PC", seq_len(k))
)

## Extract local eigenvalues (principal variances) Lambda_i ----
# Matrix of size n × k
# Each row corresponds to a spot
# Each column to a PC
EV <- as.matrix(pcagw$local.PV) # n × k
colnames(EV) <- paste0("PC", seq_len(k))
rownames(EV) <- gw_cfg$spot_ids

# Sanity check: eigenvalues should be ≥ 0
neg_mask <- EV < 0
sum(neg_mask)

# Precompute sqrt(lambda) per spot/component
sqrtEV <- sqrt(EV)

## Calculate Orthonormality diagnostics for V_i ----
# For each spot i, check t(V_i) %*% V_i ≈ I_k
orth_diag <- function(Vi) {
    G <- crossprod(Vi) # k × k Gram matrix
    diag_dev <- diag(G) - 1
    off_diag <- G - diag(diag(G))
    c(
        max_abs_diag_dev = max(abs(diag_dev)),
        max_abs_off_diag = max(abs(off_diag)),
        frobenius_I_diff = sqrt(sum((G - diag(k))^2))
    )
}

# Loop over spots ----
orth_stats <- t(vapply(seq_len(n), function(i) {
    Vi <- matrix(V_arr[i, , , drop = TRUE],
        nrow = p,
        ncol = k,
        dimnames = list(gw_cfg$genes, paste0("PC", seq_len(k)))
    )
    orth_diag(Vi)
}, numeric(3L)))
rownames(orth_stats) <- gw_cfg$spot_ids

# Summaries to print
orth_summary <- apply(orth_stats, 2, function(x) {
    c(
        max = max(x, na.rm = TRUE),
        p99 = quantile(x, 0.99, na.rm = TRUE),
        median = median(x, na.rm = TRUE)
    )
})

message("Local basis orthonormality diagnostics (t(V_i)V_i vs I_k):")
print(round(orth_summary, 6))

# Eigenvalue sanity summaries
ev_summary <- cbind(
    per_pc_min = apply(EV, 2, min, na.rm = TRUE),
    per_pc_q1  = apply(EV, 2, quantile, probs = 0.25, na.rm = TRUE),
    per_pc_med = apply(EV, 2, median, na.rm = TRUE),
    per_pc_q3  = apply(EV, 2, quantile, probs = 0.75, na.rm = TRUE),
    per_pc_max = apply(EV, 2, max, na.rm = TRUE)
)
message("Local eigenvalue (principal variance) summaries by component:")
print(round(ev_summary, 4))

# Package outputs for the next steps
gw_locals <- list(
    V_arr = V_arr, # n × p × k
    EV = EV, # n × k
    sqrtEV = sqrtEV, # n × k
    orth_stats = orth_stats
)

# Check no NA/NaN in critical structures
stopifnot(!any(!is.finite(V_arr)))
stopifnot(!any(!is.finite(EV)))

str(gw_locals, max.level = 1)

# Step 2: Build the union-of-locals reference V* (dictionary) ------------------

# Reset inputs to make sure everything is set
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k
K <- min(max(150L, k + 60L), p)
V_arr <- gw_locals$V_arr # n × p × k   (spot, gene, pc)
EV <- gw_locals$EV # n × k
sqrtEV <- gw_locals$sqrtEV # n × k

## Stack–SVD union (variance-weighted by sqrtEV)
build_dictionary_stacksvd <- function(V_arr, sqrtEV, K, weight = TRUE) {
    n <- dim(V_arr)[1]
    p <- dim(V_arr)[2]
    k <- dim(V_arr)[3]
    U <- matrix(NA_real_, nrow = p, ncol = n * k)
    col_idx <- 1L
    for (i in seq_len(n)) {
        Vi <- V_arr[i, , , drop = TRUE] # p × k
        Ui <- if (weight) Vi %*% diag(sqrtEV[i, ], k) else Vi
        U[, col_idx:(col_idx + k - 1L)] <- Ui
        col_idx <- col_idx + k
    }
    sv <- svd(U, nu = min(K, p), nv = 0)
    Vstar <- sv$u[, seq_len(min(K, ncol(sv$u))), drop = FALSE] # p × K_new
    # Enforce orthonormality again
    Vstar <- qr.Q(qr(Vstar))
    rownames(Vstar) <- gw_cfg$genes
    colnames(Vstar) <- paste0("D", seq_len(ncol(Vstar)))
    Vstar
}

# Build the dictionary
Vstar <- build_dictionary_stacksvd(V_arr, sqrtEV, K, weight = TRUE)
gw_dict$Vstar <- Vstar
gw_cfg$K <- ncol(Vstar) # update
dict_info <- list(method = "stacksvd", eig_values = dict$eig_values[seq_len(K)])

## Sanity checks on V*
# Orthonormal columns? (they should be...)
Gstar <- crossprod(Vstar)
diag_dev <- max(abs(diag(Gstar) - 1))
off_dev <- max(abs(Gstar - diag(diag(Gstar))))
message(sprintf("Dictionary V*: max |diag-1| = %.3e ; max |offdiag| = %.3e", diag_dev, off_dev))

# Dimensions and names
stopifnot(nrow(Vstar) == p, ncol(Vstar) == K)
stopifnot(identical(rownames(Vstar), gw_cfg$genes))

## Package for next step
gw_dict <- list(Vstar = Vstar, info = dict_info, method = "stacksvd")
str(gw_dict, max.level = 1)

# Step 3: Principal-angle diagnostics (local vs dictionary) --------------------
# How well does each local subspace S_i overlaps with dictionary space S_star?

# Reset input to ensure everything is present and set correctly
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k
V_arr <- gw_locals$V_arr # n × p × k
Vstar <- gw_dict$Vstar # p × k

## Helper fucntion ----
# For a given Vi (p×k), compute cosines of principal angles and proj distance
angles_one <- function(Vi, Vstar) {
    # SVD of k×k cross-basis matrix gives cosines of principal angles
    M <- crossprod(Vi, Vstar) # k × k
    sv <- svd(M, nu = 0, nv = 0)$d
    # Numerical safety: clamp to [0,1]
    cos_t <- pmin(pmax(sv, 0), 1)
    # Projection distance ||P_i - P*||_F = sqrt( k+K - 2 * sum(cos^2 theta) )
    proj_dist <- sqrt(k + ncol(Vstar) - 2 * sum(cos_t^2))
    list(cos = cos_t, proj_dist = proj_dist)
}

## Loop over spots ----
cosines_mat <- matrix(NA_real_,
    nrow = n, ncol = k,
    dimnames = list(gw_cfg$spot_ids, paste0("cos_theta_", seq_len(k)))
)
proj_dist <- numeric(n)
for (i in seq_len(n)) {
    Vi <- matrix(V_arr[i, , , drop = TRUE], nrow = p, ncol = k)
    res <- angles_one(Vi, Vstar)
    cosines_mat[i, ] <- res$cos
    proj_dist[i] <- res$proj_dist
}

## Summaries per spot ----
angle_summ_df <- data.frame(
    spot_id   = gw_cfg$spot_ids,
    cos_min   = apply(cosines_mat, 1, min, na.rm = TRUE),
    cos_q1    = apply(cosines_mat, 1, quantile, probs = 0.25, na.rm = TRUE),
    cos_med   = apply(cosines_mat, 1, median, na.rm = TRUE),
    cos_q3    = apply(cosines_mat, 1, quantile, probs = 0.75, na.rm = TRUE),
    cos_max   = apply(cosines_mat, 1, max, na.rm = TRUE),
    cos_mean  = rowMeans(cosines_mat, na.rm = TRUE),
    proj_dist = proj_dist,
    row.names = gw_cfg$spot_ids
)

## Global summaries ----
qs <- function(x) {
    quantile(x, probs = c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
}
cat("\nPrincipal-angle cosines — per-spot mean(cos θ):\n")
print(round(qs(angle_summ_df$cos_mean), 3))

cat("\nProjection distance ||P_i - P*||_F:\n")
print(round(qs(angle_summ_df$proj_dist), 3))

## Optional thresholds to flag weak overlap
thr_warn_mean_cos <- 0.6 # anything below this is dubious
thr_warn_proj <- sqrt(k + K - 2 * (k * 0.6^2)) # equivalent in proj space

n_low_cos <- sum(angle_summ_df$cos_mean < thr_warn_mean_cos)
n_high_pd <- sum(angle_summ_df$proj_dist > thr_warn_proj)

message(sprintf("Spots with mean(cos θ) < %.2f: %d / %d", thr_warn_mean_cos, n_low_cos, n))
message(sprintf("Spots with proj_dist > %.3f: %d / %d", thr_warn_proj, n_high_pd, n))

## Keep everything for later steps (alignment and masking/weighting if needed)
gw_angles <- list(
    cosines_by_spot = cosines_mat, # n × k cosines of principal angles
    per_spot        = angle_summ_df, # tidy per-spot summaries
    Vstar_k         = Vstar,
    thresholds      = list(mean_cos = thr_warn_mean_cos, proj_dist = thr_warn_proj)
)

str(gw_angles, max.level = 1)


stopifnot(identical(rownames(gw_dict$Vstar), gw_cfg$genes)) # already TRUE earlier
## Also verify local slices carry same genes:
stopifnot(identical(dimnames(gw_locals$V_arr)[[2]], gw_cfg$genes))

# ==============================================================================#
# =======================Steps for supervised tasks=============================#
# Although not used for kNN + ARI, we compute these steps to plot some QC metrics
# Step 4: Per-spot embedding (best-k match) + standardisation ------------------

## Inputs: gw_cfg (n,p,k,K), gw_locals$V_arr (n × p × k), gw_dict$Vstar (p × K),
##         pcagw$gwpca.scores (list of n matrices n × k)

# Reset input to ensure everything is present and set correctly
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k
V_arr <- gw_locals$V_arr # n × p × k (correctly sliced; spots;genes×PCs)
Vstar <- gw_dict$Vstar # p × K (K chosen large; columns orthonormal)
K_eff <- ncol(Vstar)
scores <- pcagw$gwpca.scores

## Helper: focal local scores for spot i (length k)
get_focal_scores <- function(scores_list, i, n, k) {
    Si <- scores_list[[i]]
    stopifnot(is.matrix(Si), all(dim(Si) == c(n, k)))
    Si[i, 1:k, drop = TRUE]
}

## Build Z via per-spot best-k dictionary matching ----
Z <- matrix(0,
    nrow = n, ncol = K_eff,
    dimnames = list(gw_cfg$spot_ids, colnames(Vstar))
)
subspace_energy <- numeric(n) # sum of top-k cos^2 for this spot (≤ k), QC

# Step 5: Align each local basis to the dictionary -----------------------------
# (orthogonal Procrustes)
for (i in seq_len(n)) {
    Vi <- V_arr[i, , , drop = TRUE] # p × k (DO NOT rewrap)
    Mi <- crossprod(Vi, Vstar) # k × K_eff
    sv <- svd(Mi) # Mi = U Σ V^T
    d <- pmin(pmax(sv$d, 0), 1)
    subspace_energy[i] <- sum(d[seq_len(k)]^2)

    s_i <- get_focal_scores(scores, i, n, k) # k×1
    z_local <- crossprod(sv$u[, seq_len(k), drop = FALSE], s_i) # rotate scores
    Wk <- sv$v[, seq_len(k), drop = FALSE] # selected dict axes
    Z[i, ] <- as.numeric(Wk %*% z_local) # place into K-dim frame
}

## QC: subspace energy ----
# subspace energy should be decently high now (ideally ≳12–15 of max k=20)
qs <- function(x) quantile(x, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
cat("\nSubspace energy (sum cos^2 θ over top-k dict axes):\n")
print(round(qs(subspace_energy), 3))

# Step 6: Produce per-spot, comparable coordinates -----------------------------
# (the embedding)
##  Column standardise Z (z-scores) ----
# (For fair classifiers/ARI)
Z_sd <- apply(Z, 2, sd, na.rm = TRUE)
keep_cols <- Z_sd > 0
Z_std <- base::scale(Z[, keep_cols, drop = FALSE])

## Compress to D dims by variance across spots ----
# (keeps the most useful dict axes)
D <- min(40L, ncol(Z_std)) # 20–40 is a sensible target
var_order <- order(apply(Z_std, 2, var, na.rm = TRUE), decreasing = TRUE)
sel <- var_order[seq_len(D)]
Z_compact <- Z_std[, sel, drop = FALSE]

## Save artefacts for downstream steps ----
gw_embed <- list(
    Z_raw = Z, # n × K_eff
    Z_std = Z_std, # z-scored, columns with non-zero sd
    Z_compact = Z_compact, # z-scored, top-D variance axes
    dict_axes = colnames(Z_std),
    selected_axes = colnames(Z_compact),
    subspace_energy = subspace_energy,
    mask_low = NULL,
    K_eff = K_eff,
    D = D
)

str(gw_embed, max.level = 1)

## Plotting step ---------------------------------------------------------------
## ---- Normalised subspace energy ---------------------------------------------
k <- gw_cfg$k
energy_norm <- as.numeric(subspace_energy) / k
stopifnot(length(energy_norm) == length(gw_cfg$spot_ids))
df_energy <- data.frame(
    spot_id = gw_cfg$spot_ids,
    energy_norm = energy_norm,
    stringsAsFactors = FALSE
)

## ---- 1) Histogram -----------------------------------------------------------
med <- median(df_energy$energy_norm, na.rm = TRUE)
iqr <- quantile(df_energy$energy_norm, probs = c(0.25, 0.75), na.rm = TRUE)

p_hist <- ggplot(df_energy, aes(x = energy_norm)) +
    geom_histogram(bins = 40, colour = "white") +
    geom_vline(xintercept = med, linetype = "dashed") +
    annotate("text",
        x = med, y = Inf, vjust = 1.5,
        label = sprintf("median = %.3f\nIQR = [%.3f, %.3f]", med, iqr[1], iqr[2])
    ) +
    labs(
        x = expression(paste(sum(cos^2), " / k")),
        y = "Count",
        title = "Normalised subspace energy across spots"
    ) +
    theme_classic()

print(p_hist)

## ---- 2) Spatial map (polygons) ----------------------------------------------
make_poly_map <- function(pcagw, df_energy) {
    stopifnot(!is.null(pcagw$geometry))
    sf_polys <- pcagw$geometry # sf/data.frame with 'geometry' col
    # Ensure we have a spot_id column to join on
    if (is.null(sf_polys$spot_id)) {
        sf_polys$spot_id <- rownames(sf_polys)
    }
    sf_polys %>%
        left_join(df_energy, by = "spot_id") %>%
        ggplot() +
        geom_sf(aes(fill = energy_norm), colour = NA) +
        viridis::scale_fill_viridis(direction = 1, na.value = "grey80") +
        labs(
            fill = expression(paste(sum(cos^2), " / k")),
            title = "Normalised subspace energy (spatial map)"
        ) +
        theme_void()
}

p_map <- make_poly_map(pcagw, df_energy)

print(p_map)

## ---- Save to files -------------------------------------------
ggsave("results/subspace_energy_histogram.pdf", p_hist, width = 6, height = 4)
ggsave("results/subspace_energy_spatial_map.pdf", p_map, width = 6, height = 5)



# ==============================================================================#
# ======================Steps used for unsupervised tasks=======================#
## ========= Step 4 & 6: dense shared embedding for clustering =================
## Inputs: gw_cfg (n,p,k), gw_locals$V_arr (n×p×k),
##         gw_dict$Vstar (p×K; orthonormal),
##         pcagw$gwpca.scores (list of n matrices n×k)

# Setup
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k
V_arr <- gw_locals$V_arr # n × p × k
Vstar <- gw_dict$Vstar # p × K
K_eff <- ncol(Vstar)
scores <- pcagw$gwpca.scores

# Helper: focal local scores (length k) for spot i
get_focal_scores <- function(scores_list, i, n, k) {
    Si <- scores_list[[i]]
    stopifnot(is.matrix(Si), all(dim(Si) == c(n, k)))
    Si[i, 1:k, drop = TRUE]
}

# 4) Dense shared embedding: Y_i = (V_i^T V*)^T s_i  (spots × K)
Y <- matrix(NA_real_,
    nrow = n, ncol = K_eff,
    dimnames = list(gw_cfg$spot_ids, colnames(Vstar))
)
for (i in seq_len(n)) {
    Vi <- V_arr[i, , , drop = TRUE] # p × k  (DO NOT rewrap)
    Mi <- crossprod(Vi, Vstar) # k × K
    s_i <- get_focal_scores(scores, i, n, k) # k × 1
    Y[i, ] <- as.numeric(t(Mi) %*% s_i) # K × 1
}

# 5) Row-wise L2 normalisation (good for cosine kNN graphs)
row_norms <- sqrt(rowSums(Y^2))
row_norms[row_norms == 0] <- 1
Y_l2 <- Y / row_norms

# 6) Column standardise (z-score) and pick top-D by variance (as before)
Z_sd <- apply(Y_l2, 2, sd, na.rm = TRUE)
keep_cols <- Z_sd > 0
Z_std <- base::scale(Y_l2[, keep_cols, drop = FALSE])

D <- min(40L, ncol(Z_std)) # 20–40 is sensible
var_order <- order(apply(Z_std, 2, var, na.rm = TRUE), decreasing = TRUE)
sel <- var_order[seq_len(D)]
Z_compact <- Z_std[, sel, drop = FALSE]

# Save artefacts (same names as before so downstream code needn’t change)
gw_embed <- list(
    Z_raw = Y_l2, # n × K (dense, L2-normalised rows)
    Z_std = Z_std, # column z-scored
    Z_compact = Z_compact, # top-D by across-spot variance
    dict_axes = colnames(Z_std),
    selected_axes = colnames(Z_compact),
    subspace_energy = NULL,
    mask_low = NULL,
    K_eff = K_eff,
    D = D
)
str(gw_embed, max.level = 1)


# Step 7: Fraction of local variance outside the dictionary --------------------

# Inputs:
# gw_locals$V_arr (n × p × k), gw_locals$EV (n × k), gw_dict$Vstar (p × K)
n <- gw_cfg$n
k <- gw_cfg$k
V_arr <- gw_locals$V_arr
EV <- gw_locals$EV
Vstar <- gw_dict$Vstar # columns orthonormal; K can be large

# For spot i:
#   Λ_i = diag(EV[i, ])
#   M_i = V_i^T V*   (k × K)
#   B_i = M_i M_i^T  (k × k) = V_i^T (V* V*^T) V_i
#   captured variance = tr(Λ_i B_i)
#   total local variance = tr(Λ_i) = sum(EV[i, ])
#   f_out = 1 - captured/total

f_out <- numeric(n)
captured <- numeric(n)
total <- rowSums(EV)

for (i in seq_len(n)) {
    Vi <- V_arr[i, , , drop = TRUE] # p × k
    Mi <- crossprod(Vi, Vstar) # k × K
    Bi <- Mi %*% t(Mi) # k × k
    # tr(Λ_i B_i) = sum_j λ_ij * B_i[j,j] if Λ_i is diagonal
    captured[i] <- sum(EV[i, ] * diag(Bi))
    # numerical guard: cap to [0, total]
    if (captured[i] < 0) captured[i] <- 0
    if (captured[i] > total[i]) captured[i] <- total[i]
    f_out[i] <- 1 - captured[i] / total[i]
}

# Summaries
qs <- function(x) quantile(x, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
cat("Captured fraction (1 - f_out) summary:\n")
print(round(qs(1 - f_out), 3))
cat("Outside fraction (f_out) summary:\n")
print(round(qs(f_out), 3))

# Keep for plotting / reporting
gw_outside <- data.frame(
    spot_id = gw_cfg$spot_ids,
    frac_captured = pmin(pmax(1 - f_out, 0), 1),
    frac_outside = pmin(pmax(f_out, 0), 1)
)

## Plot step -------------------------------------------------------------------
## ---------- Prep: tidy data ----------
df_out <- gw_outside
stopifnot(all(c("spot_id", "frac_outside", "frac_captured") %in% names(df_out)))

## ---------- Generic histogram helper ----------
plot_hist_frac <- function(df, var = c("frac_outside", "frac_captured"),
                           bins = 40, title = NULL, xlab = NULL) {
    var <- match.arg(var)
    xlab <- xlab %||% if (var == "frac_outside") {
        "Fraction of local variance outside dictionary (f_out)"
    } else {
        "Fraction of local variance captured (1 - f_out)"
    }
    title <- title %||% paste0(
        if (var == "frac_outside") "Outside fraction" else "Captured fraction",
        " across spots"
    )

    v <- df[[var]]
    med <- median(v, na.rm = TRUE)
    iqr <- quantile(v, probs = c(0.25, 0.75), na.rm = TRUE)

    ggplot(df, aes(x = .data[[var]])) +
        geom_histogram(bins = bins, colour = "white") +
        geom_vline(xintercept = med, linetype = "dashed") +
        annotate("text",
            x = med, y = Inf, vjust = 1.5,
            label = sprintf("median = %.3f\nIQR = [%.3f, %.3f]", med, iqr[1], iqr[2])
        ) +
        labs(x = xlab, y = "Count", title = title) +
        theme_classic()
}

## ---------- Generic spatial map helpers (polygons preferred) ----------
.make_poly_map_var <- function(pcagw, df, var, legend_lab) {
    stopifnot(!is.null(pcagw$geometry))
    sf_polys <- pcagw$geometry
    if (is.null(sf_polys$spot_id)) sf_polys$spot_id <- rownames(sf_polys)
    sf_polys %>%
        left_join(df, by = "spot_id") %>%
        ggplot() +
        geom_sf(aes(fill = .data[[var]]), colour = NA) +
        viridis::scale_fill_viridis(direction = 1, na.value = "grey80") +
        labs(fill = legend_lab, title = legend_lab) +
        theme_void()
}

plot_spatial_frac <- function(pcagw, df, var = c("frac_outside", "frac_captured")) {
    var <- match.arg(var)
    legend_lab <- if (var == "frac_outside") {
        "Fraction outside dictionary (f_out)"
    } else {
        "Fraction captured (1 - f_out)"
    }

    .make_poly_map_var(pcagw, df, var, legend_lab)
}

`%||%` <- function(a, b) if (is.null(a)) b else a # tiny util

## 1) Outside fraction (f_out) -------------------------------------------------
# A map of f_out (where GWPCA finds signal the dictionary doesn’t capture).
p_hist_out <- plot_hist_frac(df_out,
    var = "frac_outside",
    title = "Outside fraction across spots",
    xlab = "Fraction of local variance outside dictionary (f_out)"
)
p_map_out <- plot_spatial_frac(pcagw, df_out, var = "frac_outside")

print(p_hist_out)
print(p_map_out)

## 2) Captured fraction (1 - f_out) --------------------------------------------
p_hist_cap <- plot_hist_frac(df_out,
    var = "frac_captured",
    title = "Captured fraction across spots",
    xlab = "Fraction of local variance captured (1 - f_out)"
)
p_map_cap <- plot_spatial_frac(pcagw, df_out, var = "frac_captured")

print(p_hist_cap)
print(p_map_cap)

## Save to files
ggsave("results/f_out_histogram.pdf", p_hist_out, width = 6, height = 4)
ggsave("results/f_out_spatial_map.pdf", p_map_out, width = 6, height = 5)
ggsave("results/captured_histogram.pdf", p_hist_cap, width = 6, height = 4)
ggsave("results/captured_spatial_map.pdf", p_map_cap, width = 6, height = 5)


## ---- GWPCA grid: K, D, k_graph, smoothing, energy weighting -----------------
## Identify the combination of parameters that optimise GWPCA's ARI score
# Helpers
get_focal_scores <- function(scores_list, i, n, k) {
    Si <- scores_list[[i]]
    stopifnot(is.matrix(Si), all(dim(Si) == c(n, k)))
    Si[i, 1:k, drop = TRUE]
}

build_dense_Y <- function(V_arr, Vstar, scores, EV = NULL, energy = FALSE) {
    n <- dim(V_arr)[1]
    k <- dim(V_arr)[3]
    K <- ncol(Vstar)
    Y <- matrix(NA_real_, n, K, dimnames = list(rownames(V_arr), colnames(Vstar)))
    for (i in seq_len(n)) {
        Vi <- V_arr[i, , , drop = TRUE] # p×k
        Mi <- crossprod(Vi, Vstar) # k×K
        s <- get_focal_scores(scores, i, n, k)
        if (energy && !is.null(EV)) s <- as.numeric(diag(sqrt(EV[i, ]), k, k) %*% s)
        Y[i, ] <- as.numeric(t(Mi) %*% s)
    }
    Y
}

prep_Z <- function(Y, D) {
    sdv <- apply(Y, 2, sd, na.rm = TRUE)
    keep <- sdv > 0
    Z <- base::scale(Y[, keep, drop = FALSE])
    if (ncol(Z) > D) {
        v <- apply(Z, 2, var, na.rm = TRUE)
        sel <- order(v, decreasing = TRUE)[seq_len(D)]
        Z <- Z[, sel, drop = FALSE]
    }
    Z
}

smooth_spatial <- function(Z, coords, k = 12L, iters = 1L) {
    N <- nrow(Z)
    nn <- get.knn(coords, k = k)
    # Row-normalised weight matrix (sparse-ish in practice)
    Zs <- Z
    for (t in seq_len(iters)) {
        Znew <- Zs * 0
        for (i in seq_len(N)) {
            nbrs <- nn$nn.index[i, ]
            w <- nn$nn.dist[i, ]
            s <- median(w) + 1e-8
            w <- exp(-(w^2) / (2 * s^2))
            w <- w / sum(w)
            Znew[i, ] <- colSums(Zs[nbrs, , drop = FALSE] * w)
        }
        Zs <- Znew
    }
    Zs
}

knn_graph <- function(emb, k = 15L) {
    idx <- get.knn(emb, k = k)$nn.index
    edges <- cbind(rep(seq_len(nrow(emb)), each = k), as.vector(t(idx)))
    simplify(graph_from_edgelist(edges, directed = FALSE))
}

run_leiden_knn_ari <- function(emb, labels, k_graph = 15L, res = seq(0.2, 2.0, by = 0.2), seed = 2023) {
    g <- knn_graph(emb, k = k_graph)
    do.call(rbind, lapply(res, function(r) {
        set.seed(seed)
        cl <- cluster_leiden(g, resolution_parameter = r, objective_function = "modularity")
        mem <- membership(cl)
        data.frame(
            resolution = r, n_clusters = length(unique(mem)),
            ARI = adjustedRandIndex(as.integer(factor(labels)), as.integer(factor(mem)))
        )
    }))
}

# Inputs
n <- gw_cfg$n
p <- gw_cfg$p
k <- gw_cfg$k
V_arr <- gw_locals$V_arr
EV <- gw_locals$EV
scores <- pcagw$gwpca.scores
XY <- pcagw$SDF@coords[gw_cfg$spot_ids, , drop = FALSE]
labels <- labels # ground truth

# Grid
K_vals <- c(120L, 180L, 240L)
D_vals <- c(40L, 60L, 80L, 100L)
k_graph_vals <- c(10L, 15L, 30L)
iters_vals <- c(0L, 1L, 2L)
energy_vals <- c(FALSE, TRUE)

results <- list()
ctr <- 1L
for (K_new in K_vals) {
    # Rebuild dictionary Vstar for this K from current V_arr (stack-SVD approach used before)
    Vstar <- {
        # stack columns Vi * sqrtEV for all spots, then SVD → first K_new left singular vectors
        U <- matrix(NA_real_, nrow = p, ncol = n * k)
        col <- 1L
        for (i in seq_len(n)) {
            Vi <- V_arr[i, , , drop = TRUE]
            Ui <- Vi %*% diag(sqrt(EV[i, ]), k, k)
            U[, col:(col + k - 1L)] <- Ui
            col <- col + k
        }
        sv <- svd(U, nu = min(K_new, p), nv = 0)
        Q <- qr.Q(qr(sv$u[, seq_len(min(K_new, ncol(sv$u))), drop = FALSE]))
        rownames(Q) <- gw_cfg$genes
        colnames(Q) <- paste0("D", seq_len(ncol(Q)))
        Q
    }

    for (energy in energy_vals) {
        Y <- build_dense_Y(V_arr, Vstar, scores, EV, energy = energy)

        for (D in D_vals) {
            Z <- prep_Z(Y, D)
            for (it in iters_vals) {
                Zs <- if (it > 0) smooth_spatial(Z, XY, k = 12L, iters = it) else Z
                for (kg in k_graph_vals) {
                    rd <- run_leiden(Zs, labels, k_graph = kg)
                    rd$K <- K_new
                    rd$D <- D
                    rd$iters <- it
                    rd$k_graph <- kg
                    rd$energy <- energy
                    results[[ctr]] <- rd
                    ctr <- ctr + 1L
                }
            }
        }
    }
}

grid_res <- bind_rows(results)

best <- grid_res %>%
    group_by(K, D, iters, k_graph, energy) %>%
    summarise(ARI_mean = max(ARI), .groups = "drop") %>%
    arrange(desc(ARI_mean)) %>%
    slice_head(n = 10)

print(best)

#------------------------------------------------------------------------------#
# Winning combination --> these are used in the next script
# Freeze winning settings from grid
K_fixed <- 180L
D_fixed <- 100L
iters_fixed <- 2L
k_graph_fixed <- 10L
energy_fixed <- FALSE

# Rebuild dictionary V* with K_fixed, recompute dense Y, z-score columns,
# select top D_fixed, then apply 2-hop spatial smoothing (the same code used in grid).

# Rebuild dictionary Vstar for this K from current V_arr (stack-SVD approach used before)
Vstar <- {
    # stack columns Vi * sqrtEV for all spots, then SVD → first K_new left singular vectors
    U <- matrix(NA_real_, nrow = p, ncol = n * k)
    col <- 1L
    for (i in seq_len(n)) {
        Vi <- V_arr[i, , , drop = TRUE]
        Ui <- Vi %*% diag(sqrt(EV[i, ]), k, k)
        U[, col:(col + k - 1L)] <- Ui
        col <- col + k
    }
    sv <- svd(U, nu = min(K_fixed, p), nv = 0)
    Q <- qr.Q(qr(sv$u[, seq_len(min(K_fixed, ncol(sv$u))), drop = FALSE]))
    rownames(Q) <- gw_cfg$genes
    colnames(Q) <- paste0("D", seq_len(ncol(Q)))
    Q
}

Y <- build_dense_Y(V_arr, Vstar, scores, EV, energy = energy_fixed)
Z <- prep_Z(Y, D_fixed)
Z_final <- smooth_spatial(Z, XY, k = 12L, iters = iters_fixed)
rd <- run_leiden_knn_ari(Z_final, labels, k_graph = k_graph_fixed)
rd$K <- K_fixed
rd$D <- D_fixed
rd$iters <- iters_fixed
rd$k_graph <- k_graph_fixed
rd$energy <- energy_fixed

# Save the final embedding for comparison script:
gw_embed_fixed <- list(
    Z_compact = Z_final, # n × D_fixed, z-scored columns, after smoothing
    D = D_fixed
)
saveRDS(gw_embed_fixed, file = "./results/gw_embed_fixed.rds")
saveRDS(k_graph_fixed, file = "./results/k_graph_fixed.rds")
