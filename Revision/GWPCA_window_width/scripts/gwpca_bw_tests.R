#------------------------------------------------------------------------------#
## To load the packages use:
library(STExplorer)
library(SpatialFeatureExperiment)
library(tidyverse)
library(sf)
library(spdep)
library(GWmodel)
library(ggplot2)


#------------------------------------------------------------------------------#
# load and use the output from the liver analysis
load("./Liver_FGWC_DEGs/data/201124_EZ.RData")


#------------------------------------------------------------------------------#
# Helper functions

# 1) Extract PCA scores
# ref_scores: length n vector (e.g., global PCA PC1 scores ordered by spots)
# S_list: your list gwpca$gwpca.scores (length n), each an n x p matrix
# returns a length-n vector of aligned self-scores for the chosen PC
extract_aligned_pc_selfscores <- function(S_list, pc = 1, ref_scores) {
  n <- length(S_list)
  out <- numeric(n)
  for (i in seq_len(n)) {
    M <- S_list[[i]] # n x p scores from local model i
    s_vec <- M[, pc] # scores for PC across all spots
    sgn <- if (all(is.finite(s_vec))) {
      cs <- suppressWarnings(stats::cor(s_vec, ref_scores, use = "pairwise.complete.obs"))
      if (is.na(cs) || cs == 0) 1 else sign(cs)
    } else {
      1
    }
    # self-score = score of spot i under model i, sign-aligned
    out[i] <- sgn * M[i, pc]
  }
  out
}

# 2) Neighbourhood + Moran’s I
moransI_from_scores <- function(x, coords = spatialCoords(sfe), k = 6, style = "W") {
  stopifnot(is.numeric(x))
  C <- as.matrix(coords)
  # k-NN graph (robust on regular Visium grids)
  knn <- spdep::knearneigh(coords, k = k)
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)

  mt <- spdep::moran.test(x, lw, zero.policy = TRUE)
  unname(mt$estimate[["Moran I statistic"]])
}


# 3) CV error per fit
get_cv_error <- function(gwpca_obj) {
  # Try common fields; return NA_real_ if not found.
  if (!is.null(gwpca_obj$CV$CV)) {
    return(as.numeric(gwpca_obj$CV$CV))
  }
  NA_real_
}

summarise_cv <- function(df) {
  df %>%
    group_by(bw_mult) %>%
    summarise(
      n_spots = n(),
      cv_med  = median(cv_error, na.rm = TRUE),
      cv_trim = mean(cv_error, trim = 0.10, na.rm = TRUE),
      cv_q10  = quantile(cv_error, 0.10, na.rm = TRUE),
      cv_q90  = quantile(cv_error, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # keep a clean numeric for plotting
    arrange(bw_mult)
}


# 4) Local eigenvalues λ_j(x) for j=1..k as an (n_spots x k) matrix
get_local_eigenvalues <- function(gwpca_obj, k = k_pcs) {
  V <- gwpca_obj$var
  V <- as.matrix(V)
  if (ncol(V) < k) k <- ncol(V)
  V[, seq_len(k), drop = FALSE]
}


# 5) Convenience for multiple PCs (returns n x length(pcs))
extract_aligned_pc_matrix <- function(S_list, pcs = 1:4, ref_scores) {
  do.call(
    cbind,
    lapply(pcs, function(pc) extract_aligned_pc_selfscores(S_list, pc, ref_scores))
  )
}


# 6) Run one GWPCA
run_one_gwpca <- function(bw_val) {
  ## Select the sample you would like to perform a GWPCA analysis
  sfe <- getSFE(msfe, "JBO019")
  ## Get the gene names that are going to be evaluated
  vars <- top_hvgs[["JBO019"]]
  ## Set the number of components to be retained
  k_pcs <- 20
  ## Set the kernel to be used
  kernel <- "gaussian"
  ## Set the Minkowski distance power: p = 2 --> Euclidean
  p <- 2
  ## Is the bandwidth adaptive?: No because spots are fixed
  adaptive <- FALSE
  ## Cross-Validate GWPCA?
  cv <- TRUE
  ## Calculate PCA scores?
  scores <- TRUE
  ## Run a robust GWPCA?
  robust <- FALSE

  ## Make a cluster for parallel computing (otherwise GWPCA is slow!)
  my.cl <- makeClusterGWPCA(spec = 20, type = "FORK")

  # Run GWPCA
  gw <- gwpcaSTE(
    sfe = sfe,
    assay = "logcounts",
    vars = vars,
    p = p,
    k = k_pcs,
    bw = bw_val,
    kernel = kernel,
    adaptive = adaptive,
    scores = scores,
    robust = robust,
    cv = cv,
    future = future,
    strategy = "cluster",
    workers = my.cl,
    verbose = FALSE
  )
  gw
}


# 7) Plot helpers
plot_cv_curve <- function(df_sum, use = c("cv_med", "cv_trim")[1]) {
  stopifnot(is.numeric(df_sum$bw_mult), is.numeric(df_sum[[use]]))
  # choose the scalar used for selection (median by default)
  best_bw <- df_sum$bw_mult[which.min(df_sum[[use]])]

  ggplot(df_sum, aes(x = bw_mult, y = .data[[use]])) +
    geom_ribbon(aes(ymin = cv_q10, ymax = cv_q90), alpha = 0.15) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = best_bw, linetype = 2) +
    labs(
      x = "Bandwidth (× spot diameter)",
      y = if (use == "cv_med") {
        "Cross-validation error (median across spots)"
      } else {
        "Cross-validation error (10% trimmed mean)"
      },
      caption = "Ribbon: 10–90% spot-wise CV"
    ) +
    theme_classic()
}


plot_moransI <- function(df_long) {
  ggplot(df_long, aes(bw_mult, I, colour = PC, group = PC)) +
    geom_line() +
    geom_point() +
    labs(x = "Bandwidth (× spot diameter)", y = "Moran's I") +
    theme_classic()
}

plot_cumvar <- function(df) {
  ggplot(df, aes(bw_mult, med_local_cumvar)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Bandwidth (× spot diameter)",
      y = "Median local cumulative variance (PC1–k)"
    ) +
    theme_classic()
}

plot_score_map <- function(scores_vec, sfe, title = "PC1 self-scores") {
  df <- as.data.frame(spatialCoords(sfe))
  df$score <- scores_vec
  colnames(df)[1:2] <- c("x", "y")
  ggplot(df, aes(x, y, fill = score)) +
    geom_point(shape = 21, size = 4, stroke = 0.1) +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(title = title, x = NULL, y = NULL, fill = "score") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}


#------------------------------------------------------------------------------#
# Analysis
## Set a fixed bandwidth
## bw is an important parameter as it defines the neighbourhood for which the
##  PCA will be calculated. The distance is measured in ultra-high resolution
##  image pixels. The default is 3x the diameter of the Visium spot. Make sure
##  to adjust it if it is too large or too small for your setting.
spot_diam <- sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
bw_mult <- c(2, 4, 6, 8, 10)
bw_grid <- bw_mult * spot_diam

# Export coordinates for Moran's I
coords <- SpatialExperiment::spatialCoords(sfe)

## =========================
## Bandwidth sweep
## =========================
# Compute global PC1 once for sign alignment
global_pc1 <- pcagw$pca$scores[, 1]

message(
  "Running GWPCA sweep over ",
  paste(bw_mult, collapse = ", "),
  " × spot diameter"
)

gwpca_fits <- setNames(vector("list", length(bw_grid)), paste0("c", bw_mult))

for (j in seq_along(bw_grid)) {
  gwpca_fits[[j]] <- run_one_gwpca(bw_grid[j])
}

## =========================
## Collate diagnostics
## =========================
diag_df <- purrr::imap_dfr(gwpca_fits, function(gw, lab) {
  cv_err <- get_cv_error(gw)

  # Local eigenvalues & cumulative variance (PC1..k)
  eig <- get_local_eigenvalues(gw, k = k_pcs) # (n x k)
  cum_loc <- rowSums(eig[, seq_len(min(k_pcs, ncol(eig))), drop = FALSE])
  med_cum <- median(cum_loc, na.rm = TRUE)

  # Moran's I on aligned self-score fields for PC1–PC4
  pcs_use <- 1:min(4, k_pcs)
  S_list <- gw$gwpca.scores
  score_mat <- extract_aligned_pc_matrix(
    S_list,
    pcs = pcs_use,
    ref_scores = global_pc1
  )
  I_vals <- apply(score_mat, 2, function(col) moransI_from_scores(col))

  tibble(
    label = lab,
    bw = as.numeric(sub("^c", "", lab)) * spot_diam,
    bw_mult = as.numeric(sub("^c", "", lab)),
    cv_error = cv_err,
    med_local_cumvar = med_cum,
    !!!setNames(as.list(I_vals), paste0("moransI_PC", pcs_use))
  )
})

# Figure out the selected bandwidth (CV minimum)
best_idx <- which.min(diag_df$cv_error)
best_bw_mult <- diag_df$bw_mult[best_idx]

## =========================
## Figure S1: CV error vs bandwidth
## =========================
diag_df_cv_sum <- summarise_cv(diag_df[, c("bw_mult", "cv_error")])
p_s1 <- plot_cv_curve(diag_df_cv_sum, use = "cv_med")
print(p_s1)

## =========================
## Figure S3: Median local cumulative variance and Moran’s I vs bandwidth
## =========================
# Cum var line
p_cum <- plot_cumvar(diag_df)
print(p_cum)

# Moran’s I for PC1–PC4
I_long <- diag_df %>%
  select(bw_mult, starts_with("moransI_PC")) %>%
  pivot_longer(-bw_mult, names_to = "PC", values_to = "I") %>%
  mutate(PC = sub("moransI_", "", PC))
p_I <- plot_moransI(I_long)
print(p_I)

## =========================
## Figure S2: PC1 score maps at three bandwidths
## =========================
bw_choices <- c(
  undersmoothed = 2.0, # 2×
  selected      = best_bw_mult, # argmin CV
  oversmoothed  = 10.0 # 10×
)

maps <- list()
for (nm in names(bw_choices)) {
  mult <- bw_choices[[nm]]
  gw <- gwpca_fits[[paste0("c", mult)]]
  if (is.null(gw)) {
    warning("Requested map for ", nm, " (", mult, "×) not in the sweep; skipping.")
    next
  }
  pc1_vec <- extract_aligned_pc_selfscores(gw$gwpca.scores, pc = 1, ref_scores = global_pc1)
  maps[[nm]] <- plot_score_map(pc1_vec, sfe, title = paste0("PC1 (", nm, ", ", mult, "×)"))
}
# Print the three panels
for (nm in names(maps)) print(maps[[nm]])

## =========================
## Save outputs (optional)
## =========================
ggsave("S1_cv_vs_bw.pdf", p_s1, width = 5, height = 4)
ggsave("S3_cumvar_vs_bw.pdf", p_cum, width = 5, height = 4)
ggsave("S3_moransI_vs_bw.pdf", p_I, width = 5, height = 4)
lapply(names(maps), function(nm) ggsave(paste0("S2_", nm, "_pc1_map.pdf"), maps[[nm]], width = 4, height = 4))
