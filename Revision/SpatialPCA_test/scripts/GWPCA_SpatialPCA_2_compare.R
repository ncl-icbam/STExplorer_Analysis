## ================== SETUP ====================================================
library(FNN) # kNN
library(igraph) # Leiden clustering
library(mclust) # ARI
library(dplyr)
library(tidyr)
library(SpatialPCA)

## Inputs expected:
## - gw_cfg         : list with spot_ids
## - gw_embed       : list with Z_std (or Z_compact) from GWPCA step 6
## - pcagw          : your gwpca object with pca$scores
## - labels         : named factor/character vector; names == gw_cfg$spot_ids
labels <- colData(sfe)[["annotation"]]
names(labels) <- colData(sfe)[["Barcode"]]

stopifnot(is.vector(labels), length(labels) == length(gw_cfg$spot_ids))
stopifnot(identical(names(labels), gw_cfg$spot_ids))

## Controls (keep fixed across methods)
gw_embed_fixed <- readRDS("./results/gw_embed_fixed.rds")
k_graph <- readRDS("./results/k_graph_fixed.rds")
D_target <- gw_embed_fixed$D
resolutions <- seq(0.2, 2.0, by = 0.2)
seeds <- 2023
# k_graph loaded as k_graph_fixed from RDS (10)


## ================== 0) UTILITIES =============================================
## Build symmetric kNN graph from an embedding (Euclidean)
knn_graph <- function(emb, k = 15L) {
  idx <- get.knn(emb, k = k)$nn.index
  edges <- cbind(rep(seq_len(nrow(emb)), each = k), as.vector(t(idx)))
  simplify(graph_from_edgelist(edges, directed = FALSE))
}

run_leiden <- function(g, k_graph = 15L, res = seq(0.2, 2.0, by = 0.2), seed = 2023) {
  do.call(rbind, lapply(res, function(r) {
    set.seed(seed)
    cl <- cluster_leiden(g, resolution = r, objective_function = "modularity")
    mem <- membership(cl)
  }))
}

## One full sweep over resolutions & seeds for a given embedding
benchmark_embedding <- function(emb, method_name, labels, k_graph, resolutions, seeds) {
  stopifnot(nrow(emb) == length(labels))
  # kNN graph
  g <- knn_graph(emb, k = k_graph)

  res_list <- list()
  ctr <- 1L
  for (r in resolutions) {
    for (sd in seeds) {
      memb <- run_leiden(g, res = r, seed = sd)
      # ARI requires same length vectors
      ari <- adjustedRandIndex(
        as.integer(factor(labels)),
        as.integer(factor(memb))
      )
      res_list[[ctr]] <- data.frame(
        method = method_name,
        resolution = r,
        seed = sd,
        n_clusters = length(unique(factor(memb))),
        ARI = ari,
        stringsAsFactors = FALSE
      )
      ctr <- ctr + 1L
    }
  }
  bind_rows(res_list)
}


## ================== 1) GWPCA embedding =======================================
## Prefer the compact embedding from script 1
Emb_GWPCA <- gw_embed_fixed$Z_compact # already z-scored & smoothed


## ================== 2) PCA embedding (full data, same genes) =================
## Use existing global PCA scores (spots × PCs) from the GWPCA run
stopifnot(!is.null(pcagw$pca$scores))
Scores_PCA <- pcagw$pca$scores
# Align order and select top D
Scores_PCA <- Scores_PCA[gw_cfg$spot_ids, , drop = FALSE]
Emb_PCA <- Scores_PCA[gw_cfg$spot_ids, seq_len(min(D_target, ncol(Scores_PCA))), drop = FALSE]
Emb_PCA <- base::scale(Emb_PCA) # column z-score


## ================== 3) SpatialPCA embedding ==================================
## Run SpatialPCA
Emb_SpatialPCA <- NULL

message("Running SpatialPCA…")
# Expression matrix 'X' (spots × genes) and coordinates 'coords'
# Note that SpatialPCA performs a Seurat-based normalisation. For this reason we
# are providing it with the raw counts matrix, but subsetted to the cells we
# kept after QC with STExplorer.
X <- SummarizedExperiment::assay(sfe, "counts")[, gw_cfg$spot_ids]
coords <- pcagw$SDF@coords[gw_cfg$spot_ids, , drop = FALSE]

# Minimal SpatialPCA run (adjust params as needed)
sp_fit <- CreateSpatialPCAObject(
  counts = X,
  location = coords,
  project = "SpatialPCA",
  customGenelist = top_hvgs$s151673,
  min.locations = 20,
  min.features = 20
  # Leave to default since not used here
  # gene.type = "custom",
  # sparkversion = "spark",
  # numCores_spark = 1,
  # gene.number = 3000
)

sp_fit <- SpatialPCA_buildKernel(
  sp_fit,
  kerneltype = "gaussian",
  bandwidthtype = "SJ",
  bandwidth.set.by.user = pcagw$GW.arguments$bw
)

sp_fit <- SpatialPCA_EstimateLoading(
  sp_fit,
  fast = FALSE,
  SpatialPCnum = D_target
)

sp_fit <- SpatialPCA_SpatialPCs(
  sp_fit,
  fast = FALSE
)

# Extract spots × D matrix
Emb_SpatialPCA <- t(sp_fit@SpatialPCs)[gw_cfg$spot_ids, seq_len(D_target), drop = FALSE]
Emb_SpatialPCA <- base::scale(Emb_SpatialPCA)


## ================== 4) BENCHMARK: Leiden + ARI ============================
benchmarks <- list()

benchmarks[["GWPCA"]] <- benchmark_embedding(
  emb = Emb_GWPCA,
  method_name = "GWPCA",
  labels = labels,
  k_graph = k_graph,
  resolutions = resolutions,
  seeds = seeds
)

benchmarks[["PCA"]] <- benchmark_embedding(
  emb = Emb_PCA,
  method_name = "PCA",
  labels = labels,
  k_graph = k_graph,
  resolutions = resolutions,
  seeds = seeds
)

benchmarks[["SpatialPCA"]] <- benchmark_embedding(
  emb = Emb_SpatialPCA,
  method_name = "SpatialPCA",
  labels = labels,
  k_graph = k_graph,
  resolutions = resolutions,
  seeds = seeds
)

res_df <- bind_rows(benchmarks)

res_df
## ================== 5) SUMMARISE & REPORT ================================
## Two views: (A) ARI at the resolution that matches #labels; (B) best ARI over the sweep.

n_labels <- length(unique(labels))

# (A) Match-by-count: choose resolution entries whose n_clusters == n_labels
match_count <- res_df %>%
  group_by(method, resolution) %>%
  summarise(ARI_mean = mean(ARI), n_clusters = first(n_clusters), .groups = "drop") %>%
  # filter(n_clusters == n_labels) %>%
  arrange(method, resolution)

# (B) Best-over-sweep per method
best_overall <- res_df %>%
  group_by(method, resolution) %>%
  summarise(ARI_mean = mean(ARI), ARI_sd = sd(ARI), n_clusters = first(n_clusters), .groups = "drop") %>%
  group_by(method) %>%
  slice_max(order_by = ARI_mean, n = 1, with_ties = FALSE) %>%
  arrange(desc(ARI_mean))

print(match_count, n = 100)
print(best_overall)

## ================== 6) OPTIONAL: quick plots =============================
# ARI vs resolution curves (mean ± sd over seeds)
p_ari <- res_df %>%
  group_by(method, resolution) %>%
  summarise(ARI_mean = mean(ARI), ARI_sd = sd(ARI), .groups = "drop") %>%
  ggplot(aes(x = resolution, y = ARI_mean, colour = method)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ARI_mean - ARI_sd, ymax = ARI_mean + ARI_sd, fill = method), alpha = 0.15, colour = NA) +
  theme_classic() +
  labs(title = "ARI vs Leiden resolution", y = "ARI (mean over seeds)", x = "Leiden resolution")

print(p_ari)

## Save to files
ggsave("./results/ari.pdf", p_ari, width = 6, height = 4)
