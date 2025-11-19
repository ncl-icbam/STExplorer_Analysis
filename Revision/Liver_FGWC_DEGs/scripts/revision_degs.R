## ----------------------------------------------------------------------------#
## To load the packages use:
library(STExplorer)
library(tidyverse)
library(scran)
library(Matrix)
library(scran)
library(scuttle)
library(dplyr)
library(purrr)
library(readr)
library(RColorBrewer)

outdir <- file.path(getwd(), "Liver_FGWC_DEGs/results")

## ----------------------------------------------------------------------------#
## ---- helper: build hard assignments from fuzzy memberships at a cutoff ----
build_assignments <- function(fgwc_list, sid, cutoff = 0.60) {
  f <- fgwc_list[[sid]]
  M <- as.matrix(f$membership) # rows: spots; cols: clusters
  # Ensure spot IDs are present on the membership matrix:
  if (is.null(rownames(M)) && !is.null(names(f$cluster))) {
    rownames(M) <- names(f$cluster)
  }
  best_idx <- max.col(M, ties.method = "first") # 1..K
  best_score <- M[cbind(seq_len(nrow(M)), best_idx)] # max membership per spot
  data.frame(
    spot_id = rownames(M),
    sample = sid,
    best_clust = paste0("C", best_idx), # label clusters as C1..CK
    best_score = best_score,
    keep = best_score >= cutoff,
    stringsAsFactors = FALSE
  ) %>%
    mutate(best_clust = if_else(keep, best_clust, "mixed")) # nolint
}

## ---- helper: run scran::findMarkers with blocking by sample ----
run_find_markers <- function(sfe, groups,
                             test.type = "t",
                             pval.type = "any",
                             min.prop = 0.5,
                             add.summary = TRUE,
                             ...) {
  # scran expects log-expression; we'll use logNormCounts(sfe)
  stopifnot("logcounts" %in% SummarizedExperiment::assayNames(sfe))
  findMarkers(
    x = sfe,
    groups = groups,
    test.type = test.type, # "t", "wilcox", or "binom"
    pval.type = pval.type, # "any" is common for marker calling
    min.prop = min.prop, # gene should be up in at least 50% of pairwise comps
    direction = "up", # focus on up-markers for each cluster
    full.stats = TRUE,
    sorted = TRUE,
    add.summary = add.summary,
    ...
  )
}


## ----------------------------------------------------------------------------#
## ---- build SCE (skip if you already have one) ----
# counts: genes x spots (sparse Matrix recommended)
for (s in sampleNames) {
  message("Running sample: ", s)
  sfe <- getSFE(msfe, s)
  if (s == "JBO018") {
    sfe@assays@data$logcounts@seed@seed@seed@seed@seed@seed@filepath <- paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Healthy_", s, "_Results/outs/filtered_feature_bc_matrix.h5")
    sfe@assays@data$counts@seed@seed@filepath <- paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Healthy_", s, "_Results/outs/filtered_feature_bc_matrix.h5")
    sfe <- rmvImg(sfe, sample_id = s, image_id = "lowres")
    sfe <- addImg(sfe, sample_id = s, image_id = "lowres", scale_fct = 0.029013539, imageSource = paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Healthy_", s, "_Results/outs/spatial/tissue_lowres_image.png"))
  } else {
    sfe@assays@data$logcounts@seed@seed@seed@seed@seed@seed@filepath <- paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Steatotic_", s, "_Results/outs/filtered_feature_bc_matrix.h5")
    sfe@assays@data$counts@seed@seed@filepath <- paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Steatotic_", s, "_Results/outs/filtered_feature_bc_matrix.h5")
    sfe <- rmvImg(sfe, sample_id = s, image_id = "lowres")
    sfe <- addImg(sfe, sample_id = s, image_id = "lowres", scale_fct = 0.027540622, imageSource = paste0("/mnt/c/Users/Lefteris/OneDrive - Newcastle University/Projects/STExplorer_Analysis/Revision/Liver_FGWC_DEGs/data/Human_Liver_Steatotic_", s, "_Results/outs/spatial/tissue_lowres_image.png"))
  }

  ## ---- run at two cutoffs: 0.60 and 0.75 ----
  for (cutoff in c(0.60)) {
    message("Running marker detection at membership cutoff = ", cutoff)
    message("\tBuilding assignments...")
    assign_df <- build_assignments(fgwc_list, sid = s, cutoff = cutoff)

    # Groups = cluster labels
    groups <- factor(assign_df$best_clust)
    group_lbls <- unique(as.character(groups))

    # Optionally store in colData for traceability:
    colData(sfe)$FGWC_cluster <- groups
    colData(sfe)$membership <- assign_df$best_score

    # Subset the SFE object only to the assigned spots;
    # remove the spots marked as "mixed" (membership < 60%)
    selection_vS <- groups != "mixed"
    # selection_vG <- rownames(sfe) %in% top_hvgs[[s]]
    sfe_sub <- sfe[, selection_vS]
    groups_sub <- droplevels(groups[selection_vS]) # drop the "mixed" factor
    group_lbls_sub <- unique(groups_sub)
    gene_names <- as.data.frame(rowData(sfe_sub)["gene_name"])

    # Run scran::findMarkers (t-test);
    # set wilcox if robustness is preferred
    message("\tFinding markers...")
    markers <- run_find_markers(
      sfe_sub,
      groups = groups_sub,
      test.type = "wilcox", # alternatives: "wilcox" or "binom"
      pval.type = "any",
      min.prop = 0.8,
      add.summary = TRUE
    )

    # Add gene names
    for (g in group_lbls_sub) {
      markers[[g]] <- merge(
        as.data.frame(markers[[g]]),
        gene_names,
        by = "row.names", all = TRUE
      ) %>%
        rename("Row.names" = "gene_id")
    }

    ## ---- write per-cluster tables to CSV ----
    # 'markers' is a list with one data.frame per group (cluster).
    # Each data.frame has one row per gene with statistics across pairwise
    # contrasts, plus 'Top' ranking and 'summary.*' columns if add.summary=TRUE.
    out_dir <- file.path(
      outdir,
      "markers_scran",
      paste0("cutoff_", gsub("\\.", "", as.character(cutoff)))
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    lapply(names(markers), function(cl) {
      df <- markers[[cl]]
      # Keep all columns;
      # Filter genes for:
      #   FDR < 0.05;
      #   expressed in >= 80% of the spots in the cluster;
      #   log2FC > 0.8 in the comparison with the smallest p.value
      selection_params <- df$FDR < 0.05
      df_out <- df[selection_params, , drop = FALSE]
      df_out <- arrange(df_out, summary.stats)
      message("Number of genes in ", cl, ": ", nrow(df_out))
      # Write with cluster and cutoff in the filename:
      out_file <- file.path(out_dir, paste0(s, "_markers_", cl, ".csv"))
      write_csv(df_out, out_file)
      message("Wrote: ", out_file)
      invisible(NULL)
    })

    ## ---- also write a combined index of top markers per cluster ----
    top_tbl <- bind_rows(lapply(names(markers), function(cl) {
      df <- as.data.frame(markers[[cl]])
      tibble(
        cluster = cl,
        gene_name = df$gene_name,
        gene_id = df$gene_id,
        rank = df$Top,
        self_detected = df$self.detected,
        other_detected = df$other.detected,
        logFC_effectSize = df$summary.stats,
        summary_p = df$p.value,
        summary_FDR = df$FDR
      ) %>%
        arrange(cluster, rank) %>%
        filter(
          summary_FDR < 0.05
        )
    }))
    write_csv(top_tbl, file.path(out_dir, paste0(s, "_markers_top_index.csv")))

    ## ---- Get a list of top marker genes per cluster ----
    # Keep from each cluster only markers that ranked =< 100 (Top column)
    gene_list_topRank <- lapply(names(markers), function(cl) {
      df <- markers[[cl]]
      selection_params <- df$FDR < 0.05
      df <- df[selection_params, , drop = FALSE]
      df$self.diff <- df$self.average - df$other.average
      df <- arrange(df, self.diff)
      out <- df %>%
        select(gene_id) %>%
        .[[1]]
      names(out) <- df$gene_name
      message("Number of genes in ", cl, ": ", length(out))
      out
    })

    assign(paste0("gene_list_topRank_", cutoff), gene_list_topRank)

    ## ---- Create a heatmap of gene expression for the top ranking markers ----
    selection_df <- data.frame(
      gene_id = unlist(gene_list_topRank),
      gene_name = names(unlist(gene_list_topRank))
    ) %>% unique()

    message("Total number of markers left: ", length(selection_df$gene_id))

    message("\tPlotting heatmap...")
    out_file <- file.path(out_dir, paste0(s, "_heatmap.pdf"))
    scater::plotHeatmap(sfe_sub,
      features = selection_df$gene_id,
      exprs_values = "logcounts",
      center = TRUE, scale = TRUE,
      columns = colnames(sfe_sub)[order(sfe_sub$FGWC_cluster)],
      colour_columns_by = "FGWC_cluster", cluster_cols = FALSE,
      show_colnames = FALSE, cluster_rows = FALSE, filename = out_file,
      show_rownames = ifelse(length(selection_df$gene_id) > 50, FALSE, TRUE),
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100)
    )
  }
}

# PRGn
# BrBg
# YlGnBu
