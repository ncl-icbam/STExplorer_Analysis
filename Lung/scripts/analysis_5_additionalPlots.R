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
# This script is for preparing supplementary plots for the main analysis.
#
#
# ---------------------------------------------------------------------------- #
## ----Additional Plots --------------------------------------------------------
### ----Annotation -------------------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Create folder
  graphics_out <- paste0(projectFolder, "/graphics_out/Annotation/")
  dir.create(graphics_out, recursive = TRUE)

  ## Title
  title <- paste0(metadata$condition[metadata$sample_id == id],
                  " (",
                  metadata$subject_alias[metadata$sample_id == id],
                  ")")

  ## Plot Image and annotation
  p1 <- plotQC_tissueImg(msfe,
                         sample_id = id,
                         res = "lowres",
                         type = "none",
                         annotate = FALSE) +
    labs(title = title,
         subtitle = "") +
    theme(plot.title = element_text(hjust = 0.5))

  p2 <- plotQC_spotsAnnotation(msfe,
                               sample_id = id,
                               type = "hex",
                               colours = colours) +
    labs(title = "") # +
  # scale_fill_manual(limits = names(colours))

  print((p1 + p2))

  prfx <- id
  main <- "_Annotation"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 5.22,
         height = 8.23,
         units = "in",
         dpi = 300)

  ## Housekeeping
  rm(p1, p2, title)
}


### ----Gene expression of genes -----------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot Expression
  for (i in seq_along(selectedGenes_ENSG)) {
    g <- selectedGenes_ENSG[i]
    g_name <- names(g)
    message(paste0("\t... on gene: ", g_name, " (", i, "/", length(selectedGenes_ENSG), ")"))

    ## Create folder
    graphics_out <- paste0(projectFolder, "/graphics_out/Expression/", g_name, "/")
    dir.create(graphics_out, recursive = TRUE)

    plotGeneExpression(msfe,
                       sample_id = id,
                       genes = g,
                       assay = "logcounts")

    prfx <- id
    main <- "_Expression"
    sfx <- paste0("_", g_name)
    other <- ""

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           device = "png",
           width = 5.22,
           height = 8.23,
           units = "in",
           dpi = 300)
  }
  ## Housekeeping
  rm(g_name)
}


### ----Cell Abundances --------------------------------------------------------
abundance_df_lst <- vector(mode = "list", length = length(sampleNames_selected))
names(abundance_df_lst) <- sampleNames_selected

#### ----Plot: map -------------------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Extract abundances for specific cell types -------------------------------#
  abundance_df <- cellAbundance[[id]] %>%
    select(c("Barcode", paste0("q05cell_abundance_w_sf_", selectedCellTypes)))

  ## Filter out any abundance value less than 0.2 as it might be noise
  # abundance_df <- mutate(abundance_df, across(-Barcode, ~ifelse(. < 0.2, NA, .)))
  ## Fix column names
  colnames(abundance_df) <- gsub("q05cell_abundance_w_sf_", "", colnames(abundance_df))
  ## Extract geometries
  geoms <- colGeometry(msfe@sfe_data[[id]], "spotHex") %>%
    rownames_to_column(var = "Barcode")
  ## Merge abundances and geometries to plot
  abundance_df <- left_join(abundance_df, geoms, by = "Barcode")

  ## Output the list
  abundance_df_lst[[id]] <- abundance_df

  ## Plot a map per cell type -------------------------------------------------#
  for (cell in selectedCellTypes) {
    message(paste0("\t... on cell type: ", cell))

    ## Create folder
    graphics_out <- paste0(projectFolder, "/graphics_out/CellAbundance/", cell, "/")
    dir.create(graphics_out, recursive = TRUE)

    ## Prepare
    title <- gsub(".C", " c", cell)

    ## Plot
    ggplot(data = abundance_df,
           aes(fill = .data[[cell]])) +
      geom_sf(aes(geometry = geometry)) +
      scale_fill_viridis_c(option = "turbo",
                           na.value = "transparent") +
      labs(title = title,
           fill = "Cell\nabundance") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 20))


    prfx <- id
    main <- "_CellAbundance"
    sfx <- paste0("_", cell)
    other <- ""

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           device = "png",
           width = 6.64,
           height = 4.59,
           units = "in",
           dpi = 300)
  }

  ## Plot as panel ------------------------------------------------------------#
  ## Create folder
  graphics_out <- paste0(projectFolder, "/graphics_out/CellAbundance/Panel/")
  dir.create(graphics_out, recursive = TRUE)

  abundance_df <- pivot_longer(abundance_df,
                               cols = -c("Barcode", "geometry"),
                               names_to = "cellType",
                               values_to = "abundance")

  ## Plot
  ggplot(data = abundance_df,
         aes(fill = abundance)) +
    geom_sf(aes(geometry = geometry)) +
    scale_fill_viridis_c() +
    labs(fill = "Cell\nabundance") +
    facet_wrap(~cellType) +
    theme_void()

  prfx <- id
  main <- "_CellAbundance"
  sfx <- "_panel"
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 8.86,
         height = 7.43,
         units = "in",
         dpi = 300)

  ## Housekeeping
  rm(abundance_df, title, geoms)
}

#### ----Plot: heatmap ---------------------------------------------------------
graphics_out <- paste0(projectFolder, "/graphics_out/CellAbundance/Heatmaps/")
dir.create(graphics_out, recursive = TRUE)

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Extract abundances for specific cell types -------------------------------#
  abundance_df <- abundance_df_lst[[id]] %>%
    column_to_rownames(var = "Barcode") %>%
    select(!all_of("geometry")) %>%
    t() %>%
    data.frame()
  colnames(abundance_df) <- gsub("[.]", "-", colnames(abundance_df))

  ## Plot
  column_annot <- colData(msfe@sfe_data[[id]])[,c("Barcode","senMayo")] %>%
    data.frame() %>%
    select(!all_of("Barcode")) %>%
    arrange(senMayo)

  ## Order input column data by SenMayo enrichment score
  abundance_df <- abundance_df[,rownames(column_annot)]

  ## Order input row data by cell type
  abundance_df <- abundance_df[order(rownames(abundance_df)),]

  for (std in c("row", "column")) {
    ph <- pheatmap(abundance_df,
                   color = rev(cols4all::c4a(palette = "brewer.br_bg", n = 9)),
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   main = "Cell Type Density",
                   scale = std,
                   annotation_col = column_annot,
                   labels_row = names(selectedCellTypes)[order(names(selectedCellTypes))]
    )

    prfx <- id
    main <- "_CellAbundance"
    sfx <- "_heatmap"
    other <- ifelse(std == "row", "_row_std", "_col_std")

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           plot = ph,
           device = "png",
           width = 8.86,
           height = 7.43,
           units = "in",
           dpi = 300)
  }

  ## Housekeeping
  rm(abundance_df)
}

