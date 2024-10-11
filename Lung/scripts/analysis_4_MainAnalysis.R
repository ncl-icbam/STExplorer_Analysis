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
# This script is for the main analysis that also generates the pltos used in the
# manuscript figures.
#
#
# ---------------------------------------------------------------------------- #
## ----Senescence Score Analysis -----------------------------------------------

# Combine all genes across datasets
all_genes <- Reduce(intersect, lapply(msfe@sfe_data, rownames))

# Filter the 125 genes that are common across all datasets
genes_to_analyse <- intersect(all_genes, senMayoGenes_ENSG)

# Calculate presence proportions for each dataset
presence_list <- lapply(msfe@sfe_data,
                        calculate_presence_proportion,
                        genes = genes_to_analyse,
                        assay = "sct",
                        threshold = 0)

# Combine the presence proportions into a data frame
presence_df <- do.call(cbind, presence_list) %>% data.frame(check.names = FALSE)
rownames(presence_df) <- names(selectedGenes_ENSG)[match(rownames(presence_df), selectedGenes_ENSG)]

# Get the donors table
donors_df <- metadata[, c("sample_id","subject_alias", "fibrotic_extent_score_by_pathologist_0.3")]
donors <- unique(grep("IPF_", donors_df$subject_alias, value = TRUE))

# Prepare colours and annotations
column_annot <- data.frame(FS = metadata[match(colnames(presence_df), metadata$sample_id), "fibrotic_extent_score_by_pathologist_0.3"])
rownames(column_annot) <- colnames(presence_df)
column_annot <-  arrange(column_annot, FS)

presence_df <- presence_df[,rownames(column_annot)]

annotation_colours <- c("#FEF0D9", "#FDCC8A", "#FC8D59", "#D7301F")
names(annotation_colours) <- (column_annot %>% arrange(FS) %>% unique())[["FS"]]


# Plot senMayo gene proportions in each sample - row-standardised=
graphics_out <- paste0(projectFolder, "/graphics_out/Senescence/")
dir.create(graphics_out, recursive = TRUE)

for (std in c("row", "none")) {
  ph <- pheatmap(presence_df,
                 color = rev(cols4all::c4a(palette = "brewer.br_bg", n = 100)),
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 main = "Proportion of SenMayo senescence markers",
                 scale = "row",
                 annotation_col = column_annot,
                 annotation_colors = list(FS = annotation_colours),
                 angle_col = 0)


  prfx <- "senMayo_proportions_"
  main <- ""
  sfx <- "_complete"
  other <- ifelse(std == "row", "_std", "")

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         plot = ph,
         device = "png",
         width = 7.77,
         height = 14,
         units = "in",
         dpi = 300)
}


# Calculate the average presence proportion across samples with the same FS score
# Loop through each fibrosis level
for (fs in 0:3) {
  # Create dynamic variable names for vectors and columns
  av_prop_name <- paste0("AvProp_FS", fs)

  # Filter donors_df for the current fibrosis score and select sample_id
  FS <- donors_df[donors_df$fibrotic_extent_score_by_pathologist_0.3 == fs, "sample_id"]
  # if (fs == 0) { # remove sample V10T03-281-A1 because it shows not healthy patterns
  #   FS <- FS[FS != "V10T03-281-A1"]
  # }

  # Calculate the average proportion for the current fibrosis score
  presence_df[[av_prop_name]] <- rowMeans(presence_df[, FS])
}

# # Set a threshold for selecting genes (e.g., expressed in at least 30% of locations)
# threshold <- 0.3
# selected_genes <- rownames(presence_df)[presence_df$AvProp_IPF3_4 >= threshold]

# Plot senMayo gene proportions in each sample for the selected senMayo subset
column_annot <- data.frame(FS = c(0, 1, 2, 3))
rownames(column_annot) <- grep("AvProp_FS", colnames(presence_df), value = TRUE)
annotation_colours <- c("#FEF0D9", "#FDCC8A", "#FC8D59", "#D7301F")
names(annotation_colours) <- (column_annot %>% arrange(FS) %>% unique())[["FS"]]

ph <- pheatmap(presence_df[,grep("AvProp_FS", colnames(presence_df))],
               show_rownames = TRUE,
               show_colnames = TRUE,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               main = "Average proportion of SenMayo senescence markers\nper Fibrotic Score",
               scale = "row",
               annotation_col = column_annot,
               annotation_colors = list(FS = annotation_colours),
               angle_col = 0)

graphics_out <- paste0(projectFolder, "/graphics_out/Senescence/")
dir.create(graphics_out, recursive = TRUE)
prfx <- ""
main <- "senMayo_proportions"
sfx <- "_averagePerFS"
other <- "_std"

ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
       plot = ph,
       device = "png",
       width = 7.77,
       height = 16,
       units = "in",
       dpi = 300)

### ----Plot Module Scores -------------------------------------------------
#### Plot Maps ----
for (i in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[i]
  message("Working on sample: ", id, " (", i, "/", length(sampleNames_selected), ")")

  ## Fetch senescence scores names
  moduleScores <- c("senMayo", "emt", "ecm", "tgfb")

  ## Create folder
  graphics_out <- paste0(projectFolder, "/graphics_out/Expression/ModuleScores/")
  dir.create(graphics_out, recursive = TRUE)

  for (j in seq_along(moduleScores)) {
    m <- moduleScores[j]
    message(paste0("\t... on gene: ", m, " (", j, "/", length(moduleScores), ")"))

    plotModuleScores(msfe,
                     sample_id = id,
                     modules = m,
                     title_type = "") +
      cols4all::scale_fill_continuous_c4a_div(palette = "-c4a.pu_gn_div")

    prfx <- "ModuleScore"
    main <- paste0("_", m)
    sfx <- id
    other <- ""

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           device = "png",
           width = 5.22,
           height = 8.23,
           units = "in",
           dpi = 300)
  }
  ## Housekeeping
  rm(moduleScores)
}

## Plot maps as a panel
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/Expression/ModuleScores/")
## Fetch senescence scores names
moduleScores <- c("senMayo", "emt", "ecm", "tgfb")

for (j in seq_along(moduleScores)) {
  m <- moduleScores[j]
  message(paste0("\t... on gene: ", m, " (", j, "/", length(moduleScores), ")"))

  plotModuleScores(sfe_multi,
                   sample_id = NULL,
                   modules = m,
                   title_type = "",
                   facet_by = "sample_id") +
    cols4all::scale_fill_continuous_c4a_div(palette = "-c4a.pu_gn_div") +
    labs(fill = paste0(m, "\nscore"))

  prfx <- "ModuleScore"
  main <- paste0("_", m)
  sfx <- ""
  other <- "_panel"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 14.1,
         height = 11.5,
         units = "in",
         dpi = 300)
}

#### Plot rain - clouds by FS ----
## Fetch senescence scores names
moduleScores <- c("senMayo", "emt", "ecm", "tgfb")

## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/Expression/ModuleScores/")

# Export module score into a data frame
df <- colData(sfe_multi)[c("sample_id",moduleScores)] %>%
  data.frame() %>%
  mutate(FS = metadata[match(.data[["sample_id"]], metadata$sample_id), "fibrotic_extent_score_by_pathologist_0.3"]) %>%
  pivot_longer(cols = -c("sample_id", "FS"), names_to = "mScores", values_to = "scores")

# Plot
ggplot(df, aes(x = FS, y = scores, fill = FS)) +
  ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
  scale_fill_manual(values = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#D7301F")) +
  ggpubr::stat_compare_means(ref.group = "0",
                             method = "wilcox.test",
                             label = "p.signif",
                             tip.length = 0.005) +
  labs(y = "Module Score",
       x = "",
       fill = "Fibrotic\nScore") +
  facet_wrap(~mScores, nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

prfx <- "ModuleScore"
main <- "_all"
sfx <- "_rainCloud"
other <- "_byFS"

ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
       device = "png",
       width = 9.5,
       height = 6.9,
       units = "in",
       dpi = 300)

#### Plot rain - clouds by Annotation ----
# Export module score into a data frame
df <- colData(sfe_multi)[c("sample_id", "annotation", moduleScores)] %>%
  data.frame() %>%
  mutate(annotation = fct_relevel(annotation,
                                  "Normal Alveolar and Other", "Diseased",
                                  "Suspect Fibrosis", "Inflammation",
                                  "Large Airway", "Blood Vessel")) %>%
  pivot_longer(cols = -c("sample_id", "annotation"), names_to = "mScores", values_to = "scores")

## Named colour values for consistent colouring of annotation labels
colours <- c("Blood Vessel" = "#1F78C8",
             "Large Airway" = "#2C3E50",
             "Normal Alveolar and Other" = "#33A02C",
             "Inflammation" = "#FFD700",
             "Suspect Fibrosis" = "#FD8D3C",
             "Diseased" = "#CC4C00")

# Plot
ggplot(df, aes(x = annotation, y = scores, fill = annotation)) +
  ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
  scale_fill_manual(values = colours) +
  ggpubr::stat_compare_means(ref.group = "Normal Alveolar and Other",
                             method = "wilcox.test",
                             label = "p.signif",
                             tip.length = 0.005) +
  labs(y = "Module Score",
       x = "",
       fill = "Annotation") +
  facet_wrap(~mScores, nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

prfx <- "ModuleScore"
main <- "_all"
sfx <- "_rainCloud"
other <- "_byAnnotation"

ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
       device = "png",
       width = 9.5,
       height = 6.9,
       units = "in",
       dpi = 300)

#### Plot heatmap ----
## Fetch senescence scores names
moduleScores <- c("senMayo", "emt", "ecm", "tgfb")

## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/Expression/ModuleScores/")

# Export module score into a data frame
df <- colData(sfe_multi)[c("sample_id", "annotation", moduleScores)] %>%
  data.frame() %>%
  mutate(FS = metadata[match(.data[["sample_id"]], metadata$sample_id), "fibrotic_extent_score_by_pathologist_0.3"]) %>%
  rownames_to_column(var = "Barcode") %>%
  mutate(unique_id = paste(sample_id, Barcode, sep = "_"), .keep = "unused")
rownames(df) <- df$unique_id

# Prepare colours and annotations
column_annot <- df %>%
  select(all_of(c("FS", "annotation"))) %>%
  mutate(annotation = fct_relevel(annotation,
                                  "Normal Alveolar and Other", "Diseased",
                                  "Suspect Fibrosis", "Inflammation",
                                  "Large Airway", "Blood Vessel")) %>%
  arrange(FS, annotation)

df <- df[rownames(column_annot),] %>%
  select(!all_of(c("FS", "unique_id", "annotation"))) %>%
  t() %>%
  data.frame(check.names = FALSE)
# %>%
#   select(!contains("281-A1"))

annot_colours_FS <- c("#FEF0D9", "#FDCC8A", "#FC8D59", "#D7301F")
names(annot_colours_FS) <- column_annot$FS %>% unique()

annotation_colours <- list(FS = annot_colours_FS,
                           annotation = colours)

# Plot
for (module in moduleScores) {
  ph <- pheatmap(df[module,],
                 color = rev(cols4all::c4a(palette = "brewer.br_bg", n = 100)),
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 main = "Module Scores",
                 scale = "row",
                 annotation_col = column_annot,
                 annotation_colors = annotation_colours)

  prfx <- "ModuleScore"
  main <- paste0("_", module)
  sfx <- "_heatmap"
  other <- "_byFS-Annotation"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         plot = ph,
         device = "png",
         width = 14.1,
         height = 3,
         units = "in",
         dpi = 300)
}


# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
# saveRDS(msfe_gsva, file = paste0(projectFolder, "/msfe_gsva.rds"))

### ---- Spatial Autocorrelation -----------------------------------------------
## Set variables to run SA for
moduleScores <- c("senMayo", "emt", "ecm", "tgfb" )
names(moduleScores) <- c("senMayo", "emt", "ecm", "tgfb" )

#### ----Local Moran's I -------------------------------------------------------
## Run Local statistic
for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  msfe <- moranLocalIPerm(m_sfe = msfe,
                          sample_id = id,
                          genes = moduleScores,
                          mc.cores = 6)
}


#### ----Local Geary's C -------------------------------------------------------
## Run Local statistic
for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  msfe <- gearyLocalCPerm(m_sfe = msfe,
                          sample_id = id,
                          genes = moduleScores,
                          mc.cores = 6)
}


#### ----Local Getis & Ord's G -------------------------------------------------
## Run Local statistic
for (id in sampleNames_selected) {
  message("Working on sample: ", id)

  msfe <- getisLocalGPerm(m_sfe = msfe,
                          sample_id = id,
                          genes = moduleScores,
                          mc.cores = 6)
}

#### ----Plot local SA ----------------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  # ## Fetch module scores
  # senScores <- rownames(msfe_gsva@sfe_data[[id]])

  for (s in moduleScores) {
    message(paste0("\t... on score: ", s))

    for (t in c("moran", "geary", "getis")) {
      if (t == "moran") {
        message("\t\t... on Moran's I")
        graphics_out <- paste0(projectFolder, "/graphics_out/SA/SA_Local_MoranI/", s, "/")
        dir.create(graphics_out, recursive = TRUE)
        main <- "_SA_LocalI"

      } else if (t == "geary") {
        message("\t\t... on Geary's C")
        graphics_out <- paste0(projectFolder, "/graphics_out/SA/SA_Local_GearyC/", s, "/")
        dir.create(graphics_out, recursive = TRUE)
        main <- "_SA_LocalC"

      } else if (t == "getis") {
        message("\t\t... on Getis & Ord's G")
        graphics_out <- paste0(projectFolder, "/graphics_out/SA/SA_Local_GetisG/", s, "/")
        dir.create(graphics_out, recursive = TRUE)
        main <- "_SA_LocalG"
      }

      ## Plot
      for (l in c("all", "significant")) {
        p1 <- plotSA_local(m_sfe = msfe,
                           sample_id = id,
                           feature = s,
                           statistic = t,
                           test = "permutation",
                           pVal = 0.05,
                           type = "hex",
                           title = "name",
                           locations = l)

        p2 <- plotSA_localClust(m_sfe = msfe,
                                sample_id = id,
                                feature = s,
                                statistic = t,
                                test = "permutation",
                                pVal = 0.05,
                                type = "hex",
                                title = "name",
                                clust_col = NULL,
                                locations = l)

        print((p1 / p2))

        prfx <- id
        main <- main
        sfx <- paste0("_", s)
        other <- ifelse(l == "all", "_all", "_sig")

        ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
               device = "png",
               width = 5.22,
               height = 8.23,
               units = "in",
               dpi = 300)
      }
    }
  }
}
## Housekeeping
rm(p1, p2)

# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
saveRDS(msfe, file = paste0(projectFolder, "/analysis_objs/msfe.rds"))


### ----Subset - Top 25% of SenMayo scores in all samples ----------------------
## Compare different variables between areas with different senMayo scores.
## We use the top quantile (25%) of the senMayo scores from all samples to set
## a threshold for what is going to be counted as senescent and what is not.
## enrichment using violin/ rain-cloud  and other plots.
#### ----Set Cut-off ------------------------------------------------------------
## Calculate the 75th percentile (top 25%)
sen_cutoff <- quantile(sfe_multi$senMayo, probs = 0.75)

## Select cell types and genes to plot in the subsetted data
cell_types <- cellTypes
sub_genes <- c("IL6", "CDKN1A", "ACTA2", "COL1A1", "GLB1", "CCL2", "ICAM1",
               "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2")
genes <- selectedGenes_ENSG[sub_genes]
names(genes) <- sub_genes

## Prepare images list to plot
image_to_plot_lst <- vector(mode = "list", length = length(sampleNames_selected))
names(image_to_plot_lst) <- sampleNames_selected

for (id in sampleNames_selected) {
  ### Fetch Image data
  image <- imgRaster(getImg(msfe@sfe_data[[id]],
                            sample_id = id,
                            image_id = "lowres"))
  max <- image@ptr$range_max[1]
  if (is.na(max)) {
    max <- max(image[,,1])
  }
  if (max <= 1) {
    max_col_value <- 1
  } else if (max <= 255) {
    max_col_value <- 255
  } else if (max <= 65536) {
    max_col_value <- 65536
  }
  coords <- spatialCoords(msfe@sfe_data[[id]]) %>% #[gwr_dt$Barcode,] %>%
    matrix(ncol = 2)
  limits_list <- list(y_lims = c(min(coords[,2]), max(coords[,2])),
                      x_lims = c(min(coords[,1]), max(coords[,1])))

  message("\tOutputing to list...")
  ## Output to the lists
  image_to_plot_lst[[id]] <- list(image = image,
                                  max_col_value = max_col_value,
                                  limits_list = limits_list)
}

## Visualise the threshold
ggplot(as.data.frame(colData(sfe_multi)),
       aes(x = senMayo)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = diff(range(colData(sfe_multi)$senMayo))/30,
                 fill = "lightblue", colour = "black") +
  geom_density(colour = "red") +
  geom_vline(xintercept = sen_cutoff, colour = "red") +
  labs(title = "SenMayo module score distribution",
       x = "senMayo",
       y = "Density") +
  theme_classic2()

graphics_out <- paste0(projectFolder, "/graphics_out/Combo/subset_senMayo/")
dir.create(graphics_out, recursive = TRUE)
prfx <- id
main <- "_Distribution"
sfx <- "_senMayoScores"
other <- ""

ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
       device = "png",
       width = 9.4,
       height = 4.9,
       units = "in",
       dpi = 300)

#### ----Apply Cut-off ---------------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Create folder
  graphics_out <- paste0(projectFolder, "/graphics_out/Combo/subset_senMayo/")
  dir.create(graphics_out, recursive = TRUE)

  ## Extract senMayo enrichment scores
  data_senMayo <- colData(msfe@sfe_data[[id]])[c("Barcode", "senMayo")] %>%
    data.frame() %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "positive",
                                 senMayo < sen_cutoff ~ "negative",
                                 .default = "no"))
  data_senMayo_cols_to_use <- c("sen_group", "Barcode")

  ## Extract cell type abundances
  data_cellAb <- colData(msfe@sfe_data[[id]]) %>%
    data.frame() %>%
    select(starts_with("q05cell")) %>%
    select(contains(cell_types)) %>%
    mutate(Barcode = rownames(.))

  ## Fix column names
  colnames(data_cellAb) <- gsub("q05cell_abundance_w_sf_", "", colnames(data_cellAb))

  ## Merge and pivot longer for faceting
  data <- left_join(data_senMayo[data_senMayo_cols_to_use],
                    data_cellAb,
                    by = "Barcode") %>%
    select(-"Barcode") %>%
    pivot_longer(cols = -c("sen_group"),
                 names_to = "cell_type",
                 values_to = "cell_density") %>%
    filter(sen_group != "no")

  data$cell_type <- factor(data$cell_type,
                           levels = cellTypes,
                           labels = names(cellTypes))

  ## Plot SenMayo score for locations passing the cutoff only -----------------#
  geoms <- colGeometry(msfe@sfe_data[[id]], "spotHex") %>%
    data.frame() %>%
    mutate(Barcode = colnames(msfe@sfe_data[[id]]))
  rownames(geoms) <- colnames(msfe@sfe_data[[id]])

  data_senMayo <- left_join(data_senMayo, geoms, by = "Barcode") %>%
    filter(sen_group != "no")

  image <- image_to_plot_lst[[id]][["image"]]
  max_col_value <- image_to_plot_lst[[id]][["max_col_value"]]
  limits_list <- image_to_plot_lst[[id]][["limits_list"]]

  ggplot(data_senMayo,
         aes(geometry = geometry)) +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value,
                                   alpha = 0.1) +
    # ggplot2::lims(x = limits_list[[2]],
    #               y = limits_list[[1]]) +
    geom_sf(aes(fill = senMayo)) +
    scale_fill_continuous_c4a_div(reverse = TRUE) +
    theme_void()

  prfx <- id
  main <- ""
  sfx <- "_senMayoScores_"
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 7.45,
         height = 5.45,
         units = "in",
         dpi = 300)

  ## Plot cell Abundances in rain-cloud plots ---------------------------------#
  message("\tPlotting cell abundance rain-clouds...")
  ggplot(data,
         aes(x = sen_group, y = cell_density, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    # stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, colour = "black") +  # Add mean with bootstrapped confidence interval
    # stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") + # Add mean point
    ggpubr::stat_compare_means(ref.group = "negative",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "",
         y = "Absolute Cell Density",
         x = "",
         fill = "SenMayo\nenrichment\nscore") +
    facet_wrap(~cell_type, nrow = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_CellDensities"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 15,
         height = 12,
         units = "in",
         dpi = 300)

  ## Plot gene expression -----------------------------------------------------#
  ## Merge and pivot longer for expression
  data_expr <- assay(msfe@sfe_data[[id]], "sct")[genes,] %>%
    t() %>%
    data.frame() %>%
    mutate(Barcode = rownames(.))

  colnames(data_expr) <- c(names(genes), "Barcode")

  data <- data_senMayo %>%
    select(-c("senMayo", "geometry")) %>%
    left_join(., data_expr, by = "Barcode") %>%
    select(-"Barcode") %>%
    pivot_longer(cols = -c("sen_group"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  ## Plot
  message("\tPlotting gene expression rain-clouds...")
  ggplot(data,
         aes(x = sen_group, y = logcounts, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    # stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, colour = "black") +  # Add mean with bootstrapped confidence interval
    # stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") + # Add mean point
    ggpubr::stat_compare_means(ref.group = "negative",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "",
         y = "Normalised Log2-counts",
         x = "",
         fill = "SenMayo\nenrichment\nscore") +
    facet_wrap(~gene, nrow = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_GeneExpression"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 10.83,
         height = 8.83,
         units = "in",
         dpi = 300)
}


### ---- GWR -------------------------------------------------------------------
#### ----Build formulas --------------------------------------------------------
formulas_senVSMisc <- c("senMayo~ecm", "senMayo~emt", "senMayo~tgfb",
                        "senMayo~q05cell_abundance_w_sf_AT2",
                        "senMayo~q05cell_abundance_w_sf_Myofibroblasts")
names(formulas_senVSMisc) <- gsub("q05cell_abundance_w_sf_", "", formulas_senVSMisc)

#### ----Run GWR ---------------------------------------------------------------
gwr_list_senScVsMisc <- vector(mode = "list", length = length(sampleNames_selected))
gwr_runs_senScVsMisc <- vector(mode = "list", length = length(formulas_senVSMisc))
names(gwr_list_senScVsMisc) <- sampleNames_selected
names(gwr_runs_senScVsMisc) <- names(formulas_senVSMisc)

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  for (i in seq_along(formulas_senVSMisc)) {
    ## Some formulas return `Error: inv(): matrix is singular`. This is because
    ## there are too few values in the area to generate a regression model.
    ## For the sake of running everything once and tackling the erroring formulas
    ## one by one, we wrap everything inside a tryCatch function.
    f <- formulas_senVSMisc[i]
    tryCatch({
      message("\tWorking on formula: ", names(f), " (", i, "/", length(formulas_senVSMisc), ")")

      ### ----Set Bandwidth ---------------------------------------------------#
      message("\t\tSetting bandwidth...")
      spot_diameter <- msfe@sfe_data[[id]]@metadata[["spotDiameter"]][["V10T03-280-A1"]][["spot_diameter_fullres"]]
      bw <- 3*spot_diameter


      ### ----Run GWR ---------------------------------------------------------#
      message("\t\tRunning GWR...")

      gwr_runs_senScVsMisc[[names(f)]] <- gwrSTE(gwr_method = "basic",
                                                 formula = f,
                                                 m_sfe = msfe,
                                                 sample_id = id,
                                                 bw = bw,
                                                 kernel = "exponential",
                                                 assay = "sct")
    }, error = error)
  }

  ## Update output list
  message("Updating the list of GWR runs per sample.")
  gwr_list_senScVsMisc[[id]] <- gwr_runs_senScVsMisc
}

## Housekeeping
rm(spot_diameter, gwr_runs_senScVsMisc)



#### ----Plot GWR --------------------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Select GWR run for the specific sample
  gwrRuns <- gwr_list_senScVsMisc[[id]]

  for (i in seq_along(formulas_senVSMisc)) {
    f <- formulas_senVSMisc[i]
    message("\tWorking on formula: ", names(f), " (", i, "/", length(formulas_senVSMisc), ")")

    ## Create folder
    graphics_out <- paste0(projectFolder, "/graphics_out/GWR_senScore/", names(f), "/")
    dir.create(graphics_out, recursive = TRUE)

    ## Get independent variables
    independent_vars <- strsplit(f, "~")[[1]][2] %>%
      strsplit(., "\\+") %>% .[[1]]

    ## Plot GWR coefficients
    plotGWR(gwrRuns[[i]], predictorVar = independent_vars, title = names(f))

    prfx <- id
    main <- "_GWR"
    sfx <- paste0("_", names(f))
    other <- ""

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           device = "png",
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 300)

    ## Plot local R^2
    plotGWR_R2(gwrRuns[[i]])

    prfx <- id
    main <- "_GWR"
    sfx <- paste0("_", names(f))
    other <- "_R2"

    ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
           device = "png",
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 300)

  }
  ## Housekeeping
  rm(gwRuns)
}

# ---------------------------------------------------------------------------- #
## ----Checkpoint ------------------------------------------------------------ #
saveRDS(gwr_list_senScVsMisc, file = paste0(projectFolder, "/gwr_list_senScVsMisc.rds"))

### ----Subset analysis --------------------------------------------------------
#### ----Select Variables ------------------------------------------------------
sub_genes <- c("IL6", "CDKN1A", "ACTA2", "COL1A1", "GLB1", "CCL2", "ICAM1",
               "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2")
sub_genesENSG <- selectedGenes_ENSG[sub_genes]
sub_formulas <- formulas_senVSCells#[sub_formulas]
sub_samples <- sampleNames_selected
sub_cell_names <- c("AT2", "Myofibroblasts")
sub_cells <- paste0("q05cell_abundance_w_sf_", selectedCellTypes[sub_cell_names])
names(sub_cells) <- sub_cell_names

data_to_plot_lst <- vector(mode = "list", length = length(sub_samples))
names(data_to_plot_lst) <- sub_samples

image_to_plot_lst <- vector(mode = "list", length = length(sub_samples))
names(image_to_plot_lst) <- sub_samples

#### ----Prepare data ----------------------------------------------------------
for (j in seq_along(sub_samples)) {
  id <- sub_samples[j]
  message("Working on sample: ", id, " (", j, "/", length(sub_samples), ")")

  ## Select GWR run for the specific sample
  message("\tExtracting GWR data...")
  gwr_list <- gwr_list_senScVsMisc[[id]]

  gwr_dt_list <- lapply(names(sub_formulas),
                        function(x){
                          ### Get predictor variable
                          mdl <- gwr_list[[x]]$GW.arguments$formula
                          predictorVar <- gsub(".+?~", "", mdl) %>%
                            strsplit(., "[+]")
                          message(predictorVar)
                          predictorVar <- predictorVar[[1]]
                          ### Identify statistically significant areas
                          predTVcol <- paste0(predictorVar, "_TV")
                          dt <- gwr_toSF(gwr_list[[x]]) %>%
                            select(all_of(c(predictorVar, predTVcol,
                                            "Local_R2", "geometry"))) %>%
                            rownames_to_column(var = "Barcode") %>%
                            data.frame()
                          predictorVar <- ifelse(length(predictorVar) > 1,
                                                 paste(predictorVar, collapse = '+'),
                                                 predictorVar)
                          rename_with(dt, ~ paste0(predictorVar, "_R2", recycle0 = TRUE),
                                      all_of("Local_R2"))
                        })

  message("\tFiltering GWR results...")
  ## Filter for areas where senMayo is above the sen_cutoff
  senMayoScores <- colData(msfe@sfe_data[[id]])[,c("Barcode", "senMayo")] %>%
    as.data.frame()

  ### Fetch GWR results
  gwr_dt <- gwr_dt_list %>%
    reduce(left_join, by = c("Barcode", "geometry"), suffix = c("", "_ATMyo"))
  colnames(gwr_dt) <- gsub("q05cell_abundance_w_sf_", "", colnames(gwr_dt))
  gwr_dt <- gwr_dt[gwr_dt$Barcode %in% senMayoScores$Barcode, ]

  ### Fetch gene expression
  message("\tFetching gene expression data...")
  geneExpression <- SummarizedExperiment::assay(msfe@sfe_data[[id]], "sct")[sub_genesENSG, senMayoScores$Barcode] %>%
    t() %>%
    data.frame() %>%
    rownames_to_column(var = "Barcode")
  colnames(geneExpression) <- c("Barcode", names(sub_genesENSG))

  ### Fetch cell abundance
  message("\tFetching cell abundance data...")
  cell_abund <- cellAbundance[[id]][c("Barcode", sub_cells)]
  cell_abund <- cell_abund[cell_abund$Barcode %in% senMayoScores$Barcode, ]

  ### Fetch annotations
  message("\tFetching annotation data...")
  spot_annot <- colData(msfe@sfe_data[[id]])["annotation"] %>%
    data.frame() %>%
    rownames_to_column(var = "Barcode")
  spot_annot <- spot_annot[spot_annot$Barcode %in% senMayoScores$Barcode, ]

  ### Fetch Image data
  message("\tFetching image data...")
  image <- imgRaster(getImg(msfe@sfe_data[[id]],
                            sample_id = id,
                            image_id = "lowres"))
  max <- image@ptr$range_max[1]
  if (is.na(max)) {
    max <- max(image[,,1])
  }
  if (max <= 1) {
    max_col_value <- 1
  } else if (max <= 255) {
    max_col_value <- 255
  } else if (max <= 65536) {
    max_col_value <- 65536
  }
  coords <- spatialCoords(msfe@sfe_data[[id]])[gwr_dt$Barcode,] %>%
    matrix(ncol = 2)
  limits_list <- list(y_lims = c(min(coords[,2]), max(coords[,2])),
                      x_lims = c(min(coords[,1]), max(coords[,1])))

  ### Merge expression and cell abundance with GWR output
  message("\tMerging data...")
  data_to_plot <- left_join(senMayoScores, gwr_dt, by = "Barcode") %>%
    left_join(geneExpression, by = "Barcode") %>%
    left_join(cell_abund, by = "Barcode", suffix = c("", ".y")) %>%
    left_join(spot_annot, by = "Barcode")
  colnames(data_to_plot) <- gsub("q05cell_abundance_w_sf_", "abundance_", colnames(data_to_plot))

  ## Output to the lists
  message("\tOutputing to list...")
  data_to_plot_lst[[id]] <- data_to_plot
  image_to_plot_lst[[id]] <- list(image = image,
                                  max_col_value = max_col_value,
                                  limits_list = limits_list)
}

#### ---- Subset - high SenMayo senescence module score ------------------------
## Here we will plot gene expression, cell abundance and SA results only for
## areas with high SenMayo senescence module score.
##### ----Plot coefficients ----------------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/Combo/senMayo_SenPositive/")
dir.create(graphics_out, recursive = TRUE)

## Select predictor vars to plot
predictorVars <- c("ecm", "emt", "tgfb", "AT2", "Myofibroblasts")

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Export data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(senMayo > sen_cutoff)
  image <- image_to_plot_lst[[id]][["image"]]
  max_col_value <- image_to_plot_lst[[id]][["max_col_value"]]
  limits_list <- image_to_plot_lst[[id]][["limits_list"]]

  message("\tPlotting...")
  ### Master plot -------------------#
  p <- ggplot(data_to_plot,
              aes(geometry = geometry)) +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value,
                                   alpha = 0.1) +
    ggplot2::lims(x = limits_list[[2]],
                  y = limits_list[[1]]) +
    theme_void()

  ### Plot GWR coefficient ----------#
  plot_gwr <- vector(mode = "list", length = length(predictorVars))
  names(plot_gwr) <- predictorVars

  for (c in predictorVars) {
    lg <- paste0(c," ")
    plot_gwr[[c]] <- p +
      geom_sf(aes(fill = !!sym(c))) + # force evaluate aes
      scale_fill_gradient2(high = "#880A1F",
                           mid = "#FAFBF7",
                           low = "#426E92",
                           midpoint = 0,
                           n.breaks = 7) +
      labs(fill = bquote(.(lg) ~ beta[1]))
  }

  p1 <- patchwork::wrap_plots(plot_gwr, ncol = 1)

  ### Plot GWR T-values ----------#
  plot_tv <- vector(mode = "list", length = length(predictorVars))
  names(plot_tv) <- predictorVars

  for (c in predictorVars) {
    tv <- paste0(c,"_TV")
    plot_tv[[c]] <- p +
      geom_sf(aes(fill = !!sym(tv))) + # force evaluate aes
      scale_fill_gradient2(high = "#8C510A",
                           mid = "#F5F5F5",
                           low = "#01665E",
                           midpoint = 0,
                           n.breaks = 7) +
      labs(fill = tv)
  }

  p2 <- patchwork::wrap_plots(plot_tv, ncol = 1)

  ### Plot GWR local R2 ----------#
  plot_r2 <- vector(mode = "list", length = length(predictorVars))
  names(plot_r2) <- predictorVars

  for (c in predictorVars) {
    lg <- paste0(c, " ")
    r2 <- paste0(c,"_R2")
    plot_r2[[c]] <- p +
      geom_sf(aes(fill = !!sym(r2))) + # force evaluate aes
      scale_fill_viridis_c() +
      labs(fill = bquote(.(lg) ~ R^2))
  }

  p3 <- patchwork::wrap_plots(plot_r2, ncol = 1)

  ### Plot cell abundance --------#
  plot_cellAbund <- vector(mode = "list", length = length(sub_cells))
  names(plot_cellAbund) <- names(sub_cells)

  for (c in names(sub_cells)) {
    ab <- paste0("abundance_", c)
    plot_cellAbund[[c]] <- p +
      geom_sf(aes(fill = !!sym(ab))) + # force evaluate aes
      scale_fill_viridis_c(option = "turbo") +
      labs(fill = paste(c, "Density", sep = "\n"))
  }

  p4 <- patchwork::wrap_plots(plot_cellAbund, ncol = 1)

  ### Plot annotation ------------#
  # p5 <- p +
  #   geom_sf(aes(fill = annotation)) +
  #   scale_fill_manual(values = colours) + # colours as defined in the annotation plots earlier
  #   labs(fill = "Annotation")

  p1 | p2 | p3 | p4

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 19.16,
         height = 16,
         units = "in",
         dpi = 300)

  ### Plot annotation ------------#
  p +
    geom_sf(aes(fill = annotation)) +
    scale_fill_manual(values = colours) + # colours as defined in the annotation plots earlier
    labs(fill = "Annotation")

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo_annotation"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 6.09,
         height = 5.41,
         units = "in",
         dpi = 300)
}

##### ----Plot gene expression -------------------------------------------------
## Set folder
graphics_out <- paste0(projectFolder, "/graphics_out/Combo/senMayo_SenPositive/")

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Export data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(senMayo > sen_cutoff)
  image <- image_to_plot_lst[[id]][["image"]]
  max_col_value <- image_to_plot_lst[[id]][["max_col_value"]]
  limits_list <- image_to_plot_lst[[id]][["limits_list"]]

  message("\tPlotting...")
  ### Master plot
  p <- ggplot(data_to_plot,
              aes(geometry = geometry)) +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value,
                                   alpha = 0.1) +
    # ggplot2::lims(x = limits_list[[2]],
    #               y = limits_list[[1]]) +
    theme_void()

  ### Plot Gene expression
  plot_geneExpr <- vector(mode = "list", length = length(sub_genesENSG))
  names(plot_geneExpr) <- names(sub_genesENSG)

  for (g in names(sub_genesENSG)) {
    plot_geneExpr[[g]] <- p +
      geom_sf(aes(fill = !!sym(g))) + # force evaluate aes
      scale_fill_viridis_c() +
      labs(fill = paste(g, "Norm. Log2\ncounts", sep = "\n"))
  }

  p1 <- patchwork::wrap_plots(plot_geneExpr, nrow = 4, ncol = 3)

  ### Plot SenMayo senescence enrichment score
  p2 <- p +
    geom_sf(aes(fill = senMayo)) +
    scale_fill_gradient2(high = "#762A83",
                         mid = "#E8E8E8",
                         low = "#246C27",
                         midpoint = 0,
                         n.breaks = 7) +
    labs(fill = "SenMayo\nmodule\nscore")

  # c("#83488B", "#C082C7", "#E5C5EA", "#E8E8E8", "#B4D8B5", "#61A863", "#246C27")
  # c("#762A83", "#AF8DC3", "#E7D4E8", "#F7F7F7", "#D9F0D3", "#7FBF7B", "#1B7837")
  # c("#C51B7D", "#E9A3C9", "#FDE0EF", "#F7F7F7", "#E6F5D0", "#A1D76A", "#4D9221")

  p1 | p2

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo_geneExpr"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 19.16,
         height = 12,
         units = "in",
         dpi = 300)
}

##### ----Plot rain-cloud cell abundance and gene expression -------------------
## Set folder
graphics_out <- paste0(projectFolder, "/graphics_out/Combo/senMayo_SenPositive/")

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot cell abundance as rain-cloud plots ----------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    select(matches("abund|senM")) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "cell_type",
                 values_to = "cell_density") %>%
    filter(sen_group != "no") %>%
    mutate(cell_type = gsub("abundance_|.Cells", "", .data$cell_type))

  message("\tPlotting rain-cloud Cell density...")
  ggplot(data_to_plot,
         aes(x = sen_group, y = cell_density, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "Cell Densities per Group",
         y = "Absolute Cell Density",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    facet_wrap(~cell_type, nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_CellDensities"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 10.83,
         height = 6.83,
         units = "in",
         dpi = 300)
}

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot gene expression as rain-cloud plots ---------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    select(all_of(c(sub_genes, "senMayo"))) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  message("\tPlotting rain-cloud Gene Expression...")
  ggplot(data_to_plot,
         aes(x = sen_group, y = logcounts, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "",
         y = "Normalised Log2-counts",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    facet_wrap(~gene, nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_GeneExpression"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 10.83,
         height = 4.83,
         units = "in",
         dpi = 300)
}

##### ----Plot Bar graphs cell abundance and gene expression -------------------
## Set folder
graphics_out <- paste0(projectFolder, "/graphics_out/Combo/senMayo_SenPositive/")

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot cell abundance as bar-graph plots -----------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    select(matches("abund|senM")) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "cell_type",
                 values_to = "cell_density") %>%
    filter(sen_group != "no") %>%
    mutate(cell_type = gsub("abundance_|.Cells", "", .data$cell_type))

  ## Run statistical test
  ### For Positive vs Negative senescence within the AT2/Myof. groups
  s_test_senGroup <- data_to_plot %>%
    mutate(cell_type = as.factor(cell_type)) %>%
    group_by(cell_type) %>%
    rstatix::wilcox_test(cell_density ~ sen_group) %>%
    add_xy_position(fun = "mean_se", x = "cell_type", dodge = 0.8) %>%
    add_significance("p")

  group_numbers <- paste0("n = ", as.vector(as.matrix(s_test_senGroup[1, c("n1", "n2")])))

  ## Plot
  message("\tPlotting Bar-graph cell densities...")
  ggbarplot(data_to_plot,
            x = "cell_type", y = "cell_density", fill = "sen_group",
            position = position_dodge(0.9),
            add = "mean_se") +
    stat_pvalue_manual(s_test_senGroup,  label = "{p.signif}",
                       tip.length = 0.0025) +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    labs(title = "",
         y = "Mean Absolute Cell Density",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore")

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_barGraph_"
  other <- "perGroup_CellDensities"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 5.1,
         height = 6.6,
         units = "in",
         dpi = 300)
}

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot gene expression as rain-cloud plots ---------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    select(all_of(c(sub_genes, "senMayo"))) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  ## Run statistical test
  ### For Positive vs Negative senescence within the AT2/Myof. groups
  # s_test_senGroup <- data_to_plot %>%
  #   mutate(gene = as.factor(gene)) %>%
  #   group_by(gene) %>%
  #   rstatix::wilcox_test(logcounts ~ sen_group) %>%
  #   add_xy_position(fun = "mean_se", x = "gene", dodge = 0.8) %>%
  #   add_significance("p")
  #
  # group_numbers <- paste0("n = ", as.vector(as.matrix(s_test_senGroup[1, c("n1", "n2")])))

  ## Plot
  message("\tPlotting bar-graph Gene Expression...")
  ggbarplot(data_to_plot,
            x = "sen_group", y = "logcounts", fill = "sen_group",
            position = position_dodge(0.9),
            add = "mean_se",
            facet.by = "gene",
            scale = "free_y") +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005,
                               vjust = "-1") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    labs(title = "",
         y = "Mean Normalised Log-counts",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # stat_pvalue_manual(s_test_senGroup,  label = "{p.signif}",
  #                    tip.length = 0.0025, hide.ns = TRUE) +

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_barGraph_"
  other <- "perGroup_GeneExpression"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 7.5,
         height = 5.55,
         units = "in",
         dpi = 300)
}

#### ----Subset - high GWR SenMayo - ECM scores correlation --------------------
## Here we will plot gene expression, cell abundance and SA results only for
## areas with high positive correlation between SenMayo and ECM module scores.
##### ----Set cut-off ----------------------------------------------------------
gwr_cutoff <- 0.5
sen_cutoff # the 75th percentile cut-off set earlier

## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/Combo/senMayo-ecm_high_coef/")
dir.create(graphics_out, recursive = TRUE)

##### ----Plot coefficients ----------------------------------------------------
## Select predictor vars to plot
predictorVars <- c("ecm", "emt", "tgfb", "AT2", "Myofibroblasts")

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Export data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff)
  image <- image_to_plot_lst[[id]][["image"]]
  max_col_value <- image_to_plot_lst[[id]][["max_col_value"]]
  limits_list <- image_to_plot_lst[[id]][["limits_list"]]

  message("\tPlotting...")
  ### Master plot -------------------#
  p <- ggplot(data_to_plot,
              aes(geometry = geometry)) +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value,
                                   alpha = 0.1) +
    ggplot2::lims(x = limits_list[[2]],
                  y = limits_list[[1]]) +
    theme_void()

  ### Plot GWR coefficient ----------#
  plot_gwr <- vector(mode = "list", length = length(predictorVars))
  names(plot_gwr) <- predictorVars

  for (c in predictorVars) {
    lg <- paste0(c," ")
    plot_gwr[[c]] <- p +
      geom_sf(aes(fill = !!sym(c))) + # force evaluate aes
      scale_fill_gradient2(high = "#880A1F",
                           mid = "#FAFBF7",
                           low = "#426E92",
                           midpoint = 0,
                           n.breaks = 7) +
      labs(fill = bquote(.(lg) ~ beta[1]))
  }

  p1 <- patchwork::wrap_plots(plot_gwr, ncol = 1)

  ### Plot GWR T-values ----------#
  plot_tv <- vector(mode = "list", length = length(predictorVars))
  names(plot_tv) <- predictorVars

  for (c in predictorVars) {
    tv <- paste0(c,"_TV")
    plot_tv[[c]] <- p +
      geom_sf(aes(fill = !!sym(tv))) + # force evaluate aes
      scale_fill_gradient2(high = "#8C510A",
                           mid = "#F5F5F5",
                           low = "#01665E",
                           midpoint = 0,
                           n.breaks = 7) +
      labs(fill = tv)
  }

  p2 <- patchwork::wrap_plots(plot_tv, ncol = 1)

  ### Plot GWR local R2 ----------#
  plot_r2 <- vector(mode = "list", length = length(predictorVars))
  names(plot_r2) <- predictorVars

  for (c in predictorVars) {
    lg <- paste0(c, " ")
    r2 <- paste0(c,"_R2")
    plot_r2[[c]] <- p +
      geom_sf(aes(fill = !!sym(r2))) + # force evaluate aes
      scale_fill_viridis_c() +
      labs(fill = bquote(.(lg) ~ R^2))
  }

  p3 <- patchwork::wrap_plots(plot_r2, ncol = 1)

  ### Plot cell abundance --------#
  plot_cellAbund <- vector(mode = "list", length = length(sub_cells))
  names(plot_cellAbund) <- names(sub_cells)

  for (c in names(sub_cells)) {
    ab <- paste0("abundance_", c)
    plot_cellAbund[[c]] <- p +
      geom_sf(aes(fill = !!sym(ab))) + # force evaluate aes
      scale_fill_viridis_c(option = "turbo") +
      labs(fill = paste(c, "Density", sep = "\n"))
  }

  p4 <- patchwork::wrap_plots(plot_cellAbund, ncol = 1)

  ### Plot annotation ------------#
  # p5 <- p +
  #   geom_sf(aes(fill = annotation)) +
  #   scale_fill_manual(values = colours) + # colours as defined in the annotation plots earlier
  #   labs(fill = "Annotation")

  p1 | p2 | p3 | p4

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo-ECM-high_coef"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 19.16,
         height = 16,
         units = "in",
         dpi = 300)

  ### Plot annotation ------------#
  p +
    geom_sf(aes(fill = annotation)) +
    scale_fill_manual(values = colours) + # colours as defined in the annotation plots earlier
    labs(fill = "Annotation")

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo_annotation"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 6.09,
         height = 5.41,
         units = "in",
         dpi = 300)
}

##### ----Plot gene expression -------------------------------------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Export data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff)
  image <- image_to_plot_lst[[id]][["image"]]
  max_col_value <- image_to_plot_lst[[id]][["max_col_value"]]
  limits_list <- image_to_plot_lst[[id]][["limits_list"]]

  message("\tPlotting...")
  ### Master plot
  p <- ggplot(data_to_plot,
              aes(geometry = geometry)) +
    tidyterra::geom_spatraster_rgb(data = image,
                                   max_col_value = max_col_value,
                                   alpha = 0.1) +
    # ggplot2::lims(x = limits_list[[2]],
    #               y = limits_list[[1]]) +
    theme_void()

  ### Plot Gene expression
  plot_geneExpr <- vector(mode = "list", length = length(sub_genesENSG))
  names(plot_geneExpr) <- names(sub_genesENSG)

  for (g in names(sub_genesENSG)) {
    plot_geneExpr[[g]] <- p +
      geom_sf(aes(fill = !!sym(g))) + # force evaluate aes
      scale_fill_viridis_c() +
      labs(fill = paste(g, "Norm. Log2\ncounts", sep = "\n"))
  }

  p1 <- patchwork::wrap_plots(plot_geneExpr, nrow = 4, ncol = 3)

  ### Plot SenMayo senescence enrichment score
  p2 <- p +
    geom_sf(aes(fill = senMayo)) +
    scale_fill_gradient2(high = "#762A83",
                         mid = "#E8E8E8",
                         low = "#246C27",
                         midpoint = sen_cutoff,
                         n.breaks = 7) +
    labs(fill = "SenMayo\nmodule\nscore")

  # c("#83488B", "#C082C7", "#E5C5EA", "#E8E8E8", "#B4D8B5", "#61A863", "#246C27")
  # c("#762A83", "#AF8DC3", "#E7D4E8", "#F7F7F7", "#D9F0D3", "#7FBF7B", "#1B7837")
  # c("#C51B7D", "#E9A3C9", "#FDE0EF", "#F7F7F7", "#E6F5D0", "#A1D76A", "#4D9221")

  p1 | p2

  message("\tSaving ...")
  prfx <- id
  main <- "_GWR"
  sfx <- ""
  other <- "_subset_senMayo_geneExpr"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 19.16,
         height = 12,
         units = "in",
         dpi = 300)
}

##### ----Plot rain-cloud cell abundance and gene expression -------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot cell abundance as rain-cloud plots ----------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff) %>%
    select(matches("abund|senM")) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "cell_type",
                 values_to = "cell_density") %>%
    filter(sen_group != "no") %>%
    mutate(cell_type = gsub("abundance_|.Cells", "", .data$cell_type))

  message("\tPlotting rain-cloud Cell density...")
  ggplot(data_to_plot,
         aes(x = sen_group, y = cell_density, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "Cell Densities per Group",
         y = "Absolute Cell Density",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    facet_wrap(~cell_type, nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_CellDensities"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 10.83,
         height = 6.83,
         units = "in",
         dpi = 300)
}

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot gene expression as rain-cloud plots ---------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff) %>%
    select(all_of(c(sub_genes, "senMayo"))) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  message("\tPlotting rain-cloud Gene Expression...")
  ggplot(data_to_plot,
         aes(x = sen_group, y = logcounts, fill = sen_group)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.color = NA) +
    # ggdist::stat_dots(side = "left", justification = 1.1, binwidth = NA, overflow = "compress") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005) +
    labs(title = "",
         y = "Normalised Log2-counts",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    facet_wrap(~gene, nrow = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15))

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_rainCloud_"
  other <- "perGroup_GeneExpression"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 10.83,
         height = 4.83,
         units = "in",
         dpi = 300)
}

##### ----Plot Bar graphs cell abundance and gene expression -------------------
for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot cell abundance as bar-graph plots -----------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff) %>%
    select(matches("abund|senM")) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "cell_type",
                 values_to = "cell_density") %>%
    filter(sen_group != "no") %>%
    mutate(cell_type = gsub("abundance_|.Cells", "", .data$cell_type))

  if (length(levels(data_to_plot$sen_group)) != 2) {
    next
  }

  ## Run statistical test
  ### For Positive vs Negative senescence within the AT2/Myof. groups
  s_test_senGroup <- data_to_plot %>%
    mutate(cell_type = as.factor(cell_type)) %>%
    group_by(cell_type) %>%
    rstatix::wilcox_test(cell_density ~ sen_group) %>%
    add_xy_position(fun = "mean_se", x = "cell_type", dodge = 0.8) %>%
    add_significance("p")

  group_numbers <- paste0("n = ", as.vector(as.matrix(s_test_senGroup[1, c("n1", "n2")])))

  ## Plot
  message("\tPlotting Bar-graph cell densities...")
  ggbarplot(data_to_plot,
            x = "cell_type", y = "cell_density", fill = "sen_group",
            position = position_dodge(0.9),
            add = "mean_se") +
    stat_pvalue_manual(s_test_senGroup,  label = "{p.signif}",
                       tip.length = 0.0025) +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    labs(title = "",
         y = "Mean Absolute Cell Density",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore")

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_barGraph_"
  other <- "perGroup_CellDensities"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 5.1,
         height = 6.6,
         units = "in",
         dpi = 300)
}

for (j in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot gene expression as bar-graph plots ---------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    filter(ecm > gwr_cutoff) %>%
    select(all_of(c(sub_genes, "senMayo"))) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  if (length(levels(data_to_plot$sen_group)) != 2) {
    next
  }

  id <- sampleNames_selected[j]
  message("Working on sample: ", id, " (", j, "/", length(sampleNames_selected), ")")

  ## Plot gene expression as rain-cloud plots ---------------------------------#
  ## Export cell type data from lists
  message("\tExporting data from list...")
  data_to_plot <- data_to_plot_lst[[id]] %>%
    select(all_of(c(sub_genes, "senMayo"))) %>%
    mutate(sen_group = case_when(senMayo > sen_cutoff ~ "High Senescence",
                                 senMayo < sen_cutoff ~ "Low Senescence",
                                 .default = "no"),
           sen_group = fct_relevel(sen_group,
                                   "Low Senescence", "High Senescence")) %>%
    pivot_longer(cols = -c("sen_group", "senMayo"),
                 names_to = "gene",
                 values_to = "logcounts") %>%
    filter(sen_group != "no")

  ## Plot
  message("\tPlotting bar-graph Gene Expression...")
  ggbarplot(data_to_plot,
            x = "sen_group", y = "logcounts", fill = "sen_group",
            position = position_dodge(0.9),
            add = "mean_se",
            facet.by = "gene",
            scale = "free_y") +
    ggpubr::stat_compare_means(ref.group = "Low Senescence",
                               method = "wilcox.test",
                               label = "p.signif",
                               tip.length = 0.005,
                               vjust = "-1") +
    scale_fill_manual(values = c("#246C27", "#762A83")) +
    labs(title = "",
         y = "Mean Normalised Log-counts",
         x = "",
         fill = "SenMayo\nsenescence\nmodule\nscore") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  message("\tSaving ...")
  prfx <- id
  main <- ""
  sfx <- "_barGraph_"
  other <- "perGroup_GeneExpression"

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = 7.5,
         height = 5.55,
         units = "in",
         dpi = 300)
}
