## Fit mean-variance relationship
dec2 <- modelGeneVariance(msfe, sample_id = TRUE, method = "Var")

## Select top HVGs
top_hvgs2 <- getTopHighVarGenes(dec2,
                               var.field = "bio",
                               prop = 1,
                               var.threshold = 0,
                               fdr.threshold = 0.05)

## Examine overlap in HVGs
# library(ggVennDiagram)
x2 <- list(JBO019=top_hvgs2$JBO019,
           JBO018=top_hvgs2$JBO018)
ggVennDiagram(x2,label = "count") + scale_fill_gradient(low="grey90",high = "red")

x3 <- list(JBO019=top_hvgs$JBO019,
           JBO019_2=top_hvgs2$JBO019)
ggVennDiagram(x3,label = "count") + scale_fill_gradient(low="grey90",high = "red")

x4 <- list(JBO018=top_hvgs$JBO018,
           JBO018_2=top_hvgs2$JBO018)
ggVennDiagram(x4,label = "count") + scale_fill_gradient(low="grey90",high = "red")

plotGeneExpression(sfe, genes = t2g$gene[t2g$gene %in% top_hvgs2$JBO019], assay = "logcounts")

t2g_liverCellTypes <- getTerm2Gene(msig_data = msigdb, cat = "C8", subcat = "") %>%
  filter(grepl("AIZARANI_LIVER", term)) %>%
  mutate(term = gsub("AIZARANI_LIVER_C[0-9]+_", "", term))

msigdb_hepSteatosis <- msigdb[msigdb$gs_cat == "C5" & msigdb$gs_subcat == "HPO" & grepl("HEPATIC_STEATOSIS", msigdb$gs_name),]
t2g_hepSteat <- data.frame(term = msigdb_hepSteatosis$gs_name,
                           gene = msigdb_hepSteatosis$human_ensembl_gene) %>%
  mutate(term = gsub("HP_|_HEPATIC_STEATOSIS", "", term))

steatosis_gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                                 pc = 1,
                                                 genes_no = 2,
                                                 NES = 0,
                                                 minGSSize = 5,
                                                 pvalueCutoff = 0.25,
                                                 TERM2GENE = t2g_hepSteat,
                                                 pAdjustMethod = "fdr",
                                                 scoreType = "std",
                                                 nPermSimple = 10000,
                                                 mc.cores = 8)
plotGWPCA_FuncCLust(steatosis_gsea_map, count = 5, legend = "right", legend.title = "Pathways")

plotGeneExpression(sfe, genes = unique(t2g_hepSteat$gene), assay = "logcounts")



## Expression of steatosis signature genes:
steatosis_gs <- t2g$gene

sfe <- getSFE(msfe, "JBO019")

steatosis_jbo019_logCounts <- SummarizedExperiment::assay(sfe, "logcounts")
steatosis_gs_notInST <- steatosis_gs[!steatosis_gs %in% rownames(steatosis_jbo019_logCounts)]
steatosis_gs_inST <- steatosis_gs[steatosis_gs %in% rownames(steatosis_jbo019_logCounts)]

steatosis_gs_expr <- steatosis_jbo019_logCounts[steatosis_gs_inST,]

steatosis_gs_avExpr <- colSums(data.frame(steatosis_gs_expr)) %>%
  data.frame(Steatosis_Average = .) %>%
  rownames_to_column(.) %>%
  mutate(rowname = gsub("[.]", "-", rowname))

geoms_steatotic <- colGeometry(sfe, "spotHex") %>%
  rownames_to_column(.)

steatosis_gs_avExpr <- steatosis_gs_avExpr %>%
  left_join(geoms_steatotic, by = "rowname")

ggplot(data = steatosis_gs_avExpr) +
  geom_violin(aes(x = "Steatosis", y = Steatosis_Average))

ggplot(data = steatosis_gs_avExpr,
       aes(fill = Steatosis_Average)) +
  geom_sf(aes(geometry = geometry)) +
  scale_fill_viridis_c()

t2g = read.csv("./data/t2g_files/steatosis_RNA-Seq_t2g.csv")
steatosis_t2g <- t2g[t2g$gene %in% steatosis_gs_inST,]

steatosis_gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                                 pc = 1,
                                                 genes_no = 2,
                                                 NES = 0,
                                                 minGSSize = 5,
                                                 pvalueCutoff = 0.25,
                                                 TERM2GENE = steatosis_t2g,
                                                 pAdjustMethod = "fdr",
                                                 scoreType = "std",
                                                 nPermSimple = 10000,
                                                 mc.cores = 8)
plotGWPCA_FuncCLust(steatosis_gsea_map, count = 5, legend = "right", legend.title = "Pathways")


steatosis_genes_filt <- steatosis_genes %>%
  filter(adj.P.Val < 0.05)

steatosis_gs_expr_filt <- steatosis_jbo019_logCounts[rownames(steatosis_jbo019_logCounts) %in% steatosis_genes_filt$gene,]

steatosis_gs_avExpr_filt <- colSums(data.frame(steatosis_gs_expr_filt)) %>%
  data.frame(Steatosis_Average = .) %>%
  rownames_to_column(.) %>%
  mutate(rowname = gsub("[.]", "-", rowname)) %>%
  left_join(geoms_steatotic, by = "rowname")

steatosis_t2g <- t2g[t2g$gene %in% steatosis_gs_inST,]

gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 2,
                                       NES = 0,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = steatosis_t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 8)
plotGWPCA_FuncCLust(gsea_map, count = 2, legend = "right", legend.title = "Pathways")

ggplot(data = steatosis_gs_avExpr_filt) +
  geom_violin(aes(x = "Steatosis", y = Steatosis_Average))

ggplot(data = steatosis_gs_avExpr_filt,
       aes(fill = Steatosis_Average)) +
  geom_sf(aes(geometry = geometry)) +
  scale_fill_viridis_c()

steatosis_gsea_enrich <- c("ENSG00000101057", "ENSG00000162073", "ENSG00000171234")

plotGeneExpression(sfe, genes = steatosis_gsea_enrich, assay = "logcounts")

pcagw_19_genes <- colnames(pcagw_19$loadings)

sum(pcagw_19_genes %in% t2g$gene)

plotGeneExpression(sfe,
                   genes = pcagw_19_genes[pcagw_19_genes %in% t2g$gene],
                   assay = "logcounts")

## GSEA with liver cell types MSigDB gene sets
t2g_liverCellTypes <- getTerm2Gene(msig_data = msigdb, cat = "C8", subcat = "") %>%
  filter(grepl("AIZARANI_LIVER", term)) %>%
  mutate(term = gsub("AIZARANI_LIVER_C[0-9]+_", "", term))

steatosis_gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                                 pc = 1,
                                                 genes_no = 2,
                                                 NES = 0,
                                                 minGSSize = 5,
                                                 pvalueCutoff = 0.25,
                                                 TERM2GENE = t2g_liverCellTypes,
                                                 pAdjustMethod = "fdr",
                                                 scoreType = "std",
                                                 nPermSimple = 10000,
                                                 mc.cores = 8)

plotGWPCA_FuncCLust(steatosis_gsea_map, count = 1, legend = "right", legend.title = "Pathways")

## GSEA with liver cell types MSigDB gene sets merged
t2g_liverCellTypes2 <- getTerm2Gene(msig_data = msigdb, cat = "C8", subcat = "") %>%
  filter(grepl("AIZARANI_LIVER", term)) %>%
  mutate(term = gsub("AIZARANI_LIVER_C[0-9]+_", "", term)) %>%
  mutate(term = gsub("_[0-9]", "", term)) %>%
  filter(!duplicated(.))

liverCT_gsea_map <- list()
for (i in 1:3) {
  gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                         pc = i,
                                         genes_no = 2,
                                         NES = 0,
                                         minGSSize = 5,
                                         pvalueCutoff = 0.25,
                                         TERM2GENE = t2g_liverCellTypes2,
                                         pAdjustMethod = "fdr",
                                         scoreType = "std",
                                         nPermSimple = 10000,
                                         mc.cores = 8)

  liverCT_gsea_map[[i]] <- gsea_map

  print(plotGWPCA_FuncCLust(gsea_map, count = 1, legend = "right", legend.title = "Pathways"))
}

# To compute the matrix of p-value, a custom R function is used :
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

corr_list <- list()

for (s in samples) {
  pcagw <- get(ifelse(s == "JBO019", "pcagw_19", "pcagw_18"))

  clusters <- fgwc_list[[s]][["finaldata"]] %>%
    dplyr::filter(cluster == 2)
  barcodes <- rownames(clusters)

  leading_scores <- pcagw[["loadings"]][rownames(pcagw$loadings) %in% barcodes,,"PC1"]

  corrM <- cor(t(leading_scores))
  # matrix of the p-value of the correlation
  p.mat <- cor.mtest(t(leading_scores))

  col <- colorRampPalette(c("darkblue", "white", "red"))(200)

  corr_list[[s]] <- corrplot::corrplot(corr = corrM,
                                       method = "color",
                                       type = "upper",
                                       order = "hclust",
                                       p.mat = p.mat,
                                       sig.level = 0.01,
                                       insig = "blank",
                                       col = col,
                                       tl.pos = "n")

  corrPos <- corr_list[[s]][["corrPos"]] %>%
    dplyr::filter(p.value < 0.01 & corr < -0.7)

  corrPos_g1 <- corrPos %>%
    dplyr::select(xName) %>%
    unique() %>%
    dplyr::mutate(Group = 1) %>%
    dplyr::rename("Barcode" = xName)

  corrPos_g2 <- corrPos %>%
    dplyr::select(yName) %>%
    unique() %>%
    dplyr::mutate(Group = 2) %>%
    dplyr::rename("Barcode" = yName)

  corrPos_g3 <- corr_list[[s]][["corrPos"]] %>%
    dplyr::select(xName) %>%
    unique() %>%
    dplyr::filter(!(xName %in% c(corrPos_g1$Barcode, corrPos_g2$Barcode)))  %>%
    dplyr::mutate(Group = 3) %>%
    dplyr::rename("Barcode" = xName)

  corrPos_dt <- rbind(corrPos_g1, corrPos_g2, corrPos_g3)
  rownames(corrPos_dt) <- corrPos_dt$Barcode
  geoms <- data.frame(clusters["geometry"])
  rownames(geoms) <- rownames(clusters)

  corrPos_dt <- corrPos_dt %>%
    merge(., geoms, by = 0) %>%
    column_to_rownames(var = "Row.names")

  ggplot(corrPos_dt,
         aes(fill = as.factor(Group))) +
    geom_sf(aes(geometry = geometry)) +
    scale_fill_manual(values = cols_single) +
    labs(title = s,
         fill = paste0(pc,"\n", "Leading\nGene")) +
    theme_void()

}


ensgs <- c("ENSG00000105697", "ENSG00000107317", "ENSG00000135094")
names(ensgs) <- c("HAMP", "PTGDS", "SDS")

plotGWPCA_leadingScores(pcagw_19, 1, ensgs, colours = "custom", col_args = list(high = "red", mid = "white", low = "darkblue", midpoint = 0), gene_names = names(ensgs), absolute = FALSE)

ggplot(data, aes(fill = PC1_scores)) +
  geom_sf(aes(geometry = geometry)) +
  scale_fill_gradient2(high = "red",
                       mid = "white",
                       low = "darkblue",
                       midpoint = 0)

gene_data <- data.frame(gene_a = sample(1:20, 5),
                        gene_b = sample(1:20, 5),
                        gene_c = sample(1:20, 5),
                        gene_d = sample(1:20, 5))
x <- as.matrix(gene_data)
dMat <- sfe@metadata[["dMat"]][["euclidean"]][sfe@metadata[["dMat"]][["euclidean"]][,1] < bw,sfe@metadata[["dMat"]][["euclidean"]][,1] < bw]
for (i in 1:5) {
  wt <- GWmodel::gw.weight(dMat[,i], bw = bw, kernel = "gaussian", adaptive = FALSE)
  ## Centre the gene counts by subtracting a distance-weighted mean
  gene_data_centered <- t(t(x) - base::colSums(x * wt) / sum(wt))
  ## Multiply the centred matrix by the square root of the weights
  gene_data_weighted <- gene_data_centered * sqrt(wt)
  ## Perform Single Value Decomposition
  gene_data_svd <- svd(gene_data_weighted, nu = 0, nv = 2)

  # gene_data_svd$d
  print(gene_data_svd$v)
  print(gene_data_weighted)
}

res <- kegg_results %>%
  group_by(Description) %>%
  summarise(
    Frequency = n(),  # Count how many times the pathway is enriched
    Mean_pvalue = mean(p.adjust, na.rm = TRUE),
    Median_pvalue = median(p.adjust, na.rm = TRUE),
    Locations = paste(unique(location), collapse = ";")
  ) %>%
  arrange(Mean_pvalue)

