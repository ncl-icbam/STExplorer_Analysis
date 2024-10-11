import scanpy as sc

adata = sc.read_h5ad('/mnt/HDD8_16/STExplorer/ProstateCellAtlas/prostate_portal_300921.h5ad')

import scvi

import seaborn as sns

import pandas as pd

####Preprocessing####

ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"

ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)

adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

adata.var['mt'] = adata.var.index.str.startswith('MT-')

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

adata.var.sort_values('n_cells_by_counts')

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True,save='_QC.svg')

sc.pp.filter_cells(adata, min_genes=100)

import numpy as np

adata = adata[adata.obs.pct_counts_mt < 15]

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=8000, subset = True, layer = 'counts',
                            flavor = "seurat_v3", batch_key="name")

scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                             categorical_covariate_keys=["name"],
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts'])

 model = scvi.model.SCVI(adata)

 model.train()

 adata.obsm['X_scVI'] = model.get_latent_representation()

 adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)

 sc.pp.neighbors(adata, use_rep = 'X_scVI')

sc.tl.umap(adata)

####Cell type labelling####

sc.tl.leiden(adata, resolution = 0.9)

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,save='_top20_markers.svg')

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]

markers_scvi = model.differential_expression(groupby = 'leiden')

markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]

adata.uns['scvi_markers']=markers_scvi
adata.uns['markers']=markers

sc.pl.umap(adata, color = ['leiden', 'celltype','mri_grading'], legend_loc = "on data",frameon = False,save="_X_scVI_leiden_name_celltype.svg")

###Manual cell annotation (for minor adjustments compared to original manuscript)####

cell_type={"0":"T-cell",
"1":"LE-KLK3",
"2":"CE",
"3":"LE-KLK4",
"4":"Sperm",
"5":"LE-KLK3",
"6":"LE-KLK3",
"7":"LE-KLK4",
"8":"Hillock-BC",
"9":"T-cell",
"10":"Endothelial",
"11":"MNP",
"12":"NK cell",
"13":"B cell",
"14":"Fibroblast",
"15":"Mast-cell"
}

adata.obs['celltype_new']=adata.obs.leiden.map(cell_type)

###Plotting####

sc.pl.umap(adata,color = ['leiden', 'celltype','mri_grading','FABP5','LTF','TAGLN'],ncols=3,cmap='viridis',legend_loc = "on data",frameon = False,vmax=2.5,layer = 'scvi_normalized',save="_Figure_leaders_umap.svg")

sc.pl.umap(adata,color = ['leiden', 'celltype','mri_grading','RPL17','FABP5','LRRC75A'],ncols=3,cmap='viridis',legend_loc = "on data",frameon = False,vmax=3,layer = 'scvi_normalized',save="_H1_5_GG4_leaders_umap.svg")

sc.pl.umap(adata,color = ['leiden', 'celltype','mri_grading','KLK4','NPY','PDLIM5'],ncols=3,cmap='viridis',legend_loc = "on data",frameon = False,vmax=3,layer = 'scvi_normalized',save="_H2_1_GG2_leaders_umap.svg")

sc.pl.umap(adata,color = ['leiden', 'celltype','mri_grading','MPC2','SPON2','MESP1'],ncols=3,cmap='viridis',legend_loc = "on data",frameon = False,vmax=3,layer = 'scvi_normalized',save="_H2_5_stroma_leaders_umap.svg")

sc.pl.dotplot(adata, ['ACPP','AZGP1','SMS','PCA3','PCAT4','KLK2','H2AFJ','TRGC1','TFF3','MYLK','ACTA2','TPM2','RPL17','FABP5','LRRC75A','IGFBP7','LGALS1','MGP','KLK4','NPY','PDLIM5','KLK2','TMPRSS2','KLK3','TAGLN','FLNA','MYL9','MYLK','TAGLN','ACTA2','LINC01088','TMEFF2','FABP5','MPC2','SPON2','MESP1','FLNA','MYLK','ACTA2','TAGLN','MYL9','FLNA'], var_group_positions=[(0,2),(3,5),(6,8),(9,11),(12,14),(15,17),(18,20),(21,23),(24,26),(27,29),(30,32),(32,35),(36,38),(39,41)],var_group_labels=['H1_4_T10_GG2','H1_4B10_GG4','H1_4T10_stroma','H1_4B10_GG2','H1_5T10_GG4','H1_5B10_GG4','H2_1T10_GG2','H2_1T10_GG4','H2_1B10_GG2','H2_1B10_GG4','H2_5T10_GG4','H2_5T10_Stroma','V1_2B10_Benign','V1_2B10_Stroma'], groupby="highest_GLEASON_score", show=True,dendrogram=False,save="_All_sign_group_gleason_dotplot.svg")

sc.pl.dotplot(adata, ['ACPP','AZGP1','SMS','PCA3','PCAT4','KLK2','H2AFJ','TRGC1','TFF3','MYLK','ACTA2','TPM2','RPL17','FABP5','LRRC75A','IGFBP7','LGALS1','MGP','KLK4','NPY','PDLIM5','KLK2','TMPRSS2','KLK3','TAGLN','FLNA','MYL9','MYLK','TAGLN','ACTA2','LINC01088','TMEFF2','FABP5','MPC2','SPON2','MESP1','FLNA','MYLK','ACTA2','TAGLN','MYL9','FLNA'], var_group_positions=[(0,2),(3,5),(6,8),(9,11),(12,14),(15,17),(18,20),(21,23),(24,26),(27,29),(30,32),(32,35),(36,38),(39,41)],var_group_labels=['H1_4_T10_GG2','H1_4B10_GG4','H1_4T10_stroma','H1_4B10_GG2','H1_5T10_GG4','H1_5B10_GG4','H2_1T10_GG2','H2_1T10_GG4','H2_1B10_GG2','H2_1B10_GG4','H2_5T10_GG4','H2_5T10_Stroma','V1_2B10_Benign','V1_2B10_Stroma'], groupby="name", show=True,dendrogram=False,save="_All_sign_group_gleason_dotplot.svg")