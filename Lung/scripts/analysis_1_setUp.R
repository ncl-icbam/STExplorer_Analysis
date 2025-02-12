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
# This script is for loading packages and data to start the analysis.
#
#
# ---------------------------------------------------------------------------- #
# ----Load packages ------------------------------------------------------------
library(STExplorer)

## Data manipulation/ plotting
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(cols4all)
library(tidyquant)
library(ggdist)
library(ggthemes)
library(rstatix)

## Analysis
library(GSVA)
library(harmony)
library(Seurat)
library(scMiko)

# ---------------------------------------------------------------------------- #
# ----Prepare to load data -----------------------------------------------------
## Prepare vectors with the paths to the data folders
projectFolder <- "./../STExplorer_Analysis/Lung"
dataFolder <- paste0(projectFolder, "/data/hs_visium_spaceranger_output")
sampleDir <- list.dirs(dataFolder, recursive = FALSE)
sampleNames <- list.files(dataFolder)
sampleNames <- sampleNames[sampleNames != "README.txt"]
names(sampleDir) <- sampleNames


# ---------------------------------------------------------------------------- #
# ----Load Metadata ------------------------------------------------------------
## Load the metadata
metaDataFile <- paste0(projectFolder, "/data/hs_visium_metadata.tsv")
metadata <- read.table(metaDataFile, header = TRUE) %>%
  mutate(fibrotic_extent_score_by_pathologist_0.3 = as.character(fibrotic_extent_score_by_pathologist_0.3))
rownames(metadata) <- metadata$sample_id

## Load cell density-related information as generated by original authors
cellAbundanceFolder <- paste0(projectFolder, "/data/cell2location_habermann2020/")
cellGroupsFile <- paste0(projectFolder, "/data/hs_misc/habermann_cell_type_groups.csv")
cellGroups <- read.table(cellGroupsFile, header = TRUE, sep = ";")

## Load spot annotations
annotationFile <- paste0(projectFolder, "/data/hs_misc/hs_visium_merged_histo_annotations.csv")
annotation <- read.table(annotationFile, header = TRUE, sep = ",") %>%
  rename(Barcode = barcode, sample_id = sample_slide_id)

# ---------------------------------------------------------------------------- #
# ----Select samples -----------------------------------------------------------
## Select samples for which there are annotations provided by the original authors
selection_samples <- names(sampleDir) %in% unique(annotation$sample_id)
sampleDir_selected <- sampleDir[selection_samples]
sampleNames_selected <- sampleNames[selection_samples]

## Update metadata to keep only the selected samples
metadata <- metadata[metadata$sample_id %in% sampleNames_selected,]


# ---------------------------------------------------------------------------- #
# ----Load data ----------------------------------------------------------------
## Create the MSFE object
msfe <- MetaSpatialFeatureExperiment()

## Import samples
for (i in seq_along(sampleNames_selected)) {
  id <- sampleNames_selected[i]
  message("Adding sample: ", sampleNames_selected[i])
  msfe <- addSFE(msfe,
                 read10xVisiumSFE(samples = sampleDir_selected[id],
                                  sample_id = sampleNames_selected[i],
                                  type = "HDF5",
                                  data = "filtered",
                                  images = "lowres",
                                  style = "W",
                                  zero.policy = TRUE))
}

## Import annotations
gTruth_list <- list()
for (id in sampleNames_selected) {
  gTruth_list[[id]] <- annotation[annotation$sample_id == id, c("Barcode", "sample_id", "annotation")]

}

## Import cell type abundances
cellAbundance <- vector(mode = "list", length = length(sampleNames_selected))
names(cellAbundance) <- sampleNames_selected
for (s in sampleNames_selected) {
  file <- paste0(cellAbundanceFolder, s, "_spot_cell_abundances_5pc.csv")
  df <- read.csv(file) %>%
    rename(Barcode = spot_id) %>%
    mutate(Barcode = gsub(paste0(s, "_"), "", Barcode),
           sample_id = s, .after = Barcode)
  cellAbundance[[s]] <- df
}


# ---------------------------------------------------------------------------- #
# ----Select genes -------------------------------------------------------------
## Get annotations
biomart <- createBiomart()
annotation_gene <- data.frame(gene_name = rowData(msfe@sfe_data[[1]])[["symbol"]])
rownames(annotation_gene) <- rownames(msfe@sfe_data[[1]])
annotation_gene <- annotateDataFrame(annotation_gene,
                                     biomart = biomart,
                                     add = c("biotype", "chromosome", "description"))

## Senescence Markers
senescenceMarkers <- c("CDKN2A", "CDKN1A", "GLB1", "LMNB1")
senescenceMarkers_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% senescenceMarkers],
  annotation_gene$gene_name[annotation_gene$gene_name %in% senescenceMarkers]
)

## Genes Related to senescence and fibrosis
otherGenes <- c("IL6", "IL6R", "IL6ST", "GDF15", "CCL2", "TGFB1", "FAP", "VIM",
                "ACTA2", "COL1A1", "COL1A2", "FN1", "CDH1", "CDH5")
otherGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% otherGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% otherGenes]
)

## SenMayo senescence markers group
senMayoGenes <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3",
                  "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16",
                  "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5",
                  "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB",
                  "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12",
                  "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1",
                  "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2",
                  "FGF7", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1",
                  "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4",
                  "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18",
                  "IL1A", "IL1B", "IL2", "IL32", "IL7", "INHA",
                  "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF",
                  "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3",
                  "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF",
                  "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES",
                  "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4",
                  "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF",
                  "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2",
                  "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2", "CCL2", "IL6",
                  "IL6ST", "GDF15")
senMayoGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% senMayoGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% senMayoGenes]
)

## Epithelial to mesenchymal transition gene set from MSigDB HALLMARK GENESET
emtGenes <- c("ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1",
              "BDNF", "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG",
              "CD44", "CD59", "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1",
              "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
              "COL5A1", "COL5A2", "COL5A3", "COL6A2", "COL6A3", "COL7A1",
              "COL8A2", "COMP", "COPA", "CRLF1", "CCN2", "CTHRC1", "CXCL1",
              "CXCL12", "CXCL6", "CCN1", "DAB2", "DCN", "DKK1", "DPYSL3",
              "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3", "ENO2",
              "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2",
              "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1",
              "FSTL3", "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM",
              "GJA1", "GLIPR1", "COLGALT1", "GPC1", "GPX7", "GREM1", "HTRA1",
              "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6",
              "CXCL8", "INHBA", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITGB3",
              "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2",
              "P3H1", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15",
              "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5",
              "MGP", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9",
              "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "PCOLCE",
              "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1", "PLOD2",
              "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2",
              "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1", "SCG2",
              "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1",
              "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3",
              "SNAI2", "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2",
              "TGFB1", "TGFBI", "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1",
              "TIMP1", "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A",
              "TPM1", "TPM2", "TPM4", "VCAM1", "VCAN", "VEGFA", "VEGFC", "VIM",
              "WIPF1", "WNT5A")
emtGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% emtGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% emtGenes]
)


## GO Extracellular Matrix gene set from MSigDB GO::CC GENESET
ecmGenes <- c("A1BG", "A2M", "ABI3BP", "ACAN", "ACHE", "ADAM11", "ADAM19", "ADAMDEC1", "ADAMTS1",
              "ADAMTS10", "ADAMTS12", "ADAMTS13", "ADAMTS14", "ADAMTS15", "ADAMTS16", "ADAMTS17",
              "ADAMTS18", "ADAMTS19", "ADAMTS2", "ADAMTS20", "ADAMTS3", "ADAMTS4", "ADAMTS5", "ADAMTS6",
              "ADAMTS7", "ADAMTS8", "ADAMTS9", "ADAMTSL1", "ADAMTSL2", "ADAMTSL3", "ADAMTSL4", "ADAMTSL5",
              "ADIPOQ", "AEBP1", "AGRN", "AGT", "AHSG", "ALPL", "AMBP", "AMELX", "AMELY", "AMTN", "ANG",
              "ANGPT1", "ANGPT2", "ANGPT4", "ANGPTL1", "ANGPTL2", "ANGPTL3", "ANGPTL4", "ANGPTL5",
              "ANGPTL6", "ANGPTL7", "ANOS1", "ANXA1", "ANXA11", "ANXA2", "ANXA2P2", "ANXA4", "ANXA5",
              "ANXA6", "ANXA7", "ANXA8", "APCS", "APLP1", "APOA1", "APOA4", "APOC3", "APOE", "APOH",
              "ASPN", "ATRN", "ATRNL1", "AZGP1", "BCAM", "BCAN", "BGN", "BMP7", "BMPER", "C17orf58",
              "C1QA", "C1QB", "C1QC", "CALR", "CASK", "CBLN1", "CBLN4", "CCBE1", "CCDC80", "CCN1", "CCN2",
              "CCN3", "CCN4", "CCN5", "CCN6", "CD151", "CD180", "CD248", "CDH13", "CDH2", "CDON", "CFP",
              "CHAD", "CHADL", "CHI3L1", "CILP", "CLC", "CLEC14A", "CLEC3B", "CLU", "CMA1", "COCH",
              "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1",
              "COL17A1", "COL18A1", "COL19A1", "COL1A1", "COL1A2", "COL20A1", "COL21A1", "COL22A1",
              "COL23A1", "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1", "COL2A1", "COL3A1",
              "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "COL5A1", "COL5A2", "COL5A3",
              "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6", "COL7A1", "COL8A1", "COL8A2", "COL9A1",
              "COL9A2", "COL9A3", "COLEC12", "COLQ", "COMP", "CPA3", "CPN2", "CRELD1", "CRISP3",
              "CRISPLD2", "CSPG4", "CST3", "CSTB", "CTHRC1", "CTSB", "CTSC", "CTSD", "CTSF", "CTSG",
              "CTSH", "CTSL", "CTSS", "CTSZ", "CXCL12", "DAG1", "DCN", "DEFA1", "DEFA1B", "DGCR6",
              "DGCR6", "DLG1", "DMBT1", "DMP1", "DPT", "DSPP", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP1",
              "EFEMP2", "EFNA5", "EGFL6", "EGFL7", "EGFLAM", "ELANE", "ELFN1", "ELFN2", "ELN", "EMID1",
              "EMILIN1", "EMILIN2", "EMILIN3", "ENAM", "ENTPD2", "EPYC", "ERBIN", "EYS", "F12", "F13A1",
              "F2", "F3", "F7", "F9", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FBN3", "FCGBP", "FCN1",
              "FCN2", "FCN3", "FGA", "FGB", "FGF1", "FGF10", "FGF9", "FGFBP3", "FGFR2", "FGG", "FGL1",
              "FGL2", "FIBCD1", "FLG", "FLRT1", "FLRT2", "FLRT3", "FMOD", "FN1", "FRAS1", "FREM1",
              "FREM2", "FREM3", "GDF10", "GDF15", "GFOD2", "GH1", "GLG1", "GP1BA", "GPC1", "GPC2", "GPC3",
              "GPC4", "GPC5", "GPC6", "GPLD1", "GREM1", "HAPLN1", "HAPLN2", "HAPLN3", "HAPLN4", "HDGF",
              "HMCN1", "HMCN2", "HNRNPM", "HPSE", "HPSE2", "HPX", "HRG", "HRNR", "HSD17B12", "HSP90B1",
              "HSPG2", "HTRA1", "ICAM1", "IFNA2", "IGFALS", "IGFBP7", "IHH", "IL7", "IMPG1", "IMPG2",
              "INHBE", "ITGA6", "ITIH1", "ITIH2", "ITIH4", "ITIH5", "KAZALD1", "KERA", "KNG1", "KRT1",
              "L1CAM", "LAD1", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3",
              "LAMB4", "LAMC1", "LAMC2", "LAMC3", "LEFTY2", "LGALS1", "LGALS3", "LGALS3BP", "LGALS4",
              "LGALS9", "LINGO1", "LINGO2", "LINGO3", "LINGO4", "LMAN1", "LMAN1L", "LOX", "LOXL1",
              "LOXL2", "LRIG1", "LRIG2", "LRIG3", "LRRC15", "LRRC17", "LRRC24", "LRRC32", "LRRC3B",
              "LRRC3C", "LRRN1", "LRRN2", "LRRN3", "LRRTM1", "LRRTM3", "LRRTM4", "LTBP1", "LTBP2",
              "LTBP3", "LTBP4", "LUM", "MARCOL", "MATN1", "MATN2", "MATN3", "MATN4", "MBL2", "MDK",
              "MEGF9", "MEPE", "MFAP1", "MFAP2", "MFAP4", "MFAP5", "MFGE8", "MGP", "MMP1", "MMP10",
              "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP2", "MMP20",
              "MMP21", "MMP23B", "MMP24", "MMP25", "MMP26", "MMP27", "MMP28", "MMP3", "MMP7", "MMP8",
              "MMP9", "MMRN1", "MMRN2", "MST1", "MUC17", "MUC2", "MUC4", "MUC5AC", "MUC5B", "MUC6",
              "MXRA5", "MXRA7", "MYOC", "NAV2", "NCAM1", "NCAN", "NDNF", "NDP", "NID1", "NID2", "NPNT",
              "NPPA", "NTN1", "NTN3", "NTN4", "NTN5", "NTNG1", "NTNG2", "NYX", "OC90", "OGN", "OLFML2A",
              "OMD", "OPTC", "ORM1", "ORM2", "OTOG", "OTOGL", "OTOL1", "P3H1", "P3H2", "PAPLN", "PCOLCE",
              "PCSK6", "PDGFB", "PF4", "PHOSPHO1", "PI3", "PKM", "PLG", "PLOD3", "PLSCR1", "PODN",
              "PODNL1", "POMZP3", "POSTN", "PRELP", "PRG2", "PRG3", "PRG4", "PRSS1", "PRSS2", "PRTN3",
              "PSAP", "PTN", "PTPRZ1", "PTX3", "PXDN", "PZP", "RARRES2", "RBP3", "RELL2", "RELN", "RTBDN",
              "RTN4RL1", "RTN4RL2", "S100A10", "S100A4", "S100A6", "S100A7", "S100A8", "S100A9", "SBSPON",
              "SCARA3", "SDC2", "SDC3", "SEMA3B", "SEMA7A", "SERAC1", "SERPINA1", "SERPINA3", "SERPINA5",
              "SERPINB1", "SERPINB12", "SERPINB6", "SERPINB8", "SERPINB9", "SERPINC1", "SERPINE1",
              "SERPINE2", "SERPINF1", "SERPINF2", "SERPING1", "SERPINH1", "SFRP1", "SFRP2", "SHH", "SLPI",
              "SMOC1", "SMOC2", "SNORC", "SOD3", "SOST", "SPARC", "SPARCL1", "SPOCK2", "SPON1", "SPON2",
              "SPP2", "SRPX", "SRPX2", "SSC5D", "SSPOP", "SULF1", "TECTA", "TECTB", "TFIP11", "TFPI2",
              "TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBI", "TGFBR3", "TGM2", "TGM4", "THBS1", "THBS2",
              "THBS3", "THBS4", "THSD4", "TIMP1", "TIMP2", "TIMP3", "TIMP4", "TINAG", "TINAGL1", "TLR3",
              "TMEFF1", "TMEFF2", "TNC", "TNFRSF11B", "TNN", "TNR", "TNXB", "TPSAB1", "TPSB2", "TRIL",
              "UCMA", "USH2A", "VASN", "VCAN", "VEGFA", "VIT", "VTN", "VWA1", "VWA2", "VWC2", "VWF",
              "WNT11", "WNT2", "WNT2B", "WNT3", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT8A", "ZAN",
              "ZG16", "ZP1", "ZP2", "ZP3")

ecmGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% ecmGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% ecmGenes]
) # Removed "ANXA2P2" "SSPOP" because are marked as Pseudogenes


## TGF-beta gene set as found in Ma et. al., Nat Commun, 2024, DOI: https://doi.org/10.1038/s41467-023-44645-6
tgfbGenes <- c("AC003092.1", "AC022509.2", "AC083837.2", "ADAMTS14", "ADAMTS4", "ADAMTS9", "AMPD3",
               "ANGPTL4", "ANO1", "APOL1", "ATF3", "BATF2", "BCL2A1", "BDKRB1", "BHLHE40-AS1", "BIRC3",
               "C15orf48", "C3", "C5orf56", "CACNA1G", "CARD6", "CAVIN4", "CCL11", "CCL2", "CCL5", "CCL7",
               "CD274", "CD38", "CD70", "CD82", "CD83", "CFAP45", "CHI3L2", "CLDN1", "CLDN4", "CMPK2",
               "COL22A1", "CSF1", "CSF2", "CTSS", "CX3CL1", "CXCL1", "CXCL10", "CXCL2", "CXCL3", "CXCL6",
               "CXCL8", "CYP7B1", "DCLK3", "DDX58", "DDX60", "DDX60L", "EDARADD", "EGR2", "ENC1", "EREG",
               "ESM1", "FAM167A", "G0S2", "GBP1P1", "GFPT2", "GPR68", "GRIP2", "HAS1", "HBEGF", "HCK",
               "HELZ2", "HLA-F", "ICAM1", "IER3", "IFIH1", "IFIT2", "IFIT3", "IL11", "IL15", "IL15RA",
               "IL18BP", "IL1B", "IL32", "IL33", "IL6", "IRAK2", "IRF7", "ITGB3", "JUNB", "KCTD11", "KDR",
               "L1TD1", "LGALS9", "LIF", "LINC01539", "LIPG", "LRIG1", "MCC", "MED12L", "MEOX1", "MFSD2A",
               "MMP12", "MSC", "MX1", "MX2", "MYO10", "NFKB2", "NFKBIA", "NFKBIE", "NFKBIZ", "NGFR",
               "NIPAL4", "NKD1", "NKX3-1", "NPAS2", "NPTX1", "NR4A1", "NR4A2", "OAS1", "OASL", "OGFR",
               "P2RY6", "PARP12", "PARP14", "PDPN", "PHACTR1", "PKP1", "POU2F2", "PREX1", "PTGS2", "PTHLH",
               "RASD1", "RASL11B", "RCSD1", "RELB", "RFLNA", "RGS16", "RIPK2", "RIPOR3", "RND1", "RNF19B",
               "RRAD", "RSPO3", "SAA1", "SERPINA1", "SERPINB2", "SLAMF8", "SLC12A7", "SLC22A3", "SLC25A28",
               "SLC2A6", "SLC39A14", "SMOX", "SOD2", "SP6", "STAT5A", "STRA6", "STX11", "SYT12", "TAP1",
               "TCIM", "TFPI2", "TGFA", "TLR3", "TMEM51", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6",
               "TNFAIP8", "TNIP1", "TRAF1", "TRIM36", "TRPA1", "VCAM1", "XAF1", "ZC3H12A", "ZFPM2", "ZNF710")

tgfbGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% tgfbGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% tgfbGenes]
) # Removed "AC083837.2" "GBP1P1" because are marked as Pseudogenes

## Senescence-related transcription factors
transcrFact <- c("SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2")
transcrFact_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% transcrFact],
  annotation_gene$gene_name[annotation_gene$gene_name %in% transcrFact]
)

## Concatenate all but the ligand-receptor vector
selectedGenes <- c(senescenceMarkers, otherGenes, senMayoGenes, transcrFact)
selectedGenes_ENSG <- c(senescenceMarkers_ENSG, otherGenes_ENSG,
                        senMayoGenes_ENSG, transcrFact_ENSG)

# ---------------------------------------------------------------------------- #
# ----Select cell types --------------------------------------------------------
cellTypes <- colnames(cellAbundance[[1]][3:ncol(cellAbundance[[1]])]) %>%
  gsub("q05cell_abundance_w_sf_", "", .)
names(cellTypes) <- cellGroups$cell_name_special

selectedCellTypes <- c("B.Cells", "Basal", "Fibroblasts", "Endothelial.Cells",
                       "Macrophages", "T.Cells", "AT1", "AT2", "MUC5B.",
                       "KRT5..KRT17.", "Transitional.AT2", "Myofibroblasts",
                       "Ciliated", "HAS1.High.Fibroblasts",
                       "MUC5AC..High", "SCGB3A2..SCGB1A1.", "SCGB3A2.",
                       "Proliferating.Epithelial.Cells")
names(selectedCellTypes) <- c("B cells", "Basal", "Fibroblasts",
                              "Endothelial cells", "Macrophages", "T cells",
                              "AT1", "AT2", "MUC5B+", "KRT5-/KRT17+",
                              "Transitional AT2", "Myofibroblasts",
                              "Ciliated", "HAS1 High Fibroblasts",
                              "MUC5AC-high", "SCGB3A2+SCGB1A1+", "SCGB3A2+",
                              "Proliferating Epithelial Cells")

# ---------------------------------------------------------------------------- #
# ----Custom functions to use --------------------------------------------------
# Define a function to print caught errors
error = function(e){cat("ERROR: ", conditionMessage(e), "\n")}

# Define a function to calculate the proportion of spots expressing each gene
calculate_presence_proportion <- function(sfe, genes, assay, threshold = 0) {
  # Subset the dataset to include only the genes of interest
  gene_expr <- as.matrix(assay(sfe, assay))[genes, ]

  # Determine the spots where each gene is expressed above the threshold
  expressed <- gene_expr > threshold

  # Calculate the proportion of spots where each gene is expressed
  presence_proportion <- rowMeans(expressed)

  return(presence_proportion)
}

# Define a function to plot module scores
plotModuleScores <- function(m_sfe,
                             modules,
                             sample_id = NULL,
                             type = c("spot", "hex", "cntd"),
                             res = c("lowres", "hires", "fullres", "none"),
                             fill_args = list(),
                             alpha = 0.3,
                             title_type = "gene_name",
                             lab_legend = NULL,
                             facet_by = NULL,
                             ...) {

  ## Check SFE or MSFE?
  sfe <- .int_sfeORmsfe(m_sfe = m_sfe, sample_id = sample_id)

  ## Check arguments
  if (missing(type)) {
    type <- "hex"
  }
  if (missing(res)) {
    res <- "none"
  }

  ## Set geometry to use
  if (type == "hex") {
    stopifnot("spotHex" %in% names(colGeometries(sfe)))
    type <- "spotHex"
  } else if (type == "spot") {
    type <- "spotPoly"
  } else if (type == "cntd") {
    type <- "spotCntd"
  }

  ## Fetch image if needed and fix alpha
  if (res != "none") {
    ## Fetch image data and transform to raster
    image <- .int_getImgDt(sfe = sfe, sample_id = sample_id, image_id = res)
    ## Get capture area limits
    limits_list <- .int_getImgLims(sfe = sfe, res = res)
    ## Get max col value to accurately plot the image
    max_list <- lapply(image,
                       function(l){
                         max <- l@ptr$range_max[1]
                         if (is.na(max)) {
                           max <- max(l[,,1])
                         }
                         if (max <= 1) {
                           max_col_value <- 1
                         } else if (max <= 255) {
                           max_col_value <- 255
                         } else if (max <= 65536) {
                           max_col_value <- 65536
                         }
                       }
    )
  } else {
    alpha <- 1
  }

  ## Create legend title
  if (is.null(lab_legend)) {
    lab_legend <- "Module Scores"
  }

  ## Set fill arguments if not provided
  if (isEmpty(fill_args)) {
    fill_args <- list(option = "viridis",
                      na.value = "grey")
  }

  ## Fetch data
  data <- colData(sfe)[c("sample_id", modules)] %>% as.data.frame()
  geoms <- sfe@int_colData@listData[["colGeometries"]]@listData[[type]]
  rownames(geoms) <- colnames(sfe)
  data <- merge(data, geoms, by = "row.names") %>%
    pivot_longer(cols = -c("geometry", "Row.names", "sample_id"), names_to = "mScore", values_to = "score")

  ## Plot
  p <- ggplot()

  if (res != "none") {
    p <- p +
      tidyterra::geom_spatraster_rgb(data = image[[1]],
                                     max_col_value = max_list[[1]]) +
      ggplot2::lims(x = limits_list[[1]],
                    y = limits_list[[2]])
  }

  if (type != "cntd") {
    p <- p +
      ggplot2::geom_sf(data = data,
                       aes(geometry = geometry,
                           fill = score),
                       alpha = alpha,
                       colour = "transparent") +
      do.call(scale_fill_viridis_c, c(list(), fill_args)) +
      ggplot2::labs(fill = lab_legend)
  } else {
    p <- p +
      ggplot2::geom_sf(data = data,
                       aes(geometry = geometry,
                           colour = score),
                       alpha = alpha,
                       size = 0.5) +
      do.call(scale_colour_viridis_c, c(list(), fill_args)) +
      ggplot2::labs(colour = lab_legend)
  }

  if (length(unique(data$mScore)) > 1 & is.null(facet_by)) {
    p <- p +
      ggplot2::facet_wrap(~mScore, scales = "fixed", ncol = 1)
  } else if (!is.null(facet_by)) {
    p <- p +
      ggplot2::facet_wrap(facet_by, scales = "fixed", ncol = 4)
  } else {
    p <- p +
      labs(title = unique(data$mScore))
  }

  p <- p +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p
}


