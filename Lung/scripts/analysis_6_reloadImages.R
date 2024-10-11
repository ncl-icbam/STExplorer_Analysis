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
# This script is for re-loading images into the SFE objects. When saving an SFE
# (or MSFE) object and loading it again the images need to be re-loaded again.
# The below for-loop is performing this for all selected samples.
#
#
# ---------------------------------------------------------------------------- #
# ----Re-load images -----------------------------------------------------------
for (id in sampleNames_selected) {
  sfe <- getSFE(msfe, id)
  sfe@int_metadata$imgData <- NULL
  ## Get scale factors
  scaleF <- jsonlite::fromJSON(txt = file.path(sampleDir_selected[[id]],
                                               "outs/spatial",
                                               "scalefactors_json.json"))
  sfe <- SpatialFeatureExperiment::addImg(sfe,
                                          file.path(sampleDir_selected[[id]],
                                                    "outs/spatial/tissue_lowres_image.png"),
                                          sample_id = id,
                                          image_id = "lowres",
                                          scale_fct = scaleF[["tissue_lowres_scalef"]])
  sfe <- SpatialFeatureExperiment::mirrorImg(sfe,
                                             sample_id = id,
                                             image_id = "lowres")
  msfe <- addSFE(msfe, sfe, id)
  ## Housekeeping
  rm(sfe)
}
