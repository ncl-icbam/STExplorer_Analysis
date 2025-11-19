#!/usr/bin/env python3

# Run:
# python3 ./scripts/FGWC_STMiner_2_build_SpatialData_and_h5ad.py --in ./export --out ./results --name H1_5

# FILE: FGWC_STMiner_2_build_SpatialData_and_h5ad.py

import argparse, os
import pandas as pd
import numpy as np
import scipy.io
import scipy.sparse as sp
from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel

def read_mtx(path):
    m = scipy.io.mmread(path)
    return m.tocsr() if sp.issparse(m) else sp.csr_matrix(m)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir",  default="export")
    ap.add_argument("--out", dest="outdir", default="results")
    ap.add_argument("--name", dest="name",  default="sample")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    counts_path    = os.path.join(args.indir, "counts.mtx.gz")
    logcounts_path = os.path.join(args.indir, "logcounts.mtx.gz")
    features_path  = os.path.join(args.indir, "features.tsv.gz")
    barcodes_path  = os.path.join(args.indir, "barcodes.tsv.gz")
    coords_path    = os.path.join(args.indir, "coords.csv")
    labels_path    = os.path.join(args.indir, "labels.csv")

    # matrices (transpose to get spots × genes)
    X_counts = read_mtx(counts_path).T.tocsr()
    X_log    = read_mtx(logcounts_path).T.tocsr()

    var = pd.read_csv(features_path, sep="\t", header=None)
    var.columns = ["gene_id", "gene_name"]
    var.index = var["gene_id"].astype(str)

    obs = pd.read_csv(barcodes_path, sep="\t", header=None)
    obs.columns = ["barcode"]
    obs.index = obs["barcode"].astype(str)

    # sanity checks after transpose
    assert X_counts.shape == X_log.shape
    assert X_counts.shape[0] == obs.shape[0], f"{X_counts.shape[0]} rows in X vs {obs.shape[0]} obs"
    assert X_counts.shape[1] == var.shape[0], f"{X_counts.shape[1]} cols in X vs {var.shape[0]} var"

    # coords
    coords = pd.read_csv(coords_path)
    coords.index = coords["spot"].astype(str)
    coords = coords.loc[obs.index, ["x","y"]].to_numpy()

    # labels (optional)
    if os.path.exists(labels_path):
        labels = pd.read_csv(labels_path)
        labels.index = labels["spot"].astype(str)
        obs["label"] = labels.loc[obs.index, "label"].astype(str).values

    # Build AnnData: raw counts in X; logcounts in a layer
    adata = AnnData(X=X_counts, obs=obs, var=var)
    adata.layers["logcounts"] = X_log.copy()
    adata.obsm["spatial"] = coords
    adata.obs["region"] = "spots"
    adata.obs["spot_id"] = adata.obs.index.astype(str)

    # SpatialData: points + table
    df_points = pd.DataFrame(
    {"x": adata.obsm["spatial"][:, 0],
     "y": adata.obsm["spatial"][:, 1],
     "spot_id": adata.obs["spot_id"].values},
    index=adata.obs.index
)

    points = PointsModel.parse(df_points)
    
    table = TableModel.parse(
        adata,
        region="spots",
        region_key="region",
        instance_key="spot_id"
    )
    sdata = SpatialData(points={"spots": points}, table=table)

    # Write outputs
    sdata_path = os.path.join(args.outdir, f"{args.name}.sdata.zarr")
    sdata.write(sdata_path)

    h5ad_path = os.path.join(args.outdir, f"{args.name}.h5ad")
    adata.write(h5ad_path)

    print(f"Done.\nSpatialData: {sdata_path}\nAnnData for STMiner: {h5ad_path}")

if __name__ == "__main__":
    main()
