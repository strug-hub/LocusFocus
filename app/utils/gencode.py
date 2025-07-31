"""
Utility functions for gene name and gencode lookups using local files
"""

import os
from typing import List

import pandas as pd
from app.utils import parse_region_text

# /app/utils/gencode.py -> /data/
DATA_FOLDER = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../data")
)
collapsed_genes_df_hg19 = pd.read_csv(
    os.path.join(DATA_FOLDER, "collapsed_gencode_v19_hg19.gz"),
    compression="gzip",
    sep="\t",
    encoding="utf-8",
)
collapsed_genes_df_hg38 = pd.read_csv(
    os.path.join(DATA_FOLDER, "collapsed_gencode_v26_hg38.gz"),
    compression="gzip",
    sep="\t",
    encoding="utf-8",
)


def get_genes_by_location(build: str, chrom: int, startbp: int, endbp: int, gencode=False) -> List[str]:
    """
    Given a build and a region, return a list of genes that overlap that region.
    """
    collapsed_genes_df = collapsed_genes_df_hg19
    if build.lower() == "hg38":
        collapsed_genes_df = collapsed_genes_df_hg38
    regiontext = str(chrom) + ":" + str(startbp) + "-" + str(endbp)
    chrom, startbp, endbp = parse_region_text(regiontext, build)
    genes_to_draw = collapsed_genes_df.loc[
        (collapsed_genes_df["chrom"] == ("chr" + str(chrom).replace("23", "X")))
        & (
            (
                (collapsed_genes_df["txStart"] >= startbp)
                & (collapsed_genes_df["txStart"] <= endbp)
            )
            | (
                (collapsed_genes_df["txEnd"] >= startbp)
                & (collapsed_genes_df["txEnd"] <= endbp)
            )
            | (
                (collapsed_genes_df["txStart"] <= startbp)
                & (collapsed_genes_df["txEnd"] >= endbp)
            )
        )
    ]
    if gencode:
        return list(genes_to_draw["ENSG_name"])
    return list(genes_to_draw["name"])
