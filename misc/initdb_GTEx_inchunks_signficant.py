#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import significant gtex pairs only for lighter DB for local development
Files from https://www.gtexportal.org/home/downloads/adult-gtex/qtl
"""

from datetime import datetime
import gzip
import itertools
import os
import tempfile

import numpy as np
import pandas as pd
import pymongo
from pymongo import MongoClient

mongo_user = os.getenv("MONGO_USERNAME")
mongo_pass = os.getenv("MONGO_PASSWORD")
mongo_db = os.getenv("MONGO_DATABASE")
# assuming docker here
mongo_host = "mongo"


def stream_groupby_csv(path, key, agg, collection, chunk_size=1e6, pool=None, **kwargs):

    # Make sure path is a list
    if not isinstance(path, list):
        path = [path]

    # Chain the chunks
    kwargs["chunksize"] = chunk_size
    chunks = itertools.chain(*[pd.read_csv(p, **kwargs) for p in path])

    results = []
    orphans = pd.DataFrame()

    for chunk in chunks:

        # Add the previous orphans to the chunk
        chunk = pd.concat((orphans, chunk))

        # Determine which rows are orphans
        last_val = chunk[key].iloc[-1]
        is_orphan = chunk[key] == last_val

        # Put the new orphans aside
        chunk, orphans = chunk[~is_orphan], chunk[is_orphan]

    results.append(agg(chunk, collection))

    # If a pool is used then we have to wait for the results
    if pool:
        results = [r.get() for r in results]

    results.append(
        agg(orphans, collection)
    )  # ensure last chunk (gene) is pushed as well!

    return pd.concat(results)


def agg(chunk, collection):
    """lambdas can't be serialized so we need to use a function"""
    chunk.set_index("gene_id", inplace=True)
    return chunk.groupby("gene_id").apply(push_variant_dict, collection=collection)


def push_variant_dict(gene_df, collection):
    variant_id = [
        str(x).replace("chr", "").encode("utf-8") for x in list(gene_df["variant_id"])
    ]
    pval = list(gene_df["pval_nominal"])
    beta = list(gene_df["slope"])
    se = list(gene_df["slope_se"])
    ma_samples = list(gene_df["ma_samples"])
    ma_count = list(gene_df["ma_count"])
    sample_maf = list(gene_df["maf"])
    geneid = gene_df.reset_index()["gene_id"][0]
    variants_list = []
    for row in np.arange(len(variant_id)):
        variants_list.append(
            {
                "variant_id": variant_id[row].decode("utf-8"),
                "pval": float(pval[row]),
                "beta": float(beta[row]),
                "se": float(se[row]),
                "ma_samples": float(ma_samples[row]),
                "ma_count": float(ma_count[row]),
                "sample_maf": float(sample_maf[row]),
            }
        )
    gene_dict = {"gene_id": geneid, "eqtl_variants": variants_list}
    collection.insert_one(gene_dict)


##########################################
# MAIN
##########################################


def main(version: int):

    conn = (
        f"mongodb://{mongo_user}:{mongo_pass}@mongo:27017/{mongo_db}?authSource=admin"
    )
    client = MongoClient(conn)
    db_name = f"GTEx_V{version}"
    db = client[db_name]

    # we're going to start by whitelisting just a few tissues:

    tissues = ["Brain_Amygdala", "Liver"]
    files_list = [f"{t}.v{version}.signif_variant_gene_pairs.txt.gz" for t in tissues]

    for file in files_list:
        tissue_name = file.split(".")[0]
        file_path = os.path.join(os.path.dirname(__file__), "data", file)
        if tissue_name not in db.list_collection_names():
            collection = db[tissue_name]
            if file.endswith("gz") and os.path.isfile(file_path):
                print("Decompressing " + file)
                with open(file_path, "rb") as f:
                    decompressed = gzip.decompress(f.read())
                    _, file_path = tempfile.mkstemp(suffix=".txt")
                    with open(file_path, "w") as tmp:
                        tmp.write(decompressed.decode())
                # decompress(file)
            print("Parsing file " + file + " and creating tissue collection")
            stream_groupby_csv(
                path=[file_path],
                collection=collection,
                key="gene_id",
                agg=agg,
                chunk_size=1e6,
                sep="\t",
                usecols=[
                    "gene_id",
                    "variant_id",
                    "pval_nominal",
                    "slope",
                    "slope_se",
                    "ma_samples",
                    "ma_count",
                    "maf",
                ],
            )
            print(tissue_name + " collection created")
            print(datetime.now().strftime("%c"))
            print("Now indexing by gene_id")
            print(datetime.now().strftime("%c"))
            collection.create_index("gene_id")
            print("Indexing done")
            print(datetime.now().strftime("%c"))
            print("Done with tissue " + tissue_name)
            print(datetime.now().strftime("%c"))
        # break

    # Next, create the variant lookup table

    lookup_table_filename = (
        "GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz"
        if version == 7
        else "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
    )

    collection = db["variant_table"]
    tbl_chunk = pd.read_csv(
        os.path.join(
            os.path.dirname(__file__),
            "data",
            lookup_table_filename,
        ),
        sep="\t",
        chunksize=1e5,
        encoding="utf-8",
    )
    print("Pushing variant information into variant_table collection by chunks")
    for tbl in tbl_chunk:
        tbl["chr"] = [
            int(str(x).replace("chr", "").replace("X", "23")) for x in list(tbl["chr"])
        ]
        tbl["variant_id"] = [x.replace("chr", "") for x in list(tbl["variant_id"])]
        tbl_dict = tbl.to_dict(orient="records")
        collection.insert_many(tbl_dict)
    print("Variant collection created")
    print(datetime.now().strftime("%c"))

    print("Now indexing by variant_id")
    collection.create_index("variant_id")
    print("Indexing by variant_id done")
    print(datetime.now().strftime("%c"))
    print("Now indexing by chr and pos (ascending) order")
    collection.create_index(
        [("chr", pymongo.ASCENDING), ("variant_pos", pymongo.ASCENDING)]
    )
    print("Indexing by chr and pos done")
    print(datetime.now().strftime("%c"))


if __name__ == "__main__":
    [main(v) for v in [7, 8]]
