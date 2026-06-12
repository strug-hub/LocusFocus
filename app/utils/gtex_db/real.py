"""
MongoDB-backed implementation of GTExDatabase.

Wraps the GTEx_V8 and GTEx_V10 MongoDB databases.  Each tissue is stored as a
collection; each document holds all eQTL variants for one gene.  A separate
`variant_table` collection maps variant IDs to rs IDs and positional info.
"""

from typing import List

import pandas as pd
from pymongo import MongoClient

from app.utils.gtex_db.base import GTExDatabase


class RealGTExDatabase(GTExDatabase):
    """GTExDatabase backed by a live MongoDB instance."""

    def __init__(self, client: MongoClient) -> None:
        self._client = client

    # ------------------------------------------------------------------

    def list_tissues(self, version: str) -> List[str]:
        if version.upper() not in ("V8", "V10"):
            raise ValueError(f"Invalid GTEx version: {version}")
        db = self._client[f"GTEx_{version.upper()}"]
        return sorted(n for n in db.list_collection_names() if n != "variant_table")

    def get_eqtl_data(
        self, version: str, tissue: str, ensg_id_prefix: str
    ) -> pd.DataFrame:
        version = version.upper()
        if version not in ("V8", "V10"):
            raise ValueError(f"Invalid GTEx version: {version}")
        db = self._client[f"GTEx_{version}"]
        collection = db[tissue]

        results = list(
            collection.find({"gene_id": {"$regex": f"^{ensg_id_prefix}.*"}})
        )
        if not results:
            return pd.DataFrame()

        eqtl_variants = results[0].get("eqtl_variants", [])
        if not eqtl_variants:
            return pd.DataFrame()

        eqtl_df = pd.DataFrame(eqtl_variants)

        chrom = eqtl_df["variant_id"].iloc[0].split("_")[0].replace("X", "23")
        positions = eqtl_df["variant_id"].str.split("_").str[1].astype(int)
        variants_df = self.get_variants_by_region(
            int(positions.min()), int(positions.max()), chrom, version
        )

        return pd.merge(eqtl_df, variants_df, on="variant_id")

    def get_variants_by_region(
        self, start: int, end: int, chrom: str, version: str
    ) -> pd.DataFrame:
        version = version.upper()
        if version not in ("V8", "V10"):
            raise ValueError(f"Invalid GTEx version: {version}")
        if version == "V8":
            db = self._client.GTEx_V8
            rsid_col = "rs_id_dbSNP151_GRCh38p7"
            query = {"chr": int(chrom), "variant_pos": {"$gte": start, "$lte": end}}
        else:  # V10
            db = self._client.GTEx_V10
            rsid_col = "rs_id_dbSNP155_GRCh38p13"
            query = {"chr": int(chrom), "pos": {"$gte": start, "$lte": end}}

        docs = list(db["variant_table"].find(query))
        if not docs:
            return pd.DataFrame()

        df = pd.DataFrame(docs).drop(columns=["_id"]).rename(
            columns={rsid_col: "rs_id"}
        )
        # Normalise V8's "variant_pos" to "pos" so callers get a consistent column name.
        if "variant_pos" in df.columns:
            df = df.rename(columns={"variant_pos": "pos"})
        return df
