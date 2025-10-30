"""
Variant lookup functions
"""
from flask import current_app as app
from pymongo import MongoClient
import pandas as pd

from app import mongo
from app.utils.helpers import validate_chromosome
from app.utils.apis.gtex import get_variants


client = None  # type: ignore

with app.app_context():
    client: MongoClient = mongo.cx  # type: ignore


def get_variants_by_region(start: int, end: int, chrom: str, gtex_version: str) -> pd.DataFrame:
    """
    Given a region, return a dataframe of variants within that region, inclusive.

    :param start: The start position of the region
    :type start: int
    :param end: The end position of the region
    :type end: int
    :param chrom: The chromosome of the region. May be prefixed with "chr" and may use X and Y for 23.
    :type chrom: str
    :param gtex_version: The GTEx version to use. Either "V8" or "V10"
    :type gtex_version: str
    :return: A dataframe of variants within the region
    """
    if start > end:
        raise ValueError("Start must be less than or equal to end")
    if start < 0:
        raise ValueError("Start must be greater than or equal to 0")

    chrom = str(chrom).lower().replace("chr", "").replace("x", "23")
    if not validate_chromosome(chrom, prefix="", x_y_numeric=True):
        raise ValueError("Invalid chromosome format")

    if gtex_version.upper() not in ["V8", "V10"]:
        raise ValueError("gtex_version must be either 'V8' or 'V10'")
    
    if gtex_version.upper() == "V8":
        db = client.GTEx_V8
        collection = db["variant_table"]
        variants_query = collection.find(
            {
                "$and": [
                    {"chr": int(chrom)},
                    {"variant_pos": {"$gte": start, "$lte": end}},
                ]
            }
        )
        variants_list = list(variants_query)
        variants_df = pd.DataFrame(variants_list)
        if len(variants_df) == 0:
            return variants_df
        variants_df = variants_df.drop(["_id"], axis=1).rename(columns={"rs_id_dbSNP151_GRCh38p7": "rs_id"})
        return variants_df
    elif gtex_version.upper() == "V10":
        # We use the API, but convert the dataframe so that the format is the same as V8 db
        chrom = "chr" + chrom
        variant_response = get_variants(dataset_id="gtex_v10", start=start, end=end, chromosome=chrom)
        variants_df = pd.DataFrame([v.to_dict() for v in variant_response.data])
        if len(variants_df) == 0:
            return variants_df

        variants_df = variants_df.rename(columns={
                "snpId": "rs_id",
                "b37VariantId": "variant_id_b37",
                "pos": "variant_pos",
                "maf01": "num_alt_per_site",
                "variant_id": "variant_id",
                "chromosome": "chr",
            }) \
            .drop(columns=["snpIdUpper", "datasetId", "shorthand"]) \
            .astype({"num_alt_per_site": int})
        
        variants_df["chr"] = [int(x.value.replace("chr", "")) for x in variants_df["chr"]]
        variants_df["variant_id"] = variants_df["variant_id"].str.replace("chr", "")

        return variants_df

