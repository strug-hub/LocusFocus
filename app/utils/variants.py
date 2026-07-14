"""
Variant lookup functions.
"""

import pandas as pd
from flask import current_app

from app.utils.helpers import validate_chromosome


def get_variants_by_region(
    start: int, end: int, chrom: str, gtex_version: str
) -> pd.DataFrame:
    """Return variants in the region [start, end] on chrom.

    :param start: Start position (inclusive).
    :param end: End position (inclusive).
    :param chrom: Chromosome. May be prefixed with "chr" and may use X for 23.
    :param gtex_version: Either "V8" or "V10".
    :return: DataFrame with columns `variant_id`, `rs_id`, `chr`, `pos`,
             `ref`, `alt`.  Empty DataFrame when no variants are found.
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

    gtex_db = current_app.extensions["gtex_db"]
    return gtex_db.get_variants_by_region(start, end, chrom, gtex_version)
