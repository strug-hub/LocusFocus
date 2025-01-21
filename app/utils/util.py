from typing import List, Dict

from app.utils.apis.ensembl import fetch_variant_info


def format_variant(variant: Dict, suffix: str) -> str:
    ref = variant["alleles"][0]
    alt = variant["alleles"][1]
    chr = variant["seq_region_name"]
    pos = variant["start"]

    return f"{chr}_{pos}_{ref}_{alt}_{suffix}"


def format_variant_list(variants: List[Dict]) -> List[str]:
    return [format_variant(v) for v in variants]


def fetch_and_format_variants(build: str, chr: str, location: int, suffix: str):
    """Fetch variant info from ensembl API
    and return in string format
    :param build: The genome build
    :type build: str
    :param chr: The chromosome
    :type chr: str
    :param location: The location
    :type location: int
    :param suffix: The build identifier
    :type location: str
    :return: List of variant strings
    :rtype: List[str]
    """

    variants = fetch_variant_info(build, chr, location)
    return format_variant_list(variants)
