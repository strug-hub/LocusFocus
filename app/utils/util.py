from typing import List, Dict

from app.utils.apis.ensembl import fetch_ensembl_variant_info


def format_variant_gtex_style(chr: str, pos: int, ref: str, alt: str, suffix: str):
    """Format variant in gtex style"""
    return f"{chr}_{pos}_{ref}_{alt}_{suffix}"


def format_ensembl_variant(variant: Dict, suffix: str) -> str:
    """If biallelic, will add commas"""
    ref = variant["alleles"][0]
    alt = ",".join(variant["alleles"][1:])
    chr = variant["seq_region_name"]
    pos = variant["start"]

    return format_variant_gtex_style(chr, pos, ref, alt, suffix)


def format_ensembl_variant_list(variants: List[Dict]) -> List[str]:
    return [format_ensembl_variant(v) for v in variants]


def fetch_and_format_ensembl_variants(
    build: str, chr: str, suffix: str, start: int, end: int = None
):
    """Fetch variant info from ensembl API
    and return in string format
    :param build: The genome build
    :type build: str
    :param chr: The chromosome
    :type chr: str
    :param suffix: The build identifier
    :type suffix: str
    :param start: Start location
    :type start: int
    :param end: End location, optional
    :type end: int
    :return: List of variant strings
    :rtype: List[str]
    """

    variants = fetch_ensembl_variant_info(build, chr, suffix, start, end)
    return format_ensembl_variant_list(variants)
