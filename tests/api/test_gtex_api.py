from dataclasses import dataclass
from typing import Any, List

import pytest
import pandas as pd
from gtex_openapi.models.pagination_info import PaginationInfo
from gtex_openapi.exceptions import BadRequestException

from app.utils.apis.gtex import (
    fetch_all,
    get_bulk_dynamic_eqtl,
    get_dynamic_eqtl,
    get_genes,
    get_independent_eqtls,
    get_significant_single_tissue_eqtls,
    get_tissue_site_details,
    get_variants,
)


@dataclass
class DummyResponse:
    data: List[Any]
    paging_info: PaginationInfo


def test_can_fetch_v8_tissues():
    """Sanity check for v8 tissue fetch"""
    results = get_tissue_site_details(dataset_id="gtex_v8")
    assert len(results.data) == 54


def test_can_fetch_v10_tissues():
    """Sanity check for v10 tissue fetch"""
    results = get_tissue_site_details(dataset_id="gtex_v10")
    assert len(results.data) == 54


def test_fetch_all():
    """Test the recursive fetch_all wrapper on a dummy paginated API function"""

    def api_faker(page: int, some_arg: int):
        """API faker that returns 3 'pages' of data"""
        if page == 0:
            return DummyResponse(
                data=[1, some_arg],
                paging_info=PaginationInfo(
                    numberOfPages=3, page=1, maxItemsPerPage=2, totalNumberOfItems=6
                ),
            )
        elif page == 1:
            return DummyResponse(
                data=[3, 4],
                paging_info=PaginationInfo(
                    numberOfPages=3, page=2, maxItemsPerPage=2, totalNumberOfItems=6
                ),
            )
        elif page == 2:
            return DummyResponse(
                data=[5, some_arg],
                paging_info=PaginationInfo(
                    numberOfPages=3, page=3, maxItemsPerPage=2, totalNumberOfItems=6
                ),
            )
        else:
            raise ValueError(f"bad page : {page}!")

    some_arg = 100
    result = fetch_all(api_faker, some_arg=some_arg)

    assert len(result.data) == 6

    assert some_arg in result.data

    assert result.paging_info.max_items_per_page == 6

    assert result.paging_info.page == 0


def test_can_fetch_v10_variants_from_region_string():
    """Sanity check for v10 variant fetch"""
    region_string = "chr11:0-200000"
    chr, pos = region_string.split(":")
    start, end = pos.split("-")
    results = get_variants(
        dataset_id="gtex_v10", start=int(start), end=int(end), chromosome=chr
    )

    assert len(results.data) > 0
    assert results.data[0].snp_id.startswith("rs")


def test_can_fetch_v8_variants_from_region_string():
    """Sanity check for v10 variant fetch"""
    region_string = "chr11:0-200000"
    chr, pos = region_string.split(":")
    start, end = pos.split("-")

    results = get_variants(
        dataset_id="gtex_v8", start=int(start), end=int(end), chromosome=chr
    )

    assert len(results.data) > 0
    assert results.data[0].variant_id.startswith("chr11")


def test_can_fetch_eqtl():
    results = get_dynamic_eqtl(
        dataset_id="gtex_v10",
        gencode_id="ENSG00000005436.14",
        tissue_site="Liver",
        variant_id="chr7_95404491_A_T_b38",
    )

    assert results.data is not None
    assert isinstance(results.data, list)
    assert isinstance(results.p_value, float)


def test_can_fetch_eqtl_bulk():
    """Sanity check for bulk eqtl fetch"""
    body = [
        {
            "variant_id": "chr7_95404491_A_T_b38",
            "gencode_id": "ENSG00000005436.14",
            "tissue_site_detail_id": "Liver",
        },
        {
            "gencode_id": "ENSG00000225972.1",
            "tissue_site_detail_id": "Adipose_Visceral_Omentum",
            "variant_id": "chr1_629115_C_T_b38",
        },
    ]

    results = get_bulk_dynamic_eqtl(dataset_id="gtex_v10", body=body)

    assert results.data is not None
    assert len(results.errors) == 0
    assert isinstance(results.data, list)
    assert isinstance(results.data[0].to_dict()["pValue"], float)


def test_can_fetch_independent_eqtl():
    """Sanity check for independent eqtl fetch"""
    results = get_independent_eqtls(
        dataset_id="gtex_v8",
        gencode_ids=["ENSG00000005436.14", "ENSG00000225972.1"],
        tissue_sites=["Liver", "Adipose_Visceral_Omentum"],
    )
    assert results.data is not None
    assert len(results.data) > 0


def test_eqtl_fetch_pipeline():
    """
    Test fetching eQTLs using information from other API endpoints.
    """
    gtex_version = "gtex_v8"
    print("Fetching all tissues for gtex_v8...")
    tissues_response = get_tissue_site_details(dataset_id=gtex_version)
    assert tissues_response.data is not None
    assert len(tissues_response.data) > 0
    tissues = [x.tissue_site_detail_id.value for x in tissues_response.data]

    gene_symbols = ["NUCKS1", "CDK18", "RAB7L1", "SLC41A1", "PM20D1"]
    print("Fetching gene info for symbols " + str(gene_symbols) + "...")
    gene_response = get_genes(build="hg38", gene_symbols=gene_symbols)
    assert gene_response.data is not None
    assert len(gene_response.data) <= len(gene_symbols)
    gencode_ids = [x.gencode_id for x in gene_response.data]

    print("Fetching eQTLs for " + str(gencode_ids) + " in " + str(tissues) + "...")
    results = get_independent_eqtls(
        dataset_id=gtex_version,
        gencode_ids=gencode_ids,
        tissue_sites=tissues,
    )
    assert results.data is not None
    assert len(results.data) > 0


def test_can_fetch_genes():
    """Sanity check for gene fetch"""
    results = get_genes(build="hg38", gene_symbols=["NUCKS1", "CDK18"])

    assert results.data is not None
    assert len(results.data) == 2
    assert results.paging_info.number_of_pages == 1
    assert results.paging_info.page == 0


def test_eqtl_fetch_equivalence():
    """
    Check that the eQTL fetch API results are the same as bulk fetch.
    """

    variant_id = "chr1_205381100_C_T_b38"
    tissue_site = "Liver"
    gencode_id = "ENSG00000069275.12"  # "NUCKS1"
    dataset_id = "gtex_v8"

    print("Fetching single eQTL...")
    single_eqtl_results = get_dynamic_eqtl(
        dataset_id=dataset_id,
        gencode_id=gencode_id,
        tissue_site=tissue_site,
        variant_id=variant_id,
    )
    print("Done fetching single eQTL")
    assert single_eqtl_results.data is not None
    assert len(single_eqtl_results.data) > 0

    print("Fetching bulk eQTL...")
    bulk_eqtl_results = get_bulk_dynamic_eqtl(
        dataset_id=dataset_id,
        body=[
            {
                "variant_id": variant_id,
                "gencode_id": gencode_id,
                "tissue_site_detail_id": tissue_site,
            }
        ],
    )
    print("Done fetching bulk eQTL")
    assert bulk_eqtl_results.data is not None
    assert len(bulk_eqtl_results.data) > 0

    # Check equivalence
    single_dict = single_eqtl_results.to_dict()
    bulk_dict = bulk_eqtl_results.data[0].to_dict()
    assert abs(single_dict["pValue"] - bulk_dict["pValue"]) < 1e-6
    assert abs(single_dict["error"] - bulk_dict["error"]) < 1e-6
    assert abs(single_dict["nes"] - bulk_dict["nes"]) < 1e-6
    assert abs(single_dict["tStatistic"] - bulk_dict["tStatistic"]) < 1e-6
    # Missing! maf, hetCount, homoAltCount, homoRefCount


def test_significant_single_tissue_eqtl():
    """Sanity check for single-tissue eQTL fetch"""
    results = get_significant_single_tissue_eqtls(
        dataset_id="gtex_v10",
        gencode_ids=["ENSG00000005436.14", "ENSG00000225972.1"],
        tissue_site_detail_ids=["Liver", "Adipose_Visceral_Omentum"],
        variant_ids=None,
    )

    assert results.data is not None
    assert len(results.data) > 0
