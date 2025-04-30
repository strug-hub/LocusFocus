from dataclasses import dataclass
from typing import Any, List

from gtex_openapi.models.pagination_info import PaginationInfo

from app.utils.apis.gtex import (
    fetch_all,
    get_bulk_eqtl,
    get_eqtl,
    get_genes,
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
    results = get_eqtl(
        dataset_id="gtex_v10",
        gencode_id="ENSG00000005436.14",
        tissue_site="Liver",
        variant_id="chr7_95404491_A_T_b38",
    )

    assert results.data is not None
    assert isinstance(results.data, list)
    assert isinstance(results.p_value, float)


def test_can_fetch_eqtl_bulk():

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

    results = get_bulk_eqtl(dataset_id="gtex_v10", body=body)

    assert results["data"] is not None
    assert len(results["errors"]) == 0
    assert isinstance(results["data"], list)
    assert isinstance(results["data"][0]["pValue"], float)


def test_can_fetch_genes():
    """Sanity check for gene fetch"""
    results = get_genes(build="hg38", gene_symbols=["NUCKS1", "CDK18"])

    assert results.data is not None
    assert len(results.data) == 2
    assert results.paging_info.number_of_pages == 1
    assert results.paging_info.page == 0
