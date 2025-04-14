from collections.abc import Callable
from typing import Any, Dict, List

from gtex_openapi.api.datasets_endpoints_api import DatasetsEndpointsApi
from gtex_openapi.api.dynamic_association_endpoints_api import (
    DynamicAssociationEndpointsApi,
)
from gtex_openapi.api_client import ApiClient
from gtex_openapi.configuration import Configuration
from gtex_openapi.models.chromosome import Chromosome
from gtex_openapi.models.dataset_id import DatasetId
from gtex_openapi.models.dynamic_eqtl_body import DynamicEqtlBody
from gtex_openapi.models.eqtl import Eqtl
from gtex_openapi.models.post_dynamic_eqtl_result import PostDynamicEqtlResult
from gtex_openapi.models.tissue_site_detail_id import TissueSiteDetailId
from gtex_openapi.models.paginated_response_variant import PaginatedResponseVariant

from app.utils.helpers import validate_chromosome


# caller passes in string and we replace with enum value here
# this is a little backward but more convenient than callers looking up enums?
def get_chromosome_enum(chr: str):
    return Chromosome(chr)


def get_dataset_id_enum(dataset_id: str):
    return DatasetId(dataset_id)


def get_tissue_site_detail_id_enum(tissue_site_detail: str):
    return TissueSiteDetailId(tissue_site_detail)


HOSTNAME = "https://gtexportal.org"

configuration = Configuration(host=HOSTNAME)


def fetch_all(func: Callable[..., Any], page: int | None = None, **kwargs):
    page = page or 0
    results = func(page=page, **kwargs)

    if results.paging_info.number_of_pages > page + 1:
        page += 1
        results.data.extend(fetch_all(func, page=page, **kwargs).data)
    results.paging_info.number_of_pages = 1
    results.paging_info.page = 0
    results.paging_info.max_items_per_page = results.paging_info.total_number_of_items
    return results


def get_eqtl(
    dataset_id: str,
    gencode_id: str,
    tissue_site: str,
    variant_id: str,
) -> Eqtl:
    """Fetch dynamic EQTL Data

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param gencode_id: A versioned GENCODE ID of a gene
    :type gencode_id: str
    :param tissue_site: The tissue site to use in the calculation
    :type tissue_site: str, a key of gtex_openapi.models.tissue_site_detail_id.TissueSiteDetailId enum
    :param variant_id: A GTEx variant id (must begin with `chr`)
    :type variant_id: str
    :return: The calculation result
    :rtype: Eqtl
    """

    tissue_site = get_tissue_site_detail_id_enum(tissue_site)

    dataset_id = get_dataset_id_enum(dataset_id)

    if not variant_id.startswith("chr"):
        raise ValueError("Variant ID must start with 'chr' prefix!")

    with ApiClient(configuration) as api_client:
        instance = DynamicAssociationEndpointsApi(api_client)

        return instance.calculate_expression_quantitative_trait_loci_api_v2_association_dyneqtl_get(
            dataset_id=dataset_id,
            gencode_id=gencode_id,
            tissue_site_detail_id=tissue_site,
            variant_id=variant_id,
        )


def get_bulk_eqtl(dataset_id: str, body: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Fetch dynamic EQTL Data in bulk

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param body: A list of dictionaries containing the tissue site detail, variant_id, and gencode_id for each desired EQTL result
    :type body: List[Dict[Literal['gencode_id', 'tissue_site_detail_id', 'variant_id], Any]]
                see docstring for `get_eqtl` for type details.
    :return: The calculation result
    :rtype: Dict[str, Any]
    """
    dataset_id = get_dataset_id_enum(dataset_id)

    body_args = []

    for arg in body:
        body_args.append(
            DynamicEqtlBody(
                tissueSiteDetailId=get_tissue_site_detail_id_enum(
                    arg["tissue_site_detail_id"]
                ),
                variantId=arg["variant_id"],
                gencodeId=arg["gencode_id"],
            )
        )

    with ApiClient(configuration) as api_client:
        instance = DynamicAssociationEndpointsApi(api_client)

        # request is not paginated, so no need to use fetch_all
        results = instance.bulk_calculate_expression_quantitative_trait_loci_api_v2_association_dyneqtl_post(
            dataset_id=dataset_id,
            dynamic_eqtl_body=body_args,
        )

        return results.to_dict()


def get_tissue_site_details(dataset_id: str):
    """Fetch tissue details

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :return: The paginated response
    :rtype: PaginatedResponseTissueSiteDetail
    """

    dataset_id = get_dataset_id_enum(dataset_id)

    with ApiClient(configuration) as api_client:
        instance = DatasetsEndpointsApi(api_client)

        return instance.get_tissue_site_detail_api_v2_dataset_tissue_site_detail_get(
            dataset_id=dataset_id
        )


def get_variants(
    dataset_id: str, start: int, end: int, chromosome: str
) -> PaginatedResponseVariant:
    """Fetch variants from gtex api, with a return limit of 100,000

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param start: The start position
    :type start: int
    :param end: The end position
    :type end: int
    :param chromosome: The chromosome (prefixed with 'chr'), `X` and `Y` for 23
    :type chromosome: str
    :return: The paginated variant list
    :rtype: PaginatedResponseVariant
    """

    validate_chromosome(chromosome, prefix="chr", x_y_numeric=False)

    chromosome = get_chromosome_enum(chromosome)
    dataset_id = get_dataset_id_enum(dataset_id)

    with ApiClient(configuration) as api_client:
        instance = DatasetsEndpointsApi(api_client)

        return fetch_all(
            instance.get_variant_by_location_api_v2_dataset_variant_by_location_get,
            dataset_id=dataset_id,
            start=start,
            end=end,
            chromosome=chromosome,
            items_per_page=100000,
        )
