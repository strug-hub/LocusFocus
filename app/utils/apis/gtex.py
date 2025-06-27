from collections.abc import Callable
from typing import Any, Dict, List, Optional

from gtex_openapi.api.static_association_endpoints_api import StaticAssociationEndpointsApi
from gtex_openapi.api.datasets_endpoints_api import DatasetsEndpointsApi
from gtex_openapi.api.dynamic_association_endpoints_api import (
    DynamicAssociationEndpointsApi,
)
from gtex_openapi.api.reference_genome_endpoints_api import ReferenceGenomeEndpointsApi
from gtex_openapi.api_client import ApiClient
from gtex_openapi.configuration import Configuration
from gtex_openapi.models.chromosome import Chromosome
from gtex_openapi.models.dataset_id import DatasetId
from gtex_openapi.models.dynamic_eqtl_body import DynamicEqtlBody
from gtex_openapi.models.eqtl import Eqtl
from gtex_openapi.models.tissue_site_detail_id import TissueSiteDetailId
from gtex_openapi.models.paginated_response_variant import PaginatedResponseVariant
from gtex_openapi.models.paginated_response_gene import PaginatedResponseGene
from gtex_openapi.models.app_models_request_parameters_genome_build import AppModelsRequestParametersGenomeBuild
from gtex_openapi.models.post_dynamic_eqtl_result import PostDynamicEqtlResult
from gtex_openapi.models.paginated_response_independent_eqtl import PaginatedResponseIndependentEqtl
from gtex_openapi.models.tissuesitedetailid_inner import TissuesitedetailidInner
from gtex_openapi.models.paginated_response_single_tissue_eqtl import PaginatedResponseSingleTissueEqtl

from app.utils.helpers import validate_chromosome


# caller passes in string and we replace with enum value here
# this is a little backward but more convenient than callers looking up enums?
def get_chromosome_enum(chrom: str):
    return Chromosome(chrom)


def get_dataset_id_enum(dataset_id: str):
    return DatasetId(dataset_id)


def get_tissue_site_detail_id_enum(tissue_site_detail: str):
    return TissueSiteDetailId(tissue_site_detail)


def get_genome_build_enum(genome_build: str):
    return AppModelsRequestParametersGenomeBuild(genome_build)


HOSTNAME = "https://gtexportal.org"

CONFIGURATION = Configuration(host=HOSTNAME)


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


def get_genes(build: str, gene_symbols: List[str]) -> PaginatedResponseGene:
    """Fetch genes from the GTEx API.

    Useful for getting GENCODE IDs from a list of gene symbols.

    See https://gtexportal.org/api/v2/redoc#tag/Reference-Genome-Endpoints/operation/get_genes_api_v2_reference_gene_get
    for more details.

    :param build: The build to use (hg19 or hg38)
    :type build: str
    :param geneSymbols: A list of gene symbols to fetch. eg. ["NUCKS1", "CDK18"]
    :type geneSymbols: List[str]
    :return: A paginated response containing the genes
    :rtype: PaginatedResponseGene
    """
    build = build.lower()
    if build not in ["hg19", "hg38"]:
        raise ValueError("build must be 'hg19' or 'hg38'")

    if build == "hg19":
        build = "GRCh37/hg19"
    elif build == "hg38":
        build = "GRCh38/hg38"

    build = get_genome_build_enum(build)

    with ApiClient(CONFIGURATION) as api_client:
        instance = ReferenceGenomeEndpointsApi(api_client)

        return fetch_all(
            instance.get_genes_api_v2_reference_gene_get,
            gene_id=gene_symbols,
            genome_build=build,
            items_per_page=100000,
        )


def get_independent_eqtls(
    dataset_id: str,
    gencode_ids: List[str],
    tissue_sites: List[str],
) -> PaginatedResponseIndependentEqtl:
    """Fetch independent EQTL Data

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param gencode_ids: A list of versioned GENCODE IDs of genes. eg. ["ENSG00000005436.14", "ENSG00000225972.1"]
    :type gencode_ids: List[str]
    :param tissue_sites: A list of tissue sites to use in the calculation. eg. ["Liver", "Lung"]
    :type tissue_sites: List[str], a key of gtex_openapi.models.tissue_site_detail_id.TissueSiteDetailId enum.
    :return: Independent EQTL results
    :rtype: IndependentEqtl
    """
    tissue_sites = [TissuesitedetailidInner(get_tissue_site_detail_id_enum(x)) for x in tissue_sites] # type: ignore
    dataset_id = get_dataset_id_enum(dataset_id)
    with ApiClient(CONFIGURATION) as api_client:
        instance = StaticAssociationEndpointsApi(api_client)

        return fetch_all(
            instance.get_independent_eqtl_api_v2_association_independent_eqtl_get,
            gencode_id=gencode_ids,
            tissue_site_detail_id=tissue_sites,
            dataset_id=dataset_id,
            items_per_page=100000,
        )


def get_dynamic_eqtl(
    dataset_id: str,
    gencode_id: str,
    tissue_site: str,
    variant_id: str,
) -> Eqtl:
    """Fetch dynamic EQTL Data

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param gencode_id: A versioned GENCODE ID of a gene. eg. "ENSG00000005436.14"
    :type gencode_id: str
    :param tissue_site: The tissue site to use in the calculation. eg. "Liver"
    :type tissue_site: str, a key of gtex_openapi.models.tissue_site_detail_id.TissueSiteDetailId enum.
    :param variant_id: A GTEx variant id (must begin with `chr`). eg. "chr7_95404491_A_T_b38"
    :type variant_id: str
    :return: The calculation result
    :rtype: Eqtl
    """

    tissue_site = get_tissue_site_detail_id_enum(tissue_site)

    dataset_id = get_dataset_id_enum(dataset_id)

    if not variant_id.startswith("chr"):
        raise ValueError("Variant ID must start with 'chr' prefix!")

    with ApiClient(CONFIGURATION) as api_client:
        instance = DynamicAssociationEndpointsApi(api_client)

        return instance.calculate_expression_quantitative_trait_loci_api_v2_association_dyneqtl_get(
            dataset_id=dataset_id,
            gencode_id=gencode_id,
            tissue_site_detail_id=tissue_site,
            variant_id=variant_id,
        )


def get_bulk_dynamic_eqtl(dataset_id: str, body: List[Dict[str, Any]]) -> PostDynamicEqtlResult:
    """Fetch dynamic EQTL Data in bulk

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :param body: A list of dictionaries containing the tissue site detail, variant_id, and gencode_id for each desired EQTL result
    :type body: List[Dict[Literal['gencode_id', 'tissue_site_detail_id', 'variant_id], Any]]
                see docstring for `get_dynamic_eqtl` for type details.
    :return: The calculation result
    :rtype: PostDynamicEqtlResult
    :
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

    with ApiClient(CONFIGURATION) as api_client:
        instance = DynamicAssociationEndpointsApi(api_client)

        # request is not paginated, so no need to use fetch_all
        results = instance.bulk_calculate_expression_quantitative_trait_loci_api_v2_association_dyneqtl_post(
            dataset_id=dataset_id,
            dynamic_eqtl_body=body_args,
        )

        return results


def get_tissue_site_details(dataset_id: str):
    """Fetch tissue details

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`)
    :type dataset_id: str
    :return: The paginated response
    :rtype: PaginatedResponseTissueSiteDetail
    """

    dataset_id = get_dataset_id_enum(dataset_id)

    with ApiClient(CONFIGURATION) as api_client:
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

    with ApiClient(CONFIGURATION) as api_client:
        instance = DatasetsEndpointsApi(api_client)

        return fetch_all(
            instance.get_variant_by_location_api_v2_dataset_variant_by_location_get,
            dataset_id=dataset_id,
            start=start,
            end=end,
            chromosome=chromosome,
            items_per_page=100000,
        )


def get_significant_single_tissue_eqtls(dataset_id: str, gencode_ids: Optional[List[str]], variant_ids: Optional[List[str]], tissue_site_detail_ids: Optional[List[str]]) -> PaginatedResponseSingleTissueEqtl:
    """Fetch single-tissue eQTLs from GTEx API

    :param dataset_id: The identifier of the gtex dataset (`gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot`)
    :type dataset_id: str
    :param gencode_ids: A list of versioned GENCODE IDs of genes. eg. ["ENSG00000005436.14", "ENSG00000225972.1"]
                        If None, then all genes in the dataset will be fetched based on other filters.
    :type gencode_ids: Optional[List[str]]
    :param variant_ids: A list of variant IDs to fetch eQTLs for. Must be in format chr_pos_ref_alt_build, eg. chr1_205381100_C_T_b38.
                        If None, then all variants in the dataset will be fetched based on other filters.
    :type variant_ids: Optional[List[str]]
    :param tissue_site_detail_ids: A list of tissue site details to fetch eQTLs for.
                                   If None, then all tissue sites in the dataset will be fetched based on other filters.
    :type tissue_site_detail_ids: Optional[List[str]], keys of gtex_openapi.models.tissue_site_detail_id.TissueSiteDetailId enum
    :return: The paginated eqtl results
    :rtype: PaginatedResponseSingleTissueEqtl
    """
    if tissue_site_detail_ids is not None:
        tissue_site_detail_ids = [TissuesitedetailidInner(get_tissue_site_detail_id_enum(x)) for x in tissue_site_detail_ids] # type: ignore

    dataset_id = get_dataset_id_enum(dataset_id)

    with ApiClient(CONFIGURATION) as api_client:
        instance = StaticAssociationEndpointsApi(api_client)

        return fetch_all(
            instance.get_significant_single_tissue_eqtls_api_v2_association_single_tissue_eqtl_get,
            dataset_id=dataset_id,
            gencode_id=gencode_ids,
            variant_id=variant_ids,
            tissue_site_detail_id=tissue_site_detail_ids,
            items_per_page=100000,
        )
