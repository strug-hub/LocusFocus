from enum import Enum

from gtex_openapi.api.datasets_endpoints_api import DatasetsEndpointsApi
from gtex_openapi.api_client import ApiClient
from gtex_openapi.configuration import Configuration
from gtex_openapi.models.chromosome import Chromosome
from gtex_openapi.models.dataset_id import DatasetId

from app.utils.helpers import validate_chromosome


# this is circuitous but seems more convenient than passing in enum val from callers
def get_chromosome_enum(chr: str):
    return Chromosome[Chromosome(chr).name]


def get_dataset_id_enum(dataset_id: str):
    return DatasetId[DatasetId(dataset_id).name]


HOSTNAME = "https://gtexportal.org"

configuration = Configuration(host=HOSTNAME)


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


def get_variants(dataset_id: str, start: int, end: int, chromosome: str):
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

        return instance.get_variant_by_location_api_v2_dataset_variant_by_location_get(
            dataset_id=dataset_id,
            start=start,
            end=end,
            chromosome=chromosome,
            items_per_page=100000,
        )
