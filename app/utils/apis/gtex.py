from enum import Enum

from gtex_openapi.api.datasets_endpoints_api import DatasetsEndpointsApi
from gtex_openapi.api_client import ApiClient
from gtex_openapi.configuration import Configuration


# API wants these passed in as enums
class GTEX_DATASET_IDS(Enum):
    V8 = "gtex_v8"
    V10 = "gtex_v10"


HOSTNAME = "https://gtexportal.org"

configuration = Configuration(host=HOSTNAME)


def get_tissue_site_details(dataset_id: str):
    """Fetch tissue details

    :param dataset_id: The identifier of the gtex dataset (`V8`, `V10`)
    :type dataset_id: str
    :return: The paginated response
    :rtype: PaginatedResponseTissueSiteDetail
    """

    if dataset_id not in ["V8", "V10"]:
        raise ValueError("dataset_id should be `V8` or `V10`")

    with ApiClient(configuration) as api_client:
        instance = DatasetsEndpointsApi(api_client)

        return instance.get_tissue_site_detail_api_v2_dataset_tissue_site_detail_get(
            dataset_id=GTEX_DATASET_IDS[dataset_id]
        )
