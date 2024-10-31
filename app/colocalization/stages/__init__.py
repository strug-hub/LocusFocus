"""
Stages for running colocalization.

Expected order:
- create_session
- collect_user_input
- read_gwas_file
- read_secondary_datasets
- report_gtex_data
- ss_subset_gwas
- get_ld_matrix
- coloc_simple_sum
- finalize_results
"""
from .create_session import CreateSessionStage
from .collect_user_input import CollectUserInputStage
from .read_gwas_file import ReadGWASFileStage
from .read_secondary_datasets import ReadSecondaryDatasetsStage
from .report_gtex_data import ReportGTExDataStage
from .get_ld_matrix import GetLDMatrixStage
from .ss_subset_gwas import SimpleSumSubsetGWASStage
from .coloc_simple_sum import ColocSimpleSumStage
from .finalize_results import FinalizeResultsStage


__all__ = [
    'CreateSessionStage',
    'CollectUserInputStage',
    'ReadGWASFileStage',
    'ReadSecondaryDatasetsStage',
    'ReportGTExDataStage',
    'GetLDMatrixStage',
    'SimpleSumSubsetGWASStage',
    'ColocSimpleSumStage',
    'FinalizeResultsStage',
]
