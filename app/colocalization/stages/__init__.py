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

from app.colocalization.stages.create_session import CreateSessionStage
from app.colocalization.stages.collect_user_input import CollectUserInputStage
from app.colocalization.stages.read_gwas_file import ReadGWASFileStage
from app.colocalization.stages.read_secondary_datasets import ReadSecondaryDatasetsStage
from app.colocalization.stages.report_gtex_data import ReportGTExDataStage
from app.colocalization.stages.get_ld_matrix import GetLDMatrixStage
from app.colocalization.stages.ss_subset_gwas import SimpleSumSubsetGWASStage
from app.colocalization.stages.coloc_simple_sum import ColocSimpleSumStage
from app.colocalization.stages.finalize_results import FinalizeResultsStage
from app.colocalization.stages.liftover_gwas_file import LiftoverGWASFile
from app.colocalization.stages.liftover_secondary_datasets import LiftoverSecondaryDatasets


__all__ = [
    "CreateSessionStage",
    "CollectUserInputStage",
    "ReadGWASFileStage",
    "LiftoverGWASFile",
    "ReadSecondaryDatasetsStage",
    "LiftoverSecondaryDatasets",
    "ReportGTExDataStage",
    "GetLDMatrixStage",
    "SimpleSumSubsetGWASStage",
    "ColocSimpleSumStage",
    "FinalizeResultsStage",
]
