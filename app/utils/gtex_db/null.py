"""
No-op GTExDatabase implementation used when no MongoDB connection is configured.

All methods return empty results, allowing the app to start and serve non-GTEx
functionality without a live database.
"""

from typing import List, Optional

import pandas as pd

from app.utils.gtex_db.base import GTExDatabase


class NullGTExDatabase(GTExDatabase):
    """GTExDatabase that always returns empty results."""

    def list_tissues(self, version: str) -> List[str]:
        return []

    def list_genes(self, version: str) -> Optional[List[str]]:
        return []

    def get_eqtl_data(
        self, version: str, tissue: str, ensg_id_prefix: str
    ) -> pd.DataFrame:
        return pd.DataFrame()

    def get_variants_by_region(
        self, start: int, end: int, chrom: str, version: str
    ) -> pd.DataFrame:
        return pd.DataFrame()
