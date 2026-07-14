"""
Abstract interface for GTEx database access.

Both the production MongoDB-backed implementation (RealGTExDatabase) and the
in-memory fake used in tests (FakeGTExDatabase) implement this interface.
"""

from abc import ABC, abstractmethod
from typing import List, Optional

import pandas as pd


class GTExDatabase(ABC):
    """
    Provides access to GTEx eQTL and variant data.

    All methods that accept `version` expect either `"V8"` or `"V10"`.
    All methods that accept `chrom` expect a numeric string without a "chr"
    prefix (e.g. `"1"`, `"23"` for chromosome X).
    """

    @abstractmethod
    def list_tissues(self, version: str) -> List[str]:
        """Return the sorted list of tissue names available for `version`."""

    def list_genes(self, version: str) -> Optional[List[str]]:
        """Return gene names available in this database, or None if unrestricted.

        None means the caller should use its own gene source (e.g. GENCODE) without
        filtering.  A list (including an empty one) means only those genes are
        available.
        """
        return None

    @abstractmethod
    def get_eqtl_data(
        self, version: str, tissue: str, ensg_id_prefix: str
    ) -> pd.DataFrame:
        """Return eQTL data for a gene in a tissue.

        Parameters
        ----------
        version:
            GTEx version, `"V8"` or `"V10"`.
        tissue:
            Tissue name (spaces replaced with underscores). e.g. `"Artery_Aorta"`.
        ensg_id_prefix:
            ENSG gene ID **without** the version suffix, e.g. `"ENSG00000069275"`.

        Returns
        -------
        pd.DataFrame
            One row per eQTL variant.  Guaranteed columns:

            * `variant_id`   - `"{chrom}_{pos}_{ref}_{alt}_b38"`
            * `rs_id`        - dbSNP rs identifier
            * `chr`          - chromosome (int)
            * `pos`          - chromosomal position (int)
            * `ref`          - reference allele
            * `alt`          - alternate allele
            * `pval`         - association p-value (float)
            * `beta`         - effect size estimate (float)
            * `se`           - standard error (float)
            * `sample_maf`   - minor allele frequency in the GTEx sample (float)
            * `ma_count`     - minor allele count (int)

            Returns an **empty** DataFrame when no data is available for the
            given gene/tissue combination.
        """

    @abstractmethod
    def get_variants_by_region(
        self, start: int, end: int, chrom: str, version: str
    ) -> pd.DataFrame:
        """Return all variants in the chromosomal region `[start, end]`.

        Parameters
        ----------
        start:
            Start position (inclusive, 1-based).
        end:
            End position (inclusive).
        chrom:
            Chromosome as a numeric string (e.g. `"1"`, `"23"` for X).
        version:
            GTEx version, `"V8"` or `"V10"`.

        Returns
        -------
        pd.DataFrame
            One row per variant.  Guaranteed columns:
            `variant_id`, `rs_id`, `chr`, `pos`, `ref`, `alt`.

            Returns an **empty** DataFrame when the region contains no known
            variants.
        """
