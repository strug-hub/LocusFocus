"""
In-memory fake GTEx database for testing and development.

``FakeGTExDatabase`` implements the ``GTExDatabase`` interface using
pre-generated synthetic data — no running MongoDB is required.

For development purposes (e.g. running the full Flask app locally against fake
data), ``FakeGTExDatabase`` also exposes ``seed`` / ``teardown`` helpers that
write the same data into a real MongoDB instance.  Use
``setup/seed_fake_gtex_db.py`` for that workflow.

Reference data
--------------
Gene coordinates are taken from ``data/collapsed_gencode_v26_hg38.gz`` so that
gene-name lookups in ``app/utils/gtex.py`` resolve correctly.  All five genes
sit on chr1 within the default LocusFocus example region (≈205–206 Mb).
"""

import random
from typing import Any, Dict, List, Optional

import pandas as pd
from faker import Faker
from pymongo import MongoClient, ASCENDING

from app.utils.gtex_db.base import GTExDatabase


# ---------------------------------------------------------------------------
# Reference data: real ENSG IDs from collapsed_gencode_v26_hg38.gz
# ---------------------------------------------------------------------------
GENES: Dict[str, Dict[str, Any]] = {
    "NUCKS1": {
        "ensg_id": "ENSG00000069275.12",
        "chrom": 1,
        "start": 205_712_819,
        "end": 205_750_276,
    },
    "CDK18": {
        "ensg_id": "ENSG00000117266.15",
        "chrom": 1,
        "start": 205_504_595,
        "end": 205_532_793,
    },
    "ELK4": {
        "ensg_id": "ENSG00000158711.13",
        "chrom": 1,
        "start": 205_597_556,
        "end": 205_631_962,
    },
    "PM20D1": {
        "ensg_id": "ENSG00000162877.12",
        "chrom": 1,
        "start": 205_828_022,
        "end": 205_850_132,
    },
    "SLC26A9": {
        "ensg_id": "ENSG00000174502.18",
        "chrom": 1,
        "start": 205_913_048,
        "end": 205_943_460,
    },
}

TISSUES: List[str] = [
    "Liver",
    "Brain_Cortex",
    "Whole_Blood",
    "Lung",
    "Muscle_Skeletal",
]

_ALLELES = ["A", "T", "G", "C"]
_VERSIONS = ("V8", "V10")


class FakeGTExDatabase(GTExDatabase):
    """
    In-memory implementation of :class:`~app.utils.gtex_db.base.GTExDatabase`
    that uses `Faker <https://faker.readthedocs.io>`_ to generate synthetic
    eQTL and variant data.

    All data is pre-generated in ``__init__`` so that multiple calls to any
    method always return the same values (fully reproducible for a given
    ``seed``).

    Typical usage in a pytest fixture::

        @pytest.fixture()
        def fake_gtex_db(flask_app):
            fake = FakeGTExDatabase()
            flask_app.extensions["gtex_db"] = fake
            yield fake

    To populate a real MongoDB for manual UI testing::

        from pymongo import MongoClient
        client = MongoClient("mongodb://localhost:27017")
        FakeGTExDatabase().seed(client)
    """

    def __init__(
        self,
        seed: int = 42,
        genes: Dict[str, Dict[str, Any]] = GENES,
        tissues: List[str] = TISSUES,
        n_variants_per_gene: int = 50,
    ) -> None:
        self.genes = genes
        self.tissues = tissues
        self.n_variants_per_gene = n_variants_per_gene

        self._rng = random.Random(seed)
        Faker.seed(seed)
        self._fake = Faker()

        # Pre-generate all data once so the instance is immutable after init.
        self._gene_variants: Dict[str, List[Dict]] = self._build_gene_variants()
        self._tissue_eqtls: Dict[str, Dict[str, List[Dict]]] = (
            self._build_tissue_eqtls()
        )

    # ------------------------------------------------------------------
    # GTExDatabase interface
    # ------------------------------------------------------------------

    def list_tissues(self, version: str) -> List[str]:
        return sorted(self.tissues)

    def get_eqtl_data(
        self, version: str, tissue: str, ensg_id_prefix: str
    ) -> pd.DataFrame:
        gene_name = self._find_gene_by_ensg_prefix(ensg_id_prefix)
        if gene_name is None or tissue not in self._tissue_eqtls:
            return pd.DataFrame()

        eqtl_df = pd.DataFrame(self._tissue_eqtls[tissue][gene_name])
        variant_cols = pd.DataFrame(
            [
                {
                    "variant_id": v["variant_id"],
                    "rs_id": v["rs_id"],
                    "chr": v["chrom"],
                    "pos": v["pos"],
                    "ref": v["ref"],
                    "alt": v["alt"],
                }
                for v in self._gene_variants[gene_name]
            ]
        )
        return pd.merge(eqtl_df, variant_cols, on="variant_id")

    def get_variants_by_region(
        self, start: int, end: int, chrom: str, version: str
    ) -> pd.DataFrame:
        chrom_int = int(str(chrom).lower().replace("chr", "").replace("x", "23"))
        rows = [
            {
                "variant_id": v["variant_id"],
                "rs_id": v["rs_id"],
                "chr": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
            }
            for variants in self._gene_variants.values()
            for v in variants
            if v["chrom"] == chrom_int and start <= v["pos"] <= end
        ]
        return pd.DataFrame(rows) if rows else pd.DataFrame()

    # ------------------------------------------------------------------
    # Optional MongoDB seed / teardown (for development use)
    # ------------------------------------------------------------------

    def seed(
        self,
        client: MongoClient,
        versions: tuple = _VERSIONS,
    ) -> None:
        """Insert all fake data into a MongoDB instance.

        Drops and recreates every affected collection, so calling ``seed``
        multiple times is safe.
        """
        for version in versions:
            db = client[f"GTEx_{version}"]

            coll = db["variant_table"]
            coll.drop()
            coll.insert_many(self._variant_table_docs(version))
            if version == "V8":
                coll.create_index(
                    [("chr", ASCENDING), ("variant_pos", ASCENDING)]
                )
            else:
                coll.create_index([("chr", ASCENDING), ("pos", ASCENDING)])
            coll.create_index("variant_id")

            for tissue in self.tissues:
                tcoll = db[tissue]
                tcoll.drop()
                tcoll.insert_many(self._tissue_docs(tissue))
                tcoll.create_index("gene_id")

    def teardown(
        self,
        client: MongoClient,
        versions: tuple = _VERSIONS,
    ) -> None:
        """Drop all collections created by :meth:`seed`."""
        for version in versions:
            db = client[f"GTEx_{version}"]
            db["variant_table"].drop()
            for tissue in self.tissues:
                db[tissue].drop()

    # ------------------------------------------------------------------
    # Internal data builders
    # ------------------------------------------------------------------

    def _find_gene_by_ensg_prefix(self, ensg_id_prefix: str) -> Optional[str]:
        """Return the gene name whose ENSG ID starts with ``ensg_id_prefix``."""
        return next(
            (
                name
                for name, info in self.genes.items()
                if info["ensg_id"].startswith(ensg_id_prefix)
            ),
            None,
        )

    def _build_gene_variants(self) -> Dict[str, List[Dict]]:
        """Generate a shared pool of variants for every gene."""
        gene_variants: Dict[str, List[Dict]] = {}
        rs_counter = 1_000_000
        for gene_name, info in self.genes.items():
            chrom = info["chrom"]
            n = min(self.n_variants_per_gene, info["end"] - info["start"])
            positions = sorted(
                self._rng.sample(range(info["start"], info["end"]), n)
            )
            variants: List[Dict] = []
            for pos in positions:
                ref = self._rng.choice(_ALLELES)
                alt = self._rng.choice([a for a in _ALLELES if a != ref])
                maf = round(self._rng.uniform(0.01, 0.49), 4)
                variants.append(
                    {
                        "variant_id": f"{chrom}_{pos}_{ref}_{alt}_b38",
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "rs_id": f"rs{rs_counter}",
                        "maf": maf,
                        "ma_count": self._rng.randint(5, 500),
                    }
                )
                rs_counter += 1
            gene_variants[gene_name] = variants
        return gene_variants

    def _build_tissue_eqtls(self) -> Dict[str, Dict[str, List[Dict]]]:
        """Generate tissue- and gene-specific eQTL statistics for every variant."""
        tissue_eqtls: Dict[str, Dict[str, List[Dict]]] = {}
        for tissue in self.tissues:
            tissue_eqtls[tissue] = {}
            for gene_name, variants in self._gene_variants.items():
                gene_info = self.genes[gene_name]
                gene_mid = (gene_info["start"] + gene_info["end"]) // 2
                records: List[Dict] = []
                for v in variants:
                    pval = 10.0 ** (-self._rng.uniform(0.0, 15.0))
                    beta = self._rng.gauss(0.0, 0.3)
                    se = abs(self._rng.gauss(0.05, 0.02)) + 0.001
                    ma_count = v["ma_count"]
                    maf = v["maf"]
                    records.append(
                        {
                            "variant_id": v["variant_id"],
                            "pval": pval,
                            "beta": beta,
                            "se": se,
                            "sample_maf": maf,
                            "ma_count": ma_count,
                            "tss_distance": v["pos"] - gene_mid,
                            "ma_samples": max(1, int(ma_count / (2.0 * maf))),
                        }
                    )
                tissue_eqtls[tissue][gene_name] = records
        return tissue_eqtls

    def _variant_table_docs(self, version: str) -> List[Dict]:
        """Build variant_table documents for MongoDB seeding."""
        docs: List[Dict] = []
        for variants in self._gene_variants.values():
            for v in variants:
                doc: Dict[str, Any] = {
                    "variant_id": v["variant_id"],
                    "chr": v["chrom"],
                    "ref": v["ref"],
                    "alt": v["alt"],
                }
                if version.upper() == "V8":
                    doc["variant_pos"] = v["pos"]
                    doc["rs_id_dbSNP151_GRCh38p7"] = v["rs_id"]
                else:
                    doc["pos"] = v["pos"]
                    doc["rs_id_dbSNP155_GRCh38p13"] = v["rs_id"]
                docs.append(doc)
        return docs

    def _tissue_docs(self, tissue: str) -> List[Dict]:
        """Build tissue-collection documents for MongoDB seeding."""
        return [
            {
                "gene_id": self.genes[gene_name]["ensg_id"],
                "eqtl_variants": self._tissue_eqtls[tissue][gene_name],
            }
            for gene_name in self.genes
        ]
