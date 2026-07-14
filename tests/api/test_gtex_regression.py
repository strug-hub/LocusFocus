"""
Regression tests for all GTEx-related code paths.

These tests use ``FakeGTExDatabase`` so no running MongoDB is required.
Coverage:
  - FakeGTExDatabase data integrity (column types, value ranges, variant_id format)
  - get_gtex error paths (V7, unknown tissue/gene, ENSG lookup, normalization)
  - get_variants_by_region input validation
  - get_gtex_data merging (rs-ID and variant_id snp lists, bad-format errors)
  - get_gtex_data_pvalues with both ID formats
  - get_gtex_snp_matches counting
  - Flask route error responses (V7, bad version, unknown tissue/gene, variant lookup)
"""

import re

import pandas as pd
import pytest
from flask import Flask

from app.utils.errors import InvalidUsage
from tests.fake_gtex import FakeGTExDatabase, GENES, TISSUES


# ---------------------------------------------------------------------------
# FakeGTExDatabase data-integrity tests  (no Flask context needed)
# ---------------------------------------------------------------------------


class TestFakeGTExDataIntegrity:
    """Structural and value-range checks on the generated fake data."""

    def setup_method(self):
        self.db = FakeGTExDatabase()

    def test_variant_id_format(self):
        pattern = re.compile(r"^\d+_\d+_[ATGC]_[ATGC]_b38$")
        for gene, variants in self.db._gene_variants.items():
            for v in variants:
                assert pattern.match(v["variant_id"]), (
                    f"Bad variant_id for {gene}: {v['variant_id']}"
                )

    def test_pval_in_open_unit_interval(self):
        for tissue in TISSUES:
            for gene, records in self.db._tissue_eqtls[tissue].items():
                for r in records:
                    assert 0 < r["pval"] <= 1, (
                        f"pval {r['pval']} out of range for {tissue}/{gene}"
                    )

    def test_se_positive(self):
        for tissue in TISSUES:
            for gene, records in self.db._tissue_eqtls[tissue].items():
                for r in records:
                    assert r["se"] > 0, f"se must be positive, got {r['se']}"

    def test_sample_maf_in_range(self):
        for tissue in TISSUES:
            for gene, records in self.db._tissue_eqtls[tissue].items():
                for r in records:
                    assert 0 < r["sample_maf"] < 0.5, (
                        f"sample_maf {r['sample_maf']} out of (0, 0.5)"
                    )

    def test_chrom_consistent_with_variant_id(self):
        for gene, variants in self.db._gene_variants.items():
            for v in variants:
                id_chrom = int(v["variant_id"].split("_")[0])
                assert id_chrom == v["chrom"], (
                    f"chrom mismatch: variant_id says {id_chrom}, chrom field says {v['chrom']}"
                )

    def test_pos_consistent_with_variant_id(self):
        for gene, variants in self.db._gene_variants.items():
            for v in variants:
                id_pos = int(v["variant_id"].split("_")[1])
                assert id_pos == v["pos"], (
                    f"pos mismatch in {v['variant_id']}: id says {id_pos}, pos field says {v['pos']}"
                )

    def test_eqtl_variant_ids_match_gene_variant_pool(self):
        """Every eQTL record must reference a variant in the gene's pool."""
        for tissue in TISSUES:
            for gene, records in self.db._tissue_eqtls[tissue].items():
                pool_ids = {v["variant_id"] for v in self.db._gene_variants[gene]}
                for r in records:
                    assert r["variant_id"] in pool_ids, (
                        f"eQTL variant_id {r['variant_id']} not in gene pool for {gene}"
                    )

    def test_get_eqtl_data_required_columns(self):
        df = self.db.get_eqtl_data("V8", "Liver", "ENSG00000069275")
        required = {"variant_id", "rs_id", "chr", "pos", "ref", "alt",
                    "pval", "beta", "se", "sample_maf", "ma_count"}
        assert required.issubset(df.columns), f"Missing: {required - set(df.columns)}"

    def test_get_eqtl_data_column_types(self):
        df = self.db.get_eqtl_data("V8", "Liver", "ENSG00000069275")
        assert pd.api.types.is_integer_dtype(df["chr"])
        assert pd.api.types.is_integer_dtype(df["pos"])
        assert pd.api.types.is_float_dtype(df["pval"])
        assert pd.api.types.is_float_dtype(df["se"])

    def test_get_eqtl_data_empty_for_unknown_gene(self):
        df = self.db.get_eqtl_data("V8", "Liver", "ENSG99999999999")
        assert df.empty

    def test_get_eqtl_data_empty_for_unknown_tissue(self):
        df = self.db.get_eqtl_data("V8", "Kidney_Cortex", "ENSG00000069275")
        assert df.empty

    def test_get_variants_by_region_empty_outside_all_genes(self):
        df = self.db.get_variants_by_region(1, 1000, "1", "V8")
        assert df.empty

    def test_get_variants_by_region_chr_prefix_stripped(self):
        nucks1 = GENES["NUCKS1"]
        df = self.db.get_variants_by_region(nucks1["start"], nucks1["end"], "chr1", "V8")
        assert len(df) > 0

    def test_list_tissues_is_sorted(self):
        tissues = self.db.list_tissues("V8")
        assert tissues == sorted(tissues)

    def test_list_tissues_version_agnostic(self):
        assert self.db.list_tissues("V8") == self.db.list_tissues("V10")


# ---------------------------------------------------------------------------
# get_gtex error paths
# ---------------------------------------------------------------------------


class TestGetGTExErrors:
    def test_v7_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            with pytest.raises(InvalidUsage, match="V7"):
                get_gtex("V7", "Liver", "NUCKS1")

    def test_unknown_tissue_raises_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            with pytest.raises(InvalidUsage) as exc_info:
                get_gtex("V8", "Kidney_Cortex", "NUCKS1")
            assert exc_info.value.status_code == 410

    def test_unknown_gene_raises_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            with pytest.raises(InvalidUsage) as exc_info:
                get_gtex("V8", "Liver", "FAKEGENE999")
            assert exc_info.value.status_code == 410

    def test_ensg_id_with_version_suffix_resolves(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        """Passing a versioned ENSG ID (e.g. ENSG00000069275.12) resolves correctly."""
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            ensg = GENES["NUCKS1"]["ensg_id"]  # "ENSG00000069275.12"
            df = get_gtex("V8", "Liver", ensg)
            assert df.shape[0] > 0
            assert "error" not in df.columns

    def test_tissue_spaces_normalized(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        """'Whole Blood' (with space) is accepted and normalized to 'Whole_Blood'."""
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            df = get_gtex("V8", "Whole Blood", "NUCKS1")
            assert df.shape[0] > 0

    def test_version_case_insensitive(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex

            df = get_gtex("v8", "Liver", "NUCKS1")
            assert df.shape[0] > 0


# ---------------------------------------------------------------------------
# get_variants_by_region input validation
# ---------------------------------------------------------------------------


class TestGetVariantsByRegionValidation:
    def test_start_greater_than_end_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            with pytest.raises(ValueError, match="[Ss]tart"):
                get_variants_by_region(200, 100, "1", "V8")

    def test_negative_start_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            with pytest.raises(ValueError, match="[Ss]tart"):
                get_variants_by_region(-1, 100, "1", "V8")

    def test_invalid_chromosome_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            with pytest.raises(ValueError):
                get_variants_by_region(100, 200, "99", "V8")

    def test_invalid_version_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            with pytest.raises(ValueError):
                get_variants_by_region(100, 200, "1", "V7")

    def test_chr_prefix_stripped(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            nucks1 = GENES["NUCKS1"]
            df = get_variants_by_region(nucks1["start"], nucks1["end"], "chr1", "V8")
            assert len(df) > 0

    def test_x_chromosome_normalized_to_23(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        """'X' is normalized to '23'; no error raised even when no variants exist there."""
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            df = get_variants_by_region(1, 10_000_000, "X", "V8")
            assert isinstance(df, pd.DataFrame)

    def test_returns_required_columns(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.variants import get_variants_by_region

            nucks1 = GENES["NUCKS1"]
            df = get_variants_by_region(nucks1["start"], nucks1["end"], "1", "V10")
            assert {"variant_id", "rs_id", "chr", "pos", "ref", "alt"}.issubset(df.columns)


# ---------------------------------------------------------------------------
# get_gtex_data — SNP list merging
# ---------------------------------------------------------------------------


class TestGetGTExData:
    def test_variant_id_b38_snp_list(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_data

            snp_list = [v["variant_id"] for v in fake_gtex_db._gene_variants["NUCKS1"][:10]]
            df = get_gtex_data("V8", "Liver", "NUCKS1", snp_list)
            non_null = df.dropna(subset=["pval"])
            assert len(non_null) > 0

    def test_rsid_and_b38_mix_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_data

            with pytest.raises(InvalidUsage, match="[Mm]ix"):
                get_gtex_data("V8", "Liver", "NUCKS1", ["rs1000000", "1_205712820_A_T_b38"])

    def test_v7_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_data

            with pytest.raises(InvalidUsage, match="V7"):
                get_gtex_data("V7", "Liver", "NUCKS1", ["rs1000000"])

    def test_unsupported_snp_format_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_data

            with pytest.raises(InvalidUsage):
                get_gtex_data("V8", "Liver", "NUCKS1", ["notavariant"])

    def test_no_matching_snps_returns_dataframe(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        """Non-matching SNPs produce a left-join DataFrame with NaN pvals, not an error."""
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_data

            df = get_gtex_data("V8", "Liver", "NUCKS1", ["rs9999999999", "rs8888888888"])
            assert isinstance(df, pd.DataFrame)


# ---------------------------------------------------------------------------
# get_gtex_data_pvalues
# ---------------------------------------------------------------------------


class TestGetGTExDataPvalues:
    def test_rsid_format_returns_pvals(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex, get_gtex_data_pvalues

            eqtl_df = get_gtex("V8", "Liver", "NUCKS1")
            snp_list = list(eqtl_df["rs_id"].dropna()[:5])
            pvals = get_gtex_data_pvalues(eqtl_df, snp_list)
            assert len(pvals) == len(snp_list)

    def test_variant_id_format_returns_pvals(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex, get_gtex_data_pvalues

            eqtl_df = get_gtex("V8", "Liver", "NUCKS1")
            snp_list = list(eqtl_df["variant_id"][:5])
            pvals = get_gtex_data_pvalues(eqtl_df, snp_list)
            assert len(pvals) == len(snp_list)

    def test_unsupported_format_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex, get_gtex_data_pvalues

            eqtl_df = get_gtex("V8", "Liver", "NUCKS1")
            with pytest.raises(InvalidUsage):
                get_gtex_data_pvalues(eqtl_df, ["notavariant"])


# ---------------------------------------------------------------------------
# get_gtex_snp_matches
# ---------------------------------------------------------------------------


class TestGetGTExSnpMatches:
    def test_matching_snps_counted(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_snp_matches

            snp_list = [v["variant_id"] for v in fake_gtex_db._gene_variants["NUCKS1"]]
            count = get_gtex_snp_matches(snp_list, "1:205700000-205760000", "hg38", "V8")
            assert count > 0

    def test_non_matching_snps_return_zero(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_snp_matches

            count = get_gtex_snp_matches(
                ["1_1_A_T_b38", "1_2_G_C_b38"],  # positions not in fake DB
                "1:205700000-205760000",
                "hg38",
                "V8",
            )
            assert count == 0

    def test_hg19_build_raises(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        with flask_app.app_context():
            from app.utils.gtex import get_gtex_snp_matches

            with pytest.raises(InvalidUsage):
                get_gtex_snp_matches([], "1:205700000-205760000", "hg19", "V8")


# ---------------------------------------------------------------------------
# Flask route error responses
# ---------------------------------------------------------------------------


class TestGTExRoutes:
    def test_tissues_list_v7_returns_error(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V7/tissues_list")
        assert resp.status_code != 200

    def test_tissues_list_invalid_version_returns_error(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V99/tissues_list")
        assert resp.status_code != 200

    def test_gtex_route_returns_records(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V8/Liver/NUCKS1")
        assert resp.status_code == 200
        data = resp.get_json()
        assert isinstance(data, list) and len(data) > 0
        assert "pval" in data[0] and "variant_id" in data[0]

    def test_gtex_route_unknown_tissue_returns_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V8/Kidney_Cortex/NUCKS1")
        assert resp.status_code == 410

    def test_gtex_route_unknown_gene_returns_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V8/Liver/FAKEGENE999")
        assert resp.status_code == 410

    def test_gtex_variant_route_by_rs_id(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        rs_id = fake_gtex_db._gene_variants["NUCKS1"][0]["rs_id"]
        resp = flask_app.test_client().get(f"/gtex/V8/Liver/NUCKS1/{rs_id}")
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data) > 0
        assert data[0]["rs_id"] == rs_id

    def test_gtex_variant_route_by_variant_id(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        variant_id = fake_gtex_db._gene_variants["NUCKS1"][0]["variant_id"]
        resp = flask_app.test_client().get(f"/gtex/V8/Liver/NUCKS1/{variant_id}")
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data) > 0
        assert data[0]["variant_id"] == variant_id

    def test_gtex_variant_route_bad_format_returns_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V8/Liver/NUCKS1/notavariant")
        assert resp.status_code == 410

    def test_gtex_variant_route_unknown_rs_id_returns_410(self, flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
        resp = flask_app.test_client().get("/gtex/V8/Liver/NUCKS1/rs9999999999")
        assert resp.status_code == 410
