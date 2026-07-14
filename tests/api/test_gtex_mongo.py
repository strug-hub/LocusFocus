import pytest
from flask.app import Flask

from tests.fake_gtex import FakeGTExDatabase, GENES, TISSUES


def test_get_gtex(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex returns a non-empty DataFrame for a known tissue/gene pair."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex

        results = get_gtex("V8", "Liver", "NUCKS1")
        assert results.shape[0] > 0


def test_get_gtex_v10(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex works for GTEx V10."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex

        results = get_gtex("V10", "Liver", "NUCKS1")
        assert results.shape[0] > 0


def test_get_gtex_columns(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex result includes all fields required by the colocalization pipeline."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex

        df = get_gtex("V8", "Liver", "NUCKS1")
        required_cols = {"variant_id", "rs_id", "pval", "beta", "se", "sample_maf", "ma_count"}
        assert required_cols.issubset(df.columns), (
            f"Missing columns: {required_cols - set(df.columns)}"
        )


def test_get_gtex_all_tissues(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex returns data for every seeded tissue."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex

        for tissue in TISSUES:
            df = get_gtex("V8", tissue, "NUCKS1")
            assert df.shape[0] > 0, f"No rows returned for tissue '{tissue}'"


def test_get_gtex_all_genes(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex returns data for every seeded gene."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex

        for gene in GENES:
            df = get_gtex("V8", "Liver", gene)
            assert df.shape[0] > 0, f"No rows returned for gene '{gene}'"


def test_get_gtex_data(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """get_gtex_data merges correctly against a known SNP list."""
    with flask_app.app_context():
        from app.utils.gtex import get_gtex_data

        # Build an rs-ID snp_list from the first few fake variants for NUCKS1
        nucks1_variants = fake_gtex_db._gene_variants["NUCKS1"]
        snp_list = [v["rs_id"] for v in nucks1_variants[:10]]

        df = get_gtex_data("V8", "Liver", "NUCKS1", snp_list)
        # At least some SNPs should have matched
        non_null = df.dropna(subset=["pval"])
        assert len(non_null) > 0


def test_list_tissues_route(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """The /gtex/<version>/tissues_list route returns the seeded tissues."""
    client = flask_app.test_client()
    for version in ("V8", "V10"):
        resp = client.get(f"/gtex/{version}/tissues_list")
        assert resp.status_code == 200
        data = resp.get_json()
        for tissue in TISSUES:
            assert tissue in data, f"'{tissue}' missing from {version} tissue list"


def test_variant_table_v8(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """variant_table for V8 has variant_pos field (not pos)."""
    with flask_app.app_context():
        from app.utils.variants import get_variants_by_region

        nucks1 = fake_gtex_db.genes["NUCKS1"]
        df = get_variants_by_region(
            nucks1["start"], nucks1["end"], str(nucks1["chrom"]), "V8"
        )
        assert len(df) > 0
        assert "rs_id" in df.columns
        assert "variant_id" in df.columns


def test_variant_table_v10(flask_app: Flask, fake_gtex_db: FakeGTExDatabase):
    """variant_table for V10 has pos field (not variant_pos)."""
    with flask_app.app_context():
        from app.utils.variants import get_variants_by_region

        nucks1 = fake_gtex_db.genes["NUCKS1"]
        df = get_variants_by_region(
            nucks1["start"], nucks1["end"], str(nucks1["chrom"]), "V10"
        )
        assert len(df) > 0
        assert "rs_id" in df.columns


def test_fake_db_is_reproducible():
    """Two FakeGTExDatabase instances with the same seed produce identical data."""
    db1 = FakeGTExDatabase(seed=99)
    db2 = FakeGTExDatabase(seed=99)

    for gene in GENES:
        v1 = db1._gene_variants[gene]
        v2 = db2._gene_variants[gene]
        assert v1 == v2, f"Variant mismatch for gene '{gene}'"

    # Different seed → different data
    db3 = FakeGTExDatabase(seed=7)
    assert db1._gene_variants["NUCKS1"] != db3._gene_variants["NUCKS1"]


def test_fake_db_variant_count():
    """Default configuration generates the expected number of variants per gene."""
    db = FakeGTExDatabase(n_variants_per_gene=20)
    for gene, variants in db._gene_variants.items():
        assert len(variants) == 20, (
            f"Expected 20 variants for '{gene}', got {len(variants)}"
        )
