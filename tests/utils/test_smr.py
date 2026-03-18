import pytest
from flask import Flask


def test_smr(flask_app: Flask):
    with flask_app.app_context():
        from app.utils.smr import query_smr

        # contrived SNPs that are definitely in the smr dataset
        snps = [
            'chr1_751343_A_T',
            'chr1_752566_G_A',
            'chr1_752721_A_G',
            'chr1_752894_T_C',
            'chr1_753405_C_A',
            'chr1_753425_T_C',
            'chr1_753474_C_G',
            'chr1_753541_A_G',
            'chr1_754182_A_G',
            'chr1_754192_A_G'
        ]
        # 10 SNPs

        t = query_smr(
            1,
            snps,
            "LBC_BSGS_meta",
            assembly="hg38"
        )

        assert t is not None

        result = query_smr(
            1,
            snps,
            "Brain-mMeta",
            assembly="hg38"
        )

        assert result is not None
        assert len(result) == 10


@pytest.mark.skip("For dev purposes")
def test_smr_raw_query(flask_app: Flask):
    with flask_app.app_context():
        from app.utils.smr import run_smr_query

        result = run_smr_query("data/smr_mqtl/Brain-mMeta/Brain-mMeta", 1, 5.0e-8, 1, 1_000_000)

        assert len(result) > 0