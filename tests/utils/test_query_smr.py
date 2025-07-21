from unittest.mock import patch, Mock

import pandas as pd

from app.utils.smr import query_smr

mock_result = pd.DataFrame(
    [
        [
            "chr1:982513",
            1,
            982513,
            "T",
            "C",
            0.074442,
            "cg24669183",
            1,
            534242,
            pd.NA,
            "N",
            0.127188,
            0.058496,
            0.029684,
        ],
        [
            "chr1:982513",
            1,
            982513,
            "T",
            "C",
            0.074442,
            "cg12726839",
            1,
            845311,
            pd.NA,
            "N",
            0.180720,
            0.058765,
            0.002103,
        ],
        [
            "chr1:982513",
            1,
            982513,
            "A",
            "T",
            0.074442,
            "cg12726839",
            1,
            845311,
            pd.NA,
            "N",
            0.180720,
            0.058765,
            0.002103,
        ],
    ],
    columns=[
        "SNP",
        "Chr",
        "BP",
        "A1",
        "A2",
        "Freq",
        "Probe",
        "Probe_Chr",
        "Probe_bp",
        "Gene",
        "Orientation",
        "b",
        "SE",
        "p",
    ],
)


@patch("app.utils.smr.run_smr_query", return_value=mock_result)
def test_query_smr(mock: Mock):
    """Test the query function with mock data, since smr files are not committed to source,
    ensuring that filtering and query construction functions are correst.
    """
    chr = 1
    snps = ["chr1_982513_T_C"]
    dataset = "EUR"
    thresh = 1

    res = query_smr(chr, snps, dataset, thresh)

    assert len(res) == 2

    assert len(res["full_snp"].isin(snps)) == 2

    mock.assert_called_once_with(
        query_path="/code/data/smr_mqtl/EUR/EUR_chr1",
        thresh=1,
        chr=1,
        start=982,
        end=983,
    )
