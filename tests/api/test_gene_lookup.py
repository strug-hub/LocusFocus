import os
import pytest
import pandas as pd

from app.utils.apis.gtex import get_genes

try:
    GENE_DF = pd.read_csv(
        os.path.join("data", "collapsed_gencode_v26_hg38.gz"),
        compression="gzip",
        sep="\t",
        encoding="utf-8",
    )
except FileNotFoundError:
    GENE_DF = None

@pytest.mark.skipif(GENE_DF is None, reason="data/collapsed_gencode_v26_hg38.gz not found; cannot test gencode lookup parity")
def test_gencode_id_lookup_parity():
    """
    Test that gencode ID lookups using the API and the local dataframe are the same.
    """
    gene_symbols = [
        "CDK18",
        "ELK4",
        "MFSD4",
        "NUCKS1",
        "PM20D1",
        "RAB7L1"
    ]

    # Get the genes from the local dataframe
    # (columns are ENSG_name, name)
    local_gene_results = GENE_DF[GENE_DF["name"].isin(gene_symbols)]

    # Get the genes from the API
    api_gene_results = get_genes("hg38", gene_symbols)

    assert len(api_gene_results.data) == len(local_gene_results)

    # Order is expected to be different but contents are the same
    assert set(local_gene_results["name"]) == set(map(lambda x: x.gene_symbol, api_gene_results.data))

    # Check that GENCODE IDs are the same
    sorted_local_gene_results = local_gene_results.sort_values(by="name")
    api_gene_results.data.sort(key=lambda x: x.gene_symbol)

    for i in range(len(sorted_local_gene_results)):
        local_gencode = sorted_local_gene_results["ENSG_name"].iloc[i]
        api_gencode = api_gene_results.data[i].gencode_id
        assert local_gencode == api_gencode
