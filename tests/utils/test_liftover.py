import pandas as pd

from app.utils.liftover import run_liftover

test_file_hg19 = "/code/data/sample_datasets/MI_GWAS_2019_1_205500-206000kbp.tsv"
test_file_hg38 = "/code/data/sample_datasets/MI_GWAS_2019_1_205500-206000kbp_hg38.tsv"


def test_liftover_19_to_38():

    hg19_df = pd.read_csv(test_file_hg19, sep="\t")

    # rename relevant cols to app standard
    hg19_df.rename({"#CHROM": "CHROM", "BP": "POS"}, axis=1, inplace=True)

    lifted_over_df, _ = run_liftover(hg19_df, "hg19")

    assert all(lifted_over_df.columns == hg19_df.columns)
    assert len(lifted_over_df) == len(hg19_df)
    assert not all(lifted_over_df.iloc[:, 1] == hg19_df.iloc[:, 1])


def test_liftover_38_to_19():

    hg38_df = pd.read_csv(test_file_hg38, sep="\t")
    # rename relevant cols to app standard
    hg38_df.rename({"chr": "CHROM", "variant_pos": "POS"}, axis=1, inplace=True)

    lifted_over_df, dropped = run_liftover(hg38_df, "hg19")

    # we're expecting dropped rows
    assert len(dropped) > 0

    assert all(lifted_over_df.columns == hg38_df.columns)

    # assert that dropped rows plus lifted rows equal original rows
    assert len(lifted_over_df) + len(dropped) == len(hg38_df)

    # assert that the missing indexes are what we'd expect given the original input
