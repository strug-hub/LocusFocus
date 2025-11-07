import os
import subprocess
import tempfile
from typing import List, Tuple

import pandas as pd


# inspired by https://github.com/broadinstitute/liftover/blob/main/service/server.py
def run_liftover(
    original_df: pd.DataFrame, 
    source_coords: str,
    chrom_col: str = "CHROM",
    pos_col: str = "POS",
) -> Tuple[pd.DataFrame, List[int]]:
    """Lift positions from one coordinate system to the other. This script assumes that columns in the
    dataframe have standardized names (specifically "CHROM" and "POS").

    :param original_df: The GWAS dataframe whose positions neeed to be lifted
    :type original_df: pd.DataFrame
    :param source_coords: Either `hg19` or `hg38`. The positions will be lifted into the other system. If
    `source_coords` is `hg19`, then the returned DataFrame will have positions in `hg38`.
    :type source_coords: str
    :param chrom_col: The column name for the chromosome column. Defaults to "CHROM".
    :type chrom_col: str
    :param pos_col: The column name for the position column. Defaults to "POS".
    :type pos_col: str
    :return: The new dataframe with updated "CHROM" and "POS" columns and a set of indexes referring to
    rows in the input DataFrame that could not be lifted over (if they exist).
    :rtype: Tuple[pd.DataFrame, set]
    """

    if source_coords not in ["hg19", "hg38"]:
        raise ValueError("Only source coordinates in hg19 or hg38 can be lifted over!")

    # drop any existing index and replace with an index column for tracking dropped rows
    original_copy = original_df.copy().reset_index(drop=True).reset_index()

    # truncate to necessary columns only
    transformed_df = original_copy.loc[:, [chrom_col, pos_col, "index"]]
    transformed_df = transformed_df.astype({chrom_col: str})

    # make sure chr is prefixed as expected by liftOver
    transformed_df.loc[:, chrom_col] = "chr" + transformed_df.iloc[:, 0].replace(
        "chr", "", regex=True
    )

    # chrom col needs to be prefixed with #
    transformed_df.rename({chrom_col: "#CHROM"}, inplace=True, axis=1)

    # liftover expects BED-style input; add necessary columns
    transformed_df.loc[:, "end"] = transformed_df.loc[:, pos_col]
    # BED format is zero-indexed for start pos
    transformed_df.loc[:, pos_col] = transformed_df.loc[:, pos_col] - 1

    # ensure proper column order
    transformed_df = transformed_df[["#CHROM", pos_col, "end", "index"]]

    chain_root = os.path.join("/usr", "local", "share", "liftOver")

    chain_file = (
        "hg19ToHg38.over.chain.gz"
        if source_coords == "hg19"
        else "hg38ToHg19.over.chain.gz"
    )

    chain_file_path = os.path.join(chain_root, chain_file)

    with tempfile.NamedTemporaryFile(
        suffix=".bed"
    ) as input_file, tempfile.NamedTemporaryFile(
        suffix=".bed"
    ) as output_file, tempfile.NamedTemporaryFile(
        suffix=".bed"
    ) as unmapped_output_file:

        transformed_df.to_csv(input_file.name, sep="\t", index=False)

        subprocess.run(
            [
                "liftOver",
                input_file.name,
                chain_file_path,
                output_file.name,
                unmapped_output_file.name,
            ],
            check=True,
        )

        result_df = pd.read_csv(
            output_file, sep="\t", names=list(transformed_df.columns)
        )

        result_df.loc[:, pos_col] = result_df.loc[:, "end"]

        joined = original_copy.merge(
            result_df, how="inner", suffixes=("", "_lifted"), on="index"
        )

        dropped = set(original_copy["index"]).difference(joined["index"])

        # Chrom can be affected by liftover, so we'll take it from the new DataFrame
        joined = joined.drop([chrom_col, pos_col, "index", "POS_lifted"], axis=1).rename(
            {"#CHROM": chrom_col, "end": pos_col}, axis=1
        )[original_df.columns]

    return joined, list(dropped)
