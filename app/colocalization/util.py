import pandas as pd

from app.colocalization.payload import SessionPayload
from app.utils.errors import ServerError


def get_std_snp_list(payload: SessionPayload, gwas_data: pd.DataFrame) -> pd.Series:
    """
    Return standardized list of SNPs with format CHR_POS_REF_ALT_build.

    gwas_data needs to be subsetted in advance.
    """
    std_snp_list = []
    buildstr = "b37"
    if payload.get_coordinate() == "hg38":
        buildstr = "b38"

    std_snp_list = pd.Series(
        [
            f"{str(row['CHROM']).replace('23', 'X')}_{str(row['POS'])}_{str(row['REF'])}_{str(row['ALT'])}_{buildstr}"
            for _, row in gwas_data.iterrows()
        ]
    )
    # Sanity check
    try:
        assert len(std_snp_list) == len(gwas_data)
        assert len(std_snp_list) == len(payload.gwas_indices_kept)
    except AssertionError:
        raise ServerError("GWAS data and indices are not in sync")

    return std_snp_list
