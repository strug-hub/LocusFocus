from typing import List

import pandas as pd


def validate_chromosome(
    chr: str | int, prefix: str | None = "chr", x_y_numeric: bool = False
) -> bool:
    """Validate that the chromosome is in the correct format

    :param chr: The chromosome representation, Examples: chr1, 1, chrX, chrome23
    :type chr: str | int
    :param prefix: The part of `chr` that is a string, defaults to "chr"
    :type prefix: str | None, optional
    :param x_y_numeric: Whether X and Y are represented as a number, defaults to False
    :type x_y_numeric: bool, optional
    :raises ValueError: If the validation fails
    :return: The validation status (True if passing, otherwise a ValueError is raised)
    :rtype: bool
    """
    chrs: List[str] = [str(c) for c in range(1, 23)]
    chr = str(chr)
    if x_y_numeric is True:
        chrs.append("23")
    else:
        chrs.extend(["X", "Y"])
    if prefix is not None:
        if not chr.startswith(prefix):
            raise ValueError(f"Chromosome must start with {prefix}.")
        chr = chr.replace(prefix, "")
    if chr not in chrs:
        raise ValueError(f"Chromosome must be one of {', '.join(chrs)}")

    return True


def adjust_snp_column(
        snps_df: pd.DataFrame,
        target_build: str,
        snp_col: str = "SNP",
        chrom_col: str = "CHROM",
        pos_col: str = "POS",
        ref_col: str = "REF",
        alt_col: str = "ALT",
        ignore_alleles: bool = False,
    ) -> pd.DataFrame:
    """
    Adjust the SNP column to be consistent with the other columns in the dataframe.
    rsIDs are ignored. Necessary to perform after LiftOver.

    :param snps_df: The dataframe containing the SNP column
    :type snps_df: pd.DataFrame
    :param target_build: The target build, either "hg19" or "hg38"
    :type target_build: str
    :param snp_col: The name of the SNP column, defaults to "SNP"
    :type snp_col: str, optional
    :param chrom_col: The name of the chromosome column, defaults to "CHROM"
    :type chrom_col: str, optional
    :param pos_col: The name of the position column, defaults to "POS"
    :type pos_col: str, optional
    :param ref_col: The name of the reference allele column, defaults to "REF"
    :type ref_col: str, optional
    :param alt_col: The name of the alternate allele column, defaults to "ALT"
    :type alt_col: str, optional
    :param ignore_alleles: Whether alleles should be ignored. 
        If true, alleles will be pulled from the SNP column instead of the ref/alt columns. 
        Defaults to False
    :type ignore_alleles: bool, optional
    :return: The dataframe with the SNP column adjusted.
    :rtype: pd.DataFrame
    """
    if target_build.lower() not in ["hg19", "hg38"]:
        raise ValueError("Target build must be either 'hg19' or 'hg38'")

    build_suffix = "b37" if target_build.lower() == "hg19" else "b38"

    # Mask out rows with rsIDs in the SNP column
    rsid_mask = snps_df[snp_col].str.contains("rs")
    if rsid_mask.all():
        # nothing to do, rsIDs are semi-stable 
        return snps_df

    if not ignore_alleles:
        snps_df.loc[~rsid_mask, snp_col] = snps_df.loc[~rsid_mask, chrom_col].astype(str) + "_" + snps_df.loc[~rsid_mask, pos_col].astype(str) + "_" + snps_df.loc[~rsid_mask, ref_col].astype(str) + "_" + snps_df.loc[~rsid_mask, alt_col].astype(str) + "_" + build_suffix
    else:
        # get alleles from snp column
        alleles = snps_df[snp_col].str.split("_", expand=True)
        alleles.columns = ["chrom", "pos", "ref", "alt", "build"]
        snps_df.loc[~rsid_mask, snp_col] = snps_df.loc[~rsid_mask, chrom_col].astype(str) + "_" + snps_df.loc[~rsid_mask, pos_col].astype(str) + "_" + alleles.loc[~rsid_mask, "ref"].astype(str) + "_" + alleles.loc[~rsid_mask, "alt"].astype(str) + "_" + build_suffix

    return snps_df
