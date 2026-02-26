import os
import re
import subprocess
from tempfile import NamedTemporaryFile
from typing import List, Literal, TypedDict

import pandas as pd

from app.utils.liftover import run_liftover, LiftoverError

curr_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(os.path.dirname(curr_dir)), "data", "smr_mqtl")


class SMRDataset(TypedDict):
    assembly: Literal["hg19", "hg38"]
    base_filename: str
    by_chr: bool
    description: str


smr_datasets: dict[str, SMRDataset] = {
    "Brain-mMeta": {
        "assembly": "hg38",
        "by_chr": False,
        "base_filename": "Brain-mMeta",
        "description": "(estimated effective n = 1160) Qi et al. Brain-mMeta mQTL summary data",
    },
    "EAS": {
        "assembly": "hg38",
        "by_chr": True,
        "base_filename": "EAS",
        "description": "(n=2,099) mQTL summary data from a meta-analysis of samples of East Asian ancestry.",
    },
    "EUR": {
        "assembly": "hg38",
        "by_chr": True,
        "base_filename": "EUR",
        "description": "(n=3,701) mQTL summary data from a meta-analysis of samples of European ancestry.",
    },
    "Hannon et al. Blood dataset1": {
        "assembly": "hg19",
        "by_chr": False,
        "base_filename": "Aberdeen_Blood",
        "description": "(n=639) Blood dataset 1, Hannon et al. 2016 Genome Biology",
    },
    "Hannon et al. Blood dataset2": {
        "assembly": "hg19",
        "by_chr": False,
        "base_filename": "UCL_Blood",
        "description": "(n=665) Blood dataset 2, Hannon et al. 2016 Genome Biology",
    },
    "Hannon et al. FetalBrain": {
        "assembly": "hg19",
        "by_chr": False,
        "base_filename": "FB_Brain",
        "description": "(n=166) Fetal brain mQTL data (Hannon et al. 2015 Nat Neurosci)",
    },
    "LBC_BSGS_meta": {
        "assembly": "hg19",
        "by_chr": True,
        "base_filename": "bl_mqtl",
        "description": "(n=1,980) McRae et al. mQTL summary data.",
    },
    "LBC_BSGS_meta_lite": {
        "assembly": "hg19",
        "by_chr": True,
        "base_filename": "bl_mqtl_lite",
        "description": "(n=1,980) McRae et al. mQTL summary data, only SNPs with P < 1e-5 are included",
    },
    "US_mQTLS_SMR_format": {
        "assembly": "hg19",
        "by_chr": False,
        "base_filename": "US_Blood",
        "description": "(n=1,175) Whole blood mQTL data set used in Hannon et al.",
    },
}


def run_smr_query(
    query_path: str, chr: int, thresh: float, start: int, end: int
) -> pd.DataFrame | None:
    """Query the data, save to a temporary file, and return as a dataframe.

    :param query_path: The path to the data to be queried
    :type query_path: str
    :param chr: The chromosome to query
    :type chr: int
    :param thresh: The p-value threshold
    :type thresh: float
    :param start: The start bp to query
    :type start: int
    :param end: Then end bp to query
    :type end: int
    :return: The returned SNPs and pvalues, or None if none were found.
    :rtype: pd.DataFrame | None
    """
    with NamedTemporaryFile("w") as f:
        query = [
            "smr",
            "--beqtl-summary",
            query_path,
            "--query",
            str(thresh),
            "--snp-chr",
            str(chr),
            "--from-snp-kb",
            str(start),
            "--to-snp-kb",
            str(end),
            "--out",
            f.name,
        ]

        subprocess.run(query, check=True)

        try:
            df = pd.read_csv(f"{f.name}.txt", sep="\t")
            return df
        except FileNotFoundError:
            return None


def query_smr(
    chr: int,
    snps: List[str],
    dataset: str,
    thresh: float = 5.0e-8,
    assembly: Literal["hg19", "hg38"] = "hg38",
) -> pd.DataFrame | None:
    """Query mqtl data in smr format

    :param chr: The chromosome to query
    :type chr: int
    :param snps: A list of SNPS in format chr{chr}_{bp}_ref_alt
    :type snps: List[str]
    :param dataset: The dataset to query
    :type dataset: str
    :param thresh: The p-value threshold, defaults to 5e-8
    :type thresh: float
    :param assembly: The genome assembly to use, defaults to "hg38"
    :type assembly: Literal["hg19", "hg38"]
    :raises FileNotFoundError: If the dataset does not exist
    :raises ValueError: If the requested assembly does not match the dataset assembly
    :return: The SNPs and pvalues as a dataframe with the following columns:
    'SNP', 'Chr', 'BP', 'A1', 'A2', 'Freq', 'Probe', 'Probe_Chr',
    'Probe_bp', 'Gene', 'Orientation', 'b', 'SE', 'p', 'full_snp'
    Note that 'full_snp' is a combined column that takes the same format as those in ``snps``
    May be None if SMR failed to provide an output file (ie. no data for query).
    :rtype: pd.DataFrame | None
    """
    if dataset not in smr_datasets.keys():
        raise FileNotFoundError(f"Dataset {dataset} does not exist!")

    try:
        assert all(snp[:3] == "chr" for snp in snps)
    except AssertionError as e:
        raise ValueError("Some SNPs provided are in the wrong format. Please ensure that 'chr{chr}_{bp}_ref_alt' format is used") from e

    needs_liftover = False

    # Lift up smr dataset if it's hg19
    if smr_datasets[dataset]["assembly"] != assembly:
        if assembly == "hg19" and smr_datasets[dataset]["assembly"] == "hg38":
            raise ValueError(
                f"Dataset {dataset} uses {smr_datasets[dataset]['assembly']} but {assembly} was requested! Provide hg38 snps or use liftover on the provided GWAS data."
            )
        else:
            # LiftOver from hg19 -> hg38
            needs_liftover = True

    dataset_dir = os.path.join(data_dir, dataset)
    base_filepath = os.path.join(dataset_dir, dataset)
    if smr_datasets[dataset]["by_chr"]:
        base_filepath = f"{base_filepath}_chr{chr}"

    regex = r"_(\d+)_"

    if needs_liftover:
        # need to lift "down" input SNPs (hg38) for query to work properly
        # this is a lossy conversion; we might lose all the snps
        snp_df = pd.DataFrame(data=snps, columns=["SNP"])
        snp_df[["CHROM", "POS", "REF", "ALT"]] = \
            snp_df["SNP"].str.extract(r"chr(?P<CHROM>\d+)_(?P<POS>\d+)_(?P<REF>[ACTG]+)_(?P<ALT>[ACTG]+)")
        snp_df["CHROM"] = pd.to_numeric(snp_df["CHROM"])
        snp_df["POS"] = pd.to_numeric(snp_df["POS"])

        lifted, lost_snps = run_liftover(
            snp_df,
            "hg38",
            chrom_col="CHROM",
            pos_col="POS"
        )

        if len(lost_snps) == len(snp_df):
            # we lost literally all the snps
            raise LiftoverError(
                "Could not perform SMR query; "
                "Unable to liftover hg38 snps prior to query on hg19 SMR dataset. "
                "Consider deselecting this dataset or using only hg38 SMR datasets."
            )

        lifted["SNP"] = \
            lifted["CHROM"].astype(str) \
            + "_" + lifted["POS"].astype(str) \
            + "_" + lifted["REF"] \
            + "_" + lifted["ALT"]
        if not lifted["SNP"].str.contains(r"$chr").any():
            lifted["SNP"] = "chr" + lifted["SNP"]

        snps = list(lifted["SNP"])

    snp_poses = [int(re.findall(regex, snp)[0]) for snp in snps]

    start = min(snp_poses) // 1000
    end = max(snp_poses) // 1000 + 1

    query_result = run_smr_query(
        query_path=base_filepath,
        chr=chr,
        thresh=thresh,
        start=start,
        end=end,
    )

    if query_result is None:
        return None

    query_result["full_snp"] = query_result.apply(
        lambda df: f"chr{str(df['Chr'])}"
        + "_"
        + str(df["BP"])
        + "_"
        + df["A1"]
        + "_"
        + df["A2"],
        axis=1,
    )

    filtered = query_result[query_result["full_snp"].isin(snps)].drop_duplicates(["SNP"])
    
    return filtered
