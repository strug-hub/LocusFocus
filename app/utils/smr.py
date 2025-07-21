import os
import re
import subprocess
from tempfile import NamedTemporaryFile
from typing import List, TypedDict

import pandas as pd

curr_dir = os.path.dirname(__file__)
data_dir = os.path.join(os.path.dirname(os.path.dirname(curr_dir)), "data", "smr_mqtl")


class SMRDataset(TypedDict):
    base_filename: str
    by_chr: bool


datasets: dict[str, SMRDataset] = {
    "Brain-mMeta": {
        "by_chr": False,
        "base_filename": "Brain-mMeta",
    },
    "EAS": {
        "by_chr": True,
        "base_filename": "EAS",
    },
    "EUR": {
        "by_chr": True,
        "base_filename": "EUR",
    },
    "Hannon et al. Blood dataset1": {
        "by_chr": False,
        "base_filename": "Aberdeen_Blood",
    },
    "Hannon et al. Blood dataset2": {
        "by_chr": False,
        "base_filename": "UCL_Blood",
    },
    "Hannon et al. FetalBrain": {
        "by_chr": False,
        "base_filename": "FB_Brain",
    },
    # TODO: confirm
    "LBC_BSGS_meta": {
        "by_chr": True,
        "base_filename": "bl_mqtl",
    },
    "LBC_BSGS_meta_lite": {
        "by_chr": True,
        "base_filename": "bl_mqtl_lite",
    },
    "US_mQTLS_SMR_format": {
        "by_chr": False,
        "base_filename": "US_Blood",
    },
}


def run_smr_query(
    query_path: str, chr: int, thresh: float, start: int, end: int
) -> pd.DataFrame:
    """Query the data, save to a temporary file, and return as a dataframe

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
    :return: The returned SNPs and pvalues
    :rtype: pd.DataFrame
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

        return pd.read_csv(f"{f.name}.txt", sep="\t")


def query_smr(
    chr: int, snps: List[str], dataset: str, thresh: float = 5.0e-8
) -> pd.DataFrame:
    """Query mqtl data in smr format

    :param chr: The chromosome to query
    :type chr: int
    :param snps: A list of SNPS in format chr{chr}_{bp}_ref_alt
    :type snps: List[str]
    :param dataset: The dataset to query
    :type dataset: str
    :param thresh: The p-value threshold, defaults to 5e-8
    :type thresh: float
    :raises FileNotFoundError: If the dataset does not exist
    :return: The SNPs and pvalues as a dataframe with the following columns:
    'SNP', 'Chr', 'BP', 'A1', 'A2', 'Freq', 'Probe', 'Probe_Chr',
    'Probe_bp', 'Gene', 'Orientation', 'b', 'SE', 'p', 'full_snp'
    Note that 'full_snp' is a combined column that takes the same format as those in ``snps``
    :rtype: pd.DataFrame
    """
    if dataset not in datasets.keys():
        raise FileNotFoundError(f"Dataset {dataset} does not exist!")
    dataset_dir = os.path.join(data_dir, dataset)
    base_filepath = os.path.join(dataset_dir, dataset)
    if datasets[dataset]["by_chr"]:
        base_filepath = f"{base_filepath}_chr{chr}"

    regex = r"_(\d+)_"

    snp_poses = [int(re.findall(regex, snp)[0]) for snp in snps]

    start = min(snp_poses) // 1000
    end = max(snp_poses) // 1000 + 1

    query_result = run_smr_query(
        query_path=base_filepath, chr=chr, thresh=thresh, start=start, end=end
    )

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

    filtered = query_result[query_result["full_snp"].isin(snps)]

    return filtered
