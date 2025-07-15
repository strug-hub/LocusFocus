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


def query_smr(chr: int, snps: List[str], dataset: str, thresh: float = 5.0e-8):
    """Query mqtl data in smr format

    :param chr: The chromosome to query
    :type chr: int
    :param snps: A list of SNPS in format chr{chr}_{bp}_ref_alt
    :type snps: List[str]
    :param snps_rs: A list of the *same* SNPs in rsid format
    :type snps_rs: List[str]
    :param dataset: The dataset to query
    :type dataset: str
    :param thresh: The p-value threshold, defaults to 5e-8
    :type thresh: float
    :raises FileNotFoundError: If the dataset does not exist
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

    with NamedTemporaryFile("w") as f:
        query = [
            "smr",
            "--beqtl-summary",
            base_filepath,
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

        res = pd.read_csv(f"{f.name}.txt", sep="\t")

    res.sort_values(by="BP", inplace=True)
    res["full_snp"] = res.apply(
        lambda x: f"chr{str(x['Chr'])}"
        + "_"
        + str(x["BP"])
        + "_"
        + x["A1"]
        + "_"
        + x["A2"],
        axis=1,
    )

    filtered = res[res["full_snp"].isin(snps)]

    return filtered
