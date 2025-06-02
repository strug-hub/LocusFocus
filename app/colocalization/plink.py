"""
Functions for interacting with PLINK.
"""

import os
import shutil
import subprocess
from typing import List, Optional, Tuple

import pandas as pd
import numpy as np
from flask import current_app as app

from app.utils.errors import InvalidUsage, ServerError
from app.utils import write_list, x_to_23
from app.colocalization.constants import VALID_POPULATIONS


### FUNCTIONS ###


def find_plink_1kg_overlap(
    plink_filepath: str,
    snp_positions: List[int],
    snp_pvalues: Optional[List[float]] = None,
):
    """Return a Pandas dataframe containing SNP positions ("pos") and P values ("p") that were found
    in the provided 1000 Genomes dataset.

    Args:
        plink_filepath (str): Absolute path to a filename (no extension) for a given 1000Genomes dataset.
            Returned by `resolve_plink_filepath`.
        snp_positions (List[int]): List of SNP positions. Must be the same length as `snp_pvalues`.
        snp_pvalues (List[float] | None): List of SNP P values. Must be the same length as `snp_positions`. If none, then we ignore it.

    Returns:
        pd.DataFrame: A merged dataframe containing the overlap between the provided positions/pvalues, and the contents of the .bim
            file for the given 1000 Genomes population.
    """
    # Ensure lead snp is also present in 1KG; if not, choose next best lead SNP
    the1kg_snps_df = pd.read_csv(
        plink_filepath + ".bim", sep="\t", header=None
    )  # .iloc[:, 1]

    # Find lowest P-value position in snp_positions that is also in 1KG
    gwas_positions_df = pd.DataFrame({"pos": snp_positions, "p": snp_pvalues})
    # intersection
    positions_in_1kg_df = pd.merge(
        gwas_positions_df,
        the1kg_snps_df,
        how="inner",
        left_on="pos",
        right_on=the1kg_snps_df.columns[3],
    )
    return positions_in_1kg_df


def resolve_plink_filepath(build, pop, chrom):
    """
    Returns the file path of the binary plink file
    """
    if chrom == "X":
        chrom = 23
    try:
        chrom = int(chrom)
    except Exception:
        raise InvalidUsage(f"Invalid chromosome {str(chrom)}", status_code=410)
    if chrom not in np.arange(1, 24):
        raise InvalidUsage(f"Invalid chromosome {str(chrom)}", status_code=410)
    if pop not in VALID_POPULATIONS:
        raise InvalidUsage(
            f"{str(pop)} is not a recognized population", status_code=410
        )
    plink_filepath = ""
    if build.lower() in ["hg19", "grch37"]:
        if chrom == 23:
            plink_filepath = os.path.join(
                app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh37", pop, "chrX"
            )
        else:
            plink_filepath = os.path.join(
                app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh37", pop, f"chr{chrom}"
            )
    elif build.lower() in ["hg38", "grch38"]:
        if chrom == 23:
            plink_filepath = os.path.join(
                app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh38", "chrX"
            )
        else:
            plink_filepath = os.path.join(
                app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh38", f"chr{chrom}"
            )
    else:
        raise InvalidUsage(f"{str(build)} is not a recognized genome build")
    return plink_filepath


def get_plink_binary():
    """
    Return path to plink executable.
    """
    if os.name == "nt":
        return "./plink.exe"
    if shutil.which("plink") is not None:
        return "plink"
    if os.path.exists("/usr/local/bin/plink"):
        return "/usr/local/bin/plink"
    if os.path.exists("./plink"):
        return "./plink"
    raise ServerError("Could not find plink binary")


def plink_ld_pairwise(build, pop, chrom, snp_positions, snp_pvalues, outfilename):
    """
    Positions must be in hg19 coordinates.
    """
    # returns NaN for SNPs not in 1KG LD file; preserves order of input snp_positions
    plink_filepath = resolve_plink_filepath(build, pop, chrom)
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    write_list(snps, outfilename + "_snps.txt")

    positions_in_1kg_df = find_plink_1kg_overlap(
        plink_filepath, snp_positions, snp_pvalues
    )
    if len(positions_in_1kg_df) == 0:
        raise InvalidUsage(
            "No alternative lead SNP found in the 1000 Genomes. This error occurs when no provided SNPs could be found in the selected 1000 Genomes dataset. Please try a different population, or provide your own LD matrix.",
            status_code=410,
        )
    new_lead_snp_row = positions_in_1kg_df[
        positions_in_1kg_df["p"] == positions_in_1kg_df["p"].min()
    ]
    if len(new_lead_snp_row) > 1:
        app.logger.warning(
            f"Dataset has multiple lead SNPs: {new_lead_snp_row.to_json()}, taking first one..."
        )
        new_lead_snp_row = new_lead_snp_row.iloc[0]
    new_lead_snp_position = int(new_lead_snp_row["pos"])
    lead_snp = f"chr{str(int(chrom))}:{str(int(new_lead_snp_position))}"

    # plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')

    plink_binary = get_plink_binary()

    plink_args = [
        plink_binary,
        "--bfile",
        plink_filepath,
        "--chr",
        str(chrom),
        "--extract",
        outfilename + "_snps.txt",
        "--from-bp",
        str(positions_in_1kg_df["pos"].min()),
        "--to-bp",
        str(positions_in_1kg_df["pos"].max()),
        "--ld-snp",
        lead_snp,
        "--r2",
        "--ld-window-r2",
        "0",
        "--ld-window",
        "999999",
        "--ld-window-kb",
        "200000",
        "--make-bed",
        "--threads",
        "1",
        "--out",
        outfilename,
    ]

    if build.lower() in ["hg38", "grch38"]:
        if str(chrom).lower() in ["x", "23"]:
            # special case, use females only
            pop_filename = f"{pop}_female.txt"
        else:
            pop_filename = f"{pop}.txt"
        popfile = os.path.join(
            app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh38", pop_filename
        )
        plink_args.extend(["--keep", popfile])

    elif build.lower() not in ["hg19", "grch37"]:
        raise InvalidUsage(f"{str(build)} is not a recognized genome build")

    plinkrun = subprocess.run(
        plink_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    if plinkrun.returncode != 0:
        raise InvalidUsage(plinkrun.stdout.decode("utf-8"), status_code=410)
    ld_results = pd.read_csv(outfilename + ".ld", delim_whitespace=True)
    available_r2_positions = ld_results[["BP_B", "R2"]]
    pos_df = pd.DataFrame({"pos": snp_positions})
    merged_df = pd.merge(
        pos_df,
        available_r2_positions,
        how="left",
        left_on="pos",
        right_on="BP_B",
        sort=False,
    )[["pos", "R2"]]
    merged_df.fillna(-1, inplace=True)
    return merged_df, new_lead_snp_position


def plink_ldmat(
    build, pop, chrom, snp_positions, outfilename, region=None
) -> Tuple[pd.DataFrame, np.matrix]:
    """
    Generate an LD matrix using PLINK, using the provided population `pop` and the provided region information (`chrom`, `snp_positions`).
    If `region` is specified (format: (chrom, start, end)), then start and end will be used for region.

    Return a tuple containing:
    - pd.DataFrame of the generated .bim file (the SNPs used in the PLINK LD calculation).
      https://www.cog-genomics.org/plink/1.9/formats#bim
    - np.matrix representing the generated LD matrix itself
    """
    plink_filepath = resolve_plink_filepath(build, pop, chrom)
    # make snps file to extract:
    snps = [f"chr{str(int(chrom))}:{str(int(position))}" for position in snp_positions]
    write_list(snps, outfilename + "_snps.txt")

    # plink_path = subprocess.run(args=["which","plink"], stdout=subprocess.PIPE, universal_newlines=True).stdout.replace('\n','')
    if region is not None:
        from_bp = str(region[1])
        to_bp = str(region[2])
    else:
        from_bp = str(min(snp_positions))
        to_bp = str(max(snp_positions))

    plink_binary = get_plink_binary()

    plink_args = [
        plink_binary,
        "--bfile",
        plink_filepath,
        "--chr",
        str(chrom),
        "--extract",
        outfilename + "_snps.txt",
        "--from-bp",
        from_bp,
        "--to-bp",
        to_bp,
        "--r2",
        "square",
        "--make-bed",
        "--threads",
        "1",
        "--out",
        outfilename,
    ]

    if build.lower() in ["hg38", "grch38"]:
        if str(chrom).lower() in ["x", "23"]:
            # special case, females only
            pop_filename = f"{pop}_female.txt"
        else:
            pop_filename = f"{pop}.txt"
        popfile = os.path.join(
            app.config["LF_DATA_FOLDER"], "1000Genomes_GRCh38", pop_filename
        )
        plink_args.extend(["--keep", popfile])

    elif build.lower() not in ["hg19", "grch37"]:
        raise InvalidUsage(f"{str(build)} is not a recognized genome build")

    plinkrun = subprocess.run(
        args=plink_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    if plinkrun.returncode != 0:
        overlap = find_plink_1kg_overlap(plink_filepath, snp_positions, None)
        if len(overlap) == 0:
            raise InvalidUsage(
                f"No overlap found between provided SNPs and the selected 1000 Genomes dataset. Please select a different 1000 Genomes population, or provide your own LD matrix.\n\nPLINK error output:\n\n{plinkrun.stdout.decode('utf-8')}",
                status_code=410,
            )
        raise InvalidUsage(plinkrun.stdout.decode("utf-8"), status_code=410)
    # BIM file format, see https://www.cog-genomics.org/plink/1.9/formats#bim
    ld_snps_df = pd.read_csv(outfilename + ".bim", sep="\t", header=None)
    ld_snps_df.iloc[:, 0] = x_to_23(list(ld_snps_df.iloc[:, 0]))  # type: ignore
    ldmat = np.matrix(pd.read_csv(outfilename + ".ld", sep="\t", header=None))
    return ld_snps_df, ldmat
