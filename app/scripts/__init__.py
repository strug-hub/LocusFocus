"""
R scripts for LocusFocus. Functions for running R scripts and returning results.
"""

import os, subprocess
from typing import Union
import pandas as pd

# All R scripts are in app/scripts/, so we can import them here
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class ScriptError(Exception):
    """An error that occurred when running a script."""

    def __init__(self, stdout: str, stderr: str):
        self.stdout = stdout
        self.stderr = stderr
        super().__init__(
            f"Script returned non-zero exit code. Stdout: {stdout}, Stderr: {stderr}"
        )

    @property
    def message(self):
        return f"Script returned non-zero exit code.\nStdout: {self.stdout}\nStderr: {self.stderr}"


def simple_sum(
    p_values_filepath: os.PathLike,
    ldmatrix_filepath: os.PathLike,
    results_filepath: os.PathLike,
    set_based_p_threshold: Union[str, float] = 0.05,
) -> pd.DataFrame:
    """
    Run Simple Sum 2 colocalization using the getSimpleSumStats.R script.

    Does not perform any pre-processing or validation before running the script
    with the provided parameters.

    Return the DataFrame created by the R script as a result.
    Raise error if the script fails to run.
    """
    SCRIPT_PATH = os.path.join(BASE_DIR, "getSimpleSumStats.R")
    args = [
        "Rscript",
        SCRIPT_PATH,
        p_values_filepath,
        ldmatrix_filepath,
        "--set_based_p",
        str(set_based_p_threshold),
        "--outfilename",
        results_filepath,
    ]

    process_runner = subprocess.run(
        args=args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if process_runner.returncode != 0:
        raise ScriptError(process_runner.stdout, process_runner.stderr)

    simple_sum_df = pd.read_csv(results_filepath, sep="\t", encoding="utf-8")

    return simple_sum_df


def coloc2(
    coloc2_gwas_filepath: os.PathLike,
    coloc2_eqtl_filepath: os.PathLike,
    results_filepath: os.PathLike,
):
    """
    Run COLOC2 colocalization using the run_coloc2.R script.

    Does not perform any pre-processing or validation before running the script
    with the provided parameters.
    """
    SCRIPT_PATH = os.path.join(BASE_DIR, "run_coloc2.R")
    args = [
        "Rscript",
        SCRIPT_PATH,
        coloc2_gwas_filepath,
        coloc2_eqtl_filepath,
        "--outfilename",
        results_filepath,
    ]

    process_runner = subprocess.run(
        args=args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if process_runner.returncode != 0:
        raise ScriptError(process_runner.stdout, process_runner.stderr)

    coloc2_df = pd.read_csv(results_filepath, sep="\t", encoding="utf-8").fillna(-1)

    return coloc2_df


def set_based_test():
    """
    Run the set-based test using the getSimpleSumStats.R script.

    Does not perform any pre-processing or validation before running the script
    with the provided parameters.
    """
    SCRIPT_PATH = os.path.join(BASE_DIR, "getSimpleSumStats.R")
    raise NotImplementedError()
