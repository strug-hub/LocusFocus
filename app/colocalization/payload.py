from dataclasses import dataclass, field
from uuid import UUID
from typing import List, Dict, Optional, Tuple, Union
import os

import numpy as np
import pandas as pd
from werkzeug.datastructures import ImmutableMultiDict

from app.colocalization.constants import (
    ONE_SIDED_SS_WINDOW_SIZE,
    VALID_COORDINATES,
    VALID_POPULATIONS,
)
from app.utils.errors import InvalidUsage
from app.utils.gtex import get_gtex_snp_matches, gene_names
from app.utils import (
    get_session_filepath,
    parse_region_text,
)


class SessionFiles:
    """
    Files and filepaths for a given session.
    """

    def __init__(self, session_id: UUID):
        # Session files
        self.session_filepath = get_session_filepath(f"form_data-{session_id}.json")
        self.genes_session_filepath = get_session_filepath(
            f"genes_data-{session_id}.json"
        )
        self.SSPvalues_filepath = get_session_filepath(f"SSPvalues-{session_id}.json")
        self.coloc2_filepath = get_session_filepath(f"coloc2result-{session_id}.json")
        self.metadata_filepath = get_session_filepath(f"metadata-{session_id}.json")

        # Simple Sum files
        self.p_value_filepath = get_session_filepath(f"Pvalues-{session_id}.txt")
        self.ld_matrix_filepath = get_session_filepath(f"ldmat-{session_id}.txt")
        self.ld_mat_snps_filepath = get_session_filepath(f"ldmat_snps-{session_id}.txt")
        self.ld_mat_positions_filepath = get_session_filepath(
            f"ldmat_positions-{session_id}.txt"
        )
        self.simple_sum_results_filepath = get_session_filepath(
            f"SSPvalues-{session_id}.txt"
        )

        # COLOC2 files
        self.coloc2_gwas_filepath = get_session_filepath(
            f"coloc2gwas_df-{session_id}.txt"
        )
        self.coloc2_eqtl_filepath = get_session_filepath(
            f"coloc2eqtl_df-{session_id}.txt"
        )
        self.coloc2_results_filepath = get_session_filepath(
            f"coloc2result_df-{session_id}.txt"
        )

    def get_plot_template_paths(
        self, session_id: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Get the filepaths in format needed for the `plot.html` render template.

        Can be passed directly to render_template as kwargs.
        """
        paths = {
            "sessionfile": str(self.session_filepath),
            "genesfile": str(self.genes_session_filepath),
            "SSPvalues_file": str(self.SSPvalues_filepath),
            "coloc2_file": str(self.coloc2_filepath),
            "metadata_file": str(self.metadata_filepath),
        }

        adjusted_paths = {}

        # Adjust paths to be relative to the session data folder
        for key, value in paths.items():
            adjusted_paths[key] = f"session_data/{os.path.basename(value)}"

        if session_id is not None:
            adjusted_paths["sessionid"] = session_id

        return adjusted_paths


@dataclass
class SessionPayload:
    """
    Payload object for colocalization sessions.

    Contains metadata, form inputs and user data uploaded
    for this session.
    """

    # Request object
    request_form: dict
    uploaded_files: List[os.PathLike]  # Save files in advance
    # request_files: ImmutableMultiDict
    session_id: UUID

    # Form Inputs
    coordinate: Optional[str] = None
    coloc2: Optional[bool] = None
    study_type: Optional[str] = None
    num_cases: Optional[int] = None  # Used for study_type "cc" (case-control)
    ld_population: Optional[str] = None
    infer_variant: Optional[bool] = None
    gtex_tissues: Optional[List[str]] = None
    gtex_genes: Optional[List[str]] = None
    set_based_p: Optional[Union[str, float]] = None
    plot_locus: Optional[str] = None  # regionstr
    simple_sum_locus: Optional[str] = None
    lead_snp_name: Optional[str] = None

    # File data
    # GWAS data is user-uploaded, and we update gwas_indices_kept in each stage to "keep" or "discard" SNPs
    gwas_data: Optional[pd.DataFrame] = None  # Original GWAS data
    gwas_indices_kept: pd.Series = field(
        default_factory=pd.Series
    )  # Boolean Array of GWAS SNPs kept
    ld_matrix: Optional[np.matrix] = None
    secondary_datasets: Optional[Dict[str, dict]] = None
    file: SessionFiles = field(init=False)

    # LD Matrix data
    # DataFrame containing SNPs that were actually used in LD. See: https://www.cog-genomics.org/plink/1.9/formats#bim
    # has columns: "CHROM", "CHROM_POS", "POS", "ALT", "REF"
    ld_snps_bim_df: Optional[pd.DataFrame] = None

    # Other
    success: bool = False
    r2: List[float] = field(default_factory=list)

    # Simple Sum
    ss_locustext: Optional[str] = None
    ss_result_df: Optional[pd.DataFrame] = None

    # GTEx
    reported_gtex_data: dict = field(
        default_factory=dict
    )  # only used for reporting, use get_gtex_selection() instead
    gene: Optional[str] = None
    std_snp_list: pd.Series = field(default_factory=pd.Series)

    def __post_init__(self):
        # Runs after init, initializes SessionFiles object
        self.file = SessionFiles(self.session_id)

    @property
    def gwas_data_kept(self) -> pd.DataFrame:
        """
        Returns the GWAS data that was kept for Simple Sum.

        Shorthand for `payload.gwas_data.loc[payload.gwas_indices_kept]`
        """
        if self.gwas_data is None:
            raise Exception("GWAS data not loaded")
        return self.gwas_data.loc[self.gwas_indices_kept]

    def get_coordinate(self) -> str:
        """
        Get the form input for coordinate (aka. genome assembly, or 'build') for this session.

        Return either 'hg19' or 'hg38'.
        """
        if self.coordinate is None:
            if self.request_form.get("coordinate") not in VALID_COORDINATES:
                raise InvalidUsage(
                    f"Invalid coordinate: '{self.request_form.get('coordinate')}'"
                )

            self.coordinate = self.request_form.get("coordinate")
        return self.coordinate  # type: ignore

    def get_locus(self) -> str:
        """
        Gets the form input for plot locus (aka. 'regionstr' or 'regiontxt') for this session.
        """
        if self.plot_locus is None:
            self.plot_locus = self.request_form.get("locus", "1:205500000-206000000")
        return self.plot_locus

    def get_locus_tuple(self) -> Tuple[int, int, int]:
        """
        Gets the form input for plot locus but as a separated tuple: (chrom, start, end).
        """
        locus = self.get_locus()
        coordinate = self.get_coordinate()
        return parse_region_text(locus, coordinate)

    def get_infer_variant(self) -> bool:
        """
        Gets the form input for infer_variant for this session.
        True if the user would like to use dbSNP to get CHROM, POS, REF, ALT columns,
        False if they provide such columns themselves.
        """
        if self.infer_variant is None:
            self.infer_variant = bool(self.request_form.get("markerCheckbox"))
        return self.infer_variant

    def get_is_coloc2(self) -> bool:
        """
        Gets the form input for coloc2 for this session.
        True if the user wants to perform COLOC2 on this file.
        False otherwise.
        """
        if self.coloc2 is None:
            self.coloc2 = bool(self.request_form.get("coloc2check"))
        return self.coloc2

    def get_coloc2_study_type(self) -> str:
        """
        Get study type for COLOC2. Assumes that coloc2 is requested.
        """
        if self.study_type is None:
            study_type = self.request_form.get("studytype")
            if study_type not in ["quant", "cc"]:
                raise InvalidUsage(
                    f"Study type form value is invalid: {self.request_form.get('studytype')}"
                )
            self.study_type = study_type
        return self.study_type  # type: ignore

    def get_coloc2_case_control_cases(self) -> int:
        """
        Get number of case control cases specified for coloc2. Assumes coloc2 is requested, and
        selected study type is Case-Control.
        """
        if self.num_cases is None:
            try:
                num_cases = self.request_form.get("numcases", type=int)
            except ValueError:
                raise InvalidUsage(
                    "Number of cases entered must be an integer", status_code=410
                )
            self.num_cases = num_cases

        return self.num_cases  # type: ignore

    def get_ld_population(self) -> str:
        """
        Get selected LD population from user input (1000 Genomes dataset).
        If the user did not specify LD population, then use default "EUR".
        """
        if self.ld_population is None:
            pop = self.request_form.get("LD-populations", "EUR")
            if pop not in VALID_POPULATIONS:
                raise InvalidUsage(
                    f"Invalid population provided: '{pop}'. Population must be one of '{', '.join(VALID_POPULATIONS)}'"
                )
            self.ld_population = pop
        return self.ld_population

    def get_lead_snp_name(self) -> str:
        """
        Get the lead SNP name from user input, if specified.
        """

        if self.lead_snp_name is None:
            self.lead_snp_name = self.request_form.get("leadsnp", "")

        return self.lead_snp_name

    def get_lead_snp_index(self, only_kept=True) -> int:
        """
        Get the index of the lead SNP for the GWAS dataset.

        By default, only checks the kept SNPs in the GWAS dataset.

        If `only_kept` is False, then checks all SNPs in the GWAS dataset,
        including SNPs that were removed due to filtering steps (eg. bad format, duplicate SNPs, etc.).
        """
        lead_snp = self.get_lead_snp_name()
        if self.gwas_data is None:
            raise Exception("Cannot get lead SNP index when GWAS dataset is undefined.")

        gwas = self.gwas_data_kept if only_kept else self.gwas_data

        snps = [snp.split(";")[0] for snp in gwas.loc[:, "SNP"]]  # type: ignore
        if lead_snp == "":
            lead_snp = str(gwas.iloc[gwas["P"].argmin()]["SNP"]).split(";")[0]  # type: ignore
        if lead_snp not in snps:
            raise InvalidUsage("Lead SNP not found", status_code=410)
        return gwas["P"].argmin()  # type: ignore

    def get_ss_locus(self) -> str:
        """
        Get the Simple Sum region string from form input.
        If no Simple Sum region string is provided by user,
        then use default (200kbp region centered at lead SNP position).

        Prerequisite: need GWAS dataset.
        """
        if self.ss_locustext is not None:
            return self.ss_locustext

        SSlocustext = self.request_form.get("SSlocus", "")

        if SSlocustext != "":
            SSchrom, SS_start, SS_end = parse_region_text(
                SSlocustext, self.get_coordinate()
            )
        else:
            if self.gwas_data is None:
                raise Exception(
                    "Need GWAS dataset in order to find Lead SNP for SS region"
                )
            SSchrom, _, _ = self.get_locus_tuple()
            lead_snp_position_index = self.get_lead_snp_index()
            lead_snp_position = int(self.gwas_data.iloc[lead_snp_position_index, :]["POS"])  # type: ignore
            SS_start = max(int(lead_snp_position - ONE_SIDED_SS_WINDOW_SIZE), 0)
            SS_end = int(lead_snp_position + ONE_SIDED_SS_WINDOW_SIZE)
        SSlocustext = str(SSchrom) + ":" + str(SS_start) + "-" + str(SS_end)
        self.ss_locustext = SSlocustext

        return self.ss_locustext

    def get_ss_locus_tuple(self) -> Tuple[int, int, int]:
        """
        Get the Simple Sum region as a tuple: (chrom, start, end).

        Prerequisite: need GWAS dataset.
        """
        locus_text = self.get_ss_locus()
        chrom, startend = locus_text.split(":", 1)
        start, end = startend.split("-", 1)
        return int(chrom), int(start), int(end)

    def get_gtex_selection(self) -> Tuple[List[str], List[str]]:
        """
        Return a tuple containing list of selected tissues and list of selected genes.

        Format: (tissues, genes)
        """
        if self.gtex_genes is None:
            self.gtex_genes: List[str] = self.request_form.get("region-genes")  # type: ignore
        if self.gtex_tissues is None:
            self.gtex_tissues: List[str] = self.request_form.get("GTEx-tissues")  # type: ignore

        if len(self.gtex_tissues) > 0 and len(self.gtex_genes) == 0:
            raise InvalidUsage(
                "Please select one or more genes to complement your GTEx tissue(s) selection",
                status_code=410,
            )
        elif len(self.gtex_genes) > 0 and len(self.gtex_tissues) == 0:
            raise InvalidUsage(
                "Please select one or more tissues to complement your GTEx gene(s) selection",
                status_code=410,
            )

        return self.gtex_tissues, self.gtex_genes

    def get_gtex_version(self) -> str:
        """
        Get the version of GTEx needed for fetching from MongoDB.
        One of "V7" or "V8".
        """
        version = self.request_form.get("GTEx-version")
        if version is None:
            version = "V7"
        version = version.upper()
        if version not in ["V7", "V8"]:
            raise InvalidUsage(
                f"Invalid GTEx version: {version}. Must be one of 'V7' or 'V8'",
                status_code=410,
            )
        return version

    def get_p_value_threshold(self) -> Union[float, str]:
        """
        Get the user-selected P value threshold for set-based test.
        """

        if self.set_based_p is None:
            set_based_p = self.request_form["setbasedP"]  # type: ignore
            if set_based_p == "":
                set_based_p = "default"
            else:
                try:
                    set_based_p = float(set_based_p)
                    if set_based_p < 0 or set_based_p > 1:
                        raise InvalidUsage(
                            "Set-based p-value threshold given is not between 0 and 1"
                        )
                except:
                    raise InvalidUsage(
                        "Invalid value provided for the set-based p-value threshold. Value must be numeric between 0 and 1."
                    )

            self.set_based_p = set_based_p

        return self.set_based_p

    def is_ld_user_defined(self) -> bool:
        """
        Return True if the user has uploaded their own LD matrix, False otherwise.
        """
        return any(str(filepath).endswith("ld") for filepath in self.uploaded_files)

    def dump_session_data(self):
        """
        Create JSON dict of session data needed for form_data file.
        """
        # TODO: find a way to dump session data in a way that doesn't suck

        data = {}
        data["success"] = self.success
        data["sessionid"] = str(self.session_id)
        data["snps"] = list(self.gwas_data_kept["SNP"]) if self.gwas_data is not None else []
        data["inferVariant"] = self.get_infer_variant()
        data["pvalues"] = (
            list(self.gwas_data_kept["P"]) if self.gwas_data is not None else []
        )
        data["lead_snp"] = (
            self.gwas_data_kept["SNP"].iloc[self.get_lead_snp_index()]
            if self.gwas_data is not None
            else None
        )
        data["ld_values"] = self.r2
        data["positions"] = (
            list(self.gwas_data_kept["POS"]) if self.gwas_data is not None else []
        )
        data["chrom"], data["startbp"], data["endbp"] = self.get_locus_tuple()
        data["ld_populations"] = self.get_ld_population()
        data["gtex_tissues"], data["gtex_genes"] = self.get_gtex_selection()
        data["gene"] = gene_names(self.gene, self.get_coordinate())[0]
        data["gtex_version"] = self.get_gtex_version()
        data["coordinate"] = self.get_coordinate()
        data["set_based_p"] = self.set_based_p
        data["std_snp_list"] = list(self.std_snp_list)
        data["runcoloc2"] = self.coloc2 is True
        num_GTEx_matches = get_gtex_snp_matches(
            self.std_snp_list, self.get_locus(), self.get_coordinate()
        )
        data["snp_warning"] = num_GTEx_matches / len(self.std_snp_list) < 0.8
        data["thresh"] = 0.8
        data["numGTExMatches"] = num_GTEx_matches

        # secondary datasets
        if self.secondary_datasets is not None:
            data["secondary_dataset_titles"] = list(self.secondary_datasets.keys())
            for dataset_title, table in self.secondary_datasets.items():
                data[dataset_title] = table
        else:
            data["secondary_dataset_titles"] = []

        # gtex data
        if self.reported_gtex_data is not None:
            for tissue, table in self.reported_gtex_data.items():
                data[tissue] = table

        if self.coloc2:
            data["secondary_dataset_colnames"] = [
                "CHR",
                "POS",
                "SNPID",
                "PVAL",
                "BETA",
                "SE",
                "N",
                "A1",
                "A2",
                "MAF",
                "ProbeID",
            ]
        else:
            data["secondary_dataset_colnames"] = ["CHROM", "BP", "SNP", "P"]

        # simple sum
        _, ss_start, ss_end = self.get_ss_locus_tuple()
        data["SS_region"] = [ss_start, ss_end]
        data["num_SS_snps"] = len(self.std_snp_list.loc[self.gwas_indices_kept])
        if self.ss_result_df is None:
            data["first_stages"] = []
            data["first_stage_Pvalues"] = []
        else:
            data["first_stages"] = self.ss_result_df["first_stages"].tolist()
            data["first_stage_Pvalues"] = self.ss_result_df["first_stage_p"].tolist()

        return data
