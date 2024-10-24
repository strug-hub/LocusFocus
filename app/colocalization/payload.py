from dataclasses import dataclass, field
from uuid import uuid4, UUID
from typing import List, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from flask import Request

from app.colocalization.constants import (
    ONE_SIDED_SS_WINDOW_SIZE,
    VALID_COORDINATES,
    VALID_POPULATIONS,
)
from app.utils.errors import InvalidUsage, ServerError
from app.utils.gtex import get_gtex_snp_matches, gene_names
from app.utils import (
    download_file,
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

    def write_to_file(self, filename: str, data: str):
        # TODO: Implement in data-agnostic way.
        # Fail if file already exists and is not empty.
        # Write data to file. Writing method depends on data type.
        pass


@dataclass
class SessionPayload:
    """
    Payload object for colocalization sessions.

    Contains metadata, form inputs and user data uploaded
    for this session.
    """

    # Request object
    request: Request
    session_id: UUID = field(default_factory=uuid4)

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
    gwas_data: Optional[pd.DataFrame] = None
    ld_matrix: Optional[np.matrix] = None
    secondary_datasets: Optional[Dict[str, dict]] = None
    file: SessionFiles = field(init=False)

    # LD Matrix data
    # DataFrame containing SNPs that were actually used in LD. See: https://www.cog-genomics.org/plink/1.9/formats#bim
    # has columns: "CHROM", "CHROM_POS", "POS", "ALT", "REF"
    ld_snps_bim_df: Optional[pd.DataFrame] = None

    # Other
    success: bool = False
    gwas_indices_kept: List[bool] = field(default_factory=list)
    r2: List[float] = field(default_factory=list)

    # Simple Sum
    ss_locustext: Optional[str] = None
    ss_indices: List[int] = field(default_factory=list)
    ss_result_df: Optional[pd.DataFrame] = None

    # GTEx
    reported_gtex_data: dict = field(
        default_factory=dict
    )  # only used for reporting, use get_gtex_selection() instead
    gene: Optional[str] = None
    std_snp_list: List[str] = field(default_factory=list)

    def __post_init__(self):
        # Runs after init, initializes SessionFiles object
        self.file = SessionFiles(self.session_id)

    def get_coordinate(self) -> str:
        """
        Get the form input for coordinate (aka. genome assembly, or 'build') for this session.

        Return either 'hg19' or 'hg38'.
        """
        if self.coordinate is None:
            if self.request.form.get("coordinate") not in VALID_COORDINATES:
                raise InvalidUsage(
                    f"Invalid coordinate: '{self.request.form.get('coordinate')}'"
                )

            self.coordinate = self.request.form.get("coordinate")
        return self.coordinate  # type: ignore

    def get_locus(self) -> str:
        """
        Gets the form input for plot locus (aka. 'regionstr' or 'regiontxt') for this session.
        """
        if self.plot_locus is None:
            self.plot_locus = self.request.form.get("locus", "1:205500000-206000000")
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
            self.infer_variant = bool(self.request.form.get("markerCheckbox"))
        return self.infer_variant

    def get_is_coloc2(self) -> bool:
        """
        Gets the form input for coloc2 for this session.
        True if the user wants to perform COLOC2 on this file.
        False otherwise.
        """
        if self.coloc2 is None:
            self.coloc2 = bool(self.request.form.get("coloc2check"))
        return self.coloc2

    def get_coloc2_study_type(self) -> str:
        """
        Get study type for COLOC2. Assumes that coloc2 is requested.
        """
        if self.study_type is None:
            study_type = self.request.form.get("studytype")
            if study_type not in ["quant", "cc"]:
                raise InvalidUsage(
                    f"Study type form value is invalid: {self.request.form.get('studytype')}"
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
                num_cases = self.request.form.get("numcases", type=int)
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
            pop = self.request.form.get("LD-populations", "EUR")
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
            self.lead_snp_name = self.request.form.get("leadsnp", "")

        return self.lead_snp_name

    def get_lead_snp_index(self) -> int:
        """
        Get the index of the lead SNP for the GWAS dataset.
        """
        lead_snp = self.get_lead_snp_name()
        if self.gwas_data is None:
            raise Exception("Cannot get lead SNP index when GWAS dataset is undefined.")

        snps = [snp.split(";")[0] for snp in self.gwas_data.loc[:, "SNP"]]  # type: ignore
        if lead_snp == "":
            lead_snp = list(self.gwas_data.loc[self.gwas_data.loc[:, "P"] == self.gwas_data.loc[:, "P"].min()].loc[:, "SNP"])[0].split(";")[0]  # type: ignore
        if lead_snp not in snps:
            raise InvalidUsage("Lead SNP not found", status_code=410)
        return snps.index(lead_snp)

    def get_ss_locus(self) -> str:
        """
        Get the Simple Sum region string from form input.
        If no Simple Sum region string is provided by user,
        then use default (200kbp region centered at lead SNP position).

        Prerequisite: need GWAS dataset.
        """
        if self.ss_locustext is not None:
            return self.ss_locustext

        SSlocustext = self.request.form.get("SSlocus", "")

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
            self.gtex_genes = self.request.form.getlist("region-genes")
        if self.gtex_tissues is None:
            self.gtex_tissues = self.request.form.getlist("GTEx-tissues")

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
        if self.get_coordinate() == "hg38":
            return "V8"
        else:
            return "V7"

    def get_p_value_threshold(self) -> Union[float, str]:
        """
        Get the user-selected P value threshold for set-based test.
        """

        if self.set_based_p is None:
            set_based_p = self.request.form["setbasedP"]  # type: ignore
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
        return download_file(self.request, ["ld"]) is not None

    def get_current_lead_snp_index(self) -> int:
        """
        Determine and return the index of the lead SNP for the gwas dataset (ie. the SNP with the lowest P-value).

        Prerequisite: gwas_data is loaded.
        """
        try:
            assert self.gwas_data is not None
        except AssertionError:
            raise ServerError(message="Attempted to get lead SNP index before GWAS data was loaded")
        lead_snp = self.get_lead_snp_name()
        snp_list = list(self.gwas_data.loc[:,"SNP"])
        # cleaning up the SNP names a bit
        snp_list = [asnp.split(';')[0] for asnp in snp_list] # type: ignore
        if lead_snp=='': 
            lead_snp: str = list(gwas_data.loc[ gwas_data.loc[:,"P"] == gwas_data.loc[:,"P"].min() ].loc[:,"SNP"])[0].split(';')[0] # type: ignore
        if lead_snp not in snp_list:
            raise InvalidUsage(f"Lead SNP '{lead_snp}' not found in dataset", status_code=410)
        lead_snp_position_index = snp_list.index(lead_snp)

        return lead_snp_position_index

    def dump_session_data(self):
        """
        Create JSON dict of session data needed for form_data file.
        """
        # TODO: find a way to dump session data in a way that doesn't suck

        data = {}
        data["success"] = self.success
        data["sessionid"] = str(self.session_id)
        data["snps"] = list(self.gwas_data["SNP"]) if self.gwas_data is not None else []
        data["infervariant"] = self.get_infer_variant()
        data["pvalues"] = (
            list(self.gwas_data["P"]) if self.gwas_data is not None else []
        )
        data["lead_snp"] = (
            self.gwas_data["SNP"].iloc(self.get_current_lead_snp_index())
            if self.gwas_data is not None
            else None
        )
        data["ld_values"] = self.r2
        data["positions"] = (
            list(self.gwas_data["POS"]) if self.gwas_data is not None else []
        )
        data["chrom"], data["startbp"], data["endbp"] = self.get_locus_tuple()
        data["ld_populations"] = self.get_ld_population()
        data["gtex_tissues"], data["gtex_genes"] = self.get_gtex_selection()
        data["gene"] = gene_names(self.gene, self.get_coordinate())[0]
        data["gtex_version"] = self.get_gtex_version()
        data["coordinate"] = self.get_coordinate()
        data["set_based_p"] = self.set_based_p
        data["std_snp_list"] = self.std_snp_list
        data["runcoloc2"] = self.coloc2 is True
        num_GTEx_matches = get_gtex_snp_matches(
            self.std_snp_list, self.get_locus(), self.get_coordinate()
        )
        data["snp_warning"] = (
            get_gtex_snp_matches(self.std_snp_list, self.get_locus(), self.get_coordinate())
            / len(self.std_snp_list)
            < 0.8
        )
        data["thresh"] = 0.8
        data["numGTExMatches"] = num_GTEx_matches

        # secondary datasets
        if self.secondary_datasets is not None:
            data["secondary_dataset_titles"] = list(self.secondary_datasets.keys())
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
        data["num_SS_snps"] = len(self.std_snp_list)
        if self.ss_result_df is None:
            data["first_stages"] = []
            data["first_stage_Pvalues"] = []
        else:
            data["first_stages"] = self.ss_result_df["first_stages"].tolist()
            data["first_stage_Pvalues"] = self.ss_result_df["first_stage_p"].tolist()

        return data
