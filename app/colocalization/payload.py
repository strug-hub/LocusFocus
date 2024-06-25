from dataclasses import dataclass, field
from uuid import uuid4, UUID
from typing import List, Literal, Dict, Optional, Tuple

import pandas as pd
from flask import Request

from app.colocalization.utils import parse_region_text
from app.colocalization.constants import VALID_COORDINATES, VALID_POPULATIONS
from app.routes import InvalidUsage


@dataclass
class SessionPayload(object):
    """
    Payload object for colocalization sessions.

    Contains metadata, form inputs and user data uploaded
    for this session.
    """

    # Request object
    request: Request
    session_id: UUID = field(default_factory=uuid4)
    
    # Form Inputs
    coordinate: Optional[Literal['hg38', 'hg19']] = None
    coloc2: Optional[bool] = None
    study_type: Optional[Literal["quant", "cc"]] = None
    num_cases: Optional[int] = None # Used for study_type "cc" (case-control)
    ld_population: Optional[str] = None
    infer_variant: Optional[bool] = None
    gtex_tissues: Optional[List[str]] = None
    gtex_genes: Optional[List[str]] = None
    set_based_p: Optional[float] = None
    plot_locus: Optional[str] = None  # regionstr
    simple_sum_locus: Optional[str] = None
    lead_snp_name: Optional[str] = None

    # Files
    gwas_data: Optional[pd.DataFrame] = None
    ld_data: Optional[pd.DataFrame] = None
    secondary_datasets: Optional[Dict[str, pd.DataFrame]] = None

    # Other
    gwas_indices_kept: List[bool] = []
    gwas_lead_snp_index: Optional[int] = None
    r2: List[float] = []

    def get_coordinate(self) -> Literal['hg38', 'hg19']:
        """
        Gets the form input for coordinate (aka. genome assembly, or 'build') for this session.
        """
        if self.coordinate is None:
            if self.request.form.get("coordinate") not in VALID_COORDINATES:
                raise InvalidUsage(f"Invalid coordinate: '{self.request.form.get('coordinate')}'")
            
            self.coordinate = self.request.form.get("coordinate") # type: ignore
        return self.coordinate # type: ignore
    
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
            self.coloc2 = bool(self.request.form.get('coloc2check'))
        return self.coloc2
    
    def get_coloc2_study_type(self) -> Literal["quant", "cc"]:
        """
        Get study type for COLOC2. Assumes that coloc2 is requested.
        """
        if self.study_type is None:
            study_type = self.request.form.get("studytype") # type: ignore
            if study_type not in ["quant", "cc"]:
                raise InvalidUsage(f"Study type form value is invalid: {self.request.form.get('studytype')}")
            self.study_type = study_type # type: ignore
        return self.study_type # type: ignore

    def get_coloc2_case_control_cases(self) -> int:
        """
        Get number of case control cases specified for coloc2. Assumes coloc2 is requested, and
        selected study type is Case-Control.
        """
        if self.num_cases is None:
            try:
                num_cases = self.request.form.get("numcases", type=int)
            except ValueError:
                raise InvalidUsage("Number of cases entered must be an integer", status_code=410)
            self.num_cases = num_cases

        return self.num_cases # type: ignore
    
    def get_ld_population(self) -> str:
        """
        Get selected LD population from user input (1000 Genomes dataset).
        If the user did not specify LD population, then use default "EUR".
        """
        if self.ld_population is None:
            pop = self.request.form.get("LD-populations", "EUR")
            if pop not in VALID_POPULATIONS:
                raise InvalidUsage(f"Invalid population provided: '{pop}'. Population must be one of '{', '.join(VALID_POPULATIONS)}'")
            self.ld_population = pop
        return self.ld_population
    
    def get_lead_snp_name(self) -> str:
        """
        Get the lead SNP name from user input, if specified.
        """

        if self.lead_snp_name is None:
            self.lead_snp_name = self.request.form.get("leadsnp", "")

        return self.lead_snp_name
