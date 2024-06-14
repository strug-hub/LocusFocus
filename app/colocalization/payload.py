from uuid import uuid4
from typing import List, Literal, Dict, Optional, Tuple

import pandas as pd
from flask import Request

from app.colocalization.utils import parse_region_text
from app.routes import InvalidUsage


class SessionPayload(object):
    """
    Payload object for colocalization sessions.

    Contains metadata, form inputs and user data uploaded
    for this session.
    """
    # Constants
    VALID_COORDINATES = ['hg38', 'hg19']
    VALID_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR", "ASN", "NFE"]

    # Request object
    request: Request
    
    # Form Inputs
    coordinate: Optional[Literal['hg38', 'hg19']]
    coloc2: Optional[bool]
    study_type: Optional[Literal["quant", "cc"]]
    num_cases: Optional[int] # Used for study_type "cc" (case-control)
    ld_population: Optional[str]
    infer_variant: Optional[bool]
    gtex_tissues: Optional[List[str]]
    gtex_genes: Optional[List[str]]
    set_based_p: Optional[float]
    plot_locus: Optional[str]  # regionstr
    simple_sum_locus: Optional[str]
    lead_snp_name: Optional[str]

    # Files
    gwas_data: pd.DataFrame
    ld_data: pd.DataFrame
    secondary_datasets: Optional[Dict[str, pd.DataFrame]]

    # Other
    gwas_indices_kept: List[bool]
    
    def __init__(self, request: Request):
        self.session_id = uuid4()
        self.request = request

    def get_coordinate(self) -> Literal['hg38', 'hg19']:
        """
        Gets the form input for coordinate (aka. genome assembly, or 'build') for this session.
        """
        if self.coordinate is None:
            if self.request.form.get("coordinate") not in self.VALID_COORDINATES:
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