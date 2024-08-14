from dataclasses import dataclass, field
from uuid import uuid4, UUID
from typing import List, Literal, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from flask import Request

from app.colocalization.utils import parse_region_text
from app.colocalization.constants import ONE_SIDED_SS_WINDOW_SIZE, VALID_COORDINATES, VALID_POPULATIONS
from app.routes import InvalidUsage


@dataclass
class SessionPayload():
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

    # File data
    gwas_data: Optional[pd.DataFrame] = None
    ld_data: Optional[np.matrix] = None
    secondary_datasets: Optional[Dict[str, dict]] = None

    # LD Matrix data
    # DataFrame containing SNPs that were actually used in LD. See: https://www.cog-genomics.org/plink/1.9/formats#bim
    ld_snps_bim_df: Optional[pd.DataFrame] = None 

    # Other
    success: bool = False
    gwas_indices_kept: List[bool] = field(default_factory=list)
    gwas_lead_snp_index: Optional[int] = None
    r2: List[float] = field(default_factory=list)

    # Simple Sum
    ss_locustext: Optional[str] = None

    # GTEx
    gtex_data: dict = {}

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
    
    def get_lead_snp_index(self) -> int:
        """
        Get the index of the lead SNP for the GWAS dataset.
        """
        lead_snp = self.get_lead_snp_name()
        if self.gwas_data is None:
            raise Exception("Cannot get lead SNP index when GWAS dataset is undefined.")
        
        snps = [snp.split(";")[0] for snp in self.gwas_data.loc[:, "SNP"]] # type: ignore
        if lead_snp == "":
            lead_snp = list(self.gwas_data.loc[ self.gwas_data.loc[:,"P"] == self.gwas_data.loc[:,"P"].min() ].loc[:,"SNP"])[0].split(';')[0] # type: ignore
        if lead_snp not in snps:
            raise InvalidUsage('Lead SNP not found', status_code=410)
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
            SSchrom, SS_start, SS_end = parse_region_text(SSlocustext, self.get_coordinate())
        else:
            if self.gwas_data is None:
                raise Exception("Need GWAS dataset in order to find Lead SNP for SS region")
            SSchrom, _, _ = self.get_locus_tuple()
            lead_snp_position_index = self.get_lead_snp_index()
            lead_snp_position = int(self.gwas_data.iloc[lead_snp_position_index, :]["POS"]) # type: ignore
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
            raise InvalidUsage('Please select one or more genes to complement your GTEx tissue(s) selection', status_code=410)
        elif len(self.gtex_genes) > 0 and len(self.gtex_tissues) == 0:
            raise InvalidUsage('Please select one or more tissues to complement your GTEx gene(s) selection', status_code=410)

        return self.gtex_tissues, self.gtex_genes

    
    def dump_session_data(self):
        """
        Create JSON dict of session data needed for form_data file.
        """

        data = {}
        data['success'] = self.success
        data['sessionid'] = str(self.session_id)
        data['snps'] = list(self.gwas_data["SNP"]) if self.gwas_data is not None else []
        data['infervariant'] = self.get_infer_variant()
        data['pvalues'] = list(self.gwas_data["P"]) if self.gwas_data is not None else []
        data['lead_snp'] = self.gwas_data["SNP"].iloc(self.get_lead_snp_index()) if self.gwas_data is not None else None
        data['ld_values'] = self.r2
        data['positions'] = list(self.gwas_data["POS"]) if self.gwas_data is not None else []
        data['chrom'], data['startbp'], data['endbp'] = self.get_locus_tuple()
        data['ld_populations'] = self.get_ld_population()
        # TODO: Unfinished session data fields
        data['gtex_tissues'] = None
        data['gene'] = None
        data['gtex_genes'] = None
        data['gtex_version'] = None

        data['coordinate'] = self.get_coordinate()

        data['set_based_p'] = setbasedP
        SSlocustext = request.form['SSlocus'] # SSlocus defined below
        data['std_snp_list'] = std_snp_list
        data['runcoloc2'] = runcoloc2
        data['snp_warning'] = snp_warning
        data['thresh'] = thresh
        data['numGTExMatches'] = numGTExMatches

        return data