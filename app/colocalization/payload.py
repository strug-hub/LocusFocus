from uuid import uuid4
from dataclasses import dataclass
from typing import List, Literal


@dataclass
class SessionPayload(object):
    """
    Payload object for colocalization sessions.

    Contains metadata, form inputs and user data uploaded
    for this session.
    """
    
    # Form Inputs
    coordinate: Literal['hg38', 'hg19']
    coloc2: bool
    ld_population: str
    infer_variant: bool
    gtex_tissues: List[str]
    gtex_genes: List[str]
    set_based_p: float
    plot_locus: str
    simple_sum_locus: str
    lead_snp_name: str
    
    def __init__(self):
        self.session_id = uuid4()


@dataclass
class COLOC2SessionPayload(SessionPayload):
    """ 
    Payload object for colocalization sessions that also include COLOC2.
    """
