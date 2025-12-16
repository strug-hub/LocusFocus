"""
Shared constant values for colocalization.
"""

VALID_COORDINATES = ["hg38", "hg19"]
VALID_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR", "ASN", "NFE"]
LD_MAT_DIAG_CONSTANT = 1e-6
ONE_SIDED_SS_WINDOW_SIZE = 100_000  # (100 kb on either side of the lead SNP)
