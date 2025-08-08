from dataclasses import dataclass

import pandas as pd
import numpy as np
from flask import current_app as app

from app.colocalization.payload import SessionPayload
from app.utils import get_file_with_ext, standardize_snps, decompose_variant_list, x_to_23
from app.pipeline import PipelineStage
from app.utils.errors import InvalidUsage, ServerError


@dataclass
class GWASColumn:
    """
    Info about a given GWAS column.
    """

    form_id: str  # The ID of the form element where the user enters this column
    default: str  # Column default value
    coloc2: bool = False  # Required for coloc2?
    optional: bool = False  # Optional column

#                                                              (bi-allelic variants are allowed)
#                      (prefix) chrom        pos      ref       alt(s)                build
VCF_FORMAT_PATTERN = "^(?:chr)?([0-9]{1,2})_([0-9]+)_([ATCG]+)_([ATCG]+(?:,[ATCG]+)*)_b3(?:7|8)$"


class ReadGWASFileStage(PipelineStage):
    """
    Read a GWAS file into the payload as a DataFrame.

    Read the form for column names, and
    rename the columns to a standard set of column names.

    Standard column names are the .default field for each
    GWASColumn in GWAS_COLUMNS.

    Prerequisites:
    - Session is created.
    """

    VALID_GWAS_EXTENSIONS = ["txt", "tsv"]

    GWAS_COLUMNS = [
        # always needed
        GWASColumn("chrom-col", "CHROM"),
        GWASColumn("pos-col", "POS"),
        GWASColumn("pval-col", "P"),
        GWASColumn("snp-col", "SNP"),
        GWASColumn("ref-col", "REF"),
        GWASColumn("alt-col", "ALT"),
        # only needed for coloc2
        GWASColumn("beta-col", "BETA", coloc2=True),
        GWASColumn("stderr-col", "SE", coloc2=True),
        GWASColumn("numsamples-col", "N", coloc2=True),
        GWASColumn("maf-col", "MAF", coloc2=True),
        # part of coloc2 but not required, don't have form controls at the moment
        GWASColumn("type-col", "type", coloc2=True, optional=True),
        GWASColumn("ncases-col", "Ncases", coloc2=True, optional=True),
    ]

    def name(self) -> str:
        return "read-gwas-file"

    def __init__(self, enforce_one_chrom=True):
        self.enforce_one_chrom = enforce_one_chrom

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        gwas_data = self._read_gwas_file(payload)
        gwas_data = self._set_gwas_columns(payload, gwas_data)
        gwas_data = self._validate_gwas_file(payload, gwas_data)
        gwas_data = self._subset_gwas_file(payload, gwas_data)

        payload.gwas_data = gwas_data
        payload.gwas_indices_kept = pd.Series(True, index=gwas_data.index)

        # Get standardized list of SNPs
        payload.std_snp_list = self._get_std_snp_list(payload, gwas_data)

        self._snp_format_check(gwas_data)

        return payload

    def _read_gwas_file(self, payload: SessionPayload) -> pd.DataFrame:
        """
        Read any file with a valid file extension as a GWAS file.

        Save the file and return the dataframe.
        """
        gwas_filepath = get_file_with_ext(payload.uploaded_files, self.VALID_GWAS_EXTENSIONS)

        if gwas_filepath is None:
            raise InvalidUsage(
                f"GWAS file could not be found in uploaded files. Please upload a GWAS dataset in TSV format ({', '.join(self.VALID_GWAS_EXTENSIONS)})"
            )

        try:
            gwas_data = pd.read_csv(gwas_filepath, sep="\t", encoding="utf-8")
        except:
            outfile = gwas_filepath.replace(".txt", "_mod.txt").replace(
                ".tsv", "_mod.tsv"
            )
            with open(gwas_filepath) as f:
                with open(outfile, "w") as fout:
                    filestr = f.readlines()
                    for line in filestr:
                        if line[0:2] != "##":
                            fout.write(line.replace("\t\t\n", "\t\n"))
            try:
                gwas_data = pd.read_csv(outfile, sep="\t", encoding="utf-8")
            except:
                raise InvalidUsage(
                    "Failed to load primary dataset as tab-separated file. Please check formatting is adequate, and that the file is not empty.",
                    status_code=410,
                )

        return gwas_data

    def _set_gwas_columns(
        self, payload: SessionPayload, gwas_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Read form inputs for column names,
        verify that the columns match those in the DataFrame,
        and rename the columns to a standard set of column names.
        """

        # We rename the columns so that we don't have to pass variables around everywhere for each column name.
        # instead, we can just use the same column name everywhere.
        old_gwas_columns = list(gwas_data.columns)  # Columns in gwas_file

        column_mapper = {}  # maps user input OR default -> default
        column_inputs = {}  # maps default -> raw user input
        column_input_list = []  # all user input

        for column in self.GWAS_COLUMNS:
            column_input_list.append(
                payload.request_form.get(column.form_id, column.default)
            )
            column_mapper[
                payload.request_form.get(column.form_id, column.default)
            ] = column.default
            column_inputs[column.default] = payload.request_form.get(column.form_id, "")

        # Column uniqueness check
        if len(set(column_input_list)) != len(column_input_list):
            raise InvalidUsage(f"Duplicate column names provided: {column_input_list}")

        # Rename columns to match GWAS_COLUMNS.
        # All checks from this point forward should be based on GWAS_COLUMNS
        non_gwas_columns = [c for c in gwas_data.columns if c not in column_input_list]
        gwas_data = gwas_data.drop(columns=non_gwas_columns)
        gwas_data = gwas_data.rename(columns=column_mapper)

        # Get P value column (always needed)
        if "P" not in gwas_data.columns:
            raise InvalidUsage(
                f"P value column not found in provided GWAS file. Found columns: '{', '.join(old_gwas_columns)}'. Please verify that your GWAS file contains a P value column."
            )

        infer_variant = bool(payload.request_form.get("markerCheckbox"))

        # Get CHROM, POS, REF, ALT, SNP
        if infer_variant:
            gwas_data = self._infer_gwas_columns(payload, gwas_data)
        else:
            for column in ["CHROM", "POS", "REF", "ALT", "SNP"]:
                if column not in gwas_data.columns:
                    raise InvalidUsage(
                        f"'{column}' column missing where required. '{column_inputs[column]}' not in columns '{', '.join(old_gwas_columns)}'. Please update your GWAS columns to match, or type a different column name that is found in your dataset."
                    )

        # Get coloc2 if applicable
        if payload.get_is_coloc2():
            for column in [c.default for c in self.GWAS_COLUMNS if c.coloc2 and not c.optional]:
                if column not in gwas_data.columns:
                    raise InvalidUsage(
                        f"{column} column missing but is required for COLOC2. '{column_inputs[column]}' not in columns '{', '.join(old_gwas_columns)}'"
                    )

            # Extra columns that the user might provide in their GWAS file and we need for COLOC2
            studytype = payload.get_coloc2_study_type()
            if "type" not in gwas_data.columns:
                studytypedf = pd.DataFrame(
                    {"type": np.repeat(studytype, gwas_data.shape[0]).tolist()}
                )
                gwas_data = pd.concat([gwas_data, studytypedf], axis=1)
            if studytype == "cc":
                num_cases = payload.get_coloc2_case_control_cases()
                if "Ncases" not in gwas_data.columns:
                    num_cases_df = pd.DataFrame(
                        {"Ncases": np.repeat(num_cases, gwas_data.shape[0]).tolist()}
                    )
                    gwas_data = pd.concat([gwas_data, num_cases_df], axis=1)

        # Drop columns that are not used
        gwas_data = gwas_data.drop(columns=[c for c in gwas_data.columns if c not in (c.default for c in self.GWAS_COLUMNS)])

        return gwas_data

    def _infer_gwas_columns(
        self, payload: SessionPayload, gwas_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Given that the user wishes to infer the variant information using dbSNP,
        validate and rename the columns in the DataFrame accordingly, if possible.

        This assumes that the only columns provided in the GWAS file are "SNP" and "P", where "SNP" contains
        snp IDs that can be looked up using dbSNP for variant information.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """
        if "SNP" not in gwas_data.columns:
            raise InvalidUsage(
                f"SNP ID column not found when requesting to infer variant information."
            )

        coordinate = payload.get_coordinate()
        regionstr = payload.get_locus()

        variant_list = standardize_snps(list(gwas_data["SNP"]), regionstr, coordinate)
        if all(x == "." for x in variant_list):
            raise InvalidUsage(
                f"None of the variants provided could be mapped to {regionstr}!",
                status_code=410,
            )
        var_df = decompose_variant_list(variant_list)
        gwas_data = pd.concat([var_df, gwas_data], axis=1)  # add CHROM, POS, REF, ALT
        gwas_data = gwas_data.loc[
            [str(x) != "." for x in list(gwas_data["CHROM"])]
        ].copy()
        gwas_data.reset_index(drop=True, inplace=True)
        return gwas_data

    def _validate_gwas_file(
        self, payload: SessionPayload, gwas_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Validate the contents of each column in the GWAS file, for those that exist.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """

        # TODO: make helper function for error messages
        # Chromosome check
        converted_chorms = x_to_23(list(gwas_data["CHROM"]))
        if not all(isinstance(x, int) for x in converted_chorms):
            raise InvalidUsage(
                f'Chromosome column contains unrecognizable values: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(converted_chorms) if not isinstance(x, int))}',
                status_code=410,
            )

        # Position check
        if not all(isinstance(x, int) for x in list(gwas_data["POS"])):
            raise InvalidUsage(
                f'Position column has non-integer entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["POS"])) if not isinstance(x, int))}',
                status_code=410,
            )

        # P value check
        if not all(isinstance(x, float) for x in list(gwas_data["P"])):
            raise InvalidUsage(
                f'P-value column has non-numeric entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["P"])) if not isinstance(x, float))}',
                status_code=410,
            )

        # COLOC2 checks
        if payload.get_is_coloc2():
            if not all(isinstance(x, float) for x in list(gwas_data["BETA"])):
                raise InvalidUsage(
                    f'Beta column has non-numeric entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["BETA"])) if not isinstance(x, float))}',
                    status_code=410,
                )
            if not all(isinstance(x, float) for x in list(gwas_data["SE"])):
                raise InvalidUsage(
                    f'Standard error column has non-numeric entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["SE"])) if not isinstance(x, float))}',
                    status_code=410,
                )
            if not all(isinstance(x, int) for x in list(gwas_data["N"])):
                raise InvalidUsage(
                    f'Number of samples column has non-integer entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["N"])) if not isinstance(x, int))}',
                    status_code=410,
                )
            if not all(isinstance(x, float) for x in list(gwas_data["MAF"])):
                raise InvalidUsage(
                    f'MAF column has non-numeric entries: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(list(gwas_data["MAF"])) if not isinstance(x, float))}',
                    status_code=410,
                )

        # One chrom check
        if self.enforce_one_chrom and gwas_data["CHROM"].nunique() > 1:
            unique_chroms = list(gwas_data["CHROM"].unique())  # type: ignore
            raise InvalidUsage(
                f"Multiple chromosomes provided where only 1 is required: {unique_chroms}"
            )

        return gwas_data

    def _subset_gwas_file(
        self, payload: SessionPayload, gwas_data: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Subset the provided GWAS dataset based on the regionstr provided in the request.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """
        # mostly copied from subsetLocus
        chrom, start, end = payload.get_locus_tuple()

        gwas_data = gwas_data.loc[
            [str(x) != "." for x in list(gwas_data["CHROM"])]
        ].copy()
        bool1 = [x == chrom for x in x_to_23(list(gwas_data["CHROM"]))]
        bool2 = [x >= start and x <= end for x in list(gwas_data["POS"])]
        bool3 = [not x for x in list(gwas_data.isnull().any(axis=1))]
        bool4 = [str(x) != "." for x in list(gwas_data["CHROM"])]
        gwas_indices_kept = pd.Series(
            [all(conditions) for conditions in zip(bool1, bool2, bool3, bool4)]
        )
        if not all(gwas_indices_kept):
            app.logger.debug(
                f"{sum(gwas_indices_kept)} SNPs kept, {sum(~gwas_indices_kept)} SNPs removed."
            )
        gwas_data = gwas_data.loc[gwas_indices_kept].copy()
        gwas_data.sort_values(by=["POS"], inplace=True)
        gwas_data.reset_index(drop=True, inplace=True)
        gwas_data["CHROM"] = x_to_23(list(gwas_data["CHROM"]))  # type: ignore
        if gwas_data.shape[0] == 0:
            raise InvalidUsage(
                f"No data found for entered region: '{chrom}:{start}-{end}'",
                status_code=410,
            )
        # Check for invalid p=0 rows:
        zero_p = [x for x in list(gwas_data["P"]) if x == 0]
        if len(zero_p) > 0:
            raise InvalidUsage(
                f'P-values of zero detected; please replace with a non-zero p-value: {", ".join(f"Row {i+1}: {x}" for i, x in enumerate(zero_p) if x==0)}',
                status_code=410,
            )

        # payload.gwas_indices_kept = gwas_indices_kept

        return gwas_data

    def _get_std_snp_list(
        self, payload: SessionPayload, gwas_data: pd.DataFrame
    ) -> pd.Series:
        """
        Return standardized list of SNPs with format CHR_POS_REF_ALT_build.

        gwas_data needs to be subsetted in advance.
        """
        std_snp_list = []
        buildstr = "b37"
        if payload.get_coordinate() == "hg38":
            buildstr = "b38"

        std_snp_list = pd.Series(
            [
                f"{str(row['CHROM']).replace('23', 'X')}_{str(row['POS'])}_{str(row['REF'])}_{str(row['ALT'])}_{buildstr}"
                for _, row in gwas_data.iterrows()
            ]
        )
        # Sanity check
        try:
            assert len(std_snp_list) == len(gwas_data)
            assert len(std_snp_list) == len(payload.gwas_indices_kept)
        except AssertionError:
            raise ServerError("GWAS data and indices are not in sync")

        return std_snp_list
    
    def _snp_format_check(self, gwas_data: pd.DataFrame) -> None:
        """
        Perform a sanity check on the SNP column of the GWAS data.
        
        Raise an error if the SNP column contains any SNPs with format chr_pos_ref_alt_build
        have reversed alleles, or other inconsistencies within its row in the dataframe.
        """

        vcf_snps = gwas_data[gwas_data["SNP"].str.contains(VCF_FORMAT_PATTERN)]
        if len(vcf_snps) > 0:
            expanded = gwas_data["SNP"].str.extract(VCF_FORMAT_PATTERN) \
                .rename(columns={0: "CHROM", 1: "POS", 2: "REF", 3: "ALT"}) \
                .astype({"CHROM": int, "POS": int, "REF": str, "ALT": str})
                
            merged_vcf_snps = expanded.merge(vcf_snps, how="inner", on=["CHROM", "POS", "REF", "ALT"])
            if len(vcf_snps) != len(merged_vcf_snps):
                raise InvalidUsage(
                    "GWAS data contains SNPs with chrom_pos_ref_alt_build format that are not consistent (eg. reversed alleles, incorrect position or chromosome). "
                    "Please inspect your GWAS file and ensure that the SNP column is consistent with the values in other columns in the GWAS file. "
                    f"{len(vcf_snps) - len(merged_vcf_snps)} of {len(gwas_data)} SNPs are inconsistent."
                )
