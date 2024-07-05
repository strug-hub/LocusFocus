from dataclasses import dataclass

import pandas as pd
import numpy as np

from app.colocalization.payload import SessionPayload
from app.colocalization.utils import download_file, standardize_snps, decompose_variant_list, x_to_23
from app.pipeline import PipelineStage
from app.routes import InvalidUsage


@dataclass
class GWASColumn():
    """
    Info about a given GWAS column.
    """
    form_id: str  # The ID of the form element where the user enters this column
    default: str  # Column default value
    coloc2: bool = False  # Required for coloc2?


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
    VALID_GWAS_EXTENSIONS = ['txt', 'tsv']

    GWAS_COLUMNS = [
        GWASColumn("chrom-col", "CHROM"),
        GWASColumn("pos-col", "POS"),
        GWASColumn("pval-col", "P"),
        GWASColumn("snp-col", "SNP"),
        GWASColumn("ref-col", "REF"),
        GWASColumn("alt-col", "ALT"),
        GWASColumn("beta-col", "BETA", coloc2=True),
        GWASColumn("stderr-col", "SE", coloc2=True),
        GWASColumn("numsamples-col", "N", coloc2=True),
        GWASColumn("maf-col", "MAF", coloc2=True),
    ]


    def invoke(self, payload: SessionPayload) -> SessionPayload:

        gwas_data = self._read_gwas_file(payload)
        gwas_data = self._set_gwas_columns(payload, gwas_data)
        gwas_data = self._validate_gwas_file(payload, gwas_data)
        gwas_data = self._subset_gwas_file(payload, gwas_data)

        # Update lead SNP based on user input
        lead_snp_index = self._get_lead_snp(payload, gwas_data)
        payload.gwas_lead_snp_index = lead_snp_index

        payload.gwas_data = gwas_data

        return payload # TODO: finish this


    def _read_gwas_file(self, payload: SessionPayload) -> pd.DataFrame:
        """
        Read any file with a valid file extension as a GWAS file.

        Save the file and return the dataframe.
        """
        gwas_filepath = download_file(payload.request, self.VALID_GWAS_EXTENSIONS)

        if gwas_filepath is None:
            raise InvalidUsage(f"GWAS file could not be found in uploaded files. Please upload a GWAS dataset in TSV format ({', '.join(self.VALID_GWAS_EXTENSIONS)})")

        try:
            gwas_data = pd.read_csv(gwas_filepath, sep='\t', encoding='utf-8')
        except:
            outfile = gwas_filepath.replace('.txt','_mod.txt').replace('.tsv', '_mod.tsv')
            with open(gwas_filepath) as f:
                with open(outfile, 'w') as fout:
                    filestr = f.readlines()
                    for line in filestr:
                        if line[0:2] != "##":
                            fout.write(line.replace('\t\t\n','\t\n'))
            try:
                gwas_data = pd.read_csv(outfile, sep="\t", encoding='utf-8')
            except:
                raise InvalidUsage('Failed to load primary dataset. Please check formatting is adequate.', status_code=410)

        return gwas_data


    def _set_gwas_columns(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> pd.DataFrame:
        """
        Read form inputs for column names,
        verify that the columns match those in the DataFrame,
        and rename the columns to a standard set of column names.
        """

        # We rename the columns so that we don't have to pass variables around everywhere for each column name.
        # instead, we can just use the same column name everywhere.

        # TODO: In the future, implement a way to infer columns based on a sample of the column's content

        # Get all form inputs, rename columns accordingly
        # Extra labels do not cause an error, so missing columns need to be checked
        old_gwas_columns = list(gwas_data.columns)

        column_mapper = {}  # maps user input -> default
        column_inputs = {}  # maps default -> raw user input
        column_input_list = []

        for column in self.GWAS_COLUMNS:
            column_input_list.append(payload.request.form.get(column.form_id, column.default))
            column_mapper[payload.request.form.get(column.form_id, column.default)] = column.default
            column_inputs[column.default] = payload.request.form.get(column.form_id, "")

        # Column uniqueness check
        if len(set(column_input_list)) != len(column_input_list):
            raise InvalidUsage(f'Duplicate column names provided: {column_input_list}')

        gwas_data = gwas_data.rename(columns=column_mapper)

        # Get P value column (always needed)
        if "P" not in gwas_data.columns:
            raise InvalidUsage(f"P value column not found in provided GWAS file.")

        infer_variant = bool(payload.request.form.get('markerCheckbox'))

        # Get CHROM, POS, REF, ALT
        if infer_variant:
            gwas_data = self._infer_gwas_columns(payload, gwas_data)
        else:
            for column in ["CHROM", "POS", "REF", "ALT"]:
                if column not in gwas_data.columns:
                    raise InvalidUsage(f"{column} column missing where required. '{column_inputs[column]}' not in columns '{', '.join(old_gwas_columns)}'")

        # Get coloc2 if applicable
        if payload.get_is_coloc2():
            for column in [c for c in self.GWAS_COLUMNS if c.coloc2]:
                if column not in gwas_data.columns:
                    raise InvalidUsage(f"{column} column missing but is required for COLOC2. '{column_inputs[column]}' not in columns '{', '.join(old_gwas_columns)}'")

            # Extra columns that the user might provide in their GWAS file and we need for COLOC2
            studytype = payload.get_coloc2_study_type()
            if 'type' not in gwas_data.columns:
                studytypedf = pd.DataFrame({'type': np.repeat(studytype, gwas_data.shape[0]).tolist()})
                gwas_data = pd.concat([gwas_data, studytypedf], axis=1)
            if studytype == "cc":
                num_cases = payload.get_coloc2_case_control_cases()
                if 'Ncases' not in gwas_data.columns:
                    num_cases_df = pd.DataFrame({'Ncases': np.repeat(num_cases, gwas_data.shape[0]).tolist()})
                    gwas_data = pd.concat([gwas_data, num_cases_df], axis=1)

        return gwas_data


    def _infer_gwas_columns(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> pd.DataFrame:
        """
        Given that the user wishes to infer the variant information using dbSNP,
        validate and rename the columns in the DataFrame accordingly, if possible.

        This assumes that the only columns provided in the GWAS file are "SNP" and "P", where "SNP" contains
        snp IDs that can be looked up using dbSNP for variant information.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """
        if "SNP" not in gwas_data.columns:
            raise InvalidUsage(f"SNP ID column not found when requesting to infer variant information.")

        coordinate = payload.get_coordinate()
        regionstr = payload.get_locus()

        variant_list = standardize_snps(list(gwas_data["SNP"]), regionstr, coordinate)
        if all(x=="." for x in variant_list):
            raise InvalidUsage(f"None of the variants provided could be mapped to {regionstr}!", status_code=410)
        var_df = decompose_variant_list(variant_list)
        gwas_data = pd.concat([var_df, gwas_data], axis=1) # add CHROM, POS, REF, ALT
        gwas_data = gwas_data.loc[ [str(x) != '.' for x in list(gwas_data["CHROM"])] ].copy()
        gwas_data.reset_index(drop=True, inplace=True)
        return gwas_data


    def _validate_gwas_file(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> pd.DataFrame:
        """
        Validate the contents of each column in the GWAS file, for those that exist.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """

        # Chromosome check
        if not all(isinstance(x, int) for x in x_to_23(list(gwas_data["CHROM"]))):
            raise InvalidUsage(f'Chromosome column contains unrecognizable values', status_code=410)

        # Position check
        if not all(isinstance(x, int) for x in list(gwas_data["POS"])):
            raise InvalidUsage(f'Position column has non-integer entries', status_code=410)

        # P value check
        if not all(isinstance(x, float) for x in list(gwas_data["P"])):
            raise InvalidUsage(f'P-value column has non-numeric entries', status_code=410)

        # COLOC2 checks
        if payload.get_is_coloc2():
            if not all(isinstance(x, float) for x in list(gwas_data["BETA"])):
                raise InvalidUsage(f'Beta column has non-numeric entries')
            if not all(isinstance(x, float) for x in list(gwas_data["SE"])):
                raise InvalidUsage(f'Standard error column has non-numeric entries')
            if not all(isinstance(x, int) for x in list(gwas_data["N"])):
                raise InvalidUsage(f'Number of samples column has non-integer entries')
            if not all(isinstance(x, float) for x in list(gwas_data["MAF"])):
                raise InvalidUsage(f'MAF column has non-numeric entries')

        return gwas_data

    def _subset_gwas_file(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> pd.DataFrame:
        """
        Subset the provided GWAS dataset based on the regionstr provided in the request.

        Prerequisite: gwas_data columns are set as default values (eg. "SNP", "P", etc.)
        """
        # mostly copied from subsetLocus
        chrom, start, end = payload.get_locus_tuple()

        gwas_data = gwas_data.loc[ [str(x) != '.' for x in list(gwas_data["CHROM"])] ].copy()
        bool1 = [x == chrom for x in x_to_23(list(gwas_data["CHROM"]))]
        bool2 = [x>=start and x<=end for x in list(gwas_data["POS"])]
        bool3 = [not x for x in list(gwas_data.isnull().any(axis=1))]
        bool4 = [str(x) != '.' for x in list(gwas_data["CHROM"])]
        gwas_indices_kept = [ all(conditions) for conditions in zip(bool1,bool2,bool3,bool4)]
        gwas_data = gwas_data.loc[ gwas_indices_kept ].copy()
        gwas_data.sort_values(by=[ "POS" ], inplace=True)
        chromcolnum = list(gwas_data.columns).index("CHROM")
        gwas_data.reset_index(drop=True, inplace=True)
        gwas_data.iloc[:,chromcolnum] = x_to_23(list(gwas_data["CHROM"])) # type: ignore
        if gwas_data.shape[0] == 0:
            raise InvalidUsage('No data found for entered region', status_code=410)
        # Check for invalid p=0 rows:
        zero_p = [x for x in list(gwas_data["P"]) if x==0]
        if len(zero_p)>0:
            raise InvalidUsage('P-values of zero detected; please replace with a non-zero p-value')

        payload.gwas_indices_kept = gwas_indices_kept

        return gwas_data
    
    def _get_lead_snp(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> int:
        """
        Determine the lead SNP index for this gwas dataset. Might not be used but is handy to have stored ahead of time.
        """

        lead_snp = payload.get_lead_snp_name()
        snp_list = list(gwas_data.loc[:,"SNP"])
        # cleaning up the SNP names a bit 
        snp_list = [asnp.split(';')[0] for asnp in snp_list] # type: ignore
        if lead_snp=='': lead_snp = list(gwas_data.loc[ gwas_data.loc[:,"P"] == gwas_data.loc[:,"P"].min() ].loc[:,"SNP"])[0].split(';')[0] # type: ignore
        if lead_snp not in snp_list:
            raise InvalidUsage('Lead SNP not found', status_code=410)
        lead_snp_position_index = snp_list.index(lead_snp)

        return lead_snp_position_index