import os
from typing import Tuple
import pandas as pd
from flask import current_app as app
from app.colocalization.payload import SessionPayload
from app.colocalization.utils import clean_snps
from app.pipeline.pipeline_stage import PipelineStage
from app.routes import InvalidUsage


class SimpleSumSubsetGWASStage(PipelineStage):
    """
    Stage responsible for subsetting the SNPs in the GWAS dataset such that
    Simple Sum colocalization may be performed at a later stage.

    Must take place after a GWAS dataset is added to the payload.
    Must also take place before any stages that rely on GWAS dataset (eg. LD matrix).
    """

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        if payload.gwas_data is None:
            raise Exception("GWAS dataset not found")

        SS_gwas_data, ss_indices = self._subset_gwas(payload, payload.gwas_data)
        payload.gwas_data = SS_gwas_data
        self._write_gwas_to_file(payload, payload.gwas_data, ss_indices) # type: ignore
        self._check_pos_duplicates(SS_gwas_data)

        return payload


    def _subset_gwas(self, payload: SessionPayload, gwas_data: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Given a GWAS dataset, return a new GWAS dataset that is ready for Simple Sum.

        SNPs are removed if they fall outside of the specified locus for Simple Sum.
        """
        chrom, SS_start, SS_end = payload.get_ss_locus_tuple()

        chromList = [('chr' + str(chrom).replace('23','X')), str(chrom).replace('23','X')]
        if 'X' in chromList:
            chromList.extend(['chr23','23'])
        gwas_chrom_col = pd.Series([str(x) for x in list(gwas_data["CHROM"])])
        SS_chrom_bool = [str(x).replace('23','X') for x in gwas_chrom_col.isin(chromList) if x == True]
        SS_indices = SS_chrom_bool & (gwas_data["POS"] >= SS_start) & (gwas_data["POS"] <= SS_end)
        SS_gwas_data = gwas_data.loc[ SS_indices ]

        # TODO: Include this step somewhere in COLOC2 handling stages
        # if runcoloc2:
        #     coloc2_gwasdf = SS_gwas_data.rename(columns={
        #         chromcol: 'CHR'
        #         ,poscol: 'POS'
        #         ,snpcol: 'SNPID'
        #         ,pcol: 'PVAL'
        #         ,refcol: REF
        #         ,altcol: ALT
        #         ,betacol: BETA
        #         ,stderrcol: SE
        #         ,mafcol: MAF
        #         ,numsamplescol: 'N'
        #     })
        #     coloc2_gwasdf = coloc2_gwasdf.reindex(columns=coloc2gwascolnames)

        if SS_gwas_data.shape[0] == 0:
            raise InvalidUsage('No data points found for entered Simple Sum region', status_code=410)

        return SS_gwas_data, SS_indices
    
    
    def _write_gwas_to_file(self, payload: SessionPayload, gwas_data: pd.DataFrame, ss_indices: pd.Series):
        """
        Writes data from subsetted GWAS dataset to file for reporting purposes.
        """
        regionstr = payload.get_locus()
        coordinate = payload.get_coordinate()

        std_snp_list = clean_snps(list(gwas_data["SNP"]), regionstr, coordinate)
        ss_std_snp_list = [e for i,e in enumerate(std_snp_list) if ss_indices[i]]  # TODO: Is this not redundant?

        gwas_df = pd.DataFrame({
            'Position': list(gwas_data["POS"]),
            'SNP': std_snp_list,
            'variant_id': ss_std_snp_list,
            'P': list(gwas_data["P"])
        })
        with app.app_context():
            gwas_df.to_csv(
                os.path.join(app.config["SESSION_FOLDER"], f'gwas_df-{payload.session_id}.txt'), 
                index=False, 
                encoding='utf-8', 
                sep="\t"
            )


    def _check_pos_duplicates(self, subsetted_gwas_data: pd.DataFrame):
        """
        Raise error if there are duplicate positions in the now-subsetted GWAS data.

        TODO: positions with different alt alleles count as duplicates if same position occurs. 
        Determine if this is okay, or if there's a better way to check here.
        """

        positions = subsetted_gwas_data["POS"]
        if len(positions) != len(set(positions)):
            # collect duplicates for error message
            dups = set([x for x in positions if positions.count(x) > 1])
            dup_counts = [(x, positions.count(x)) for x in dups]
            raise InvalidUsage(f'Duplicate chromosome basepair positions detected: {[f"bp: {dup[0]}, num. duplicates: {dup[1]}" for dup in dup_counts]}')
        return None
