from tqdm import tqdm
from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.gtex import get_gtex_data


class GetGTExDataStage(PipelineStage):
    """
    A stage to get selected GTEx genes and tissues.

    Prerequisite: GWAS data is loaded in session (full region, not just SS region).
    """
    def invoke(self, payload: SessionPayload) -> SessionPayload:

        gtex_version = "V7"
        if payload.get_coordinate() == "hg38":
            gtex_version = "V8"

        gtex_tissues, gtex_genes = payload.get_gtex_selection()
        if len(gtex_tissues)>0:
            for tissue in tqdm(gtex_tissues):
                # for the full region (not just the SS region)
                eqtl_df = get_gtex_data(
                    gtex_version, 
                    tissue,
                    gene,
                    snp_list,
                    raiseErrors=True
                )
                if len(eqtl_df) > 0:
                    eqtl_df.fillna(-1, inplace=True)
                gtex_data[tissue] = eqtl_df.to_dict(orient='records')

        return payload
    