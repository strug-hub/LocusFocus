from tqdm import tqdm
from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.gtex import get_gtex_data


class ReportGTExDataStage(PipelineStage):
    """
    A stage to report on selected GTEx data.

    Prerequisite: GWAS data is loaded in session *before* subsetting (original uploaded data).
    """
    def invoke(self, payload: SessionPayload) -> SessionPayload:

        if payload.gwas_data is None:
            raise Exception("GWAS data not loaded; needed for GTEx data selection stage")

        gtex_data = {}
        coordinate = payload.get_coordinate()

        gtex_version = "V7"
        if coordinate == "hg38":
            gtex_version = "V8"

        gtex_tissues, gtex_genes = payload.get_gtex_selection()

        if len(gtex_genes) > 0:
            gene = gtex_genes[0]
        elif coordinate == 'hg19':
            gene = 'ENSG00000174502.14'
        elif coordinate == 'hg38':
            gene = 'ENSG00000174502.18'

        snp_list = [asnp.split(";")[0] for asnp in payload.gwas_data["SNP"]]

        if len(gtex_tissues) > 0:
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

        payload.gtex_data = gtex_data

        return payload
    