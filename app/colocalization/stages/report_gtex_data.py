import json
from tqdm import tqdm
from app.colocalization.payload import SessionPayload
from app.pipeline.pipeline_stage import PipelineStage
from app.utils.gtex import (
    get_gtex_data,
    collapsed_genes_df_hg19,
    collapsed_genes_df_hg38,
)

import numpy as np


class ReportGTExDataStage(PipelineStage):
    """
    A stage to report on selected GTEx data.

    Prerequisite:
    - GWAS data is loaded in session *before* subsetting (original uploaded data).
    - GWAS dataset contains only one chromosome.
    """

    def name(self):
        return "report-gtex-data"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        if payload.gwas_data is None:
            raise Exception(
                "GWAS data not loaded; needed for GTEx data selection stage"
            )

        gtex_data = {}
        gtex_version = payload.get_gtex_version()

        gtex_tissues, gtex_genes = payload.get_gtex_selection()

        if len(gtex_genes) > 0:
            gene = gtex_genes[0]
        elif gtex_version == "V7":
            gene = "ENSG00000174502.14"
        elif gtex_version == "V8":
            gene = "ENSG00000174502.18"

        snp_list = [asnp.split(";")[0] for asnp in payload.gwas_data["SNP"]]

        if len(gtex_tissues) > 0:
            for tissue in tqdm(gtex_tissues):
                # for the full region (not just the SS region)
                eqtl_df = get_gtex_data(
                    gtex_version, tissue, gene, snp_list, raiseErrors=True
                )
                if len(eqtl_df) > 0:
                    eqtl_df.fillna(-1, inplace=True)
                gtex_data[tissue] = eqtl_df.to_dict(orient="records")

        payload.reported_gtex_data = gtex_data
        payload.gene = gene

        self._report_genes_to_plot(payload)

        return payload

    def _report_genes_to_plot(self, payload: SessionPayload):
        """
        Saves a JSON file of genes in the region of interest for later plotting.

        The file is saved as genes_data-<session_id>.json and contains a list of
        dictionaries, where each dictionary represents a gene and contains the
        following keys:
            - name: the name of the gene
            - txStart: the start of the gene's transcript
            - txEnd: the end of the gene's transcript
            - exonStarts: a list of the starts of the gene's exons
            - exonEnds: a list of the ends of the gene's exons

        This data is used to draw the genes in the region of interest.
        """
        chrom, startbp, endbp = payload.get_locus_tuple()

        gtex_version = payload.get_gtex_version()
        if gtex_version == "V7":
            collapsed_genes_df = collapsed_genes_df_hg19
        else:
            collapsed_genes_df = collapsed_genes_df_hg38

        genes_to_draw = collapsed_genes_df.loc[
            (collapsed_genes_df["chrom"] == ("chr" + str(chrom).replace("23", "X")))
            & (
                (
                    (collapsed_genes_df["txStart"] >= startbp)
                    & (collapsed_genes_df["txStart"] <= endbp)
                )
                | (
                    (collapsed_genes_df["txEnd"] >= startbp)
                    & (collapsed_genes_df["txEnd"] <= endbp)
                )
                | (
                    (collapsed_genes_df["txStart"] <= startbp)
                    & (collapsed_genes_df["txEnd"] >= endbp)
                )
            )
        ]
        genes_data = []
        for i in np.arange(genes_to_draw.shape[0]):
            genes_data.append(
                {
                    "name": list(genes_to_draw["name"])[i],
                    "txStart": list(genes_to_draw["txStart"])[i],
                    "txEnd": list(genes_to_draw["txEnd"])[i],
                    "exonStarts": [
                        int(bp)
                        for bp in list(genes_to_draw["exonStarts"])[i].split(",")
                    ],
                    "exonEnds": [
                        int(bp) for bp in list(genes_to_draw["exonEnds"])[i].split(",")
                    ],
                }
            )

        json.dump(genes_data, open(payload.file.genes_session_filepath, "w"))
