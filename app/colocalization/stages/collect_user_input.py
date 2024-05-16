from flask import Request

from app.colocalization.payload import SessionPayload
from app.pipeline import PipelineStage
from app.routes import InvalidUsage


class CollectUserInputStage(PipelineStage):
    """
    Stage 1 of the Colocalization pipeline.

    Reads and validates the provided user input (files, form inputs, etc.),
    and creates a new payload object that will be used for the rest of the
    pipeline chain.
    """
    VALID_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR", "ASN", "NFE"]
    VALID_COORDINATES = ["hg19", "hg38"]
    
    def invoke(self, request: Request) -> SessionPayload:
        
        new_payload = SessionPayload()

        new_payload = self._read_form_inputs(request, new_payload)
        new_payload = self._read_files(request, new_payload)
        new_payload = self._fix_gwas_columns(request, new_payload)
        
        return new_payload

    def _read_form_inputs(self, request: Request, new_payload: SessionPayload) -> SessionPayload:
        """
        Populate the payload with form inputs read from the given request.
        """
        # Collect all form errors, check length at end, and then raise exception if > 0
        errors = []

        if request.form.get("coordinate") not in self.VALID_COORDINATES:
            errors.append(f"Invalid coordinate: '{request.form.get('coordinate')}'")
        new_payload.coordinate = request.form.get("coordinate", "hg19") # type: ignore
        
        new_payload.coloc2 = bool(request.form.get("coloc2check"))

        if request.form.get("LD-populations") not in self.VALID_POPULATIONS:
            errors.append(f"Invalid 1000 Genomes population: '{request.form.get('LD-populations')}'")
        new_payload.ld_population = request.form.get("LD-populations", "EUR")

        new_payload.infer_variant = bool(request.form.get("markerCheckbox"))
        new_payload.plot_locus = request.form.get("locus", "1:205500000-206000000")
        new_payload.simple_sum_locus = request.form.get("SSlocus", "")

        # GTEx
        new_payload.gtex_tissues = request.form.getlist("GTEx-tissues")
        new_payload.gtex_genes = request.form.getlist("region-genes")
        if len(new_payload.gtex_tissues) > 0 and len(new_payload.gtex_genes) == 0:
            errors.append('Please select one or more genes to complement your GTEx tissue(s) selection')
        elif len(new_payload.gtex_genes) > 0 and len(new_payload.gtex_tissues) == 0:
            errors.append('Please select one or more tissues to complement your GTEx gene(s) selection')

        # First stage set-based test P value threshold
        p_threshold = request.form.get("setbasedP", "")
        if p_threshold == "":
            p_threshold = "default"
        else:
            try:
                p_threshold = float(p_threshold)
                if p_threshold < 0 or p_threshold > 1:
                    errors.append('Set-based p-value threshold given is not between 0 and 1')
                else:
                    new_payload.set_based_p = p_threshold
            except:
                errors.append('Invalid value provided for the set-based p-value threshold. Value must be numeric between 0 and 1.')

        new_payload.lead_snp_name = request.form.get("leadsnp", "")

        return new_payload


    def _read_files(self, request: Request, new_payload: SessionPayload) -> SessionPayload:
        """
        Populate the payload with dataframes for the uploaded files.
        Gets file information from the request.
        """
        return new_payload

    def _fix_gwas_columns(self, request: Request, new_payload: SessionPayload) -> SessionPayload:
        """
        Using info in the provided request, check and rename the columns in the GWAS file.
        """
        return new_payload