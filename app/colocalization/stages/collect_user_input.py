from flask import Request

from app.colocalization.payload import SessionPayload
from app.pipeline import PipelineStage
from app.utils.errors import InvalidUsage


class CollectUserInputStage(PipelineStage):
    """
    Stage 1 of the Colocalization pipeline.

    Reads and validates the provided user input (files, form inputs, etc.),
    and creates a new payload object that will be used for the rest of the
    pipeline chain.
    """

    VALID_POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR", "ASN", "NFE"]
    VALID_COORDINATES = ["hg19", "hg38"]

    def name(self) -> str:
        """
        The name of the stage.
        """
        return "collect-user-input"

    def invoke(self, payload: SessionPayload) -> SessionPayload:

        new_payload = self._read_form_inputs(payload)

        return new_payload

    def _read_form_inputs(self, payload: SessionPayload) -> SessionPayload:
        """
        Populate the payload with form inputs read from the given request.
        """
        # Collect all form errors, check length at end, and then raise exception if > 0
        errors = []

        if payload.request.form.get("coordinate") not in self.VALID_COORDINATES:
            errors.append(
                f"Invalid coordinate: '{payload.request.form.get('coordinate')}'"
            )
        payload.coordinate = payload.request.form.get("coordinate", "hg19")  # type: ignore

        payload.coloc2 = bool(payload.request.form.get("coloc2check"))

        if payload.request.form.get("LD-populations") not in self.VALID_POPULATIONS:
            errors.append(
                f"Invalid 1000 Genomes population: '{payload.request.form.get('LD-populations')}'"
            )
        payload.ld_population = payload.request.form.get("LD-populations", "EUR")

        payload.infer_variant = bool(payload.request.form.get("markerCheckbox"))
        payload.plot_locus = payload.request.form.get("locus", "1:205500000-206000000")
        payload.simple_sum_locus = payload.request.form.get("SSlocus", "")

        # GTEx
        payload.gtex_tissues = payload.request.form.getlist("GTEx-tissues")
        payload.gtex_genes = payload.request.form.getlist("region-genes")
        if len(payload.gtex_tissues) > 0 and len(payload.gtex_genes) == 0:
            errors.append(
                "Please select one or more genes to complement your GTEx tissue(s) selection"
            )
        elif len(payload.gtex_genes) > 0 and len(payload.gtex_tissues) == 0:
            errors.append(
                "Please select one or more tissues to complement your GTEx gene(s) selection"
            )

        # First stage set-based test P value threshold
        p_threshold = payload.request.form.get("setbasedP", "")
        if p_threshold == "":
            p_threshold = "default"
        else:
            try:
                p_threshold = float(p_threshold)
                if p_threshold < 0 or p_threshold > 1:
                    errors.append(
                        "Set-based p-value threshold given is not between 0 and 1"
                    )
                else:
                    payload.set_based_p = p_threshold
            except:
                errors.append(
                    "Invalid value provided for the set-based p-value threshold. Value must be numeric between 0 and 1."
                )

        payload.lead_snp_name = payload.request.form.get("leadsnp", "")

        if len(errors) > 0:
            raise InvalidUsage(
                message=f"Error(s) found in uploaded form",
                payload={f"error_{i+1}": error for i, error in enumerate(errors)},
            )

        return payload
