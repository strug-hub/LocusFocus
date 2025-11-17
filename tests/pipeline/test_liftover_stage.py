from unittest.mock import Mock

import numpy as np
import pandas as pd


def test_update_indices_kept(app):
    """Test that the stage is accurately adding unlifted indices to payload.gwas_indices_kept"""
    with app.app_context():
        from app.colocalization.payload import SessionPayload
        from app.colocalization.stages.liftover_gwas_file import LiftoverGWASFile

        mock_payload = Mock(spec=SessionPayload)
        stage = LiftoverGWASFile()

        mock_payload.gwas_data_original = pd.DataFrame(np.arange(20).reshape((10, 2)))
        dropped_indices = [1, 3, 8]
        mock_payload.gwas_indices_kept = pd.Series(
            np.ones(len(mock_payload.gwas_data_original), dtype=bool)
        )
        mock_payload.gwas_indices_kept[dropped_indices] = False
        unlifted_over_indices = [1, 2, 9]
        indices_kept = stage.update_indices_kept(
            payload=mock_payload, unlifted_over=unlifted_over_indices
        )

        # Assert that positions corresponding to `dropped_indices` and `unlifted_over` indices are False
        assert list(mock_payload.gwas_data_original[~indices_kept].index) == list(
            set(dropped_indices + unlifted_over_indices)
        )
