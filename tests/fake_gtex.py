"""Re-export FakeGTExDatabase from its canonical location in the app package."""

from app.utils.gtex_db.fake import FakeGTExDatabase, GENES, TISSUES

__all__ = ["FakeGTExDatabase", "GENES", "TISSUES"]
