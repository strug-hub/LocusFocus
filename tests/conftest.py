import pytest
import os

from app import create_app
from app.config import DevConfig


@pytest.fixture()
def flask_app():
    """Create a Flask app for testing."""
    config = DevConfig()
    config.SECRET_KEY = "test"
    config.DISABLE_CACHE = True
    config.DISABLE_CELERY = True
    config.CACHE_TYPE = "NullCache"
    config.MONGO_URI = None
    app = create_app(config=config)
    yield app


@pytest.fixture()
def fake_gtex_db(flask_app):
    """Inject a `FakeGTExDatabase` into the Flask app for the duration of a test.

    No running MongoDB is required.  The fake is installed as
    `app.extensions["gtex_db"]`, replacing the `RealGTExDatabase` that
    `create_app` registers by default.

    Yields the `FakeGTExDatabase` instance so tests can inspect generated
    variants, tissues, and eQTL records.
    """
    from tests.fake_gtex import FakeGTExDatabase

    fake = FakeGTExDatabase()
    flask_app.extensions["gtex_db"] = fake
    print(f"FakeGTExDatabase injected into app.extensions['gtex_db']")
    yield fake
