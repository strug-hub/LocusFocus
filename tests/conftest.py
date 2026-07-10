import pytest
import os

from app import create_app
from app.config import DevConfig


@pytest.fixture(scope="session")
def flask_app():
    """Create a single Flask app shared across all tests.

    Session scope is required because routes.py is cached by Python's module
    system after first import, so routes only register on the first app
    instance created per process.
    """
    config = DevConfig()
    config.SECRET_KEY = "test"
    config.DISABLE_CACHE = True
    config.DISABLE_CELERY = True
    config.CACHE_TYPE = "NullCache"
    config.MONGO_URI = None
    app = create_app(config=config)
    app.debug = True  # disables Talisman's force_https redirect in tests
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

    original = flask_app.extensions.get("gtex_db")
    fake = FakeGTExDatabase()
    flask_app.extensions["gtex_db"] = fake
    print(f"FakeGTExDatabase injected into app.extensions['gtex_db']")
    yield fake
    flask_app.extensions["gtex_db"] = original
