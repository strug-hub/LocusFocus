import pytest

from app import create_app


@pytest.fixture()
def flask_app():
    app = create_app()
    app.config.update(
        {
            "TESTING": True,
            "DISABLE_CACHE": True,
            "DISABLE_CELERY": True,
            "CACHE_TYPE": "NullCache",
            "MONGO_CONNECTION_STRING": "mongodb://localhost:27017",
        }
    )
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
    yield fake
