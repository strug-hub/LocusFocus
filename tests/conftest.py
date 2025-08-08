import pytest
from app import create_app

@pytest.fixture()
def flask_app():
    app = create_app()
    app.config.update({
        "TESTING": True,
        "DISABLE_CACHE": True,
        "DISABLE_CELERY": True,
        "CACHE_TYPE": "NullCache",
        "MONGO_CONNECTION_STRING": "mongodb://localhost:27017",
    })

    # other setup can go here

    yield app

    # clean up / reset resources here
