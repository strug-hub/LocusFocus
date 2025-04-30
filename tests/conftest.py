import pytest
from app import create_app

@pytest.fixture()
def app():
    app = create_app()
    app.config.update({
        "TESTING": True,
        "DISABLE_CACHE": True,
        "DISABLE_CELERY": True,
        "CACHE_TYPE": "NullCache",
    })

    # other setup can go here

    yield app

    # clean up / reset resources here
