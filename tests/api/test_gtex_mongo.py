import pytest
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError
from flask.app import Flask


def skip_if_no_mongo(flask_app):
    """Skip the test if mongo is not available"""
    with flask_app.app_context():
        try:
            client = MongoClient(flask_app.config["MONGO_CONNECTION_STRING"], serverSelectionTimeoutMS=10000)
            client.server_info()
        except ServerSelectionTimeoutError:
            pytest.skip("MongoDB not available")

def test_get_gtex(flask_app: Flask):
    """Sanity check for get_gtex"""
    with flask_app.app_context():
        print("checking for mongo")
        skip_if_no_mongo(flask_app)

        from app.utils.gtex import get_gtex

        version = "V8"
        tissue = "Liver"
        gene = "NUCKS1"
        print("testing get_gtex")
        results = get_gtex(version, tissue, gene)
        print("got results")
        assert results.shape[0] > 0
