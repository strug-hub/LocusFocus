"""
Config classes for LocusFocus
"""

import os
from dotenv import load_dotenv


load_dotenv(os.path.dirname(__file__))


class BaseConfig:
    """
    Common configuration options for all environments.
    """

    # Content Security Policy to use with Flask-Talisman
    CSP_POLICY = {
        "default-src": "'self'",
        "img-src": [
            "*",
            "data:",
            "https://*.google-analytics.com",
            "https://*.googletagmanager.com",
        ],
        "script-src": [
            "'self'",
            "https://*.googletagmanager.com",
            "https://cdnjs.cloudflare.com",
            "https://cdn.plot.ly",
            "'unsafe-hashes'",
            "'sha256-hphOYdb9WX9dW4pYcQdXa8E450mGtzl7k4kSIg1GOIo='", # d3.js
            "'sha256-47DEQpj8HBSa+/TImW+5JCeuQeRkm5NMpJWZG3hSuFU='", # mdb.min.js
            "'sha256-OTeu7NEHDo6qutIWo0F2TmYrDhsKWCzrUgGoxxHGJ8o='", # mdb.min.js
        ],
        "style-src": [
            "'self'",
            "https://use.fontawesome.com",
            "https://cdnjs.cloudflare.com",
            "https://stackpath.bootstrapcdn.com",
            "'unsafe-hashes'",
            "'sha256-hphOYdb9WX9dW4pYcQdXa8E450mGtzl7k4kSIg1GOIo='", # d3.js
            "'sha256-47DEQpj8HBSa+/TImW+5JCeuQeRkm5NMpJWZG3hSuFU='", # mdb.min.js
            "'sha256-OTeu7NEHDo6qutIWo0F2TmYrDhsKWCzrUgGoxxHGJ8o='", # mdb.min.js
        ],
        "font-src": [
            "'self'",
            "https://use.fontawesome.com",
            "https://cdnjs.cloudflare.com",
            "https://stackpath.bootstrapcdn.com",
        ],
        "connect-src": [
            "'self'",
            "https://*.google-analytics.com",
            "https://*.analytics.google.com",
            "https://*.googletagmanager.com",
        ],
    }
    NONCE_LIST = ["script-src", "style-src"]

    APP_ENV = os.environ.get("APP_ENV", "production")

    FLASK_APP_DEBUG = os.environ.get("FLASK_APP_DEBUG", "False").lower() == "true"

    UPLOAD_FOLDER = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "static", "upload")
    )
    SESSION_FOLDER = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "static", "session_data")
    )
    # Flask-Uploads
    UPLOADED_FILES_DEST = UPLOAD_FOLDER
    UPLOADED_FILES_ALLOW = set(["txt", "tsv", "ld", "html"])
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 500MB limit

    SEND_FILE_MAX_AGE_DEFAULT = 300  # 5 min cache
    SECRET_KEY = os.environ["FLASK_SECRET_KEY"]

    LF_DATA_FOLDER = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),  # project root below app folder
        "data",
    )

    CELERY = {
        "timezone": "America/Toronto",
        "broker_url": "redis://localhost:6379/0",
        "result_backend": f"file://{os.path.join(SESSION_FOLDER, 'celery_results')}",
    }

    DISABLE_CELERY = os.environ.get("DISABLE_CELERY", "False").lower() == "true"


class DevConfig(BaseConfig):
    """
    Configuration options for development environment only.
    """

    MONGO_URI = os.environ.get("MONGO_CONNECTION_STRING")


class ProdConfig(BaseConfig):
    """
    Configuration options for production environment only.
    """

    MONGO_URI = os.environ.get("MONGO_CONNECTION_STRING")
