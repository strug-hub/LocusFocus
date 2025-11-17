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
            "cdnjs.cloudflare.com",
            "cdn.datatables.net",
            "cdn.jsdelivr.net",
            "cdn.plot.ly",
            "stackpath.bootstrapcdn.com",
        ],
        "style-src": [
            "'self'",
            "'sha256-OTeu7NEHDo6qutIWo0F2TmYrDhsKWCzrUgGoxxHGJ8o='",  # mdb
            "'sha256-47DEQpj8HBSa+/TImW+5JCeuQeRkm5NMpJWZG3hSuFU='",  # plotly
            "use.fontawesome.com",
            "cdn.jsdelivr.net",
            "cdnjs.cloudflare.com",
            "cdn.plot.ly",
            "stackpath.bootstrapcdn.com",
            "cdn.datatables.net",
        ],
        "font-src": [
            "'self'",
            "use.fontawesome.com",
            "cdnjs.cloudflare.com",
            "stackpath.bootstrapcdn.com",
        ],
        "connect-src": [
            "'self'",
            "cdnjs.cloudflare.com",
            "stackpath.bootstrapcdn.com",
            "https://*.google-analytics.com",
            "https://*.analytics.google.com",
            "https://*.googletagmanager.com",
        ],
    }

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

    REDIS_PASSWORD = os.environ.get("REDIS_PASSWORD", None)

    CELERY = {
        "timezone": "America/Toronto",
        "broker_url": (
            "redis://localhost:6379/0"
            if REDIS_PASSWORD is None
            else f"redis://:{os.environ['REDIS_PASSWORD']}@localhost:6379/0"
        ),
        "result_backend": f"file://{os.path.join(SESSION_FOLDER, 'celery_results')}",
        "broker_connection_retry_on_startup": True,
    }

    DISABLE_CELERY = os.environ.get("DISABLE_CELERY", "False").lower() == "true"

    # Caching
    DISABLE_CACHE = os.environ.get("DISABLE_CACHE", "False").lower() == "true"
    CACHE_TYPE = "NullCache" if DISABLE_CACHE else "FileSystemCache"
    CACHE_DEFAULT_TIMEOUT = 60 * 60 * 24 * 7  # 1 week
    CACHE_DIR = os.path.join(LF_DATA_FOLDER, "cache")
    CACHE_KEY_PREFIX = "locusfocus-"


class DevConfig(BaseConfig):
    """
    Configuration options for development environment only.
    """

    DISABLE_CACHE = os.environ.get("DISABLE_CACHE", "False").lower() == "true"
    MONGO_URI = os.environ.get("MONGO_CONNECTION_STRING")
    CACHE_TYPE = "NullCache" if DISABLE_CACHE else "FileSystemCache"


class ProdConfig(BaseConfig):
    """
    Configuration options for production environment only.
    """

    MONGO_URI = os.environ.get("MONGO_CONNECTION_STRING")
    CACHE_TYPE = "FileSystemCache"
