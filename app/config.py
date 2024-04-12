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
            "'unsafe-inline'",  # TODO
            "'unsafe-eval'",  # TODO
            "https://*.googletagmanager.com",
            "cdnjs.cloudflare.com",
            "cdn.plot.ly",
        ],
        "style-src": [
            "'self'",
            "'unsafe-inline'",  # TODO
            "use.fontawesome.com",
            "cdnjs.cloudflare.com",
            "stackpath.bootstrapcdn.com",
        ],
        "font-src": [
            "'self'",
            "use.fontawesome.com",
            "cdnjs.cloudflare.com",
            "stackpath.bootstrapcdn.com",
        ],
        "connect-src": [
            "'self'",
            "https://*.google-analytics.com",
            "https://*.analytics.google.com",
            "https://*.googletagmanager.com",
        ],
    }
    UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), "static", "upload")
    # Flask-Uploads
    UPLOADED_FILES_DEST = UPLOAD_FOLDER
    UPLOADED_FILES_ALLOW = set(['txt', 'tsv', 'ld', 'html'])
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 500MB limit

    SEND_FILE_MAX_AGE_DEFAULT = 300  # 5 min cache
    SECRET_KEY = os.environ.get("FLASK_SECRET_KEY", "")

    LF_DATA_FOLDER = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),  # project root below app folder
            "data"
        )


class DevConfig(BaseConfig):
    """
    Configuration options for development environment only.
    """

    pass


class ProdConfig(BaseConfig):
    """
    Configuration options for production environment only.
    """

    pass
