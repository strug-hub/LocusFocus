"""
Flask application factory.
"""

from flask import Flask
from flask_sitemap import Sitemap
from flask_talisman import Talisman
from flask_uploads import UploadSet, DATA, configure_uploads
from .config import ProdConfig


ext = Sitemap()
talisman = Talisman()
files = UploadSet("files", DATA)


def create_app(config_class=ProdConfig):
    """
    Create an instance of a Flask app for LocusFocus.
    """
    app = Flask(__name__, instance_relative_config=False)

    app.config.from_object(config_class())
    if app.config["SECRET_KEY"] is None:
        raise Exception("SECRET_KEY is not set! Add FLASK_SECRET_KEY to environment!")

    ext.init_app(app)
    talisman.init_app(app, content_security_policy=app.config["CSP_POLICY"])
    configure_uploads(app, files)

    with app.app_context():
        from . import routes

        return app
