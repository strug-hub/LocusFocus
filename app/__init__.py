"""
Flask application factory.
"""

from celery import Celery, Task
from flask import Flask
from flask_pymongo import PyMongo
from flask_sitemap import Sitemap
from flask_talisman import Talisman
from .config import ProdConfig


ext = Sitemap()
talisman = Talisman()
mongo = PyMongo()


def celery_init_app(app: Flask) -> Celery:
    """
    Create a Celery app instance for LocusFocus.
    """
    # Copied from https://flask.palletsprojects.com/en/stable/patterns/celery/
    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app


def create_app(config_class=ProdConfig):
    """
    Create an instance of a Flask app for LocusFocus.
    """
    app = Flask(__name__, instance_relative_config=False)

    app.config.from_object(config_class())
    if app.config["SECRET_KEY"] is None or app.config["SECRET_KEY"] == "":
        raise Exception("SECRET_KEY is not set! Add FLASK_SECRET_KEY to environment!")

    ext.init_app(app)
    talisman.init_app(app, content_security_policy=app.config["CSP_POLICY"])
    mongo.init_app(app)
    celery_init_app(app)

    with app.app_context():
        from . import routes

        return app
