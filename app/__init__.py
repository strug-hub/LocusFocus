"""
Flask application factory.
"""

from celery.app import Celery
from celery.app.task import Task
from flask import Flask
from flask_pymongo import PyMongo
from flask_sitemap import Sitemap
from flask_talisman import Talisman

from app.config import BaseConfig, DevConfig, ProdConfig
from app.cache import cache

ConfigClass = ProdConfig if BaseConfig.APP_ENV == "production" else DevConfig

ext = Sitemap()
talisman = Talisman()
mongo = PyMongo()


def celery_init_app(app: Flask) -> Celery:
    """
    Create a Celery app instance for LocusFocus.
    """

    # Copied from https://flask.palletsprojects.com/en/stable/patterns/celery/
    # ensures that tasks are executed in the Flask application context
    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app


def create_app(config_class=ConfigClass):
    """
    Create an instance of a Flask app for LocusFocus.
    """
    app = Flask(__name__, instance_relative_config=False)

    #    Initialize debugger (for attaching to later with vscode).
    #    Note that hot reloading creates a new server process and will detach debugger if active.
    #    Loading the app a second time while the server is running (e.g., scripts)
    #    will call this line again and raise address conflict, which we swallow
    #    Note also that this will signifcantly slow down the app
    if ConfigClass.FLASK_APP_DEBUG:
        import debugpy

        try:
            debugpy.listen(("0.0.0.0", 5678))
        except RuntimeError:
            pass

    app.config.from_object(config_class())
    if app.config["SECRET_KEY"] is None or app.config["SECRET_KEY"] == "":
        raise Exception("SECRET_KEY is not set! Add FLASK_SECRET_KEY to environment!")

    ext.init_app(app)
    talisman.init_app(app, content_security_policy=app.config["CSP_POLICY"])
    mongo.init_app(app)
    celery_init_app(app)
    cache.init_app(app)

    with app.app_context():
        from app import routes
        from app.jobs import routes as jobs_routes

        app.register_blueprint(jobs_routes.jobs_bp)

        return app
