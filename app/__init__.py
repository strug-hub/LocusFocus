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

DEFAULT_CONFIG = ProdConfig() if BaseConfig.APP_ENV == "production" else DevConfig()

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


def create_app(config=DEFAULT_CONFIG):
    """
    Create an instance of a Flask app for LocusFocus.
    """
    app = Flask(__name__, instance_relative_config=False)

    #    Initialize debugger (for attaching to later with vscode).
    #    Note that hot reloading creates a new server process and will detach debugger if active.
    #    Loading the app a second time while the server is running (e.g., scripts)
    #    will call this line again and raise address conflict, which we swallow
    #    Note also that this will signifcantly slow down the app
    if config.FLASK_APP_DEBUG:
        import debugpy

        try:
            debugpy.listen(("0.0.0.0", 5678))
        except RuntimeError:
            pass

    app.config.from_object(config)
    if app.config["SECRET_KEY"] is None or app.config["SECRET_KEY"] == "":
        raise RuntimeError("SECRET_KEY is not set! Add FLASK_SECRET_KEY to environment!")

    ext.init_app(app)
    talisman.init_app(app, content_security_policy=app.config["CSP_POLICY"])

    from app.utils.gtex_db import FakeGTExDatabase, NullGTExDatabase, RealGTExDatabase

    is_production = app.config.get("APP_ENV") == "production"

    if app.config.get("MONGO_URI"):
        mongo.init_app(app)
        try:
            app.logger.debug("MongoDB connection test")
            _version = mongo.cx.server_info().get("version")
            app.logger.debug(f"Connected to MongoDB {_version}")
            app.extensions["gtex_db"] = RealGTExDatabase(mongo.cx)
        except Exception as e:
            if is_production:
                raise RuntimeError("MongoDB connection failed in production") from e
            app.logger.error(f"MongoDB connection failed: {e}")
            app.logger.warning("Falling back to FakeGTExDatabase (synthetic GTEx data)")
            app.extensions["gtex_db"] = FakeGTExDatabase()
    else:
        if is_production:
            raise RuntimeError("MONGO_CONNECTION_STRING is required in production")
        app.logger.warning("No MONGO_CONNECTION_STRING set; using FakeGTExDatabase (synthetic GTEx data)")
        app.extensions["gtex_db"] = FakeGTExDatabase()

    celery_init_app(app)
    cache.init_app(app)

    if app.config["CACHE_TYPE"] == "NullCache":
        app.logger.debug("Cache is disabled")
    else:
        app.logger.debug(f"Cache is enabled: {app.config['CACHE_TYPE']}")

    with app.app_context():
        from app import routes  # noqa: F401
        from app.jobs import routes as jobs_routes

        app.register_blueprint(jobs_routes.jobs_bp)

        return app
