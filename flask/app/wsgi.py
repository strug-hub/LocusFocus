from os import getenv

from app import create_app

app = create_app()

if __name__ == "__main__":
    app.run(port=getenv("FLASK_PORT"), host="0.0.0.0")
