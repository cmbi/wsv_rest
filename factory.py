

from flask import Flask


def create_app():
    app = Flask(__name__)
    app.config.from_object('wiws.settings')

    from wiws.interface import bp as api_bp
    app.register_blueprint(api_bp)

    return app
