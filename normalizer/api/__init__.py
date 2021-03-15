from flask import Flask
from flask_cors import CORS

from .api import blueprint

app = Flask(__name__)

app.register_blueprint(blueprint, url_prefix="/api")

CORS(app, resources={r"/*": {"origins": "*"}})
