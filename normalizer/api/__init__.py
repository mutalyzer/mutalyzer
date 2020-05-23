from flask import Flask, render_template
from flask_cors import CORS
from .api import blueprint


app = Flask(__name__)
app.register_blueprint(blueprint, url_prefix='/api')

# enable CORS
CORS(app, resources={r'/*': {'origins': '*'}})


@app.route('/')
def index():
    return render_template('index.html')
