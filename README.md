# Normalizer
## Development server
### Installation

Install the backend.

    mkvirtualenv -p /usr/bin/python3 normalizer
    git clone git@github.com:mutalyzer/normalizer.git
    cd normalizer
    git checkout refactor
    pip install -e .

Install the frontend.

    git clone git@git.lumc.nl:mlefter/mutalyzer-visualization-vuetify.git
    cd mutalyzer-visualization-vuetify
    npm install

### Running

Start the backend.

    normalizer-api

Navigate to `http://localhost:5000/api` to interact with the API.

Start the frontend.

    npm run serve

Navigate to `http://localhost:8080` to interact with the website.
