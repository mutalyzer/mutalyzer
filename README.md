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

### Enable the cache

Create a cache folder and a configuration file.

    mkdir cache
    echo MUTALYZER_CACHE_DIR = \'$(pwd)/cache\' > config.txt

Populate the cache.

    retriever --reference_id NC_000022.11 --parse > cache/NC_000022.11

### Running

Start the backend.

    MUTALYZER_SETTINGS="$(pwd)/config.txt" normalizer-api

Navigate to `http://localhost:5000/api` to interact with the API.

Start the frontend.

    npm run serve

Navigate to `http://localhost:8080` to interact with the website.
