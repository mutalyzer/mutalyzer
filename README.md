# Normalizer
## Development server
### Installation

Install the backend.

    mkvirtualenv -p /usr/bin/python3 normalizer
    git clone https://github.com/mutalyzer/normalizer.git
    cd normalizer
    pip install -e .

Install the frontend.

    git clone https://github.com/mutalyzer/website.git
    cd website
    npm install

### Enable the cache

Create a cache folder and a configuration file.

    mkdir cache
    echo MUTALYZER_CACHE_DIR = \'$(pwd)/cache\' > config.txt

Populate the cache.

    mutalyzer_retriever --id NC_000022.11 --parse > cache/NC_000022.11

### Running

Start the backend.

    MUTALYZER_SETTINGS="$(pwd)/config.txt" normalizer-api

Navigate to `http://localhost:5000/api` to interact with the API.

Start the frontend.

    npm run serve

Navigate to `http://localhost:8080` to interact with the website.
