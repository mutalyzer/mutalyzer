Usage
=====

This package provides a :doc:`command line interface <cli>`.

To run the Normalizer on a description variant:

.. code-block:: console

    $ mutalyzer_normalizer "NG_012337.1(NM_003002.2):c.274G>T"

Enable the file based cache
---------------------------

Create a cache directory and a configuration file:

.. code-block:: console

    $ mkdir cache
    $ echo MUTALYZER_CACHE_DIR = $(pwd)/cache > config.txt

To ensure that any uncached references encountered during a Mutalyzer run will be added to the cache directory:

.. code-block:: console

    $ echo MUTALYZER_FILE_CACHE_ADD = true >> config.txt

Setup the email address used to communicate with the NCBI:

.. code-block:: console

    $ echo EMAIL = your.email@address.com >> config.txt

Optionally, setup the NCBI API key_:

.. code-block:: console

    $ echo NCBI_API_KEY = your_NCBI_key >> config.txt

Populate the cache:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_retriever --id NG_012337.3 --parse --split --output cache

For GRCh37 and GRCh38 chromosomal references (``NC_``):

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_retriever ncbi_assemblies --ref_id_start NC_ --assembly_id_start GRCh --output cache --include_sequence

Now the tool can be run with the cache:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_normalizer "NC_000022.11(NM_182984.5):c.95del"


.. _key: https://support.nlm.nih.gov/knowledgebase/article/KA-05316/en-us
