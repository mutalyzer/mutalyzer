Usage
=====

This package provides a :doc:`command line interface <cli>`.

The Normalizer
--------------

To normalize a description variant:

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

The Mapper
----------

To map `NM_003002.4:c.274G>T` to the `GRCh38` assembly chromosome:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_mapper "NM_003002.4:c.274G>T" --reference-id GRCh38

To map `NM_003002.2:c.274G>T` to the `GRCh38` assembly chromosome using the `NM_003002.4` transcript:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_mapper "NM_003002.4:c.274G>T" --reference-id GRCh38 --selector-id NM_003002.4

For the above note that the sequences are sliced to the transcript's exons and that the variants introduced by the sequences differences are filtered out.

To map `NG_012337.1(NM_003002.2):c.274G>T` to `NG_012337.3(NM_003002.4)`:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_mapper "NG_012337.1(NM_003002.2):c.274G>T" --reference-id NG_012337.3 --selector-id NM_003002.4

To map `NG_012337.1(NM_003002.2):c.274G>T` to `NG_012337.3(NM_003002.4)` with variant filtering:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_mapper "NG_012337.1(NM_003002.2):c.274G>T" --reference-id NG_012337.3 --selector-id NM_003002.4 --filter


.. _key: https://support.nlm.nih.gov/knowledgebase/article/KA-05316/en-us
