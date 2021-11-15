Usage
=====

This package provides a :doc:`command line interface <cli>`.

To run the Name Checker on a description variant:

.. code-block:: console

    $ mutalyzer_name_checker "NG_012337.1(NM_003002.2):c.274G>T"

Enable the file based cache
---------------------------

Create a cache directory and a configuration file:

.. code-block:: console

    $ mkdir cache
    $ echo MUTALYZER_CACHE_DIR = \'$(pwd)/cache\' > config.txt

Populate the cache:

.. code-block:: console

    $ mutalyzer_retriever --id NC_000022.11 --parse > cache/NC_000022.11


Now the tool can be run with the cache:

.. code-block:: console

    $ MUTALYZER_SETTINGS="$(pwd)/config.txt" mutalyzer_name_checker "NC_000022.11(NM_182984.5):c.95del"