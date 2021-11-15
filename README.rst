Mutalyzer HGVS Parser
=====================

.. image:: https://img.shields.io/github/last-commit/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser/graphs/commit-activity
.. image:: https://readthedocs.org/projects/mutalyzer-hgvs-parser/badge/?version=latest
   :target: https://mutalyzer-hgvs-parser.readthedocs.io/en/latest
.. image:: https://img.shields.io/github/release-date/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser/releases
.. image:: https://img.shields.io/github/release/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser/releases
.. image:: https://img.shields.io/pypi/v/mutalyzer-hgvs-parser.svg
   :target: https://pypi.org/project/mutalyzer-hgvs-parser/
.. image:: https://img.shields.io/github/languages/code-size/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser
.. image:: https://img.shields.io/github/languages/count/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser
.. image:: https://img.shields.io/github/languages/top/mutalyzer/hgvs-parser.svg
   :target: https://github.com/mutalyzer/hgvs-parser
.. image:: https://img.shields.io/github/license/mutalyzer/hgvs-parser.svg
   :target: https://raw.githubusercontent.com/mutalyzer/hgvs-parser/master/LICENSE.md

----

Package to syntax check and convert Mutalyzer HGVS variant descriptions into
a dictionary model to easily access descriptions information in a
programmatically manner.

**Features:**

- Accepts HGVS descriptions with multiple variants (one HGVS allele).
- Any description sub-part can be parsed and converted as well.
- Supports common deviations to the HGVS guidelines.
- Command line and library interfaces available.


Quick start
-----------

Parse and convert a description from the command line:

.. code-block:: console

    $ mutalyzer_hgvs_parser -c "NG_012337.1:c.20del"
    {
      "reference": {
        "id": "NG_012337.1"
      },
      "coordinate_system": "c",
      "variants": [
        {
          "location": {
            "type": "point",
            "position": 20
          },
          "type": "deletion",
          "source": "reference"
        }
      ]
    }


The ``to_model()`` function can be used for the same purpose:

.. code:: python

    >>> from mutalyzer_hgvs_parser import to_model
    >>> model = to_model("NG_012337.1:c.20del")
    >>> model['reference']
    {'id': 'NG_012337.1'}


Please see ReadTheDocs_ for the latest documentation.

.. _ReadTheDocs: https://mutalyzer-hgvs-parser.readthedocs.io/en/latest/
