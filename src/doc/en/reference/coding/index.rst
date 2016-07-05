.. _sage.coding:

Coding Theory
=============

Coding Theory objects
---------------------


.. toctree::
   :maxdepth: 2

   Linear code and its abstract base class <sage/coding/linear_code>
   sage/coding/channel_constructions
   sage/coding/decoder
   sage/coding/encoder

Catalogs
--------

.. toctree::
   :maxdepth: 2

   sage/coding/channels_catalog
   sage/coding/codes_catalog
   sage/coding/decoders_catalog
   sage/coding/encoders_catalog

Code constructions
------------------


.. toctree::
   :maxdepth: 2

   Constructing a linear code given a generator matrix <sage/coding/linear_code>


The named code families below are represented in Sage by their own classes,
allowing specialised implementations of e.g. decoding or computation of properties:

.. toctree::
   :maxdepth: 2

   sage/coding/extended_code
   sage/coding/grs
   sage/coding/hamming_code
   sage/coding/punctured_code
   sage/coding/reed_muller_code

In contrast, for some code families Sage can only construct their generator
matrix and has no other a priori knowledge on them:

.. toctree::
   :maxdepth: 2

   Code families whose generator matrix can be constructed <sage/coding/code_constructions>
   Generator matrices constructed using GAP-Guava <sage/coding/guava>
   sage/coding/two_weight_db
   sage/coding/sd_codes
   sage/coding/binary_code

Methods and Operations related to Linear Codes
----------------------------------------------


.. toctree::
   :maxdepth: 2

   Bounds on linear codes <sage/coding/bounds_catalog>
   sage/coding/delsarte_bounds
   sage/coding/codecan/codecan
   sage/coding/codecan/autgroup_can_label

Source coding
-------------

.. toctree::
   :maxdepth: 1

   sage/coding/source_coding/huffman

.. include:: ../footer.txt
