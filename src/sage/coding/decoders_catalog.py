r"""
Index of decoders

The ``codes.decoders`` object may be used to access the decoders that Sage can build.

**Generic decoders**

- :func:`linear_code.LinearCodeSyndromeDecoder <sage.coding.linear_code.LinearCodeSyndromeDecoder>`
- :func:`linear_code.LinearCodeNearestNeighborDecoder <sage.coding.linear_code.LinearCodeNearestNeighborDecoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.decoders_catalog import *
"""

from linear_code import (LinearCodeSyndromeDecoder, LinearCodeNearestNeighborDecoder)
