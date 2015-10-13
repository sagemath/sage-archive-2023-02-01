r"""
Index of encoders

The ``codes.encoders`` object may be used to access the encoders that Sage can build.

:class:`linear_code.LinearCodeGeneratorMatrixEncoder <sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.encoders_catalog import *
"""

from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import('sage.coding.linear_code', 'LinearCodeGeneratorMatrixEncoder')
