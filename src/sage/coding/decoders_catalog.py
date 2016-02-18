r"""
Index of decoders

The ``codes.decoders`` object may be used to access the decoders that Sage can build.

**Generic decoders**

- :func:`linear_code.LinearCodeSyndromeDecoder <sage.coding.linear_code.LinearCodeSyndromeDecoder>`
- :func:`linear_code.LinearCodeNearestNeighborDecoder <sage.coding.linear_code.LinearCodeNearestNeighborDecoder>`

**Subfield subcode decoder**

- :class:`subfield_subcode.SubfieldSubcodeOriginalCodeDecoder <sage.coding.subfield_subcode.SubfieldSubcodeOriginalCodeDecoder>`

**Generalized Reed-Solomon code decoders**

- :func:`grs.GRSBerlekampWelchDecoder <sage.coding.grs.GRSBerlekampWelchDecoder>`
- :func:`grs.GRSErrorErasureDecoder <sage.coding.grs.GRSErrorErasureDecoder>`
- :func:`grs.GRSGaoDecoder <sage.coding.grs.GRSGaoDecoder>`
- :func:`grs.GRSKeyEquationSyndromeDecoder <sage.coding.grs.GRSKeyEquationDecoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.decoders_catalog import *
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#                     2015 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from linear_code import (LinearCodeSyndromeDecoder, LinearCodeNearestNeighborDecoder)
from subfield_subcode import SubfieldSubcodeOriginalCodeDecoder
from grs import (GRSBerlekampWelchDecoder,
                 GRSGaoDecoder,
                 GRSKeyEquationSyndromeDecoder,
                 GRSErrorErasureDecoder)
