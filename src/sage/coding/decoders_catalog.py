r"""
Index of decoders

The ``codes.decoders`` object may be used to access the decoders that Sage can build.

**Generic decoders**

- :class:`linear_code.LinearCodeSyndromeDecoder <sage.coding.linear_code.LinearCodeSyndromeDecoder>`
- :class:`linear_code.LinearCodeNearestNeighborDecoder <sage.coding.linear_code.LinearCodeNearestNeighborDecoder>`

**Generalized Reed-Solomon code decoders**

- :class:`grs.GRSBerlekampWelchDecoder <sage.coding.grs.GRSBerlekampWelchDecoder>`
- :class:`grs.GRSErrorErasureDecoder <sage.coding.grs.GRSErrorErasureDecoder>`
- :class:`grs.GRSGaoDecoder <sage.coding.grs.GRSGaoDecoder>`
- :class:`grs.GRSKeyEquationSyndromeDecoder <sage.coding.grs.GRSKeyEquationSyndromeDecoder>`
- :class:`guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder <sage.coding.guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder>`

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
from guruswami_sudan.gs_decoder import GRSGuruswamiSudanDecoder
from grs import (GRSBerlekampWelchDecoder,
                 GRSGaoDecoder,
                 GRSKeyEquationSyndromeDecoder,
                 GRSErrorErasureDecoder)
