r"""
Index of decoders

The ``codes.decoders`` object may be used to access the decoders that Sage can build.

**Generic decoders**

- :func:`linear_code.LinearCodeSyndromeDecoder <sage.coding.linear_code.LinearCodeSyndromeDecoder>`
- :func:`linear_code.LinearCodeNearestNeighborDecoder <sage.coding.linear_code.LinearCodeNearestNeighborDecoder>`

**Generalized Reed-Solomon code decoders**

- :func:`guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder <sage.coding.guruswami_sudan.gs_decoder.GRSGuruswamiSudanDecoder>`

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
