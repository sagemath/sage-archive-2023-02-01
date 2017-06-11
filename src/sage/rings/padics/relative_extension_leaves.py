"""
Relative extensions of p-adic fields.
"""

#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .generic_nodes import pAdicFixedModRingGeneric
from .eisenstein_extension_generic import EisensteinExtensionGeneric

class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, shift_seed, names):
        pass
