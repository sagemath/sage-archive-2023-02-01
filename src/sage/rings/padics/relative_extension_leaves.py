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
from .relative_ramified_FM import RelativeRamifiedFixedModElement
from .pow_computer_relative import PowComputer_relative_maker

class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, shift_seed, names):
        unram_prec = (prec + poly.degree() - 1) // poly.degree()
        self.prime_pow = PowComputer_relative_maker(poly.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, poly, 'fixed-mod')
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, RelativeRamifiedFixedModElement)
