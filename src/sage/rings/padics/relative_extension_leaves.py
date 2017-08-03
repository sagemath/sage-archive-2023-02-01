"""
Relative extensions of `p`-adic rings
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

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from .generic_nodes import pAdicFixedModRingGeneric
from .eisenstein_extension_generic import EisensteinExtensionGeneric
from .relative_ramified_FM import RelativeRamifiedFixedModElement
from .pow_computer_relative import PowComputer_relative_maker

class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'fixed-mod')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFixedModElement)
        from .relative_ramified_FM import pAdicCoercion_ZZ_FM, pAdicConvert_QQ_FM
        self.register_coercion(pAdicCoercion_ZZ_FM(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FM(self))

class pAdicRelativeBaseringInjection(Morphism):
    def __init__(self, R, S):
        if not R.is_field() or S.is_field():
            Morphism.__init__(self, Hom(R, S, R.category()))
        else:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            Morphism.__init__(self, Hom(R, S, SetsWithPartialMaps()))

    def _call_(self, x):
        return self.codomain()([x])

    def _call_with_extra_args(self, x, args=(), kwds={}):
        return self.codomain()([x], *args, **kwds)
