from pow_computer_ext import PowComputer_ext_maker
from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX

from unramified_extension_generic import UnramifiedExtensionGeneric
from eisenstein_extension_generic import EisensteinExtensionGeneric
from padic_general_extension_generic import pAdicGeneralExtensionGeneric

from padic_capped_relative_ring_generic import pAdicCappedRelativeRingGeneric
from padic_capped_relative_field_generic import pAdicCappedRelativeFieldGeneric
from padic_lazy_ring_generic import pAdicLazyRingGeneric
from padic_lazy_field_generic import pAdicLazyFieldGeneric
from padic_capped_absolute_ring_generic import pAdicCappedAbsoluteRingGeneric
from padic_fixed_mod_ring_generic import pAdicFixedModRingGeneric

from padic_ring_fixed_mod import pAdicRingFixedMod

from unramified_extension_absolute_element import UnramifiedExtensionAbsoluteElement
from unramified_extension_capped_relative_element import UnramifiedExtensionCappedRelativeElement
from unramified_extension_lazy_element import UnramifiedExtensionLazyElement
from eisenstein_extension_absolute_element import EisensteinExtensionAbsoluteElement
from eisenstein_extension_capped_relative_element import EisensteinExtensionCappedRelativeElement
from eisenstein_extension_lazy_element import EisensteinExtensionLazyElement
from padic_general_extension_absolute_element import pAdicGeneralExtensionAbsoluteElement
from padic_general_extension_capped_relative_element import pAdicGeneralExtensionCappedRelativeElement
from padic_general_extension_lazy_element import pAdicGeneralExtensionLazyElement

from unramified_fixed_mod_element import UnramifiedFixedModElement
from eisenstein_fixed_mod_element import EisensteinFixedModElement

from sage.rings.integer_ring import ZZ

class UnramifiedExtensionRingCappedRelative(UnramifiedExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedExtensionCappedRelativeElement)

class UnramifiedExtensionFieldCappedRelative(UnramifiedExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedExtensionCappedRelativeElement)

class UnramifiedExtensionRingLazy(UnramifiedExtensionGeneric, pAdicLazyRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedExtensionLazyElement)
        pAdicLazyRingGeneric.__init__(self, poly.base_ring().prime(), prec, print_mode, names, halt)

class UnramifiedExtensionFieldLazy(UnramifiedExtensionGeneric, pAdicLazyFieldGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedExtensionLazyElement)
        pAdicLazyFieldGeneric.__init__(self, poly.base_ring().prime(), prec, print_mode, names, halt)

class UnramifiedExtensionRingCappedAbsolute(UnramifiedExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedExtensionAbsoluteElement)

class UnramifiedExtensionRingFixedMod(UnramifiedExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), max(min(prec - 1, 30), 1), prec, False, ntl_poly, "FM", "u")
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, UnramifiedFixedModElement)

    def coerce_map_explicit(self, S):
        from sage.rings.padics.morphism import Morphism_ZZ_UnrFM, Morphism_ZpFM_UnrFM
        if S is ZZ:
            return Morphism_ZZ_UnrFM(self)
        elif isinstance(S, pAdicRingFixedMod) and S.prime() == self.prime():
            return Morphism_ZpFM_UnrFM(S, self)
        return None

class EisensteinExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinExtensionCappedRelativeElement)

class EisensteinExtensionFieldCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinExtensionCappedRelativeElement)

class EisensteinExtensionRingLazy(EisensteinExtensionGeneric, pAdicLazyRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinExtensionLazyElement)
        pAdicLazyRingGeneric.__init__(self, poly.base_ring().prime(), prec, print_mode, names, halt)

class EisensteinExtensionFieldLazy(EisensteinExtensionGeneric, pAdicLazyFieldGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinExtensionLazyElement)
        pAdicLazyFieldGeneric.__init__(self, poly.base_ring().prime(), prec, print_mode, names, halt)

class EisensteinExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinExtensionAbsoluteElement)

class EisensteinExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, poly, prec, halt, print_mode, names):
        prec = prec // poly.degree()
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        #print poly.base_ring().prime(), prec, poly.degree(), ntl_poly
        # deal with prec not a multiple of e better.
        self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), max(min(prec - 1, 30), 1), prec, False, ntl_poly, "FM", "e")
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, EisensteinFixedModElement)

    def coerce_map_explicit(self, S):
        from sage.rings.padics.morphism import Morphism_ZZ_EisFM, Morphism_ZpFM_EisFM
        if S is ZZ:
            return Morphism_ZZ_EisFM(self)
        elif isinstance(S, pAdicRingFixedMod) and S.prime() == self.prime():
            return Morphism_ZpFM_EisFM(S, self)
        return None

    # Temp fix
    def gen(self, n=0):
        if n != 0:
            raise IndexError, "extension has only one generator"
        #this won't work for padic_general extensions
        return self([0,1])


class pAdicGeneralExtensionRingCappedRelative(pAdicGeneralExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionCappedRelativeElement)

class pAdicGeneralExtensionFieldCappedRelative(pAdicGeneralExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionCappedRelativeElement)

class pAdicGeneralExtensionRingLazy(pAdicGeneralExtensionGeneric, pAdicLazyRingGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionLazyElement)
        pAdicLazyRingGeneric.__init__(self, upoly.base_ring().prime(), prec, print_mode, names, halt)

class pAdicGeneralExtensionFieldLazy(pAdicGeneralExtensionGeneric, pAdicLazyFieldGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionLazyElement)
        pAdicLazyFieldGeneric.__init__(self, upoly.base_ring().prime(), prec, print_mode, names, halt)

class pAdicGeneralExtensionRingCappedAbsolute(pAdicGeneralExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionAbsoluteElement)

class pAdicGeneralExtensionRingFixedMod(pAdicGeneralExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, upoly, epoly, poly, prec, halt, print_mode, names):
        pAdicGeneralExtensionGeneric.__init__(self, upoly, epoly, poly, prec, print_mode, names, pAdicGeneralExtensionAbsoluteElement)
