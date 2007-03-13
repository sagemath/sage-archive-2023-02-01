import sage.rings.polynomial_element
import sage.rings.padics.padic_ring_capped_relative
import sage.rings.padics.eisenstein_extension_generic
import sage.rings.padics.padic_ring_extension_generic
import sage.rings.padics.eisenstein_extension_capped_relative_element
import sage.rings.padics.padic_ring_generic
import sage.rings.infinity

infinity = sage.rings.infinity.infinity
Polynomial = sage.rings.polynomial_element.Polynomial
pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
EisensteinExtensionGeneric = sage.rings.padics.eisenstein_extension_generic.EisensteinExtensionGeneric
pAdicRingExtensionGeneric = sage.rings.padics.padic_ring_extension_generic.pAdicRingExtensionGeneric
EisensteinExtensionCappedRelativeElement = sage.rings.padics.eisenstein_extension_capped_relative_element.EisensteinExtensionCappedRelativeElement
pAdicRingGeneric = sage.rings.padics.padic_ring_generic.pAdicRingGeneric

class EisensteinRingExtensionCappedRelative(EisensteinExtensionGeneric, pAdicRingExtensionGeneric):
    r"""
    Eisenstein Extension of a p-adic ring with capped relative precision

    You should not create this class directly unless you know what you're doing.  Use ExtensionFactory.
    """

    def __init__(self, base, poly, names, prec, halt, print_mode):
        EisensteinExtensionGeneric.__init__(self, poly, names, prec, print_mode)

    def __call__(self, x, absprec = infinity, relprec = infinity):
        return EisensteinExtensionCappedRelativeElement(self, x, absprec, relprec)

    def gen(self, n = 0):
        if n == 0:
            try:
                return self._gen
            except AttributeError:
                self._gen = EisensteinExtensionCappedRelativeElement(self, self.polynomial_ring().gen(), check = False, construct = True)
                return self._gen
