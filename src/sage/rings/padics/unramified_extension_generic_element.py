import sage.rings.finite_field
import sage.rings.polynomial.polynomial_quotient_ring_element
from sage.rings.padics.misc import min
import sage.rings.integer
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.finite_field_element
import sage.rings.polynomial.polynomial_ring_constructor
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_ring_lazy
import sage.rings.polynomial.polynomial_element
import sage.rings.infinity
import sys
import sage.rings.polynomial.polynomial_ring_constructor
import sage.rings.padics.padic_extension_generic_element

GF = sage.rings.finite_field.GF
PolynomialRing = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing
infinity = sage.rings.infinity.infinity
PQRElement = sage.rings.polynomial.polynomial_quotient_ring_element.PolynomialQuotientRingElement
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
is_IntegerMod = sage.rings.integer_mod.is_IntegerMod
is_FiniteFieldElement = sage.rings.finite_field_element.is_FiniteFieldElement
pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
Polynomial = sage.rings.polynomial.polynomial_element.Polynomial
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy
pAdicExtensionGenericElement = sage.rings.padics.padic_extension_generic_element.pAdicExtensionGenericElement

class UnramifiedExtensionGenericElement(pAdicExtensionGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        if construct:
            pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check = False, construct = True)
            return
        if check:
            from sage.rings.padics.unramified_extension_generic import UnramifiedExtensionGeneric
            if not isinstance(parent, UnramifiedExtensionGeneric):
                raise TypeError, "parent must be an UnramifiedExtensionGeneric"
        if is_FiniteFieldElement(x):
            #This currently only works for relative extensions of Zp
            if x.parent().order().is_prime():
                x = parent._PQR.polynomial_ring()(x)
            else:
                x = x.polynomial()
                x = parent._PQR.polynomial_ring()(x)
                #x = x.parent().change_ring(parent.ground_ring())(x)
        pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check, construct)

    def __rshift__(self, shift):
        poly = self._value._polynomial
        return self.parent()._element_class(self.parent(), self.parent()._PQR.polynomial_ring()([poly[i].__rshift__(shift) for i in range(poly.degree() + 1)]), construct = True)

    def __lshift__(self, shift):
        poly = self._value._polynomial
        return self.parent()._element_class(self.parent(), self.parent()._PQR.polynomial_ring()([poly[i].__lshift__(shift) for i in range(poly.degree() + 1)]), construct = True)

    def list(self):
        #Need to change this to allow for base_rings that are not Zp.  Also, this is slow.
        answer = []
        me = self
        while me != 0:
            #print type(me)
            answer.append(self.parent()(me.residue(1), self.parent().precision_cap()))
            me = me >> 1
        zero = self.parent()(0)
        answer = answer + [zero]*(self.precision_absolute() - len(answer))
        return answer

    def precision_absolute(self):
        if self._value == 0:
            if self._is_exact_zero():
                return infinity
            else:
                return self.parent().precision_cap()
        return min([c.precision_absolute() for c in self._value._polynomial.list()])

    def residue(self, n):
        if n == 1:
            # will currently only work over Zp and Qp, need to implement conway reduction model for finite fields.
            K = self.parent().residue_class_field()
            a = K.gen()
            polylist = self._value._polynomial.list()
            apow = K(1)
            ans = K(0)
            for i in range(len(polylist)):
                ans += polylist[i].residue(1) * apow
                apow *= a
            return ans
        else:
            raise NotImplementedError, "Need to think about extensions of finite rings first."

    def valuation(self):
        clist = self._value._polynomial.list()
        if len(clist) == 0:
            if self._is_exact_zero():
                return infinity
            else:
                return self.parent().precision_cap()
        return min([a.valuation() for a in clist])
