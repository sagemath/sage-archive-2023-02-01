import sage.rings.finite_field
import sage.rings.polynomial_quotient_ring_element
from sage.rings.padics.misc import min
import sage.rings.integer
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.finite_field_element
import sage.rings.polynomial_ring_constructor
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_ring_lazy
import sage.rings.polynomial_element
import sage.rings.infinity
import sys
import sage.rings.polynomial_ring_constructor
import sage.rings.padics.padic_extension_generic_element

GF = sage.rings.finite_field.GF
PolynomialRing = sage.rings.polynomial_ring_constructor.PolynomialRing
infinity = sage.rings.infinity.infinity
PQRElement = sage.rings.polynomial_quotient_ring_element.PolynomialQuotientRingElement
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
is_IntegerMod = sage.rings.integer_mod.is_IntegerMod
is_FiniteFieldElement = sage.rings.finite_field_element.is_FiniteFieldElement
pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
Polynomial = sage.rings.polynomial_element.Polynomial
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy
pAdicExtensionGenericElement = sage.rings.padics.padic_extension_generic_element.pAdicExtensionGenericElement

class EisensteinExtensionGenericElement(pAdicExtensionGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        if construct:
            pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check = False, construct = True)
            return
        if check:
            from sage.rings.padics.eisenstein_extension_generic import EisensteinExtensionGeneric
            if not isinstance(parent, EisensteinExtensionGeneric):
                raise TypeError, "parent must be an EisensteinExtensionGeneric"
        pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check, construct)

    def __rshift__(self, shift, mode = 'simple'):
        shift = Integer(shift)
        if shift < 0:
            return self.__lshift__(-shift)
        #Should speed this up later
        if shift == 0:
            return self
        if self.parent().is_field():
            return self / self.parent()._element_class(self.parent(), self.parent()._PQR.polynomial_ring().gen()**shift, construct = True)
        else:
            #a = self._polynomial.list()
            #a = a + [self.base_ring()(0)] * (self.parent().degree() - len(a))
            #alpha = a.pop(0) >> 1
            #returnval = [c + alpha*d for (c, d) in zip(a, self.parent()._shift1list)]
            #returnval.append(alpha)


            #PQR = self.parent()._PQR
            #poly_ring_field = PQR.polynomial_ring().base_extend(PQR.polynomial_ring().base_ring().fraction_field())
            #field_PQR = PolynomialQuotientRing(poly_ring_field, poly_ring_field(PQR.modulus()), PQR.variable_name())
            #negnoconst = poly_ring_field(self._value._polynomial) - self.parent().ground_ring().fraction_field()(self._value._polynomial[0].residue(1), self.parent().ground_ring().precision_cap())
            #fieldPQRquot = field_PQR(negnoconst) / field_PQR.gen()
            #returnval = PQR.polynomial_ring()(fieldPQRquot._polynomial)
            #return self.parent()._element_class(self.parent(), returnval, construct = True).__rshift__(shift - 1)

            f = self.parent().defining_polynomial()
            a = self._value._polynomial
            f0inv = self.parent().ground_ring()(~(f[0].unit_part()))
            a0sh = a[0].__rshift__(1)
            return self.parent()._element_class(self.parent(), self.parent().polynomial_ring()([a[i + 1] - f[i + 1] * f0inv * a0sh for i in range(f.degree())]), construct = True).__rshift__(shift - 1)

    def __lshift__(self, shift, mode = 'simple'):
        shift = Integer(shift)
        if shift < 0:
            return self.__rshift__(-shift)
        return self * self.parent()(self.parent().polynomial_ring().gen()**shift)

    def list(self):
        #Note: if you change this implementation you may need to fix eisenstein_extension_capped_relative_element's list function.
        if self.valuation() < 0:
            return self.unit_part().list()
        answer = []
        a = self
        while a != 0:
            answer.append(self.parent()(a.residue(1), self.parent().precision_cap()))
            a = a >> 1
        if self.parent().is_field():
            return answer[min([i for i in range(len(answer)) if answer[i] != 0]):]
        else:
            return answer

    def residue(self, prec):
        if prec == 1:
            return self._value._polynomial[0].residue(1)
        else:
            raise NotImplementedError

    def precision_absolute(self):
        if self._value._polynomial == 0:
            if self._is_exact_zero():
                return infinity #is this really what we want to do?  Should we eliminate leading zeros from p-adic polynomials
            else:
                return self.parent().precision_cap()
        clist = self._value._polynomial.list()
        e = self.parent().e()
        return min([clist[i].precision_absolute()*e + i for i in range(self._value._polynomial.degree() + 1)])

    def valuation(self):
        clist = self._value._polynomial.list()
        e = self.parent().e()
        if len(clist) == 0:
            if self._is_exact_zero():
                return infinity #see comment above
            else:
                return self.parent().precision_cap()
        return min([clist[i].valuation()*e + i for i in range(len(clist))])
