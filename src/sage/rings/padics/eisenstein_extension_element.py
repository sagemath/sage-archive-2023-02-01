from sage.rings.finite_field import GF
import sage.rings.polynomial_quotient_ring_element
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
import sage.rings.padics.padic_extension_element

PolynomialRing = sage.rings.polynomial_ring_constructor.PolynomialRing
infinity = sage.rings.infinity.infinity
PQRElement = sage.rings.polynomial_quotient_ring_element.PolynomialQuotientRingElement
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
is_IntegerMod = sage.rings.integer_mod.is_IntegerMod
is_FiniteFieldElement = sage.rings.finite_field_element.is_FiniteFieldElement
pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
Polynomial = sage.rings.polynomial_element.Polynomial
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy
pAdicExtensionGenericElement = sage.rings.padics.padic_extension_generic_element.pAdicExtensionGenericElement

class EisensteinExtensionElement(pAdicExtensionGenericElement):
    def __init__(self, parent, x, absprec = None, relprec = None, check = True, construct = False):
        if construct:
            PQRElement.__init__(self, parent, x, check=False)
            return
        if check:
            from sage.rings.padics.eisenstein_extension_generic import EisensteinExtensionGeneric
            if not isinstance(parent, EisensteinExtensionGeneric):
                raise TypeError, "parent must be an EisensteinExtensionGeneric"
        pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check, construct)

    def _add_(self, right):
        return EisensteinExtensionElement(self.parent(), self._polynomial + right._polynomial, check = False, construct = True)

    def __rshift__(self, shift):
        shift = Integer(shift)
        if shift < 0:
            return self.__lshift__(-shift)
        #Should speed this up later
        if shift == 0:
            return self
        else:
            a = self.list()
            alpha = a.pop(0) >> 1
            returnval = [c + alpha*d for (c, d) in zip(a, self.parent()._shift1list)]
            returnval.append(-alpha)
            return EisensteniExtensionElement(self.parent(), self.parent().polynomial_ring()(returnval), construct = True).__rshift__(shift - 1)

            #f = self.parent().defining_polynomial()
            #a = self._polynomial
            #f0inv = self.base_ring()(~(f[0].unit_part()))
            #a0sh = a[0].__rshift__(1)
            #return EisensteinExtensionElement(self.parent(), self.parent().polynomial_ring()([a[i + 1] - f[i + 1] * f0inv * a0sh for i in range(f.degree())]), construct = True)

    def __lshift__(self, shift):
        shift = Integer(shift)
        if shift < 0:
            return self.__rshift__(-shift)
        return self * self.parent()(self.parent().polynomial_ring().gen()**shift)

    def slice(self, i, j, k = 1):
        raise NotImplementedError

    def _mul_(self, right):
        return EisensteinExtensionElement(self.parent(), self._polynomial * right._polynomial, check = False, construct = True)

    def _neg_(self):
        return EisensteinExtensionElement(self.parent(), -self._polynomial, check = False, construct = True)

    def _repr_(self, mode = None, do_latex = False, caprel = False):
        if mode is None:
            mode = self.parent().print_mode()
        elif not ((mode == 'val-unit') or (mode == 'series') or (mode == 'val-unit-p') or (mode == 'series-p') or (mode == 'integer') or (mode == 'integer-p')):
            raise TypeError, "printing-mode must be one of 'val-unit', 'series', 'integer', 'val-unit-p', 'series-p' or 'integer-p'"
        if self._polynomial == 0:
            if mode == 'val-unit' or mode == 'series' or mode == 'integer':
                if do_latex:
                    return "O(%s^{%s})"%(self.parent().prime(), self.precision_absolute())
                else:
                    return "O(%s^%s)"%(self.parent().prime(), self.precision_absolute())
            else:
                if do_latex:
                    return "O(p^{%s})"%(self.precision_absolute())
                else:
                    return "O(p^%s)"%(self.precision_absolute())
        if mode == 'val-unit':
            return "set the print mode to series"
        elif mode == 'val-unit-p':
            return "set the print mode to series"
        elif mode == 'integer':
            return "set the print mode to series"
        elif mode == 'integer-p':
            return "set the print mode to series"
        else:
            if mode == 'series':
                p = self.parent()._uniformizer_sym(do_latex)
            else:
                p = "p"
            selflist = self.list()
            triples = [(selflist[i], p, i) for i in range(2, self.precision_absolute())]
            #print triples
            triples = [a for a in triples if a[0] != 0] #need to change this later to account for __cmp__ throwing error on lazies
            def addparen(c):
                c = "%s"%c
                if c.find("+") != -1:
                    c = "(%s)"%c
                return c
            triples = [(addparen(a[0]), a[1], a[2]) for a in triples]
            #return "bad"
            if do_latex:
                s = " + ".join(["%s\\cdot%s^{%s}"%(a) for a in triples]) + " + O(%s^{%s})"%(p, self.precision_absolute())
                if self.precision_absolute() > 1 and selflist[1] != 0:
                    s = "%s\\cdot%s + "%(addparen(selflist[1]), p) + s
                s = s.replace(" 1\\cdot", " ")
                if s.startswith("1\\cdot"):
                    s = s.replace("1\\cdot", "", 1)
                if self.precision_absolute() > 0 and selflist[0] != 0:
                    s = "%s + "%(addparen(selflist[0])) + s
            else:
                s = " + ".join(["%s*%s^%s"%(a) for a in triples]) + " + O(%s^%s)"%(p, self.precision_absolute())
                if self.precision_absolute() > 1 and selflist[1] != 0:
                    s = "%s*%s + "%(addparen(selflist[1]), p) + s
                s = s.replace(" 1*", " ")
                if s.startswith("1*"):
                    s = s.replace("1*", "", 1)
                if self.precision_absolute() > 0 and selflist[0] != 0:
                    s = "%s + "%(addparen(selflist[0])) + s
            s = s.replace(" +  + ", " + ")
            return s

    def _sub_(self, right):
        return EisensteinExtensionElement(self.parent(), self._polynomial - right._polynomial, check = False, construct = True)

    def copy(self):
        return EisensteinExtensionElement(self.parent(), self._polynomial, check = False, construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def list(self):
        answer = []
        a = self
        while a != 0:
            answer.append(a.residue(1))
            a = a >> 1
        if self.parent().is_field():
            return answer[min([i for i in range(len(answer)) if answer[i] != 0]):]
        else:
            return answer

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def precision_absolute(self):
        if self._polynomial == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity #is this really what we want to do?  Should we eliminate leading zeros from p-adic polynomials
            else:
                return self.parent().precision_cap()
        clist = self._polynomial.list()
        e = self.parent().e()
        return min(clist[i].precision_absolute()*e + i for i in range(self._polynomial.degree()))

    def valuation(self):
        clist = self._polynomial.list()
        if len(clist) == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity #see comment above
            else:
                return self.parent().precision_cap()
        return min(clist[i].valuation()*e + i for i in range(self._polynomial.degree()))

    def _val_unit(self):
        return self.valuation(), self.unit_part()
