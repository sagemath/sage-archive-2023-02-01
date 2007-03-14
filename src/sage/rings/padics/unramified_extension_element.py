import sage.rings.finite_field
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
import sage.rings.padics.padic_extension_generic_element

GF = sage.rings.finite_field.GF
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

class UnramifiedExtensionElement(pAdicExtensionGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        if construct:
            PQRElement.__init__(self, parent, x, check=False)
            return
        if check:
            from sage.rings.padics.unramified_extension_generic import UnramifiedExtensionGeneric
            if not isinstance(parent, UnramifiedExtensionGeneric):
                raise TypeError, "parent must be an UnramifiedExtensionGeneric"
        pAdicExtensionGenericElement.__init__(self, parent, x, absprec, relprec, check, construct)

    def _add_(self, right):
        return UnramifiedExtensionElement(self.parent(), self._polynomial + right._polynomial, check = False, construct = True)

    def __rshift__(self, shift):
        return UnramifiedExtensionElement(self.parent(), self.parent().polynomial_ring()([self._polynomial[i].__rshift__(shift) for i in range(self._polynomial.degree() + 1)]), construct = True)

    def __lshift__(self, shift):
        return UnramifiedExtensionElement(self.parent(), self.parent().polynomial_ring()([self._polynomial[i].__lshift__(shift) for i in range(self._polynomial.degree() + 1)]), construct = True)

    def slice(self, i, j, k = 1):
        raise NotImplementedError

    def _mul_(self, right):
        return UnramifiedExtensionElement(self.parent(), self._polynomial * right._polynomial, check = False, construct = True)

    def _neg_(self):
        return UnramifiedExtensionElement(self.parent(), -self._polynomial, check = False, construct = True)

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
            pprint = self.base_ring()._uniformizer_print()
            selflist = self.list()
            triples = [(selflist[i], pprint, i) for i in range(2, self.precision_absolute())]
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
                s = " + ".join(["%s\\cdot%s^{%s}"%(a) for a in triples]) + " + O(%s^{%s})"%(pprint, self.precision_absolute())
                if self.precision_absolute() > 1 and selflist[1] != 0:
                    s = "%s\\cdot%s + "%(addparen(selflist[1]), pprint) + s
                s = s.replace(" 1\\cdot", " ")
                if s.startswith("1\\cdot"):
                    s = s.replace("1\\cdot", "", 1)
                if self.precision_absolute() > 0 and selflist[0] != 0:
                    s = "%s + "%(addparen(selflist[0])) + s
            else:
                s = " + ".join(["%s*%s^%s"%(a) for a in triples]) + " + O(%s^%s)"%(pprint, self.precision_absolute())
                if self.precision_absolute() > 1 and selflist[1] != 0:
                    s = "%s*%s + "%(addparen(selflist[1]), pprint) + s
                s = s.replace(" 1*", " ")
                if s.startswith("1*"):
                    s = s.replace("1*", "", 1)
                if self.precision_absolute() > 0 and selflist[0] != 0:
                    s = "%s + "%(addparen(selflist[0])) + s
            s = s.replace(" +  + ", " + ")
            return s

    def _sub_(self, right):
        return UnramifiedExtensionElement(self.parent(), self._polynomial - right._polynomial, check = False, construct = True)

    def copy(self):
        return UnramifiedExtensionElement(self.parent(), self._polynomial, check = False, construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def list(self):
        #Need to change this to allow for base_rings that are not Zp.  Also, this is slow.
        answer = []
        me = self
        while me != 0:
            answer.append(me.residue(1))
            me = me >> 1
        zero = self.parent().residue_class_field()(0)
        answer = answer + [zero]*(self.precision_absolute() - len(answer))
        return answer

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def precision_absolute(self):
        if self._polynomial == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity
            else:
                return self.parent().precision_cap()
        return min(c.precision_absolute() for c in self._polynomial.list())

    def residue(self, n):
        if n == 1:
            # will currently only work over Zp and Qp, need to implement conway reduction model for finite fields.
            K = self.parent().residue_class_field()
            a = K.gen()
            polylist = self._polynomial.list()
            apow = K(1)
            ans = K(0)
            for i in range(len(polylist)):
                ans += polylist[i].residue(1) * apow
                apow *= a
            return ans
        else:
            raise NotImplementedError, "Need to think about extensions of finite rings first."

    def valuation(self):
        clist = self._polynomial.list()
        if len(clist) == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity
            else:
                return self.parent().precision_cap()
        return min([a.valuation() for a in clist])

    def _val_unit(self):
        return self.valuation(), self.unit_part()
