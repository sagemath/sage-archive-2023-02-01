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

class UnramifiedRingExtensionElement(PQRElement, pAdicRingGenericElement):
    def __init__(self, parent, x, absprec = None, relprec = None, check = True, construct = False):
        if construct:
            PQRElement.__init__(self, parent, x, check=False)
            return
        if check:
            from sage.rings.padics.unramified_ring_extension import UnramifiedRingExtension
            if not isinstance(parent, UnramifiedRingExtension):
                raise TypeError, "parent must be an UnramifiedRingExtension"
        if isinstance(x, PQRElement):
            if check:
                if not (x.parent() is parent):
                    raise NotImplementedError, "we have not yet implemented coercion between different unramified extensions"
            PQRElement.__init__(self, parent, x._polynomial, check = False)
            return
        if isinstance(x, (int, long, Integer, Rational)) or is_IntegerMod(x) or is_FiniteFieldElement(x) or (isinstance(x, sage.structure.element.Element) and (isinstance(x.parent(), pAdicRingBaseGeneric) or x.parent() is parent.base_ring())):
            x = parent.base_ring()(x, absprec = absprec)
            PQRElement.__init__(self, parent, sage.rings.polynomial_ring_constructor.PolynomialRing(parent.base_ring(), parent.polynomial_ring().variable_name())(x), check=False)
            return
        if isinstance(x, list):
            x = [parent.base_ring()(v, absprec = absprec) for v in x]
            PQRElement.__init__(self, parent, sage.rings.polynomial_ring_constructor.PolynomialRing(parent.base_ring(), parent.polynomial_ring().variable_name())(x), check=False)
        if isinstance(x, Polynomial):
            if x in parent.polynomial_ring():
                PQRElement.__init__(self, parent, x, check = False)
                return
            raise NotImplementedError, "changing polynomial rings not yet supported"
        raise TypeError, "Cannot form an unramified ring extension element from %s"%(x)

    def _add_(self, right):
        return UnramifiedRingExtensionElement(self.parent(), self._polynomial + right._polynomial, check = False, construct = True)

    def _div_(self, right):
        return self * right.__invert__()

    def __floordiv__(self, right):
        raise NotImplementedError

    def __getitem__(self, n):
        if isinstance(n, slice):
            if n.start == 0:
                raise ValueError, "due to limitations in Python 2.5, you must call the slice() function rather than using the [:] syntax in this case"
            if n.stop == sys.MAXINT:
                return self.slice(n.start, None, n.step)
            return self.slice(n.start, n.stop, n.step)
        if n < 0:
            return self.parent().residue_class_field()(0)
        return self.list()[n]

    def __len__(self):
        return 0 #this is so that __getitem__ with negative inputs works.  It may also have the consequence of having p-adics evaluate to False in a Boolean context.

    def slice(self, i, j, k = 1):
        raise NotImplementedError

    def __invert__(self):
        #The following is a hack to get inverses to work until we fix the inheritance structure
        if self == 0:
            raise ZeroDivisionError
        g = self._polynomial._xgcd(self.parent().defining_polynomial())
        return self.parent()((~g[0].leading_coefficient())*g[1])
        #this is what we should be returning
        return self.parent().fraction_field()(self).__invert__()

    def _latex_(self, mode = None):
        return self._repr_(mode, do_latex = True)

    def __mod__(self, right):
        raise NotImplementedError

    def _mul_(self, right):
        return UnramifiedRingExtensionElement(self.parent(), self._polynomial * right._polynomial, check = False, construct = True)

    def _neg_(self):
        return UnramifiedRingExtensionElement(self.parent(), -self._polynomial, check = False, construct = True)

    def _pari_init_(self):
        raise TypeError, "Pari does not support p-adic extension rings"

    #def __pow__(self, right):  #Using PQRElement's

    def _repr_(self, mode = None, do_latex = False, caprel = False):
        if mode is None:
            mode = self.parent().get_print_mode()
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
            return "not-done"
        elif mode == 'val-unit-p':
            return "not-done"
        elif mode == 'integer':
            return "not-done"
        elif mode == 'integer-p':
            return "not-done"
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
        return UnramifiedRingExtensionElement(self.parent(), self._polynomial - right._polynomial, check = False, construct = True)

    def copy(self):
        return UnramifiedRingExtensionElement(self.parent(), self._polynomial, check = False, construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_square(self):
        if self._polynomial == 0:
            return True
        if self.parent().prime() != 2:
            return self.valuation() % 2 == 0 and self.residue(1).is_square()
        raise NotImplementedError

    def is_zero(self, prec):
        return self.valuation() >= prec

    def is_equal_to(self, right, prec):
        return (self - right).is_zero(prec)

    def list(self):
        K = self.parent().residue_class_field()
        a = K.gen()
        polylist = self._polynomial.list()
        return [sum(polylist[i][j] * (a ** i) for i in range(len(polylist))) for j in range(self.precision_absolute())]

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def minimal_polynomial(self, name):
        return self.minpoly()

    def ordp(self):
        return self.valuation()

    def padded_list(self, n):
        return self.list()[:n] + [0 for w in range(n, self.precision_absolute())]

    def precision_absolute(self):
        if self._polynomial == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity
            else:
                return self.parent().precision_cap()
        return min(c.precision_absolute() for c in self._polynomial.list())

    def precision_relative(self):
        val = self.valuation()
        if val is infinity:
            return Integer(0)
        return self.precision_absolute() - self.valuation()

    def rational_reconstruction(self):
        raise NotImplementedError #should this be a type error instead?

    def residue(self, prec):
        if prec == 1:
            return self.list()[0]
        else:
            raise NotImplementedError

    def square_root(self):
        raise NotImplementedError

    def unit_part(self):
        raise NotImplementedError

    def _unit_part(self):
        raise TypeError, "Extension elements don't implement this function"

    def valuation(self):
        clist = self._polynomial.list()
        if len(clist) == 0:
            if isinstance(self.parent().ground_ring_of_tower(), (pAdicRingCappedRelative, pAdicRingLazy)):
                return infinity
            else:
                return self.parent().precision_cap()
        return min([a.valuation() for a in clist])

    def norm(self, K = None):
        if K is None:
            return PQRElement.norm(self)
        else:
            raise NotImplementedError

    def trace(self, K = None):
        if K is None:
            return PQRElement.trace(self)
        else:
            raise NotImplementedError
