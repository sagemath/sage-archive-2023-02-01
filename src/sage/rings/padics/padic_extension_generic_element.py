from __future__ import with_statement
from sage.rings.finite_field import GF
from sage.rings.padics.misc import min
import sage.rings.polynomial_quotient_ring_element
import sage.rings.integer
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.finite_field_element
import sage.rings.polynomial_ring_constructor
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_ring_lazy
import sage.rings.padics.padic_ring_capped_relative
import sage.rings.polynomial_element
import sage.rings.infinity
import sys
import padic_generic

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
local_print_mode = sage.rings.padics.padic_generic.local_print_mode

class pAdicExtensionGenericElement(PQRElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        if construct:
            PQRElement.__init__(self, parent, x, check=False)
            return
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if isinstance(x, PQRElement):
            if x.parent() is parent:
                if absprec is infinity and relprec is infinity:
                    PQRElement.__init__(self, parent, x._polynomial, check = False)
                elif absprec is infinity:
                    PQRElement.__init__(self, parent, x.add_bigoh(relprec + x.valuation())._polynomial, check = False) #not the fastest method
                else:
                    PQRElement.__init__(self, parent, x.add_bigoh(absprec)._polynomial, check = False) #not the fastest method
                return
            elif parent.base_ring().has_coerce_map_from(x.parent().base_ring()) and parent._coerce_(x.parent().defining_polynomial()) == parent.defining_polynomial(): #has_coerce_map should maybe be is_isomorphic? or fraction_field?
                if absprec is infinity and relprec is infinity:
                    PQRElement.__init__(self, parent, x.parent().polynomial_ring().base_extend(parent.base_ring())(x._polynomial), check = False)
                elif absprec is infinity:
                    PQRElement.__init__(self, parent, x.parent().polynomial_ring().base_extend(parent.base_ring())(x.add_bigoh(relprec + x.valuation())._polynomial), check = False)
                else:
                    PQRElement.__init__(self, parent, x.parent().polynomial_ring().base_extend(parent.base_ring())(x.add_bigoh(absprec)._polynomial), check = False)
                return
        elif isinstance(x, Polynomial):
            if x.parent() is parent.polynomial_ring():
                PQRElement.__init__(self, parent, x, check = False)
                return
        if isinstance(x, list):
            if len(x) == parent.degree(): #don't do this with length checking.  I don't remember what William and I decided
                #should probably be wrapped in a try block, but I'm not sure what types of exceptions I'm worried about.  TypeError?
                absprec = min([v.precision_absolute() for v in x] + [absprec])
                val = min([v.valuation() for v in x])
                absprec = min(absprec, relprec + val)
                x = [parent.base_ring()(v, absprec = absprec) for v in x]
                PQRElement.__init__(self, parent, parent.polynomial_ring()(x), check=False)
                return
        #We now try casting into the base ring
        try:
            x = parent.base_ring()(x, absprec = absprec, relprec = relprec)
            PQRElement.__init__(self, parent, parent.polynomial_ring()(x), check=False)
        except TypeError:
            raise TypeError, "cannot create element of %s out of %s of type %s"%(parent, x, type(x))

    def _repr_(self, mode = None, do_latex = False):
        if self._is_exact_zero():
            return "0"
        if mode is None:
            mode = self.parent().print_mode()
        elif not ((mode == 'val-unit') or (mode == 'series') or (mode == 'terse')):
            raise TypeError, "printing-mode must be one of 'val-unit', 'series', or 'terse'"
        pprint = self.parent().variable_name()
        if self._polynomial == 0:
            if mode == 'val-unit' or mode == 'series':
                if do_latex:
                    return "O(%s^{%s})"%(pprint, self.precision_absolute())
                else:
                    return "O(%s^%s)"%(pprint, self.precision_absolute())
            elif mode == 'terse':
                if do_latex:
                    return "0 + O(%s^{%s})"%(pprint, self.precision_absolute())
                else:
                    return "0 + O(%s^%s)"%(pprint, self.precision_absolute())
        if mode == 'val-unit':
            if self.valuation() == 0:
                return self._repr_(mode = 'terse', do_latex = do_latex)
            elif self.valution() == 1:
                return "%s + %s"%(self._uniformizer_print(), self._repr_(mode = 'terse', do_latex = do_latex))
            else:
                if do_latex:
                    return "%s^{%s} + %s"%(self._uniformizer_print(), self.valuation(), self._repr_(mode = 'terse', do_latex = do_latex))
                else:
                    return "%s^%s + %s"%(self._uniformizer_print(), self.valuation(), self._repr_(mode = 'terse', do_latex = do_latex))
        elif mode == 'terse':
            with local_print_mode(self.base_ring(), 'terse'):
                if do_latex:
                    return self._polynomial._latex_(name = self.parent().variable_name())
                else:
                    return self._polynomial._repr(name = self.parent().variable_name())
        else:
            pprint = self.parent()._uniformizer_print()
            selflist = [a.residue(1) for a in self.padded_list(self.precision_absolute())]
            triples = [(selflist[i], pprint, i) for i in range(2, self.precision_absolute())]
            #print triples
            triples = [a for a in triples if a[0] != 0] #need to change this later to account for __cmp__ throwing error on lazies
            def addparen(c):
                c = "%s"%c
                if c.find("+") != -1:
                    c = "(%s)"%c
                return c
            triples = [(addparen(a[0]), a[1], a[2]) for a in triples]
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

    def _div_(self, right):
        return self * right.__invert__()

    def __floordiv__(self, right):
        if self.parent().is_field():
            return self._div_(right)
        right = self.parent()(right)
        return self.parent()(self / right.unit_part()).__rshift__(right.valuation())

    def __mod__(self, right):
        if self.parent().is_field():
            return self.parent()(0)
        return self - right * (self // right)

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

    def __invert__(self):
        if self == 0:
            raise ZeroDivisionError
        if self.parent().is_field():
            g = self._polynomial._xgcd(self.parent().defining_polynomial())
            return self.parent()((~g[0].leading_coefficient())*g[1])
        else:
            return self.parent().fraction_field()(self).__invert__()

    def _latex_(self, mode = None):
        return self._repr_(mode, do_latex = True)

    def _pari_init_(self):
        raise TypeError, "Pari does not support p-adic extension rings"

    #def __pow__(self, right):  #Using PQRElement's

    def is_square(self): #is this implementation still correct in the eisenstein case?
        if self._polynomial == 0:
            return True
        if self.parent().prime() != 2:
            return self.valuation() % 2 == 0 and self.residue(1).is_square()
        raise NotImplementedError

    def _is_exact_zero(self):
        return False

    def is_zero(self, prec):
        return self.valuation() >= prec

    def is_equal_to(self, right, prec):
        return (self - right).is_zero(prec)

    def minimal_polynomial(self, name):
        return self.minpoly()

    def ordp(self):
        return self.valuation()

    def padded_list(self, absprec = infinity, relprec = infinity):
        if absprec is infinity and relprec is infinity:
            raise ValueError, "must specify at least one of absprec and relprec"
        if self.parent().is_field():
            if self.valuation() is infinity:
                if relprec < infinity:
                    return [self.parent().residue_class_field()(0)]*relprec
                else:
                    return []
            relprec = min(relprec, absprec - self.valuation())
            retlist = self.list()
            return retlist[:relprec] + [self.parent()(0)]*(relprec - len(retlist))
        if self.valuation() is infinity:
            if absprec < infinity:
                return [self.parent()(0)]*absprec
            else:
                raise ValueError, "would return infinite list"
        absprec = min(absprec, relprec + self.valuation())
        retlist = self.list()
        return retlist[:absprec] + [self.parent()(0)]*(absprec - len(retlist))

    def precision_relative(self):
        val = self.valuation()
        if val is infinity:
            return Integer(0)
        return self.precision_absolute() - self.valuation()

    def rational_reconstruction(self):
        raise NotImplementedError #should this be a type error instead?

    def unit_part(self):
        return self.__rshift__(self.valuation())

    def _unit_part(self):
        raise TypeError, "Extension elements don't implement this function"

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
