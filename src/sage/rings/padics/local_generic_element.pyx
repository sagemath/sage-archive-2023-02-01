import sage.rings.commutative_ring_element
import sage.rings.integer
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.finite_field_element
from sage.rings.infinity import infinity
import sage.rings.arith
import sage.rings.padics.precision_error
import sage.libs.pari.gen

from sage.structure.element cimport ModuleElement, RingElement, CommutativeRingElement

cdef class LocalGenericElement(CommutativeRingElement):
    #cdef ModuleElement _add_c_impl(self, ModuleElement right):
    #    raise NotImplementedError

    cdef RingElement _div_c_impl(self, RingElement right):
        r"""
        Returns the quotient of self by right.

        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); b = R(5); a / b
            2 + 4*7 + 5*7^2 + 2*7^3 + O(7^4)
        """
        return self * right.__invert__()

    #def __getitem__(self, n):
    #    raise NotImplementedError

    def slice(self, i, j, k = 1):
        if i is None:
            i = self.valuation()
        if j is None:
            j = self.precision_absolute()
        if i is infinity:
            return self
        pk = self.parent().uniformizer_pow(k)
        ppow = self.parent().uniformizer_pow(i)
        ans = self.parent()(0)
        selflist = self.list()
        if self.parent().is_field():
            i -= self.valuation()
            j -= self.valuation()
            i = max(i, 0)
            j = min(j, self.precision_relative())
        else:
            i = max(i, 0)
            j = min(j, self.precision_absolute())
        for c in selflist[i:j:k]:
            ans += ppow * c
            ppow *= pk
        ans += self.parent()(0, absprec = j)
        return ans

    def _latex_(self):
        return self._repr(do_latex = True)

    #def __mod__(self, right):
    #    raise NotImplementedError

    #cdef RingElement _mul_c_impl(self, RingElement right):
    #    raise NotImplementedError

    #cdef _neg_c_impl(self):
    #    raise NotImplementedError

    def _pari_init_(self):
        return self._repr(mode = 'series')

    #def __pow__(self, right):
    #    raise NotImplementedError

    def _repr_(self):
        return self._repr()

    #def _repr(self, mode = None, do_latex = False):
    #    raise NotImplementedError

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        r"""
            Returns the difference between self and right.

            EXAMPLE:
                sage: R = Zp(7,4,'capped-rel','series'); a = R(12); b = R(5); a - b
                7 + O(7^4)
        """
        return self + (-right)

    def add_bigoh(self, prec):
        return self.slice(None, prec)

    def additive_order(self, prec):
        r"""
        Returns the additive order of self, where self is considered to be zero if it is zero modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the additive order of self
        """
        if self.is_zero(prec):
            return Integer(1)
        else:
            return infinity

    #def copy(self):
    #    raise NotImplementedError

    #def exp(self):
    #    raise NotImplementedError

    def is_integral(self):
        """
        Returns whether self is an integral element.

        INPUT:
            self -- a local ring element
        OUTPUT:
            boolean -- whether self is an integral element.

        EXAMPLES:
            sage: R = Qp(3,20)
            sage: a = R(7/3); a.is_integral()
            False
            sage: b = R(7/5); b.is_integral()
            True
        """
        return bool(self.valuation() >= 0)

    #def is_square(self):
    #    raise NotImplementedError

    def is_unit(self):
        """
        Returns whether self is a unit

        INPUT:
            self -- a local ring element
        OUTPUT:
            boolean -- whether self is a unit
        EXAMPLES:
            sage: R = Zp(3,20,'capped-rel'); K = Qp(3,20,'capped-rel')
            sage: R(0).is_unit()
            False
            sage: R(1).is_unit()
            True
            sage: R(2).is_unit()
            True
            sage: R(3).is_unit()
            False
            sage: R(4).is_unit()
            True
            sage: R(6).is_unit()
            False
            sage: R(9).is_unit()
            False
            sage: K(0).is_unit()
            False
            sage: K(1).is_unit()
            True
            sage: K(2).is_unit()
            True
            sage: K(3).is_unit()
            False
            sage: K(4).is_unit()
            True
            sage: K(6).is_unit()
            False
            sage: K(9).is_unit()
            False
            sage: K(1/3).is_unit()
            False
            sage: K(1/9).is_unit()
            False
        """
        return bool(self.valuation() == 0)

    #def is_zero(self, prec):
    #    raise NotImplementedError

    #def is_equal_to(self, right, prec):
    #    raise NotImplementedError

    #def lift(self):
    #    raise NotImplementedError

    #def list(self):
    #    raise NotImplementedError

    #def log(self):
    #    raise NotImplementedError

    #def multiplicative_order(self, prec):
    #    raise NotImplementedError

    #def padded_list(self):
    #    raise NotImplementedError

    #def precision_absolute(self):
    #    raise NotImplementedError

    #def precision_relative(self):
    #    raise NotImplementedError

    #def residue(self, prec):
    #    raise NotImplementedError

    def sqrt(self, extend = True, all = False):
        r"""
        Returns the square root of this local ring element

        INPUT:
            self -- a local ring element
        OUTPUT:
            local ring element -- the square root of self
        """
        return self.square_root(extend, all)

    #def square_root(self):
    #    raise NotImplementedError

    #def unit_part(self):
    #    raise NotImplementedError

    #def valuation(self):
    #    raise NotImplementedError

    def _min_valuation(self):
        return self.valuation()
