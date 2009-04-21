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
    #cpdef ModuleElement _add_(self, ModuleElement right):
    #    raise NotImplementedError

    cpdef RingElement _div_(self, RingElement right):
        r"""
        Returns the quotient of self by right.

        INPUT:
            self -- a p-adic element.
            right -- a p-adic element distinguishable from zero.
                     In a fixed-modulus ring, this element must be a unit.

        EXAMPLES:
            sage: R = Zp(7, 4, 'capped-rel', 'series'); R(3)/R(5)
            2 + 4*7 + 5*7^2 + 2*7^3 + O(7^4)
            sage: R(2/3) / R(1/3)
            2 + O(7^4)
            sage: R(49) / R(7)
            7 + O(7^5)
            sage: R = Zp(7, 4, 'capped-abs', 'series'); 1/R(7)
            7^-1 + O(7^2)
            sage: R = Zp(7, 4, 'fixed-mod'); 1/R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        return self * right.__invert__()

    #def __getitem__(self, n):
    #    raise NotImplementedError

    def slice(self, i, j, k = 1):
        r"""
        Returns the sum of the $p^{i + FOO * k}$ terms of the series expansion
        of self, for $i + FOO*k$ between i and j-1 inclusive, and FOO an arbitrary
        integer. Behaves analogously to the slice function for lists.

        INPUT:
            self -- a p-adic element
            i -- an integer
            j -- an integer
            k -- a positive integer, default value 1

        EXAMPLES:
            sage: R = Zp(5, 6, 'capped-rel')
            sage: a = R(1/2); a
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + O(5^6)
            sage: a.slice(2, 4)
            2*5^2 + 2*5^3 + O(5^4)
            sage: a.slice(1, 6, 2)
            2*5 + 2*5^3 + 2*5^5 + O(5^6)
            sage: a.slice(5, 4)
            O(5^4)
        """
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
        return self._repr_(do_latex = True)

    #def __mod__(self, right):
    #    raise NotImplementedError

    #cpdef RingElement _mul_(self, RingElement right):
    #    raise NotImplementedError

    #cdef _neg_(self):
    #    raise NotImplementedError

    def _pari_init_(self):
        return self._repr_(mode = 'series')

    #def __pow__(self, right):
    #    raise NotImplementedError

    cpdef ModuleElement _sub_(self, ModuleElement right):
        r"""
        Returns the difference between self and right.

        EXAMPLE:
            sage: R = Zp(7, 4, 'capped-rel', 'series'); a = R(12); b = R(5); a - b
            7 + O(7^4)
            sage: R(4/3) - R(1/3)
            1 + O(7^4)
        """
        return self + (-right)

    def add_bigoh(self, prec):
        return self.slice(None, prec)

    def additive_order(self, prec):
        r"""
        Returns the additive order of self, where self is considered to be zero
        if it is zero modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer

        OUTPUT:
            integer -- the additive order of self

        EXAMPLES:
                sage: R = Zp(7, 4, 'capped-rel', 'series'); a = R(7^3); a.additive_order(3)
                1
                sage: a.additive_order(4)
                +Infinity
                sage: R = Zp(7, 4, 'fixed-mod', 'series'); a = R(7^5); a.additive_order(6)
                1
        """
        if self.is_zero(prec):
            return sage.rings.integer.Integer(1)
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
        return self.valuation() >= 0

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
        return self.valuation() == 0

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
        Returns the square root of this local ring element.

        TODO: doc what "extend" and "all" do

        INPUT:
            self -- a local ring element

        OUTPUT:
            local ring element -- the square root of self

        EXAMPLES:
            sage: R = Zp(13, 10, 'capped-rel', 'series')
            sage: a = sqrt(R(-1)); a * a
            12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + 12*13^7 + 12*13^8 + 12*13^9 + O(13^10)
            sage: sqrt(R(4))
            2 + O(13^10)
            sage: sqrt(R(4/9)) * 3
            2 + O(13^10)

        """
        return self.square_root(extend, all)

    #def square_root(self):
    #    raise NotImplementedError

    #def unit_part(self):
    #    raise NotImplementedError

    #def valuation(self):
    #    raise NotImplementedError

    def normalized_valuation(self):
        r"""
        Returns the normalized valuation of this local ring element,
        i.e., the valuation divided by the absolute ramification index.

        INPUT:
            self -- a local ring element.

        OUTPUT:
            rational -- the normalized valuation of self.

        EXAMPLES:
            sage: Q7 = Qp(7)
            sage: R.<x> = Q7[]
            sage: F.<z> = Q7.ext(x^3+7*x+7)
            sage: z.normalized_valuation()
            1/3
        """
        F = self.parent()
        return self.valuation()/F.ramification_index()

    def _min_valuation(self):
        r"""
        Returns the valuation of this local ring element.

        TODO: this docstring can't be right. Surely this is different from "valuation".

        INPUT:
            self -- a local ring element.

        OUTPUT:
            integer -- the valuation of self.

        EXAMPLES:
            sage: R = Qp(7, 4, 'capped-rel', 'series')
            sage: R(7)._min_valuation()
            1
            sage: R(1/7)._min_valuation()
            -1
        """
        return self.valuation()
