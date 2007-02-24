import sage.rings.commutative_ring_element
import sage.rings.integer
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.finite_field_element
import sage.rings.infinity
import sage.rings.arith
import sage.rings.padics.precision_error
import sage.libs.pari.gen

infinity = sage.rings.infinity.infinity

class LocalGenericElement(sage.rings.commutative_ring_element.CommutativeRingElement):
    def _add_(self, right):
        raise NotImplementedError

    def _div_(self, right):
        r"""
        Returns the quotient of self by right.

        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); b = R(5); a / b
            2 + 4*7 + 5*7^2 + 2*7^3 + O(7^4)
        """
        return self * right.__invert__()

    def __getitem__(self, n):
        raise NotImplementedError

    def __getslice__(self, i, j):
        raise NotImplementedError

    def _latex_(self):
        return self._repr_(do_latex = True)

    def __mod__(self, right):
        raise NotImplementedError

    def _mul_(self, right):
        raise NotImplementedError

    def _neg_(self):
        raise NotImplementedError

    def _pari_init_(self):
        return self._repr_(mode = 'series')

    def __pow__(self, right):
        raise NotImplementedError

    def _repr_(self, mode = None, do_latex = False):
        raise NotImplementedError

    def _sub_(self, right):
        r"""
            Returns the difference between self and right.

            EXAMPLE:
                sage: R = Zp(7,4,'capped-rel','series'); a = R(12); b = R(5); a - b
                7 + O(7^4)
        """
        return self + right._neg_()

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

    def copy(self):
        raise NotImplementedError

    def exp(self):
        raise NotImplementedError

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

    def is_square(self):
        raise NotImplementedError

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

    def is_zero(self, prec):
        raise NotImplementedError

    def is_equal_to(self, right, prec):
        raise NotImplementedError

    def lift(self):
        raise NotImplementedError

    def list(self):
        raise NotImplementedError

    def log(self):
        raise NotImplementedError

    def multiplicative_order(self, prec):
        raise NotImplementedError

    def padded_list(self):
        raise NotImplementedError

    def precision_absolute(self):
        raise NotImplementedError

    def precision_relative(self):
        raise NotImplementedError

    def residue(self, prec):
        raise NotImplementedError

    def set_precision_absolute(self, prec):
        raise NotImplementedError

    def set_precision_relative(self, prec):
        raise NotImplementedError

    def sqrt(self):
        r"""
        Returns the square root of this local ring element

        INPUT:
            self -- a local ring element
        OUTPUT:
            local ring element -- the square root of self
        """
        return self.square_root()

    def square_root(self):
        raise NotImplementedError

    def unit_part(self):
        raise NotImplementedError

    def valuation(self):
        raise NotImplementedError

    def _min_valuation(self):
        return self.valuation()
