"""
Dense univariate polynomials over Z/nZ, implemented using NTL.

AUTHORS:
    -- Robert Bradshaw: split off from polynomial_element_generic.py (2007-09)

"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

from sage.rings.polynomial.polynomial_element import is_Polynomial, Polynomial_generic_dense

from sage.libs.all import pari, pari_gen

from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZX, zero_ZZX, ZZ_p, ZZ_pX, set_modulus
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
import sage.rings.integer as integer

from sage.rings.fraction_field_element import FractionFieldElement
import sage.rings.polynomial.polynomial_ring

import polynomial_singular_interface
from sage.interfaces.all import singular as singular_default

def make_element(parent, args):
    return parent(*args)


cdef class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite.

    Much of the underlying arithmetic is done using NTL.

    EXAMPLES:

        sage: R.<x> = PolynomialRing(Integers(16))
        sage: f = x^3 - x + 17
        sage: f^2
        x^6 + 14*x^4 + 2*x^3 + x^2 + 14*x + 1

        sage: loads(f.dumps()) == f
        True

        sage: R.<x> = Integers(100)[]
        sage: p = 3*x
        sage: q = 7*x
        sage: p+q
        10*x
        sage: R.<x> = Integers(8)[]
        sage: parent(p)
        Univariate Polynomial Ring in x over Ring of integers modulo 100
        sage: p + q
        10*x

    """
    def __init__(self, parent, x=None, check=True,
                 is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

        if construct:
            if isinstance(x, ZZ_pX):
                self.__poly = x
                return
            self.__poly = ZZ_pX(x, parent.modulus())
            return

        self.__poly = ZZ_pX([], parent.modulus())

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = x.__poly.__copy__()
                return
            else:
                R = parent.base_ring()
                x = [ZZ(R(a)) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            x = self._dict_to_list(x, R(0))


        elif isinstance(x, ZZX):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.Vecrev()]
            check = False

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator() == 1:
                x = x.numerator().__poly
                check = False

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            R = parent.base_ring()
            x = [ZZ(R(a)) for a in x]

        self.__poly = ZZ_pX(x, parent.modulus())

    def _ntl_set_modulus(self):
        self.parent()._ntl_set_modulus()

    def __reduce__(self):
        return make_element, (self.parent(), (self.list(), False, self.is_gen()))

    def int_list(self):
        return eval(str(self.__poly).replace(' ',','))

    def _pari_(self, variable=None):
        """
        EXAMPLES:
            sage: t = PolynomialRing(IntegerModRing(17),"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            Mod(1, 17)*t^3 + Mod(3, 17)*t
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.int_list()).Polrev(variable) * \
               pari(1).Mod(self.parent().base_ring().order())

    def ntl_ZZ_pX(self):
        r"""
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality is available through this
        function.

        WARNING:
        You must call \code{ntl.set_modulus(ntl.ZZ(n))} before doing
        arithmetic with this object!
        """
        return self.__poly

    def __getitem__(self, n):
        return self.parent().base_ring()(self.__poly[n]._sage_())

    def __getslice__(self, i, j):
        R = self.base_ring()
        if i < 0:
            i = 0
        if j > self.__poly.degree()+1:
            j = self.__poly.degree()+1
        v = [R(self.__poly[k]._sage_()) for k in range(i,j)]
        return self.parent()([0]*int(i) + v)

    def _unsafe_mutate(self, n, value):
        n = int(n)
        if n < 0:
            raise IndexError, "n must be >= 0"
        self._ntl_set_modulus()
        self.__poly[n] = int(value)

    def _pow(self, n):
        n = int(n)
        self._ntl_set_modulus()
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        self._ntl_set_modulus()
        return self.parent()(self.__poly + (<Polynomial_dense_mod_n>right).__poly, construct=True)

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(Integers(100), 'x').0
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
        self._ntl_set_modulus()
        return self.parent()(self.__poly * (<Polynomial_dense_mod_n>right).__poly, construct=True)

    def _rmul_(self, c):
        try:
            self._ntl_set_modulus()
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError, msg: # should this realy be a TypeError
            raise TypeError, msg

    def _lmul_(self, c):
        try:
            self._ntl_set_modulus()
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError, msg: # should this realy be a TypeError
            raise TypeError, msg

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_dense_mod_n):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        self._ntl_set_modulus()
        v = self.__poly.quo_rem(right.__poly)
        P = self.parent()
        return (P(v[0], construct=True), P(v[1], construct=True) )

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power $x^n$. If $n$ is negative,
        terms below $x^n$ will be discarded. Does not change this polynomial.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(Integers(12345678901234567890))
            sage: p = x^2 + 2*x + 4
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(-5)
             0
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        AUTHOR:
            -- David Harvey (2006-08-06)
        """
        if n == 0:
            return self
        self._ntl_set_modulus()
        return self.parent()(self.__poly.left_shift(n),
                             construct=True)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        self._ntl_set_modulus()
        return self.parent()(self.__poly - (<Polynomial_dense_mod_n>right).__poly, construct=True)

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and \
                         right[0] in self.parent().base_ring():
            d = right[0]
        elif (right in self.parent().base_ring()):
            d = right
        else:
            return Polynomial.__floordiv__(self, right)
        return self.parent()([c // d for c in self.list()], construct=True)

##     def __copy__(self):
##         self.parent()._ntl_set_modulus()
##         f = self.parent()()
##         f.__poly = self.__poly.copy()
##         return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def is_irreducible(self):
        return bool(self._pari_().polisirreducible())

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES:
            sage: _.<x> = Integers(100)[]
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [83, 3, 0, 1]
        """
        R = self.base_ring()
        return [R(x) for x in self.int_list()]

    def ntl_set_directly(self, v):
        r"""
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers modulo n are stored internally
        using NTL's ZZ_pX class.  Use this function to set the value
        of this polynomial using the NTL constructor, which is
        potentially \emph{very} fast.  The input v is either a vector
        of ints or a string of the form '[ n1 n2 n3 ... ]' where the
        ni are integers and there are no commas between them.  The
        optimal input format is the string format, since that's what
        NTL uses by default.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(Integers(100))
            sage: R([1,-2,3])
            3*x^2 + 98*x + 1
            sage: f = R(0)
            sage: f.ntl_set_directly([1,-2,3])
            sage: f
            3*x^2 + 98*x + 1
            sage: f.ntl_set_directly('[1 -2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 98*x + 1
        """
        if self.is_gen():
            raise TypeError, "Cannot change the value of the generator."
        self._ntl_set_modulus()
        self.__poly = ZZ_pX(v, self.parent().modulus())
        try:
            del self.__list
        except AttributeError:
            pass

    # Polynomial_singular_repr

    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        return polynomial_singular_interface._singular_func(self, singular, have_ring, force)
    def _singular_init_func(self, singular=singular_default, have_ring=False, force=False):
        return polynomial_singular_interface._singular_init_func(self, singular, have_ring, force)
    def lcm(self, singular=singular_default, have_ring=False):
        return polynomial_singular_interface.lcm_func(self, singular, have_ring)
    def diff(self, variable, have_ring=False):
        return self.derivative()
    def resultant(self, other, variable=None):
        return polynomial_singular_interface.resultant_func(self, other, variable)



cdef class Polynomial_dense_mod_p(Polynomial_dense_mod_n):
    """
    A dense polynomial over the integers modulo p, where p is prime.
    """
    def _gcd(self, right):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_dense_mod_p):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.ntl_ZZ_pX().gcd(right.ntl_ZZ_pX())
        return self.parent()(g, construct=True)

    def _xgcd(self, right):
        """
        Return $g, u, v$ such that \code{g = u*self + v*right}.
        """
        r, s, t = self.ntl_ZZ_pX().xgcd(right.ntl_ZZ_pX())
        return self.parent()(r, construct=True), self.parent()(s, construct=True), \
               self.parent()(t, construct=True)

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        EXAMPLES:
            sage: R.<x> = GF(19)['x']
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True
        """
        self.parent()._ntl_set_modulus()
        other = self.parent()._coerce_(other)
        return self.base_ring()(str(self.ntl_ZZ_pX().resultant(other.ntl_ZZ_pX())))

    def discriminant(self):
        """
        EXAMPLES:
            sage: _.<x> = PolynomialRing(GF(19))
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            12
        """
        self.parent()._ntl_set_modulus()
        return self.base_ring()(str(self.ntl_ZZ_pX().discriminant()))

    # PARI is way better than NTL for poly factor for certain degrees, and is called
    # by default in the base class.
    #def factor(self, verbose=False):
    #    M = self.monic()
    #    self.parent()._ntl_set_modulus()
    #    F = [(self.parent()(f, construct=True), n) for f, n in M.ntl_ZZ_pX().factor(verbose)]
    #    return factorization.Factorization(F)