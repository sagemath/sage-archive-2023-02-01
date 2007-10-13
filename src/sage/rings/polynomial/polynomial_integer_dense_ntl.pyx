"""
Dense univariate polynomials over Z, implemented using NTL.

AUTHORS:
    -- David Harvey: split off from polynomial_element_generic.py (2007-09)

"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################


# todo: should probably be cimporting some of these:

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport ModuleElement, RingElement


from sage.rings.polynomial.polynomial_element import is_Polynomial

from sage.libs.ntl.all import ZZX, zero_ZZX

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer

from sage.libs.all import pari, pari_gen
from sage.structure.factorization import Factorization

from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.arith import lcm
import sage.rings.polynomial.polynomial_ring


cdef class Polynomial_integer_dense_ntl(Polynomial):
    """
    A dense polynomial over the integers.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if construct:
            if isinstance(x, ZZX):
                self.__poly = x
                return
            self.__poly = ZZX(x)
            return

        self.__poly = zero_ZZX.copy()

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = (<Polynomial_integer_dense_ntl>x).__poly.copy()
                return
            else:
                x = [ZZ(a) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            x = self._dict_to_list(x, ZZ(0))


        elif isinstance(x, ZZX):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.Vecrev()]
            check = False

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_integer_dense_ntl):
            if x.denominator() == 1:
                x = (<Polynomial_integer_dense_ntl>x.numerator()).__poly
                check = False

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            x = [ZZ(z) for z in x]

        self.__poly = ZZX(x)

    def content(self):
        """
        Return the greatest common divisor of the coefficients of this
        polynomial.
        """
        return ZZ(self.__poly.content())

    def ntl_ZZX(self):
        """
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality may be available through
        this function.
        """
        return self.__poly

    def __reduce__(self):
        return Polynomial_integer_dense_ntl, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        return self.__poly[n].get_as_sage_int()

    def __getslice__(self, i, j):
        i = max(0,i)
        j = min(j, self.__poly.degree()+1)
        v = [self.__poly[k].get_as_sage_int() for k in range(i,j)]
        P = self.parent()
        return P([0] * int(i) + v)

    def _pow(self, n):
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        return self.parent()(self.__poly + (<Polynomial_integer_dense_ntl>right).__poly, construct=True)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        return self.parent()(self.__poly - (<Polynomial_integer_dense_ntl>right).__poly, construct=True)

    # todo: write _neg_c_impl

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_integer_dense_ntl):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        v = self.__poly.quo_rem((<Polynomial_integer_dense_ntl>right).__poly)
        return self.parent()(v[0], construct=True), \
               self.parent()(v[1], construct=True)

    def gcd(self, right):
        """
        Return the GCD of self and other.  The leading
        coefficient need not be 1.
        """
        if not isinstance(right, Polynomial_integer_dense_ntl):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.gcd((<Polynomial_integer_dense_ntl>right).__poly)
        return self.parent()(g, construct=True)

    def lcm(self, right):
        """
        Return the LCM of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_integer_dense_ntl):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.lcm((<Polynomial_integer_dense_ntl>right).__poly)
        return self.parent()(g, construct=True)

    def xgcd(self, right):
        """
        Return $g, u, v$ such that \code{g = u*self + v*right}.

        If self and right are coprime as polynomials over the
        rationals, then $g$ is guaranteed to be the resultant of self
        and right, as a constant polynomial.

        EXAMPLES:
            sage: P.<x> = PolynomialRing(ZZ)
            sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (27*x^2 + 54, 1, -x^2 - 3*x - 9)
            sage: u*F + v*G
            27*x^2 + 54
            sage: x.xgcd(P(0))
            (1, 0, x)
            sage: f = P(0)
            sage: f.xgcd(x)
            (x, 0, 1)
            sage: F = (x-3)^3; G = (x-15)^2
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (2985984, -432*x + 8208, 432*x^2 + 864*x + 14256)
            sage: u*F + v*G
            2985984
        """
        r, s, t = self.ntl_ZZX().xgcd(right.ntl_ZZX())
        K = self.base_ring()
        rr = K(str(r))   # optimize in future
        if rr == 0:
            f = self.base_extend(QQ)
            g, u, v = f.xgcd(right.base_extend(QQ))
            d = lcm([g.denominator(), u.denominator(), v.denominator()])
            R = self.parent()
            return R(d*g), R(d*u), R(d*v)
        else:
            S = self.parent()
            return S(rr), S(s, construct=True), \
                   S(t, construct=True)


    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: _.<x> = PolynomialRing(ZZ)
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 - 10*x^2 + 32*x - 32
        """
        return self.parent()(self.__poly * (<Polynomial_integer_dense_ntl>right).__poly, construct=True)

    def _lmul_(self, right):
        # todo: does this get called in the cython version?
        try:
            return self.parent()(ZZX([right]) * self.__poly, construct=True)
        except RuntimeError, msg:
            raise TypeError, msg

    def _rmul_(self, left):
        # todo: does this get called in the cython version?
        try:
            return self.parent()(ZZX([left]) * self.__poly, construct=True)
        except RuntimeError, msg:
            raise TypeError, msg

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and right[0] in ZZ:
            d = ZZ(right[0])
        elif (right in self.parent().base_ring()):
            d = ZZ(right)
        else:
            return Polynomial.__floordiv__(self, right)
        return self.parent()([c // d for c in self.list()], construct=True)

    def _unsafe_mutate(self, n, value):
        n = int(n)
        if n < 0:
            raise IndexError, "n must be >= 0"
        self.__poly[n] = int(value)

    def complex_roots(self, flag=0):
        """
        Returns the complex roots of this polynomial.
        INPUT:
            flag -- optional, and can be
                    0: (default), uses Schonhage's method modified by Gourdon,
                    1: uses a modified Newton method.
        OUTPUT:
            list of complex roots of this polynomial, counted with multiplicities.

        NOTE: Calls the pari function polroots.

        EXAMPLE:
        We compute the roots of the characteristic polynomial of some Salem numbers:
            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: alpha = f.complex_roots()[0]; alpha
            0.713639173536901
            sage: f(alpha)
            0
        """
        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(QQ, 'x')
        return R(self.list()).complex_roots()

    def real_root_intervals(self):
        """
        Returns isolating intervals for the real roots of this polynomial.

        EXAMPLE:
        We compute the roots of the characteristic polynomial of some Salem numbers:
            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.real_root_intervals()
            [((1/2, 3/4), 1), ((1, 3/2), 1)]
        """

        from sage.rings.polynomial.real_roots import real_roots

        return real_roots(self)

##     def __copy__(self):
##         f = Polynomial_integer_dense(self.parent())
##         f.__poly = self.__poly.copy()
##         return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def discriminant(self):
        """
        EXAMPLES:
            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            -7911
        """
        return ZZ(str(self.__poly.discriminant()))

    def _pari_(self, variable=None):
        """
        EXAMPLES:
            sage: t = PolynomialRing(ZZ,"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            t^3 + 3*t - 17
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.list()).Polrev(variable)

    def squarefree_decomposition(self):
        """
        Return the square-free decomposition of self.  This is
        a partial factorization of self into square-free, relatively
        prime polynomials.

        This is a wrapper for the NTL function SquareFreeDecomp.

        EXAMPLES:
            sage: x = polygen(ZZ)
            sage: p = 37 * (x-1)^2 * (x-2)^2 * (x-3)^3 * (x-4)
            sage: p.squarefree_decomposition()
            (37) * (x - 4) * (x^2 - 3*x + 2)^2 * (x - 3)^3
        """
        p = self.__poly
        c = p.content()
        if c != 1:
            p = p / ZZX([c])
        decomp = p.square_free_decomposition()
        pdecomp = [(self.parent()(f, construct=True), exp) for (f, exp) in decomp]
        return Factorization(pdecomp, unit=c, sort=False)

    def factor_mod(self, p):
        """
        Return the factorization of self modulo the prime p.

        INPUT:
            p -- prime

        OUTPUT:
            factorization of self reduced modulo p.

        EXAMPLES:
            sage: R.<x> = ZZ['x']
            sage: f = -3*x*(x-2)*(x-9) + x
            sage: f.factor_mod(3)
            x
            sage: f = -3*x*(x-2)*(x-9)
            sage: f.factor_mod(3)
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

            sage: f = 2*x*(x-2)*(x-9)
            sage: f.factor_mod(7)
            (2) * x * (x + 5)^2
        """
        from sage.rings.finite_field import FiniteField
        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        f = self._pari_()
        if f * pari('Mod(1,%s)'%p) == pari(0):
            raise ValueError, "factorization of 0 not defined"
        G = f.factormod(p)
        k = FiniteField(p)
        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(k, names=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, unit=R(self).leading_coefficient())


    def factor_padic(self, p, prec=10):
        """
        Return p-adic factorization of self to given precision.

        INPUT:
            p -- prime
            prec -- integer; the precision

        OUTPUT:
            factorization of self reduced modulo p.
        """
        from sage.rings.padics.factory import Zp
        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        prec = Integer(prec)
        if prec <= 0:
            raise ValueError, "prec must be positive"
        K = Zp(p, prec, type='capped-abs')
        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(K, names=self.parent().variable_name())
        return R(self).factor()

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES:
            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [-17, 3, 0, 1]
            sage: f = PolynomialRing(ZZ,'x')(0)
            sage: f.list()
            []
        """
        return [self[i] for i in range(self.degree()+1)]

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        NOTES:
            Implemented using NTL's polresultant function.

        EXAMPLES:
            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            -8
            sage: r.parent() is ZZ
            True
        """
        other = self.parent()._coerce_(other)
        return ZZ(str(self.__poly.resultant((<Polynomial_integer_dense_ntl>other).__poly, 0)))

    def ntl_set_directly(self, v):
        """
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers are stored internally using NTL's ZZX
        class.  Use this function to set the value of this polynomial using
        the NTL constructor, which is potentially quicker.   The input v
        is either a vector of ints or a string of the form '[ n1 n2 n3 ... ]'
        where the \code{n_i} are integers and there are no commas between them.
        The optimal input format is the string format, since that's what NTL uses.

        EXAMPLES:
            sage: R.<w> = PolynomialRing(ZZ)
            sage: R([1,2,3])
            3*w^2 + 2*w + 1
            sage: f = R(0)
            sage: f.ntl_set_directly([1,2,3])
            sage: f
            3*w^2 + 2*w + 1
            sage: f.ntl_set_directly('[1 2 3 4]')
            sage: f
            4*w^3 + 3*w^2 + 2*w + 1
        """
        if self.is_gen():
            raise TypeError, "Cannot change the value of the generator."
        self.__poly = ZZX(v)
        try:
            del self.__list
        except AttributeError:
            pass

