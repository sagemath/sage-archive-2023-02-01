"""
Univariate Polynomial Base Class

AUTHORS:
    -- William Stein: first version
    -- Martin Albrecht: Added singular coercion.
"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

import operator

import copy

import sage.rings.rational
import integer
import sage.rings.polynomial_ring
import arith
import sage.rings.ring_element as ring_element
import integer_ring
import rational_field
import integer_mod_ring
import polynomial_pyx
import rational_field
import complex_field
#import padic_field
from infinity import infinity
import sage.misc.misc as misc
from sage.misc.sage_eval import sage_eval
from sage.misc.latex import latex
from sage.structure.factorization import Factorization

from sage.interfaces.all import singular as singular_default, is_SingularElement
from sage.libs.all import pari, pari_gen

from real_mpfr import RealField, is_RealNumber, is_RealField
RR = RealField()

from sage.structure.element import RingElement
from sage.structure.element cimport Element

from rational_field import QQ
from integer_ring import ZZ

from integral_domain import is_IntegralDomain

import polynomial_fateman

def is_Polynomial(f):
    return PY_TYPE_CHECK(f, Polynomial)

cdef class Polynomial(CommutativeAlgebraElement):
    """
    A polynomial.

    EXAMPLE:
        sage: R.<y> = QQ['y']
        sage: S.<x> = R['x']
        sage: f = x*y; f
        y*x
        sage: type(f)
        <class 'sage.rings.polynomial_element_generic.Polynomial_generic_dense'>
    """
    def __init__(self, parent, is_gen = False, construct=False):
        """
        The following examples illustrate creation of elements of
        polynomial rings, and some basic arithmetic.

        First we make a polynomial over the integers and do some arithmetic:
            sage: R.<x> = ZZ[]
            sage: f = x^5 + 2*x^2 + (-1); f
            x^5 + 2*x^2 - 1
            sage: f^2
            x^10 + 4*x^7 - 2*x^5 + 4*x^4 - 4*x^2 + 1

        Next we do arithmetic in a sparse polynomial ring over the integers:
            sage: R.<x> = ZZ[ ]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: S.<Z> = R[ ]; S
            Univariate Polynomial Ring in Z over Univariate Polynomial Ring in x over Integer Ring
            sage: f = Z^3 + (x^2-2*x+1)*Z - 3; f
            Z^3 + (x^2 - 2*x + 1)*Z + -3
            sage: f*f
            Z^6 + (2*x^2 - 4*x + 2)*Z^4 + (-6)*Z^3 + (x^4 - 4*x^3 + 6*x^2 - 4*x + 1)*Z^2 + (-6*x^2 + 12*x - 6)*Z + 9
            sage: f^3 == f*f*f
            True
        """
        CommutativeAlgebraElement.__init__(self, parent)
        self._is_gen = is_gen

    def _add_(self, right):
        if self.degree() >= right.degree():
            x = list(self.list())
            y = right.list()
        else:
            x = list(right.list())
            y = self.list()

        for i in xrange(len(y)):
            x[i] += y[i]

        return self.polynomial(x)

    def _lmul_(self, left):
        """
        Multiply self on the left by a scalar.

        EXAMPLE:
            sage: R.<x> = ZZ[]
            sage: f = (x^3 + x + 5)
            sage: f._lmul_(7)
            7*x^3 + 7*x + 35
            sage: 7*f
            7*x^3 + 7*x + 35
        """
        # todo -- should multiply individual coefficients??
        #         that could be in derived class.
        #         Note that we are guaranteed that right is in the base ring, so this could be fast.
        return self.parent()(left) * self

    def _rmul_(self, right):
        """
        Multiply self on the right by a scalar.

        EXAMPLE:
            sage: R.<x> = ZZ[]
            sage: f = (x^3 + x + 5)
            sage: f._rmul_(7)
            7*x^3 + 7*x + 35
            sage: f*7
            7*x^3 + 7*x + 35
        """
        # todo -- Should multiply individual coefficients??
        #         that could be in derived class.
        #         Note that we are guaranteed that right is in the base ring, so this could be fast.
        return self * self.parent()(right)

    def __call__(self, *a):
        """
        Evaluate polynomial at x=a using Horner's rule

        INPUT:
            a -- ring element a; need not be in the coefficient
                 ring of the polynomial.

        OUTPUT:
            the value of f at a.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = x/2 - 5
            sage: f(3)
            -7/2
            sage: R.<x> = ZZ[]
            sage: f = (x-1)^5
            sage: f(2/3)
            -1/243

        We evaluate a polynomial over a quaternion algebra:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: R.<w> = PolynomialRing(A,sparse=True)
            sage: f = i*j*w^5 - 13*i*w^2 + (i+j)*w + i
            sage: f(i+j+1)
            24 + 26*i - 10*j - 25*k
            sage: w = i+j+1; i*j*w^5 - 13*i*w^2 + (i+j)*w + i
            24 + 26*i - 10*j - 25*k

        AUTHORS:
            -- David Joyner, 2005-04-10
            -- William Stein, 2006-01-22; change so parent
               is determined by the arithmetic
        """
        a = a[0]
        if isinstance(a, tuple):
            a = a[0]
        d = self.degree()
        result = self[d]
        i = d - 1
        while i >= 0:
            result = result * a + self[i]
            i -= 1
        return result

    cdef int _cmp_c_impl(self, Element other) except -2:
        """
        Compare the two polynomials self and other.

        We order polynomials first by degree, then in dictionary order
        starting with the coefficient of largest degree.

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: 3*x^3  + 5 > 10*x^2 + 19
            True
            sage: x^2 - 2*x - 1 < x^2 - 1
            True
            sage: x^2 - 2*x - 1 > x^2 - 1
            False
            sage: R(-1) < 0
            False
            sage: x^3 - 3 > 393939393
            True
        """
        d1 = self.degree(); d2 = other.degree()
        c = cmp(d1, d2)
        if c: return c
        for i in reversed(xrange(d1+1)):
            c = cmp(self[i], other[i])
            if c: return c
        return 0

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: P = PolynomialRing(ZZ,'x')(0)
            sage: bool(P)
            False
            sage: P = PolynomialRing(ZZ, 'x')([1,2,3])
            sage: bool(P)
            True
        """
        return self.degree() >= 0

    def __getitem__(self, n):
        raise NotImplementedError

    def __iter__(self):
        return iter(self.list())

    def __hash__(self):
        if self.degree() >= 1:
            return hash(tuple(self.list()))
        else:
            return hash(self[0])

    def __float__(self):
         if self.degree() > 0:
             raise TypeError, "cannot coerce nonconstant polynomial to float"
         return float(self[0])

    def __int__(self):
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial to int"
        return int(self[0])

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: H = Hom(R, QQ); H
            Set of Homomorphisms from Univariate Polynomial Ring in x over Integer Ring to Rational Field
            sage: f = H([5]); f
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Rational Field
              Defn: x |--> 5
            sage: f(x)
            5
            sage: f(x^2 + 3)
            28
        """
        a = im_gens[0]
        P = a.parent()
        d = self.degree()
        result = P._coerce_(self[d])
        i = d - 1
        while i >= 0:
            result = result * a + P._coerce_(self[i])
            i -= 1
        return result

    def _integer_(self):
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial"
        return integer.Integer(self[0])

    def __invert__(self):
        """
        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = x - 90283
            sage: f.__invert__()
            1/(x - 90283)
            sage: ~f
            1/(x - 90283)
        """
        return self.parent()(1)/self

    def inverse_of_unit(self):
        """
        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = x - 90283
            sage: f.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ValueError: self is not a unit.
            sage: f = R(-90283); g = f.inverse_of_unit(); g
            -1/90283
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field
        """
        if self.degree() > 0:
            raise ValueError, "self is not a unit."
        return self.parent()(~(self[0]))

    def __long__(self):
        """
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: f = x - 902384
            sage: long(f)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
            sage: long(R(939392920202))
            939392920202L
        """
        if self.degree() > 0:
            raise TypeError, "cannot coerce nonconstant polynomial to long"
        return long(self[0])

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: (x - 4)*(x^2 - 8*x + 16)
            x^3 - 12*x^2 + 48*x - 64
        """
        if right == 0 or self == 0:
            return self.polynomial(0)
        return self._mul_karatsuba(right)

    def square(self):
        """
        Returns the square of this polynomial.

        TODO:
          -- This is just a placeholder; for now it just uses ordinary
          multiplication. But generally speaking, squaring is faster than
          ordinary multiplication, and it's frequently used, so subclasses
          may choose to provide a specialised squaring routine.

          -- Perhaps this even belongs at a lower level? ring_element
          or something?

        AUTHOR:
          -- David Harvey (2006-09-09)

        """
        return self * self

    def __div__(self, right):
        """
        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = (x^3 + 5)/3; f
            1/3*x^3 + 5/3
            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field

        If we do the same over $\ZZ$ the result is in the polynomial
        ring over $\QQ$.

            sage: x  = ZZ['x'].0
            sage: f = (x^3 + 5)/3; f
            1/3*x^3 + 5/3
            sage: f.parent()
            Univariate Polynomial Ring in x over Rational Field

        Divides can make elements of the fraction field:

            sage: R.<x> = QQ['x']
            sage: f = x^3 + 5
            sage: g = R(3)
            sage: h = f/g; h
            1/3*x^3 + 5/3
            sage: h.parent()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field

        This is another example over a non-prime finite field
        (submited by a student of Jon Hanke).  It illustrates
        cancellation between the numerator and denominator
        over a non-prime finite field.
            sage: R.<x> = PolynomialRing(GF(5^2, 'a'), 'x')
            sage: f = x^3 + 4*x
            sage: f / (x - 1)
            x^2 + x
        """
        try:
            if not isinstance(right, Element) or right.parent() != self.parent():
                R = self.parent().base_ring()
                x = R(right)
                return ~x * self
        except (TypeError, ValueError, ZeroDivisionError):
            pass
        return RingElement.__div__(self, right)


    def _pow(self, right):
        """
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: f = x - 1
            sage: f._pow(3)
            x^3 - 3*x^2 + 3*x - 1
            sage: f^3
            x^3 - 3*x^2 + 3*x - 1
        """
        if self.degree() <= 0:
            return self.parent()(self[0]**right)
        if right < 0:
            return (~self)**(-right)
        if self._is_gen:   # special case x**n should be faster!
            v = [0]*right + [1]
            return self.parent()(v, check=True)
        return arith.generic_power(self, right, self.parent()(1))

    def _repr(self, name=None):
        s = " "
        m = self.degree() + 1
        r = reversed(xrange(m))
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x != 0:
                if n != m-1:
                    s += " + "
                x = str(x)
                if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if s==" ":
            return "0"
        return s[1:]

    def _repr_(self):
        r"""
        EXAMPLES:
            sage: f = x^3+2/3*x^2 - 5/3
            sage: f._repr_()
            'x^3 + 2/3*x^2 - 5/3'
            sage: f.rename('vaughn')
            sage: f
            vaughn
        """
        return self._repr()

    def _latex_(self, name=None):
        r"""
        EXAMPLES:
            sage: latex(x^3+2/3*x^2 - 5/3)
             x^{3} + \frac{2}{3}x^{2} - \frac{5}{3}
        """
        s = " "
        m = self.degree() + 1
        r = reversed(xrange(m))
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = self.list()
        for n in reversed(xrange(m)):
            x = coeffs[n]
            if x != 0:
                if n != m-1:
                    s += " + "
                x = latex(x)
                if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "|%s^{%s}"%(name,n)
                elif n==1:
                    var = "|%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")
        if s==" ":
            return "0"
        return s[1:]


    def __setitem__(self, n, value):
        raise IndexError, "polynomials are immutable"


    def __floordiv__(self,right):
        """
        Quotient of division of self by other.  This is denoted //.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def div(self,right):
        """
        Quotient of division of self by other.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def __mod__(self, other):
        """
        Remainder of division of self by other.
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: x % (x+1)
            -1
            sage: (x^3 + x - 1) % (x^2 - 1)
            2*x - 1
        """
        _, R = self.quo_rem(other)
        return R

    def _is_atomic(self):
        return PyBool_FromLong(self.degree() == self.valuation())

    def _mul_generic(self, right):
        d1 = self.degree()
        d2 = right.degree()
        d = d1 + d2
        w = [sum([self[i]*right[k-i] for i in range(0,min(d1,k)+1) if \
                  i <= d1 and k-i <= d2 and self[i]!=0 and right[k-i]!=0]) \
                for k in range(d+1)]
        return self.parent()(w)

    def _mul_fateman(self, right):
        r"""
        Returns the product of two polynomials using Kronecker's trick
        to do the multiplication.  This could be used used over a
        generic base ring.

        NOTES:
        \begin{itemize}
          \item Since this is implemented in interpreted Python, it
                could be hugely sped up by reimplementing it in Pyrex.
          \item Over the reals there is precision loss, at least in
                the current implementation.
        \end{itemize}

        INPUT:
           self -- Polynomial
           right -- Polynomial (over same base ring as self)

        OUTPUT: Polynomial
           The product self*right.

        ALGORITHM:
        Based on a paper by R. Fateman

          {\tt http://www.cs.berkeley.edu/~fateman/papers/polysbyGMP.pdf}

        The idea is to encode dense univariate polynomials as big
        integers, instead of sequences of coefficients. The paper
        argues that because integer multiplication is so cheap, that
        encoding 2 polynomials to big numbers and then decoding the
        result might be faster than popular multiplication algorithms.
        This seems true when the degree is larger than 200.

        EXAMPLES:
            sage: S.<y> = PolynomialRing(RR)
            sage: f = y^10 - 1.393493*y + 0.3
            sage: f._mul_karatsuba(f)
            1.00000000000000*y^20 - 2.78698600000000*y^11 + 0.600000000000000*y^10 + 0.000000000000000111022302462516*y^8 - 0.000000000000000111022302462516*y^6 - 0.000000000000000111022302462516*y^3 + 1.94182274104900*y^2 - 0.836095800000000*y + 0.0900000000000000
            sage: f._mul_fateman(f)
            1.00000000000000*y^20 - 2.78698600000000*y^11 + 0.600000000000000*y^10 + 1.94182274104900*y^2 - 0.836095800000000*y + 0.0900000000000000

        Advantages:

        \begin{itemize}

        \item Faster than Karatsuba over $\Q$ and $\Z$
             (but much slower still than calling NTL's
             optimized C++ implementation, which is the
             default over $\Z$)

        \item Potentially less complicated.

        \end{itemize}

        Drawbacks:
        \begin{itemize}
        \item Slower over R when the degree of both of polynomials is less
              than 250 (roughly).
        \item Over R, results may not be as accurate as the Karatsuba
              case. This is because we represent coefficients of
              polynomials over R as fractions, then convert them back to
              floating-point numbers.
        \end{itemize}

        AUTHOR:
           -- Didier Deshommes (2006-05-25)
        """
        return self.parent()(polynomial_fateman._mul_fateman_mul(self,right))

    def _mul_karatsuba(self, right):
        r"""
        Returns the product of two polynomials using the Karatsuba
        divide and conquer multiplication algorithm.  This is only
        used over a generic base ring.  (Special libraries like NTL
        are used, e.g., for the integers and rationals, which are much
        faster.)

        INPUT:
           self: Polynomial
           right: Polynomial (over same base ring as self)

        OUTPUT: Polynomial
           The product self*right.

        ALGORITHM:
           The basic idea is to use that
           $$
               (aX + b) (cX + d) = acX^2 + ((a+b)(c+d)-ac-bd)X + bd
           $$
           where ac=a*c and bd=b*d, which requires three
           multiplications instead of the naive four.  (In my examples,
           strangely just doing the above with four multiplications
           does tend to speed things up noticeably.)
           Given f and g of arbitrary degree bigger than one, let e
           be min(deg(f),deg(g))/2.  Write
           $$
                  f = a X^e + b   \text{ and }   g = c X^e + d
           $$
           and use the identity
           $$
                 (aX^e + b) (cX^e + d) = ac X^{2e} +((a+b)(c+d) - ac - bd)X^e + bd
           $$
           to recursively compute $fg$.

        TIMINGS:
        On a Pentium M 1.8Ghz laptop:
           f=R.random(1000,bound=100)
           g=R.random(1000,bound=100)
           time h=f._mul_karatsuba(g)
           Time: 0.42 seconds
           The naive multiplication algorithm takes 14.58 seconds.
           In contrast, MAGMA does this sort of product almost
           instantly, and can easily deal with degree 5000.  Basically
           MAGMA is 100 times faster at polynomial multiplication.

           Over Z using NTL, multiplying two polynomials constructed
           using R.random(10000,bound=100) takes 0.10 seconds.  Using
           MAGMA V2.11-10 the same takes 0.14 seconds.  So in this
           case NTL is somewhat faster than MAGMA.

           Over Q using PARI, multiplying two polynomials constructed
           using R.random(10000,bound=100) takes 1.23 seconds.  Not
           good!  TODO: use NTL polynomials over Z with a denominator
           instead of PARI.

        NOTES:
         * Karatsuba multiplication of polynomials is also implemented in PARI in
                src/basemath/polarit3.c
         * The MAGMA documentation appears to give no information about how
           polynomial multiplication is implemented.
        """
        return self.parent()(do_karatsuba(self.list(), right.list()))

    def base_ring(self):
        """
        Return the base ring of the parent of self.

        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: x.base_ring()
            Integer Ring
            sage: (2*x+3).base_ring()
            Integer Ring
        """
        return self.parent().base_ring()

    def base_extend(self, R):
        """
        Return a copy of this polynomial but with coefficients in R, if there
        is a natural map from coefficient ring of self to R.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = x^3 - 17*x + 3
            sage: f.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no such base extension
            sage: f.change_ring(GF(7))
            x^3 + 4*x + 3
        """
        S = self.parent().base_extend(R)
        return S(self)

    def change_ring(self, R):
        """
        Return a copy of this polynomial but with coefficients in R, if at
        all possible.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: f = K.defining_polynomial()
            sage: f.change_ring(GF(7))
            x^2 + x + 1
        """
        S = self.parent().change_ring(R)
        return S(self)

    def __copy__(self):
        """
        Return a "copy" of self.  This is just self, since in SAGE polynomials are
        immutable this just returns self again.

        EXAMPLES:
        We create the polynomial $f=x+3$, then note that the copy is just
        the same polynomial again, which is fine since polynomials are immutable.

            sage: x = ZZ['x'].0
            sage: f = x + 3
            sage: g = copy(f)
            sage: g is f
            True
        """
        return self

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.

        EXAMPLES:
            sage: x = ZZ['x'].0
            sage: f = x^93 + 2*x + 1
            sage: f.degree()
            93
            sage: x = PolynomialRing(QQ, 'x', sparse=True).0
            sage: f = x^100000
            sage: f.degree()
            100000

            sage: x = QQ['x'].0
            sage: f = 2006*x^2006 - x^2 + 3
            sage: f.degree()
            2006
            sage: f = 0*x
            sage: f.degree()
            -1
            sage: f = x + 33
            sage: f.degree()
            1

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        raise NotImplementedError

    def denominator(self):
        """
        Return the least common multiple of the denominators of
        the entries of self, when this makes sense, i.e., when the
        coefficients have a denominator function.

        WARNING: This is not the denominator of the rational function
        defined by self, which would always be 1 since self is a polynomial.

        EXAMPLES:
        First we compute the denominator of a polynomial with integer
        coefficients, which is of course 1.
            sage: R.<x> = ZZ[]
            sage: f = x^3 + 17*x + 1
            sage: f.denominator()
            1

        Next we compute the denominator of a polynomial with rational coefficients.
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = (1/17)*x^19 - (2/3)*x + 1/3; f
            1/17*x^19 - 2/3*x + 1/3
            sage: f.denominator()
            51

        Finally, we try to compute the denominator of a polynomial with
        coefficients in the real numbers, which is a ring whose elements
        do not have a denominator method.
            sage: R.<x> = RR[]
            sage: f = x + RR('0.3'); f
            1.00000000000000*x + 0.300000000000000
            sage: f.denominator()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.real_mpfr.RealNumber' object has no attribute 'denominator'
        """
        if self.degree() == -1:
            return 1
        R = self.base_ring()
        x = self.list()
        d = x[0].denominator()
        for y in x:
            d = d.lcm(y.denominator())
        return d

    def derivative(self):
        return self.polynomial([self[n]*n for n in xrange(1,self.degree()+1)])

    def integral(self):
        try:
            return self.polynomial([0] + [self[n]/(n+1) for n in xrange(0,self.degree()+1)])
        except TypeError:
            raise ArithmeticError, "coefficients of integral cannot be coerced into the base ring"


    def dict(self):
        X = {}
        Y = self.list()
        for i in xrange(len(Y)):
            X[i] = Y[i]
        return X

    def factor(self):
        r"""
        Return the factorization of self.

        INPUT:
            a polynomial

        OUTPUT:
            Factorization -- the factorization of self, which is
            a product of a unit with a product of powers of irreducible
            factors.

        Over a field the irreducible factors are all monic.

        EXAMPLES:
        We factor some polynomials over $\Q$.
            sage: x = QQ['x'].0
            sage: f = (x^3 - 1)^2
            sage: f.factor()
            (x - 1)^2 * (x^2 + x + 1)^2

        Notice that over the field $\Q$ the irreducible factors are monic.
            sage: f = 10*x^5 - 1
            sage: f.factor()
            (10) * (x^5 - 1/10)
            sage: f = 10*x^5 - 10
            sage: f.factor()
            (10) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)

        Over $\Z$ the irreducible factors need not be monic:
            sage: x = ZZ['x'].0
            sage: f = 10*x^5 - 1
            sage: f.factor()
            10*x^5 - 1


        We factor a non-monic polynomial over the finite field $F_{25}$.
            sage: k.<a> = GF(25)
            sage: R.<x> = k[]
            sage: f = 2*x^10 + 2*x + 2*a
            sage: F = f.factor(); F
            (2) * (x + a + 2) * (x^2 + 3*x + 4*a + 4) * (x^2 + (a + 1)*x + a + 2) * (x^5 + (3*a + 4)*x^4 + (3*a + 3)*x^3 + 2*a*x^2 + (3*a + 1)*x + 3*a + 1)

        Notice that the unit factor is included when we multiply $F$ back out.
            sage: F.mul()
            2*x^10 + 2*x + 2*a

        Factorization also works even if the variable of the finite field is nefariously
        labeled "x".
            sage: x = GF(3^2, 'a')['x'].0
            sage: f = x^10 +7*x -13
            sage: G = f.factor(); G
            (x + a) * (x + 2*a + 1) * (x^4 + (a + 2)*x^3 + (2*a + 2)*x + 2) * (x^4 + 2*a*x^3 + (a + 1)*x + 2)
            sage: prod(G) == f
            True

            sage: f.parent().base_ring()._assign_names(['a'])
            sage: f.factor()
            (x + a) * (x + 2*a + 1) * (x^4 + (a + 2)*x^3 + (2*a + 2)*x + 2) * (x^4 + 2*a*x^3 + (a + 1)*x + 2)

            sage: k = GF(9,'x')    # purposely calling it x to test robustness
            sage: x = PolynomialRing(k,'x0').gen()
            sage: f = x^3 + x + 1
            sage: f.factor()
            (x0 + 2) * (x0 + x) * (x0 + 2*x + 1)
            sage: f = 0*x
            sage: f.factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

            sage: f = x^0
            sage: f.factor()
            1

        Arbitrary precision real and complex factorization:
            sage: R.<x> = RealField(100)[]
            sage: F = factor(x^2-3); F
            (1.0000000000000000000000000000*x + 1.7320508075688772935274463415) * (1.0000000000000000000000000000*x - 1.7320508075688772935274463415)
            sage: F.mul()
            1.0000000000000000000000000000*x^2 - 3.0000000000000000000000000000
            sage: factor(x^2 + 1)
            1.0000000000000000000000000000*x^2 + 1.0000000000000000000000000000
            sage: C = ComplexField(100)
            sage: R.<x> = C[]
            sage: F = factor(x^2+3); F
            (1.0000000000000000000000000000*x + -1.7320508075688772935274463415*I) * (1.0000000000000000000000000000*x + 1.7320508075688772935274463415*I)
            sage: F.mul()
            1.0000000000000000000000000000*x^2 + 3.0000000000000000000000000000
            sage: factor(x^2+1)
            (1.0000000000000000000000000000*x + -1.0000000000000000000000000000*I) * (1.0000000000000000000000000000*x + 1.0000000000000000000000000000*I)
            sage: f = C.0 * (x^2 + 1) ; f
            1.0000000000000000000000000000*I*x^2 + 1.0000000000000000000000000000*I
            sage: F=factor(f); F
            (1.0000000000000000000000000000*I) * (1.0000000000000000000000000000*x + -1.0000000000000000000000000000*I) * (1.0000000000000000000000000000*x + 1.0000000000000000000000000000*I)
            sage: F.mul()
            1.0000000000000000000000000000*I*x^2 + 1.0000000000000000000000000000*I
        """

        # PERFORMANCE NOTE:
        #     In many tests with SMALL degree PARI is substantially
        #     better than NTL.  (And magma is better yet.)  And the
        #     timing difference has nothing to do with moving Python
        #     data to NTL and back.
        #     For large degree ( > 1500) in the one test I tried, NTL was
        #     *much* better than MAGMA, and far better than PARI.  So probably
        #     NTL's implementation is asymptotically better.  I could use
        #     PARI for smaller degree over other rings besides Z, and use
        #     NTL in general.

        R = self.parent().base_ring()
        if self.degree() < 0:
            raise ValueError, "factorization of 0 not defined"
        G = None

        from sage.rings.number_field.all import is_NumberField
        from sage.rings.finite_field import is_FiniteField

        n = None
        if integer_mod_ring.is_IntegerModRing(R) or \
              integer_ring.is_IntegerRing(R) or \
              rational_field.is_RationalField(R):

            G = list(self._pari_('x').factor())

        elif is_NumberField(R) or is_FiniteField(R):

            v = [x._pari_("a") for x in self.list()]
            f = pari(v).Polrev()
            G = list(f.factor())

        elif is_RealField(R):
            n = pari.set_real_precision(int(3.5*R.prec()) + 1)
            G = list(self._pari_('x').factor())

        elif complex_field.is_ComplexField(R):
            # This is a hack to make the polynomial have complex coefficients, since
            # otherwise PARI will factor over RR.
            n = pari.set_real_precision(int(3.5*R.prec()) + 1)
            if self.leading_coefficient() != R.gen():
                G = list((pari(R.gen())*self._pari_('x')).factor())
            else:
                G = self._pari_('x').factor()

        #elif padic_field.is_pAdicField(R):
        #    G = list(self._pari_('x').factorpadic(R.prime(), R.prec()))

        if G is None:
            raise NotImplementedError

        return self._factor_pari_helper(G, n)

    def _factor_pari_helper(self, G, n=None, unit=None):
        pols = G[0]
        exps = G[1]
        F = []
        R = self.parent()
        c = R.base_ring()(1)
        for i in xrange(len(pols)):
            f = R(pols[i])
            e = int(exps[i])
            if unit is None:
                c *= f.leading_coefficient()
            F.append((f,e))

        if unit is None:

            unit = R.base_ring()(self.leading_coefficient()/c)

        if not unit.is_unit():

            F.append((R(unit), 1))
            unit = R.base_ring()(1)

        elif R.base_ring().is_field():
            # When the base ring is a field we normalize
            # the irreducible factors so they have leading
            # coefficient 1.
            one = R.base_ring()(1)
            for i in range(len(F)):
                c = F[i][0].leading_coefficient()
                if c != 1:
                    unit *= c
                    F[i] = (F[i][0].monic(), F[i][1])

        if not n is None:
            pari.set_real_precision(n)  # restore precision
        return Factorization(F, unit)

    def _lcm(self, other):
        """
        Let f and g be two polynomials.  Then this function
        returns the monic least common multiple of f and g.
        """
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q.leading_coefficient())*q  # make monic  (~ is inverse in python)

    def is_constant(self):
        return self.degree() <= 0

    def root_field(self, names, check_irreducible=True):
        """
        Return the field generated by the roots of self.  The output
        is either a number field, relative number field, a quotient of
        a polynomial ring over a field, or the fraction field of the
        base ring.

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('a')
            Number Field in a with defining polynomial x^3 + x + 17

            sage: R.<x> = QQ['x']
            sage: f = x - 3
            sage: f.root_field('b')
            Rational Field

            sage: R.<x> = ZZ['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('b')
            Number Field in b with defining polynomial x^3 + x + 17

            sage: y = QQ['x'].0
            sage: L.<a> = NumberField(y^3-2)
            sage: R.<x> = L['x']
            sage: f = x^3 + x + 17
            sage: f.root_field('c')
            Extension by x^3 + x + 17 of the Number Field in a with defining polynomial x^3 - 2

            sage: R.<x> = PolynomialRing(GF(9,'a'))
            sage: f = x^3 + x^2 + 8
            sage: K.<alpha> = f.root_field(); K
            Univariate Quotient Polynomial Ring in alpha over Finite Field in a of size 3^2 with modulus x^3 + x^2 + 2
            sage: alpha^2 + 1
            alpha^2 + 1
            sage: alpha^3 + alpha^2
            1
        """
        from sage.rings.number_field.number_field import is_NumberField, NumberField

        R = self.base_ring()
        if not is_IntegralDomain(R):
            raise ValueError, "the base ring must be a domain"

        if self.degree() <= 1:
            return R.fraction_field()

        if integer_ring.is_IntegerRing(R):
            return NumberField(self, names)


        if rational_field.is_RationalField(R) or is_NumberField(R):
            return NumberField(self, names)

        if check_irreducible and not self.is_irreducible():
            raise ValueError, "polynomial must be irreducible"

        return sage.rings.polynomial_ring.PolynomialRing(R.fraction_field(),
                              self.parent().variable_name()).quotient(self, names)


    def constant_coefficient(self):
        return self[0]

    def is_monic(self):
        """
        Returns True if this polynomial is monic.  The zero
        polynomial is by definition not monic.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = x + 33
            sage: f.is_monic()
            True
            sage: f = 0*x
            sage: f.is_monic()
            False
            sage: f = 3*x^3 + x^4 + x^2
            sage: f.is_monic()
            True
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.is_monic()
            False

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return bool(not self.is_zero() and self[self.degree()] == 1)

    def is_unit(self):
        r"""
        Return True if this polynomial is a unit.

        EXAMPLES:
            sage: a = Integers(90384098234^3)
            sage: b = a(2*191*236607587)
            sage: b.is_nilpotent()
            True
            sage: R.<x> = a[]
            sage: f = 3 + b*x + b^2*x^2
            sage: f.is_unit()
            True
            sage: f = 3 + b*x + b^2*x^2 + 17*x^3
            sage: f.is_unit()
            False

        EXERCISE (Atiyah-McDonald, Ch 1): Let $A[x]$ be a polynomial
        ring in one variable.  Then $f=\sum a_i x^i \in A[x]$ is a
        unit if and only if $a_0$ is a unit and $a_1,\ldots, a_n$ are
        nilpotent.
        """
        if self.degree() > 0:
            for i in range(1,self.degree()+1):
                if not self[i].is_nilpotent():
                    return False
            return True
        return self[0].is_unit()

    def is_nilpotent(self):
        r"""
        Return True if this polynomial is nilpotent.

        EXAMPLES:
            sage: R = Integers(12)
            sage: S.<x> = R[]
            sage: f = 5 + 6*x
            sage: f.is_nilpotent()
            False
            sage: f = 6 + 6*x^2
            sage: f.is_nilpotent()
            True
            sage: f^2
            0

        EXERCISE (Atiyah-McDonald, Ch 1): Let $A[x]$ be a polynomial
        ring in one variable.  Then $f=\sum a_i x^i \in A[x]$ is
        nilpotent if and only if every $a_i$ is nilpotent.
        """
        for i in range(self.degree()+1):
            if not self[i].is_nilpotent():
                return False
        return True

    def is_gen(self):
        return bool(self._is_gen)

    def is_zero(self):
        return bool(self.degree() == -1)

    def leading_coefficient(self):
        return self[self.degree()]

    def monic(self):
        """
        Return this polynomial divided by its leading coefficient.
        Does not change this polynomial.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = 2*x^2 + x^3 + 56*x^5
            sage: f.monic()
            x^5 + 1/56*x^3 + 1/28*x^2
            sage: f = (1/4)*x^2 + 3*x + 1
            sage: f.monic()
            x^2 + 12*x + 4

    The following happens because $f = 0$ cannot be made into a monic polynomial
            sage: f = 0*x
            sage: f.monic()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        Notice that the monic version of a polynomial over the
        integers is defined over the rationals.
            sage: x = ZZ['x'].0
            sage: f = 3*x^19 + x^2 - 37
            sage: g = f.monic(); g
            x^19 + 1/3*x^2 - 37/3
            sage: g.parent()
            Univariate Polynomial Ring in x over Rational Field


        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        if self.is_monic():
            return self
        a = ~self.leading_coefficient()
        R = self.parent()
        if a.parent() != R.base_ring():
            S = R.base_extend(a.parent())
            return a*S(self)
        else:
            return a*self


    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.
        """
        raise NotImplementedError

    def coeffs(self):
        r"""
        Returns \code{self.list()}.

        (It potentially slightly faster better to use
        \code{self.list()} directly.)

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = 10*x^3 + 5*x + 2/17
            sage: f.coeffs()
            [2/17, 5, 0, 10]
        """
        return self.list()

    def newton_raphson(self, n, x0):
        """
        Return a list of n iterative approximations to a root of this
        polynomial, computed using the Newton-Raphson method.

        The Newton-Raphson method is an iterative root-finding algorithm.
        For f(x) a polynomial, as is the case here, this is essentially
        the same as Horner's method.

        INPUT:
           n -- an integer (=the number of iterations),
           x0 -- an initial guess x0.

        OUTPUT:
           A list of numbers hopefully approximating a root of f(x)=0.

           ** If one of the iterates is a critical point of f then
              a ZeroDivisionError exception is raised.

        EXAMPLES:
            sage: x = PolynomialRing(RealField(), 'x').gen()
            sage: f = x^2 - 2
            sage: f.newton_raphson(4, 1)
            [1.50000000000000, 1.41666666666667, 1.41421568627451, 1.41421356237469]

        AUTHORS: David Joyner and William Stein (2005-11-28)
        """
        n = integer.Integer(n)
        df = self.derivative()
        K = self.parent().base_ring()
        a = K(x0)
        L = []
        for i in range(n):
            a -= self(a) / df(a)
            L.append(a)
        return L

    def polynomial(self, *args, **kwds):
        return self.parent()(*args, **kwds)

    def newton_slopes(self, p):
        """
        Return the $p$-adic slopes of the Newton polygon of self,
        when this makes sense.

        OUTPUT:
            -- list of rational numbers

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = x^3 + 2
            sage: f.newton_slopes(2)
            [1/3, 1/3, 1/3]

        ALGORITHM: Uses PARI.
        """
        f = self._pari_()
        v = list(f.newtonpoly(p))
        return [sage.rings.rational.Rational(x) for x in v]


    #####################################################################
    # Conversions to other systems
    #####################################################################
    def _pari_(self, variable=None):
        """
        Return polynomial as a PARI object.

        EXAMPLES:
            sage: f = PolynomialRing(QQ, 'X')([0,1,2/3,3])
            sage: pari(f)
            3*X^3 + 2/3*X^2 + X
        """
        try:
            return self.__pari
        except AttributeError:
            K = self.base_ring()
            n = None
            if is_RealField(K) or complex_field.is_ComplexField(K):
                n = pari.get_real_precision()
                pari.set_real_precision(int(K.prec()*3.5)+1)
            v = self.list()
            try:
                v = [x._pari_() for x in v]
            except AttributeError:
                pass
            if variable is None:
                variable = self.parent().variable_name()
            self.__pari = pari(v).Polrev(variable)
            if not n is None:
                pari.set_real_precision(n)
            return self.__pari

    def _pari_init_(self):
        return str(self._pari_())

    def _magma_init_(self):
        """
        Return a string that evaluates in Magma to this polynomial.

        EXAMPLES:
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: f._magma_init_()
            'Polynomial(IntegerRing(), [5,-17,0,1])'
        """
        return 'Polynomial(%s, [%s])'%(self.base_ring()._magma_init_(), ','.join([a._magma_init_() for a in self.list()]))

    def _magma_(self, G=None):
        """
        Return the Magma version of this polynomial.

        EXAMPLES:
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: g = magma(f); g              # optional -- requires Magma
            y^3 - 17*y + 5

        Note that in Magma there is only one polynomial ring over each base,
        so if we make the polynomial ring over ZZ with variable $z$, then
        this changes the variable name of the polynomial we already defined:
            sage: R.<z> = ZZ[]
            sage: magma(R)                     # optional -- requires Magma
            Univariate Polynomial Ring in z over Integer Ring
            sage: g                            # optional -- requires Magma
            z^3 - 17*z + 5

        In SAGE the variable name does not change:
            sage: f
            y^3 - 17*y + 5
        """
        if G is None:
            import sage.interfaces.magma
            G = sage.interfaces.magma.magma
        self.parent()._magma_(G)  # defines the variable name
        f = G(self._magma_init_())
        return f

    def _gap_init_(self):
        return str(self)

    def _gap_(self, G):
        """
        EXAMPLES:
            sage: R.<y> = ZZ[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g
            y^3-17*y+5
            sage: f._gap_init_()
            'y^3 - 17*y + 5'
            sage: R.<z> = ZZ[]
            sage: gap(R)
            PolynomialRing(..., [ z ])
            sage: g
            y^3-17*y+5
            sage: gap(z^2 + z)
            z^2+z

        We coerce a polynomial with coefficients in a finite field:

            sage: R.<y> = GF(7)[]
            sage: f = y^3 - 17*y + 5
            sage: g = gap(f); g
            y^3+Z(7)^4*y+Z(7)^5
            sage: g.Factors()
            [ y+Z(7)^0, y+Z(7)^0, y+Z(7)^5 ]
            sage: f.factor()
            (y + 5) * (y + 1)^2
        """
        if G is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
        self.parent()._gap_(G)
        return G(self._gap_init_())

    ######################################################################


    def resultant(self, other, flag=0):
        raise NotImplementedError

        ## This should be switched to use NTL, which can apparently compute
        ## resultants!
        ##        void XGCD(ZZ& r, ZZX& s, ZZX& t, const ZZX& a, const ZZX& b,
        ##          long deterministic=0);
        ##// r = resultant of a and b; if r != 0, then computes s and t such
        ##// that: a*s + b*t = r; otherwise s and t not affected.  if
        ##// !deterministic, then resultant computation may use a randomized
        ##// strategy that errs with probability no more than 2^{-80}.
        #m = magma.Magma()
        #cmd = "R<%s> := PolynomialRing(RationalField()); "%self.parent().variable_name() + \
        #      "Resultant(%s, %s);"%(self,other)
        #s = m.cmd(cmd)
        #i = s.find("\r")
        #return eval(s[:i])

    def reverse(self):
        v = list(self.list())
        v.reverse()
        return self.parent()(v)

    def roots(self, multiplicities=True):
        """
        Return all roots of this polynomial.

        INPUT:
            multiplicities -- bool (default: True, except over RR or CC)
                              if True return list of pairs (r, n), where r is
                              the root and n is the multiplicity.

        If the polynomial is over RR or CC returns all roots in CC
        with multiplicities all set to 1.

        Over all other rings it just returns the roots that lie in the
        base ring.

        EXAMPLES:
            sage: x = QQ['x'].0
            sage: f = x^3 - 1
            sage: f.roots()
            [(1, 1)]
            sage: f = (x^3 - 1)^2
            sage: f.roots()
            [(1, 2)]

            sage: f = -19*x + 884736
            sage: f.roots()
            [(884736/19, 1)]
            sage: (f^20).roots()
            [(884736/19, 20)]

            sage: K.<z> = CyclotomicField(3)
            sage: f = K.defining_polynomial()
            sage: g = f.change_ring(GF(7))
            sage: g.roots()
            [(4, 1), (2, 1)]
            sage: g.roots(multiplicities=False)
            [4, 2]

            sage: x = RR['x'].0
            sage: f = x^2 - 1e100
            sage: f.roots()
            [-1.00000000000000e50, 1.00000000000000e50]
            sage: f = x^10 - 2*(5*x-1)^2
            sage: f.roots()
            [-1.67726703399418, 0.199954796285057, 0.200045306115242, 1.57630351618444, 1.10042307611716 + 1.15629902427493*I, 1.10042307611716 - 1.15629902427493*I, -1.20040047425969 + 1.15535958549432*I, -1.20040047425969 - 1.15535958549432*I, -0.0495408941527470 + 1.63445468021367*I, -0.0495408941527470 - 1.63445468021367*I]

            sage: x = CC['x'].0
            sage: f = (x-1)*(x-I)
            sage: f.roots()
            [1.00000000000000*I, 1.00000000000000]
        """
        seq = []

        K = self.parent().base_ring()

        if is_RealField(K) or complex_field.is_ComplexField(K):
            if is_RealField(K):
                K = K.complex_field()
            n = pari.get_real_precision()
            pari.set_real_precision(int(K.prec()/3.2)+1)
            r = pari(self).polroots()
            seq = [K(root) for root in r]
            pari.set_real_precision(n)
            return seq

        try:
            rts = self.factor()
        except NotImplementedError:
            raise NotImplementedError, "root finding for this polynomial not implemented"
        for fac in rts:
            g = fac[0]
            if g.degree() == 1:
                if multiplicities:
                    seq.append((-g[0]/g[1],fac[1]))
                else:
                    seq.append(-g[0]/g[1])
        return seq

    def valuation(self, p=None):
        r"""
        If $f = a_r x^r + a_{r+1}x^{r+1} + \cdots$, with $a_r$ nonzero,
        then the valuation of $f$ is $r$.  The valuation of the zero
        polynomial is $\infty$.

        If a prime (or non-prime) $p$ is given, then the valuation is
        the largest power of $p$ which divides self.

        The valuation at $\infty$ is -self.degree().

        EXAMPLES:
        sage: P,x=PolynomialRing(ZZ,'x').objgen()
        sage: (x^2+x).valuation()
        1
        sage: (x^2+x).valuation(x+1)
        1
        sage: (x^2+1).valuation()
        0
        sage: (x^3+1).valuation(infinity)
        -3
        sage: P(0).valuation()
        +Infinity
        """
        cdef int k
        if self.is_zero():
            return infinity
        if p == infinity:
            return -self.degree()
        if p is None:
            p = self.parent().gen()
        if not isinstance(p, Polynomial) or not p.parent() is self.parent():
            raise TypeError, "The polynomial, p, must have the same parent as self."
        if p is None or p == self.parent().gen():
            for i in xrange(self.degree()+1):
                if self[i] != 0:
                    return ZZ(i)
        else:
            if p.degree() == 0:
                raise ArithmeticError, "The polynomial, p, must have positive degree."
            k = 0
            while self % p == 0:
                k = k + 1
                self = self.__floordiv__(p)
            return integer.Integer(k)
        raise RuntimeError, "bug in computing valuation of polynomial"

    def ord(self, p=None):
        """Synonym for valuation

        EXAMPLES:
        sage: P,x=PolynomialRing(ZZ,'x').objgen()
        sage: (x^2+x).ord(x+1)
        1
        """
        return self.valuation(p)

    def name(self):
        return self.parent().variable_name()

    def _xgcd(self, other):
        r"""
        Extended gcd of self and polynomial other.

        Returns g, u, and v such that
              \code{g = u*self + v*other.}

        EXAMPLES:
            sage: P.<x> = QQ[]
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
        """
        if other.is_zero():
            R = self.parent()
            return R(1), R(0), self
        # Algorithm 3.2.2 of Cohen, GTM 138
        R = self.parent()
        A = self
        B = other
        U = R(1)
        G = A
        V1 = R(0)
        V3 = B
        while not V3.is_zero():
            Q, R = G.quo_rem(V3)
            T = U - V1*Q
            U = V1
            G = V3
            V1 = T
            V3 = R
        V = (G-A*U)//B
        return G, U, V

    def is_irreducible(self):
        """
        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: (x^3 + 1).is_irreducible()
            True
            sage: (x^2 - 1).is_irreducible()
            False
            sage: R(0).is_irreducible()
            ..?

        $4$ is irreducible as a polynomial, since as a polynomial
        it doesn't factor:
            sage: R(4).is_irreducible()
            True
        """
        if self.is_zero():
            raise ValueError, "self must be nonzero"
        if self.degree() == 0:
            return True

        F = self.factor()
        if len(F) > 1 or F[0][1] > 1:
            return False
        return True

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power $x^n$. If $n$ is negative,
        terms below $x^n$ will be discarded. Does not change this polynomial (since
        polynomials are immutable).

        EXAMPLES:
            sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'w'),'x')
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
            return self   # safe because immutable.
        if n > 0:
            output = [self.base_ring()(0)] * n
            output.extend(self.coeffs())
            return self.polynomial(output, check=False)
        if n < 0:
            if n > self.degree():
                return self.polynomial([])
            else:
                return self.polynomial(self.coeffs()[-int(n):], check=False)

    def truncate(self, n):
        r"""
        Returns the polynomial of degree $ < n$ which is equivalent to self
        modulo $x^n$.
        """
        return self.parent()(self[:int(n)], check=False)

    def radical(self):
        """
        Returns the radical of self; over a field, this is the product of the
        distinct irreducible factors of self. (This is also sometimes called the
        "square-free part" of self, but that term is ambiguous; it is sometimes used
        to mean the quotient of self by its maximal square factor.)

        EXAMPLES:
            sage: P.<x> = ZZ[]
            sage: t = (x^2-x+1)^3 * (3*x-1)^2
            sage: t.radical()
            3*x^3 - 4*x^2 + 4*x - 1
        """
        return self // self.gcd(self.derivative())


# ----------------- inner functions -------------
# Sagex can't handle function definitions inside other function


def _karatsuba_sum(v,w):
    if len(v)>=len(w):
        x = list(v)
        y = w
    else:
        x = list(w)
        y = v
    for i in range(len(y)):
        x[i] = x[i] + y[i]
    return x

def _karatsuba_dif(v,w):
    if len(v)>=len(w):
        x = list(v)
        y = w
    else:
        x = list(w)
        y = v
    for i in range(len(y)):
        x[i] -= y[i]
    return x

def do_karatsuba(left, right):
    if len(left) == 0 or len(right) == 0:
        return []
    if len(left) == 1:
        return [left[0]*a for a in right]
    if len(right) == 1:
        return [right[0]*a for a in left]
    if len(left) == 2 and len(right) == 2:
        b = left[0]
        a = left[1]
        d = right[0]
        c = right[1]
        ac = a*c
        bd = b*d
        return [bd,(a+b)*(c+d)-ac-bd,ac]
    e = min(len(left), len(right))/2
    assert e>=1, "bug in karatsuba"
    a, b = left[e:], left[:e]
    c, d = right[e:], right[:e]
    ac = do_karatsuba(a,c)
    bd = do_karatsuba(b,d)
    zeros = [0] * e
    t2 = zeros + zeros + ac
    t1 = zeros + _karatsuba_dif(do_karatsuba(_karatsuba_sum(a,b),_karatsuba_sum(c,d)),_karatsuba_sum(ac,bd))
    t0 = bd
    return _karatsuba_sum(t0,_karatsuba_sum(t1,t2))
