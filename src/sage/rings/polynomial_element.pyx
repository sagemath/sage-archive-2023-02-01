"""
Univariate Polynomials

AUTHORS:
    -- William Stein: first version
    -- Martin Albrecht: Added singular coercion.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

import copy

#from sage.structure.element import Element, IntegralDomainElement, CommutativeAlgebraElement
from sage.structure.element import EuclideanDomainElement, PrincipalIdealDomainElement, RingElement
import sage.rings.rational_field
import sage.rings.integer_ring
import sage.rings.rational
import integer
import finite_field
import padic_field
import sage.rings.polynomial_ring
import arith
import sage.rings.ring_element as ring_element
import integer_ring
import integer_mod_ring
import polynomial_pyx
import rational_field
import complex_field
import fraction_field_element
import fraction_field
from infinity import infinity
import sage.misc.misc as misc
from sage.misc.sage_eval import sage_eval
from sage.libs.all import pari, pari_gen
from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZX, ZZX_class, ZZ_p, ZZ_pX, ZZ_pX_class, set_modulus
import sage.misc.latex as latex
from sage.structure.factorization import Factorization

from sage.interfaces.all import singular as singular_default, is_SingularElement

from sage.rings.polynomial_singular_interface import Polynomial_singular_repr

from real_field import RealField, is_RealNumber, is_RealField
RR = RealField()

# Faster than SAGE's
from math import log as pylog
from math import ceil as pyceil
from math import floor as pyfloor
from arith import gcd

QQ = sage.rings.rational_field.RationalField()

ZZ = sage.rings.integer_ring.IntegerRing()

def is_Polynomial(f):
    return isinstance(f, Polynomial)

cdef class Polynomial(CommutativeAlgebraElement):
    """
    A polynomial.

    EXAMPLE:
        sage: R.<y> = QQ['y']
        sage: S.<x> = R['x']
        sage: f = x*y; f
        y*x
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_dense'>
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

    def __cmp__(self, other):
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
            sage: R(-1) < R(0)
            False
        """
        d1 = self.degree(); d2 = other.degree()
        c = cmp(d1, d2)
        if c: return c
        for i in reversed(xrange(d1+1)):
            c = cmp(self[i], other[i])
            if c: return c
        return 0

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
        return self.parent()(1)/self

    def inverse_of_unit(self):
        if self.degree() > 0:
            raise ValueError, "self is not a unit."
        return self.parent()(~(self[0]))

    def __long__(self):
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
        return self._repr()

    def _latex_(self, name=None):
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
                x = latex.latex(x)
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
        return self.degree() == self.valuation()

    def _mul_generic(self, right):
        d1 = self.degree()
        d2 = right.degree()
        d = d1 + d2
        w = [sum([self[i]*right[k-i] for i in range(0,min(d1,k)+1) if \
                  i <= d1 and k-i <= d2 and self[i]!=0 and right[k-i]!=0]) \
                for k in range(d+1)]
        return Polynomial(self.parent(), w)

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
            1.00000000000000*y^20 - 2.78698600000000*y^11 + 0.599999999999999*y^10 - 0.000000000000000222044604925031*y^8 + 0.000000000000000222044604925031*y^6 + 0.000000000000000222044604925031*y^3 + 1.94182274104900*y^2 - 0.836095799999999*y + 0.0899999999999999
            sage: f._mul_fateman(f)
            1.00000000000000*y^20 - 2.78698600000000*y^11 + 0.599999999999999*y^10 + 1.94182274104899*y^2 - 0.836095800000000*y + 0.0899999999999999

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
        return self.parent()(_mul_fateman_mul(self,right))

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
        Return a copy of this polynomial but with coefficients in R.
        """
        S = sage.rings.polynomial_ring.PolynomialRing(R,
                                      names = self.parent().variable_name())
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
            1.00000000000000*x + 0.299999999999999
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
        from sage.rings.integer_ring import IntegerRing
        from sage.rings.rational_field import RationalField
        from sage.rings.finite_field import is_FiniteField

        if integer_mod_ring.is_IntegerModRing(R) or \
               isinstance(R, (IntegerRing, RationalField)):

            G = list(self._pari_('x').factor())

        elif is_NumberField(R) or is_FiniteField(R):

            v = [x._pari_("a") for x in self.list()]
            f = pari(v).Polrev()
            G = list(f.factor())


        if G is None:
            raise NotImplementedError

        return self._factor_pari_helper(G)

    def _factor_pari_helper(self, G, unit=None):
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

        return Factorization(F, unit)

    def _lcm(self, other):
        """
        Let f and g be two polynomials.  Then this function
        returns the monic least common multiple of f and g.
        """
        f = self*other
        g = self.gcd(other)
        q = f//g
        return ~(q[q.degree()])*q  # make monic  (~ is inverse in python)

    def is_constant(self):
        return self.degree() <= 0

    def root_field(self, names):
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
        from integral_domain import is_IntegralDomain
        from rational_field import is_RationalField
        from sage.rings.number_field.number_field import is_NumberField, NumberField

        R = self.base_ring()
        if not is_IntegralDomain(R):
            raise ValueError, "the base ring must be a domain"

        if self.degree() <= 1:
            return R.fraction_field()

        if isinstance(R, sage.rings.integer_ring.IntegerRing):
            return NumberField(self, names)


        if is_RationalField(R) or is_NumberField(R):
            return NumberField(self, names)

        if not self.is_irreducible():
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
            [1.50000000000000, 1.41666666666666, 1.41421568627450, 1.41421356237468]

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
            sage: g = f.base_extend(GF(7))
            sage: g.roots()
            [(4, 1), (2, 1)]
            sage: g.roots(multiplicities=False)
            [4, 2]
        """
        seq = []

        K = self.parent().base_ring()

        if is_RealField(K) or complex_field.is_ComplexField(K):
            if is_RealField(K):
                K = K.complex_field()
            n = pari.get_real_precision()
            pari.set_real_precision(int(K.prec()/3.2)+1)
            r = pari(self).polroots()
            r = str(r).rstrip('~')
            seq = sage_eval(r, locals={'I':K.gen(), 'RealNumber':K._real_field()})
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

    def valuation(self):
        r"""
        If $f = a_r x^r + a_{r+1}x^{r+1} + \cdots$, with $a_r$ nonzero,
        then the valuation of $f$ is $r$.  The valuation of the zero
        polynomial is $\infty$.
        """
        if self.is_zero():
            return infinity
        for i in xrange(self.degree()+1):
            if self[i] != 0:
                return i
        raise RuntimeError, "bug in computing valuation of polynomial"

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

class Polynomial_generic_dense(Polynomial):
    """
    A generic dense polynomial.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'y'))
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_dense'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = []
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = list(x.list())
            elif x.parent() == R:
                x = [x]
            else:
                x = [R(a) for a in x.list()]
                check = False
        elif isinstance(x, dict):
            zero = R(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.Vecrev()]
            check = True
        elif not isinstance(x, list):
            x = [x]   # constant polynomials
        if check:
            self.__coeffs = [R(z) for z in x]
        else:
            self.__coeffs = x
        if check:
            self.__normalize()

    def __normalize(self):
        x = self.__coeffs
        n = len(x)-1
        while n>=0 and x[n] == 0:
            del x[n]
            n -= 1

    def __getitem__(self,n):
        if n < 0 or n >= len(self.__coeffs):
            return self.base_ring()(0)
        return self.__coeffs[n]

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        return self.__coeffs[i:j]

    def _unsafe_mutate(self, n, value):
        if (<Polynomial>self)._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        value = self.base_ring()(value)
        if n >= 0 and n < len(self.__coeffs):
            self.__coeffs[n] = value
            if n == len(self.__coeffs) and value == 0:
                self.__normalize()
        elif n < 0:
            raise IndexError, "polynomial coefficient index must be nonnegative"
        elif value != 0:
            zero = self.base_ring()(0)
            for _ in xrange(len(self.__coeffs), n):
                self.__coeffs.append(zero)
            self.__coeffs.append(value)

    def __floordiv__(self, right):
        if right.parent() == self.parent():
            return Polynomial.__floordiv__(self, right)
        d = self.parent().base_ring()(right)
        return self.polynomial([c // d for c in self.__coeffs], check=false)

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.
        """
        return list(self.__coeffs)

    def degree(self):
        return len(self.__coeffs) - 1

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power $x^n$. If $n$ is negative,
        terms below $x^n$ will be discarded. Does not change this polynomial.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(PolynomialRing(QQ,'y'), 'x')
            sage: p = x^2 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial_element.Polynomial_generic_dense'>
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        AUTHOR:
            -- David Harvey (2006-08-06)
        """
        if n == 0:
            return self
        if n > 0:
            output = [self.base_ring()(0)] * n
            output.extend(self.__coeffs)
            return self.polynomial(output, check=False)
        if n < 0:
            if n > self.degree():
                return self.polynomial([])
            else:
                return self.polynomial(self.__coeffs[-int(n):], check=False)


class Polynomial_generic_sparse(Polynomial):
    """
    A generic sparse polynomial.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(PolynomialRing(QQ, 'y'), sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_sparse'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = {}
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = dict(x.dict())
            elif x.parent() == R:
                x = {0:x}
            else:
                w = {}
                for n, c in x.dict().iteritems():
                    w[n] = R(c)
                #raise TypeError, "Cannot coerce %s into %s."%(x, parent)
        elif isinstance(x, list):
            y = {}
            for i in xrange(len(x)):
                if x[i] != 0:
                    y[i] = x[i]
            x = y
        elif not isinstance(x, dict):
            x = {0:x}   # constant polynomials
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.Vecrev()]
            check = True
        if check:
            self.__coeffs = {}
            for i, z in x.iteritems():
                self.__coeffs[i] = R(z)
        else:
            self.__coeffs = x
        if check:
            self.__normalize()


    def _repr(self, name=None):
        r"""
        AUTHOR:
            -- David Harvey (2006-08-05), based on Polynomial._repr()
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        coeffs = list(self.__coeffs.iteritems())
        coeffs.sort()
        for (n, x) in reversed(coeffs):
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

    def __normalize(self):
        x = self.__coeffs
        zero = self.base_ring()(0)
        D = [n for n, z in x.iteritems() if z == 0]
        for n in D:
            del x[n]

    def __getitem__(self,n):
        if not self.__coeffs.has_key(n):
            return self.base_ring()(0)
        return self.__coeffs[n]

    def __getslice__(self, i, j):
        if i < 0:
            i = 0
        zero = self.base_ring()(0)
        v = [zero for _ in xrange(i,j)]
        x = self.__coeffs
        for k in set(x.keys()).intersection(set(xrange(i,j))):
            v[k] = x[k]
        return v

    def _unsafe_mutate(self, n, value):
        if (<Polynomial>self)._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        value = self.base_ring()(value)
        x = self.__coeffs
        if n < 0:
            raise IndexError, "polynomial coefficient index must be nonnegative"
        if value == 0:
            if x.has_key(n):
                del x[n]
        else:
            x[n] = value

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.
        """
        zero = self.base_ring()(0)
        v = [zero for _ in xrange(self.degree()+1)]
        for n, x in self.__coeffs.iteritems():
            v[n] = x
        return v

    #def _pari_(self, variable=None):
    #    if variable is None:
    #        return self.__pari
    #    else:
    #        return self.__pari.subst('x',variable)

    def degree(self):
        v = self.__coeffs.keys()
        if len(v) == 0:
            return -1
        return max(v)

    def _add_(self, right):
        r"""
        EXAMPLES:
            sage: R.<x> = PolynomialRing(Integers(), sparse=True)
            sage: (x^100000 + 2*x^50000) + (4*x^75000 - 2*x^50000 + 3*x)
             x^100000 + 4*x^75000 + 3*x

        AUTHOR:
            -- David Harvey (2006-08-05)
        """
        output = copy.copy(self.__coeffs)

        for (index, coeff) in right.__coeffs.iteritems():
            if index in output:
                output[index] += coeff
            else:
                output[index] = coeff

        output = self.polynomial(output, check=False)
        output.__normalize()
        return output

    def _mul_(self, right):
        r"""
        EXAMPLES:
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
             x^200000 - x^100000
            sage: (x^100000 - x^50000) * R(0)
             0

        AUTHOR:
            -- David Harvey (2006-08-05)
        """
        output = {}

        for (index1, coeff1) in self.__coeffs.iteritems():
            for (index2, coeff2) in right.__coeffs.iteritems():
                product = coeff1 * coeff2
                index = index1 + index2
                if index in output:
                    output[index] += product
                else:
                    output[index] = product

        output = self.polynomial(output, check=False)
        output.__normalize()
        return output

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power $x^n$. If $n$ is negative,
        terms below $x^n$ will be discarded. Does not change this polynomial.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^100000 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial_element.Polynomial_generic_sparse'>
            sage: p.shift(0)
             x^100000 + 2*x + 4
            sage: p.shift(-1)
             x^99999 + 2
            sage: p.shift(-100002)
             0
            sage: p.shift(2)
             x^100002 + 2*x^3 + 4*x^2

        AUTHOR:
            -- David Harvey (2006-08-06)
        """
        n = int(n)
        if n == 0:
            return self
        if n > 0:
            output = {}
            for (index, coeff) in self.__coeffs.iteritems():
                output[index + n] = coeff
            return self.polynomial(output, check=False)
        if n < 0:
            output = {}
            for (index, coeff) in self.__coeffs.iteritems():
                if index + n >= 0:
                    output[index + n] = coeff
            return self.polynomial(output, check=False)


class Polynomial_generic_domain(Polynomial, IntegralDomainElement):
    def __init__(self, parent, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

    def is_unit(self):
        r"""
        Return True if this polynomial is a unit.

        EXERCISE (Atiyah-McDonald, Ch 1): Let $A[x]$ be a polynomial
        ring in one variable.  Then $f=\sum a_i x^i \in A[x]$ is a
        unit if and only if $a_0$ is a unit and $a_1,\ldots, a_n$ are
        nilpotent.

        EXAMPLES:
        """
        if self.degree() > 0:
            return False
        return self[0].is_unit()

class Polynomial_generic_field(Polynomial_generic_domain,
                               EuclideanDomainElement,
                               Polynomial_singular_repr):
    def quo_rem(self, other):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.

        EXAMPLES:
            sage: R.<y> = PolynomialRing(QQ)
            sage: K.<t> = NumberField(y^2 - 2)
            sage: P.<x> = PolynomialRing(K)
            sage: x.quo_rem(K(1))
            (x, 0)
            sage: x.xgcd(K(1))
            (1, 0, 1)
        """
        other = self.parent()(other)
        if other.is_zero():
            raise ZeroDivisionError, "other must be nonzero"

        # This is algorithm 3.1.1 in Cohen GTM 138
        A = self
        B = other
        R = A
        Q = self.polynomial(0)
        X = self.parent().gen()
        while R.degree() >= B.degree():
            S =  (R.leading_coefficient()/B.leading_coefficient()) * X**(R.degree()-B.degree())
            Q += S
            R -= S*B
        return (Q, R)

    def _gcd(self, other):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        g = EuclideanDomainElement._gcd(self, other)
        c = g.leading_coefficient()
        if c.is_unit():
            return (1/c)*g
        return g


class Polynomial_generic_sparse_field(Polynomial_generic_sparse, Polynomial_generic_field):
    """
    EXAMPLES:
        sage: R.<x> = PolynomialRing(RR, sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial_element.Polynomial_generic_sparse_field'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_sparse.__init__(self, parent, x, check, is_gen)


class Polynomial_generic_dense_field(Polynomial_generic_dense, Polynomial_generic_field):
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_dense.__init__(self, parent, x, check, is_gen)


class Polynomial_rational_dense(Polynomial_generic_field):
    """
    A dense polynomial over the rational numbers.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

        if construct:
            self.__poly = x
            return

        self.__poly = pari([]).Polrev()

        if x is None:
            return         # leave initialized to 0 polynomial.


        if fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError, "denominator must be 1"
            else:
                x = x.numerator()

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = x.__poly.copy()
                return
            else:
                x = [QQ(a) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = QQ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, pari_gen):
            f = x.Polrev()
            self.__poly = f
            assert self.__poly.type() == "t_POL"
            return

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            x = [QQ(z) for z in x]

        self.__list = list(x)
        while len(self.__list) > 0 and self.__list[-1] == 0:
            del self.__list[-1]

        # NOTE: It is much faster to convert to string and let pari's parser at it,
        # which is why we pass str(x) in.
        self.__poly = pari(str(x)).Polrev()
        assert self.__poly.type() == "t_POL"

    def _repr(self, name=None):
        if name is None:
            name = self.parent().variable_name()
        return str(self.__poly).replace("x", name)

    def _repr_(self):
        return self._repr()

    def __reduce__(self):
        return Polynomial_rational_dense, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        return QQ(self.__poly[n])

    def __getslice__(self, i, j):
        return [QQ(x) for x in self.__poly[i:j]]

    def _pow(self, n):
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return Polynomial_rational_dense(self.parent(), self.__poly**n, construct=True)

    def _add_(self, right):
        return Polynomial_rational_dense(self.parent(),
                                         self.__poly + right.__poly, construct=True)

    def is_irreducible(self):
        try:
            return self.__poly.polisirreducible()
        except NotImplementedError:
            F = self.__poly.factor()
            if len(F) > 1 or F[0][1] > 1:
                return False
            return True

    def galois_group(self, pari_group=False, use_kash=False):
        r"""
        Return the Galois group of f as a permutation group.

        INPUT:
            self -- an irreducible polynomial

            pari_group -- bool (default: False); if True instead return
                          the Galois group as a PARI group.  This has
                          a useful label in it, and may be slightly faster
                          since it doesn't require looking up a group in
                          Gap.  To get a permutation group from a PARI
                          group P, type PermutationGroup(P).

            use_kash --   bool (default: False); if True use KASH's Galois
                          command instead of using the PARI C library.
                          An attempt is always made to use KASH if the
                          degree of the polynomial is >= 12.

        ALGORITHM: The Galois group is computed using PARI in C
        library mode, or possibly kash if available.

        \note{ The PARI documentation contains the following warning:
        The method used is that of resolvent polynomials and is
        sensitive to the current precision. The precision is updated
        internally but, in very rare cases, a wrong result may be
        returned if the initial precision was not sufficient.}

        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(); G            # uses optional database_gap package
            Transitive group number 5 of degree 4
            sage: G.gens()                           # uses optional database_gap package
            ((1,2,3,4), (1,2))
            sage: G.order()                          # uses optional database_gap package
            24

        It is potentially useful to instead obtain the corresponding
        PARI group, which is little more than a $4$-tuple.  See the
        PARI manual for the exact details.  (Note that the third
        entry in the tuple is in the new standard ordering.)
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(pari_group=True); G
            PARI group [24, -1, 5, "S4"] of degree 4
            sage: PermutationGroup(G)                # uses optional database_gap package
            Transitive group number 5 of degree 4

        You can use KASH to compute Galois groups as well.  The
        avantage is that KASH can compute Galois groups of fields up
        to degree 23, whereas PARI only goes to degree 11.  (In my
        not-so-thorough experiments PARI is faster than KASH.)

            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(use_kash=true)      # requires optional KASH
            Transitive group number 5 of degree 4

        """
        from sage.groups.all import PariGroup, PermutationGroup, TransitiveGroup
        if not self.is_irreducible():
            raise ValueError, "polynomial must be irreducible"
        if self.degree() > 11 or use_kash:
            # TODO -- maybe use KASH if available or print message that user should install KASH?
            try:
                from sage.interfaces.all import kash
                kash.eval('X := PolynomialRing(RationalField()).1')
                s = self._repr(name='X')
                G = kash('Galois(%s)'%s)
                d = int(kash.eval('%s.ext1'%G.name()))
                n = int(kash.eval('%s.ext2'%G.name()))
                return TransitiveGroup(d, n)
            except RuntimeError:
                raise NotImplementedError, "Sorry, computation of Galois groups of fields of degree bigger than 11 is not yet implemented.  Try installing the optional free (closed source) KASH package, which supports up to degree $23$."
        G = self.__poly.polgalois()
        H = PariGroup(G, self.degree())
        if pari_group:
            return H
        else:
            return PermutationGroup(H)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_rational_dense):
            right = self.parent()(right)
        if right.parent() != self.parent():
            raise TypeError
        v = self.__poly.divrem(right.__poly)
        return Polynomial_rational_dense(self.parent(), v[0], construct=True), \
               Polynomial_rational_dense(self.parent(), v[1], construct=True)


    def _mul_(self, right):
        """
        EXAMPLES:
            sage: R.<x> = PolynomialRing(QQ)
            sage: (x - QQ('2/3'))*(x^2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def _unsafe_mutate(self, n, value):
        try:
            del self.__list
        except AttributeError:
            pass
        if (<Polynomial>self)._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        if n < 0:
            raise IndexError, "n must be >= 0"
        if n <= self.__poly.poldegree():
            self.__poly[n] = QQ(value)
        else:
            self.__poly = self.__poly + pari('(%s)*x^%s'%(QQ(value),n))
        if hasattr(self, "__list"):
            del self.__list

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
            sage: R.<x> = PolynomialRing(QQ)
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.complex_roots()[0]
            0.713639173536900
        """
        R = self.__poly.polroots(flag)
        C = complex_field.CC
        return [C(a) for a in R]

    def copy(self):
        f = Polynomial_rational_dense(self.parent())
        f.__poly = self.__poly.copy()
        return f

    def degree(self):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.poldegree(), -1)

    def discriminant(self):
        """
        EXAMPLES:
            sage: _.<x> = PolynomialRing(QQ)
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            -7911
        """
        return QQ(self.__poly.poldisc())

    def disc(self):
        """
        Same as discriminant().
        """
        return self.discriminant()

    def factor_mod(self, p):
        """
        Return the factorization of self modulo the prime p.

        INPUT:
            p -- prime

        OUTPUT:
            factorization of self reduced modulo p.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        G = self._pari_().factormod(p)
        K = finite_field.FiniteField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, names=self.parent().variable_name())
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
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        prec = integer.Integer(prec)
        if prec <= 0:
            raise ValueError, "prec must be positive"
        G = self._pari_().factorpadic(p, prec)
        K = padic_field.pAdicField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, names=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, K(self.leading_coefficient()))

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES:
            sage: _.<x> = PolynomialRing(QQ)
            sage: f = x^3 + 3*x - 17/13; f
            x^3 + 3*x - 17/13
            sage: v = f.list(); v
            [-17/13, 3, 0, 1]
            sage: v[0] = 0
            sage: f
            x^3 + 3*x - 17/13
            sage: f.list()
            [-17/13, 3, 0, 1]
        """
        return [QQ(x) for x in self.__poly.Vecrev()]

##     def partial_fraction(self, g):
##         """
##         Return partial fraction decomposition of self/g, where g
##         has the same parent as self.
##         """
##         g = self.parent()(g)
##         from sage.interfaces.maxima import maxima
##         h = maxima(self)/maxima(g)
##         k = h.partfrac(self.parent().variable())

    def rescale(self, a):
        """
        Return f(a*X).
        """
        b = 1
        c = []
        for i in range(self.degree()+1):
            c.append(b*self[i])
            b *= a
        return self.parent()(c)

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial
        OUTPUT:
            an element of the base ring of the polynomial ring

        NOTES:
            Implemented using pari's polresultant function.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: f.resultant(g)
            -8
        """
        if not isinstance(other, Polynomial):
            other = self.polynomial(other)
        if self.parent() != other.parent():
            raise TypeError
        return QQ(str(self.__poly.polresultant(other.__poly, 0)))

    def hensel_lift(self, p, e):
        """
        Assuming that self factors modulo $p$ into distinct factors,
        computes the Hensel lifts of these factors modulo $p^e$.  We
        assume that $p$ has integer coefficients.
        """
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        e = integer.Integer(e)
        if e < 1:
            raise ValueError, "e must be at least 1"
        F = self.factor_mod(p)
        y = []
        for g, n in F:
            if n > 1:
                raise ArithmeticError, "The polynomial must be square free modulo p."
            y.append(g)
        H = self._pari_().polhensellift(y, p, e)
        R = integer_mod_ring.IntegerModRing(p**e)
        S = sage.rings.polynomial_ring.PolynomialRing(R, self.parent().variable_name())
        return [S(eval(str(m.Vec().Polrev().Vec()))) for m in H]

class Polynomial_integer_dense(Polynomial_generic_domain,
                               IntegralDomainElement):
    """
    A dense polynomial over the integers.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if construct:
            if isinstance(x, ZZX_class):
                self.__poly = x
                return
            self.__poly = ZZX(x)
            return

        self.__poly = ZZX([])

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = x.__poly.copy()
                return
            else:
                x = [ZZ(a) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = ZZ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, ZZX_class):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.Vecrev()]
            check = False

        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_integer_dense):
            if x.denominator() == 1:
                x = x.numerator().__poly
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
        return Polynomial_integer_dense, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        return ZZ(self.__poly[n])

    def __getslice__(self, i, j):
        i = max(0,i)
        j = min(j, self.__poly.degree()+1)
        return [ZZ(self.__poly[k]) for k in range(i,j)]

    def _pow(self, n):
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    def _add_(self, right):
        return self.parent()(self.__poly + right.__poly, construct=True)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        v = self.__poly.quo_rem(right.__poly)
        return self.parent()(v[0], construct=True), \
               self.parent()(v[1], construct=True)

    def gcd(self, right):
        """
        Return the GCD of self and other.  The leading
        coefficient need not be 1.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.gcd(right.__poly)
        return self.parent()(g, construct=True)

    def lcm(self, right):
        """
        Return the LCM of self and other, as a monic polynomial.
        """
        if not isinstance(right, Polynomial_integer_dense):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        g = self.__poly.lcm(right.__poly)
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
            QQ = sage.rings.rational_field.QQ
            f = self.base_extend(QQ)
            g, u, v = f.xgcd(right.base_extend(QQ))
            d = arith.lcm([g.denominator(), u.denominator(), v.denominator()])
            R = self.parent()
            return R(d*g), R(d*u), R(d*v)
        else:
            S = self.parent()
            return S(rr), S(s, construct=True), \
                   S(t, construct=True)


    def _mul_(self, right):
        """
        EXAMPLES:
            sage: _.<x> = PolynomialRing(ZZ)
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 - 10*x^2 + 32*x - 32
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and right[0] in ZZ:
            d = ZZ(right[0])
        elif (right in self.parent().base_ring()):
            d = ZZ(right)
        else:
            return Polynomial.__floordiv__(self, right)
        return self.parent()([c // d for c in self.list()], construct=True)

    def _unsafe_mutate(self, n, value):
        if (<Polynomial>self)._is_gen:
            raise ValueError, "the generator cannot be changed"
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
            0.713639173536900
            sage: f(alpha)
            -0.000000000000000222044604925031
        """
        QQ = sage.rings.rational_field.RationalField()
        R = sage.rings.polynomial_ring.PolynomialRing(QQ, 'x')
        return R(self.list()).complex_roots()

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
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        f = self._pari_()
        if f * pari('Mod(1,%s)'%p) == pari(0):
            raise ValueError, "factorization of 0 not defined"
        G = f.factormod(p)
        k = finite_field.FiniteField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(k, names=self.parent().variable_name())
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
        p = integer.Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        prec = integer.Integer(prec)
        if prec <= 0:
            raise ValueError, "prec must be positive"
        G = self._pari_().factorpadic(p, prec)
        K = padic_field.pAdicField(p)
        R = sage.rings.polynomial_ring.PolynomialRing(K, names=self.parent().variable_name())
        return R(1)._factor_pari_helper(G, K(self.leading_coefficient()))

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES:
            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [-17, 3, 0, 1]
        """
        return [ZZ(str(self.__poly[i])) for i in xrange(self.degree()+1)]

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
            sage: f.resultant(g)
            -8
        """
        if not isinstance(other, Polynomial) or self.parent() != other.parent():
            other = self.polynomial(other)
        return ZZ(str(self.__poly.resultant(other.__poly, 0)))

    def ntl_set_directly(self, v):
        """
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers are stored internally using NTL's ZZX
        class.  Use this function to set the value of this polynomial using
        the NTL constructor, which is potentially quicker.   The input v
        is either a vector of ints or a string of the form '[ n1 n2 n3 ... ]'
        where the ni are integers and there are no commas between them.
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
        if (<Polynomial>self)._is_gen:
            raise TypeError, "Cannot change the value of the generator."
        self.__poly = ZZX(v)
        try:
            del self.__list
        except AttributeError:
            pass



class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(Integers(16))
        sage: f = x^3 - x + 17
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True,
                 is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

        if construct:
            if isinstance(x, ZZ_pX_class):
                self.__poly = x
                return
            parent._ntl_set_modulus()
            self.__poly = ZZ_pX(x)
            return

        self.__poly = ZZ_pX([])

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                parent._ntl_set_modulus()
                self.__poly = x.__poly.__copy__()
                return
            else:
                R = parent.base_ring()
                x = [ZZ(R(a)) for a in x.list()]
                check = False

        if isinstance(x, dict):
            zero = ZZ(0)
            n = max(x.keys())
            v = [zero for _ in xrange(n+1)]
            for i, z in x.iteritems():
                v[i] = z
            x = v

        elif isinstance(x, ZZX_class):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.Vecrev()]
            check = False

        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator() == 1:
                x = x.numerator().__poly
                check = False

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        if check:
            R = parent.base_ring()
            x = [ZZ(R(a)) for a in x]

        parent._ntl_set_modulus()
        self.__poly = ZZ_pX(x)

    def __reduce__(self):
        return Polynomial_dense_mod_n, \
               (self.parent(), self.list(), False, self.is_gen())

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
        """
        Return underlying NTL representation of this polynomial.
        Additional ``bonus'' functionality is available through this
        function.
        """
        return self.__poly

    def __getitem__(self, n):
        return self.base_ring()(self.__poly[n])

    def __getslice__(self, i, j):
        R = self.base_ring()
        if i < 0:
            i = 0
        if j > self.__poly.degree()+1:
            j = self.__poly.degree()+1
        return [R(self.__poly[k]) for k in range(i,j)]

    def _unsafe_mutate(self, n, value):
        if self._is_gen:
            raise ValueError, "the generator cannot be changed"
        n = int(n)
        if n < 0:
            raise IndexError, "n must be >= 0"
        self.parent()._ntl_set_modulus()
        self.__poly[n] = int(value)

    def _pow(self, n):
        n = int(n)
        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    def _add_(self, right):
        return self.parent()(self.__poly + right.__poly, construct=True)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: x = PolynomialRing(Integers(100), 'x').0
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
        return self.parent()(self.__poly * right.__poly, construct=True)

    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.
        """
        if not isinstance(right, Polynomial_dense_mod_n):
            right = self.parent()(right)
        elif self.parent() != right.parent():
            raise TypeError
        v = self.__poly.quo_rem(right.__poly)
        P = self.parent()
        return P(v[0], construct=True), P(v[1], construct=True)

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
        return self.parent()(self.__poly.left_shift(n), construct=True)

    def _sub_(self, right):
        return self.parent()(self.__poly - right.__poly, construct=True)

    def __floordiv__(self, right):
        if is_Polynomial(right) and right.is_constant() and right[0] in self.parent().base_ring():
            d = right[0]
        elif (right in self.parent().base_ring()):
            d = right
        else:
            return Polynomial.__floordiv__(self, right)
        return self.parent()([c // d for c in self.list()], construct=True)

##     def __setitem__(self, n, value):
##         if self._is_gen:
##             raise ValueError, "the generator cannot be changed"
##         n = int(n)
##         if n < 0:
##             raise IndexError, "n must be >= 0"
##         self.parent()._ntl_set_modulus()
##         self.__poly[n] = int(value)

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
        if (<Polynomial>self)._is_gen:
            raise TypeError, "Cannot change the value of the generator."
        self.parent()._ntl_set_modulus()
        self.__poly = ZZ_pX(v)
        try:
            del self.__list
        except AttributeError:
            pass

class Polynomial_dense_mod_p(Polynomial_dense_mod_n,
                             Polynomial_singular_repr,
                             PrincipalIdealDomainElement):
    """
    A dense polynomial over the integers modulo p, where p is prime.
    """
    def __reduce__(self):
        return Polynomial_dense_mod_p, \
               (self.parent(), self.list(), False, self.is_gen())

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
            sage: _.<x> = PolynomialRing(GF(19))
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: f.resultant(g)
            11
        """
        if not isinstance(other, Polynomial) or self.parent() != other.parent():
            other = self.polynomial(other)
        self.parent()._ntl_set_modulus()
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


# ----------------- inner functions -------------
# Pyrex can't handle function definitions inside other function

def _mul_fateman_to_int2(f_list,g_list):
    """
    Convert an polynomial to an integer by evaluating it
    INPUT: p, a list of integers
    OUTPUT: padding
    """
    max_coeff_f = max([abs(i) for i in f_list])
    max_coeff_g = max([abs(i) for i in g_list])
    b = (1+min(len(f_list),len(g_list)))*max_coeff_f*max_coeff_g
    return int(pyceil(pylog(b,2)))

def _mul_fateman_to_poly(number,padding):
    """
    Converts a number to a polynomial, according
    to a padding
    OUTPUT: a list containing the coefficient of
    a polynomial of degree len(list)

    """
    coeffs = []
    flag=0
    append = coeffs.append
    if number < 0:
        number = -number
        flag=1

    while number > 0:
        r =  number%(1<<padding)
        number = (number-r) >> padding
        if r > (1<<(padding-1)):
            r -= 1<<padding
            number+=1
        append(r)

    if flag==1:
        return [-c for c in coeffs]
    return coeffs

def _mul_fateman_mul(f,g):
    """
    Multiply 2 polynomials
    """

    f=f.base_extend(QQ)
    g=g.base_extend(QQ)

    f_list = f.list()
    g_list = g.list()

    # If these polynomials have real
    # coefficients, convert them to
    # rational coeficients.
    # Note: no precision is lost in this
    # direction

    fgcd = gcd(f_list)
    ggcd = gcd(g_list)

    # Need to change ring to ZZ
    z_poly_f=(f*fgcd.denominator()).base_extend(ZZ)
    z_poly_g=(g*ggcd.denominator()).base_extend(ZZ)

    div = 1/(fgcd.denominator()*ggcd.denominator())

    z_poly_f_list = z_poly_f.coeffs()
    z_poly_g_list = z_poly_g.coeffs()
    padding = _mul_fateman_to_int2(z_poly_f_list,z_poly_g_list)

    n_f = z_poly_f(1<<padding)
    n_g = z_poly_g(1<<padding)

    if div == 1: return _mul_fateman_to_poly(n_f*n_g,padding)
    #return to_poly(n_f*n_g,padding)
    else:
        l=_mul_fateman_to_poly(n_f*n_g,padding)
        return [QQ(i*div) for i in l]



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
