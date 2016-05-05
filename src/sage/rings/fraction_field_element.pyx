"""
Fraction Field Elements

AUTHORS:

- William Stein (input from David Joyner, David Kohel, and Joe Wetherell)

- Sebastian Pancratz (2010-01-06): Rewrite of addition, multiplication and
  derivative to use Henrici's algorithms [Ho72]

REFERENCES:

.. [Ho72] \E. Horowitz, "Algorithms for Rational Function Arithmetic
   Operations", Annual ACM Symposium on Theory of Computing, Proceedings of
   the Fourth Annual ACM Symposium on Theory of Computing, pp. 108--118, 1972

"""

#*****************************************************************************
#
#   Sage: System for Algebra and Geometry Experimentation
#
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

from sage.structure.element cimport FieldElement, ModuleElement, RingElement, \
        Element
from sage.structure.element import parent

import integer_ring
from integer_ring import ZZ
from rational_field import QQ

import sage.misc.latex as latex

def is_FractionFieldElement(x):
    """
    Returns whether or not ``x`` is a :class`FractionFieldElement`.

    EXAMPLES::

        sage: from sage.rings.fraction_field_element import is_FractionFieldElement
        sage: R.<x> = ZZ[]
        sage: is_FractionFieldElement(x/2)
        False
        sage: is_FractionFieldElement(2/x)
        True
        sage: is_FractionFieldElement(1/3)
        False
    """
    return isinstance(x, FractionFieldElement)

cdef class FractionFieldElement(FieldElement):
    """
    EXAMPLES::

        sage: K = FractionField(PolynomialRing(QQ, 'x'))
        sage: K
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: loads(K.dumps()) == K
        True
        sage: x = K.gen()
        sage: f = (x^3 + x)/(17 - x^19); f
        (x^3 + x)/(-x^19 + 17)
        sage: loads(f.dumps()) == f
        True

    TESTS:

    Test if :trac:`5451` is fixed::

        sage: A = FiniteField(9,'theta')['t']
        sage: K.<t> = FractionField(A)
        sage: f= 2/(t^2+2*t); g =t^9/(t^18 + t^10 + t^2);f+g
        (2*t^15 + 2*t^14 + 2*t^13 + 2*t^12 + 2*t^11 + 2*t^10 + 2*t^9 + t^7 + t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)/(t^17 + t^9 + t)

    Test if :trac:`8671` is fixed::

        sage: P.<n> = QQ[]
        sage: F = P.fraction_field()
        sage: P.one()//F.one()
        1
        sage: F.one().quo_rem(F.one())
        (1, 0)
    """
    cdef object __numerator
    cdef object __denominator

    def __init__(self, parent, numerator, denominator=1,
                 coerce=True, reduce=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: K.<x> = Frac(ZZ['x'])
            sage: FractionFieldElement(K, x, 4)
            x/4
            sage: FractionFieldElement(K, x, x, reduce=False)
            x/x
            sage: f = FractionFieldElement(K, 'hi', 1, coerce=False, reduce=False)
            sage: f.numerator()
            'hi'

            sage: x = var('x')
            sage: K((x + 1)/(x^2 + x + 1))
            (x + 1)/(x^2 + x + 1)
            sage: K(355/113)
            355/113

        """
        FieldElement.__init__(self, parent)
        if coerce:
            self.__numerator   = parent.ring()(numerator)
            self.__denominator = parent.ring()(denominator)
        else:
            self.__numerator   = numerator
            self.__denominator = denominator
        if reduce and parent.is_exact():
            try:
                self.reduce()
            except ArithmeticError:
                pass
        if self.__denominator.is_zero():
            raise ZeroDivisionError, "fraction field element division by zero"

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: K = GF(7)['a,b'].fraction_field()
            sage: a,b = K.gens()

        ::

            sage: phi = F.hom([a+b, a*b], K)
            sage: phi(x+y) # indirect doctest
            a*b + a + b

        ::

            sage: (x^2/y)._im_gens_(K, [a+b, a*b])
            (a^2 + 2*a*b + b^2)/(a*b)
            sage: (x^2/y)._im_gens_(K, [a, a*b])
            a/b
        """
        nnum = codomain.coerce(self.__numerator._im_gens_(codomain, im_gens))
        nden = codomain.coerce(self.__denominator._im_gens_(codomain, im_gens))
        return codomain.coerce(nnum/nden)

    def reduce(self):
        """
        Divides out the gcd of the numerator and denominator.

        Automatically called for exact rings, but because it may be
        numerically unstable for inexact rings it must be called manually
        in that case.

        EXAMPLES::

            sage: R.<x> = RealField(10)[]
            sage: f = (x^2+2*x+1)/(x+1); f
            (x^2 + 2.0*x + 1.0)/(x + 1.0)
            sage: f.reduce(); f
            x + 1.0
        """
        try:
            g = self.__numerator.gcd(self.__denominator)
            if not g.is_unit():
                num, _ = self.__numerator.quo_rem(g)
                den, _ = self.__denominator.quo_rem(g)
            else:
                num = self.__numerator
                den = self.__denominator
            if not den.is_one() and den.is_unit():
                try:
                    num *= den.inverse_of_unit()
                    den  = den.parent().one()
                except Exception:
                    pass
            self.__numerator   = num
            self.__denominator = den
        except AttributeError:
            raise ArithmeticError, "unable to reduce because lack of gcd or quo_rem algorithm"
        except TypeError:
            raise ArithmeticError, "unable to reduce because gcd algorithm doesn't work on input"
        except NotImplementedError:
            raise ArithmeticError("unable to reduce because gcd algorithm not implemented on input")

    def __copy__(self):
        """
        Make a copy of ``self``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: copy(f)
            (x + y)/y
        """
        return self.__class__(self._parent, self.__numerator,
                self.__denominator, coerce=False, reduce=False)

    def numerator(self):
        """
        Return the numerator of ``self``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.numerator()
            x + y
        """
        return self.__numerator

    def denominator(self):
        """
        Return the denominator of ``self``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x/y+1; f
            (x + y)/y
            sage: f.denominator()
            y
        """
        return self.__denominator


    def is_square(self,root=False):
        """
        Returns whether or not ``self`` is a perfect square. If the optional
        argument ``root`` is ``True``, then also returns a square root (or
        ``None``, if the fraction field element is not square).

        INPUT:

        -  ``root`` -- whether or not to also return a square
           root (default: ``False``)

        OUTPUT:

        -  ``bool`` - whether or not a square

        -  ``object`` - (optional) an actual square root if
           found, and None otherwise.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: (1/t).is_square()
            False
            sage: (1/t^6).is_square()
            True
            sage: ((1+t)^4/t^6).is_square()
            True
            sage: (4*(1+t)^4/t^6).is_square()
            True
            sage: (2*(1+t)^4/t^6).is_square()
            False
            sage: ((1+t)/t^6).is_square()
            False

            sage: (4*(1+t)^4/t^6).is_square(root=True)
            (True, (2*t^2 + 4*t + 2)/t^3)
            sage: (2*(1+t)^4/t^6).is_square(root=True)
            (False, None)

            sage: R.<x> = QQ[]
            sage: a = 2*(x+1)^2 / (2*(x-1)^2); a
            (2*x^2 + 4*x + 2)/(2*x^2 - 4*x + 2)
            sage: a.numerator().is_square()
            False
            sage: a.is_square()
            True
            sage: (0/x).is_square()
            True
        """
        a = self.numerator()
        b = self.denominator()
        if not root:
            return  (a*b).is_square( root = False )
        is_sqr, sq_rt = (a*b).is_square( root = True )
        if is_sqr:
            return True, self._parent( sq_rt/b )
        return False, None


    def __hash__(self):
        """
        This function hashes in a special way to ensure that generators of
        a ring `R` and generators of a fraction field of `R` have the same
        hash. This enables them to be used as keys interchangeably in a
        dictionary (since ``==`` will claim them equal). This is particularly
        useful for methods like ``subs`` on ``ParentWithGens`` if you are
        passing a dictionary of substitutions.

        EXAMPLES::

            sage: R.<x>=ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: ((x+1)/(x^2+1)).subs({x:1})
            1
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x>=QQ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x,y,z>=ZZ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: R.<x,y,z>=QQ[] # this probably has a separate implementation from ZZ[]
            sage: hash(R.0)==hash(FractionField(R).0)
            True
            sage: ((x+1)/(x^2+1)).subs({x:1})
            1
            sage: d={x:1}
            sage: d[FractionField(R).0]
            1
            sage: hash(R(1)/R(2))==hash(1/2)
            True
        """
        # This is same algorithm as used for members of QQ
        #cdef long n, d
        n = hash(self.__numerator)
        d = hash(self.__denominator)
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

    def __call__(self, *x, **kwds):
        """
        Evaluate the fraction at the given arguments. This assumes that a
        call function is defined for the numerator and denominator.

        EXAMPLES::

            sage: x = PolynomialRing(RationalField(),'x',3).gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            -2*x1*x2 + x0 + x1
            sage: f(1,2,5)
            -17
            sage: h = f /(x[1] + x[2])
            sage: h
            (-2*x1*x2 + x0 + x1)/(x1 + x2)
            sage: h(1,2,5)
            -17/7
            sage: h(x0=1)
            (-2*x1*x2 + x1 + 1)/(x1 + x2)
        """
        return self.__numerator(*x, **kwds) / self.__denominator(*x, **kwds)

    def _is_atomic(self):
        """
        EXAMPLES::

            sage: K.<x> = Frac(ZZ['x'])
            sage: x._is_atomic()
            True
            sage: f = 1/(x+1)
            sage: f._is_atomic()
            False
        """
        return self.__numerator._is_atomic() and self.__denominator._is_atomic()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: K.<x> = Frac(ZZ['x'])
            sage: repr(x+1) # indirect doctest
            'x + 1'
            sage: repr((x+1)/(x-1))
            '(x + 1)/(x - 1)'
            sage: repr(1/(x-1))
            '1/(x - 1)'
            sage: repr(1/x)
            '1/x'
        """
        if self.is_zero():
            return "0"
        s = "%s"%self.__numerator
        if self.__denominator != 1:
            denom_string = str( self.__denominator )
            if self.__denominator._is_atomic() and not ('*' in denom_string or '/' in denom_string):
                s = "%s/%s"%(self.__numerator._coeff_repr(no_space=False),denom_string)
            else:
                s = "%s/(%s)"%(self.__numerator._coeff_repr(no_space=False),denom_string)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of this fraction field element.

        EXAMPLES::

            sage: R = PolynomialRing(QQ, 'x')
            sage: F = R.fraction_field()
            sage: x = F.gen()
            sage: a = x^2 / 1
            sage: latex(a) # indirect doctest
            x^{2}
            sage: latex(x^2/(x^2+1))
            \frac{x^{2}}{x^{2} + 1}
            sage: a = 1/x
            sage: latex(a)
            \frac{1}{x}

        TESTS::

            sage: R = RR['x']     # Inexact, so no reduction.
            sage: F = Frac(R)
            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: z = FractionFieldElement(F, 0, R.gen(), coerce=False)
            sage: z.numerator() == 0
            True
            sage: z.denominator() == R.gen()
            True
            sage: latex(z) # indirect doctest
            0
        """
        if self.is_zero():
            return "0"
        if self.__denominator == 1:
            return latex.latex(self.__numerator)
        return "\\frac{%s}{%s}"%(latex.latex(self.__numerator),
                                 latex.latex(self.__denominator))

    def _magma_init_(self, magma):
        """
        Return a string representation of ``self`` Magma can understand.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: magma((x^2 + x + 1)/(x + 1)) # optional - magma # indirect doctest
            (x^2 + x + 1)/(x + 1)

        ::

            sage: R.<x,y> = QQ[]
            sage: magma((x+y)/x)                        # optional - magma
            (x + y)/x
        """
        pgens = magma(self._parent).gens()

        s = self._repr_()
        for i, j in zip(self._parent.variable_names(), pgens):
            s = s.replace(i, j.name())

        return s

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Computes the sum of ``self`` and ``right``.

        INPUT:

        - ``right`` -- ``ModuleElement`` to add to ``self``

        OUTPUT:

        - Sum of ``self`` and ``right``

        EXAMPLES::

            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: x+y # indirect doctest
            x + y
            sage: 1/x + 1/y
            (x + y)/(x*y)
            sage: 1/x + 1/(x*y)
            (y + 1)/(x*y)
            sage: Frac(CDF['x']).gen() + 3
            x + 3.0
        """
        rnum = self.__numerator
        rden = self.__denominator
        snum = (<FractionFieldElement> right).__numerator
        sden = (<FractionFieldElement> right).__denominator

        if (rnum.is_zero()):
            return <FractionFieldElement> right
        if (snum.is_zero()):
            return self

        if self._parent.is_exact():
            try:
                d = rden.gcd(sden)
                if d.is_unit():
                    return self.__class__(self._parent, rnum*sden + rden*snum,
                        rden*sden, coerce=False, reduce=False)
                else:
                    rden = rden // d
                    sden = sden // d
                    tnum = rnum * sden + rden * snum
                    if tnum.is_zero():
                        return self.__class__(self._parent, tnum,
                            self._parent.ring().one(), coerce=False,
                            reduce=False)
                    else:
                        tden = self.__denominator * sden
                        e    = tnum.gcd(d)
                        if not e.is_unit():
                            tnum = tnum // e
                            tden = tden // e
                        if not tden.is_one() and tden.is_unit():
                            try:
                                tnum = tnum * tden.inverse_of_unit()
                                tden = self._parent.ring().one()
                            except AttributeError:
                                pass
                            except NotImplementedError:
                                pass
                        return self.__class__(self._parent, tnum, tden,
                            coerce=False, reduce=False)
            except AttributeError:
                pass
            except NotImplementedError:
                pass
            except TypeError:
                pass

        rnum = self.__numerator
        rden = self.__denominator
        snum = (<FractionFieldElement> right).__numerator
        sden = (<FractionFieldElement> right).__denominator

        return self.__class__(self._parent, rnum*sden + rden*snum, rden*sden,
            coerce=False, reduce=False)

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Computes the difference of ``self`` and ``right``.

        INPUT:

        - ``right`` - ``ModuleElement`` to subtract from ``self``

        OUTPUT:

        - Difference of ``self`` and ``right``

        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: t - 1/t # indirect doctest
            (t^2 + 6)/t
        """
        return self._add_(-right)

    cpdef RingElement _mul_(self, RingElement right):
        """
        Computes the product of ``self`` and ``right``.

        INPUT:

        - ``right`` - ``RingElement`` to multiply with ``self``

        OUTPUT:

        - Product of ``self`` and ``right``

        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: a = t/(1+t)
            sage: b = 3/t
            sage: a*b # indirect doctest
            3/(t + 1)
        """
        rnum = self.__numerator
        rden = self.__denominator
        snum = (<FractionFieldElement> right).__numerator
        sden = (<FractionFieldElement> right).__denominator

        if (rnum.is_zero() or snum.is_zero()):
            return self._parent.zero()

        if self._parent.is_exact():
            try:
                d1 = rnum.gcd(sden)
                d2 = snum.gcd(rden)
                if not d1.is_unit():
                    rnum = rnum // d1
                    sden = sden // d1
                if not d2.is_unit():
                    rden = rden // d2
                    snum = snum // d2
                tnum = rnum * snum
                tden = rden * sden
                if not tden.is_one() and tden.is_unit():
                    try:
                        tnum = tnum * tden.inverse_of_unit()
                        tden = self._parent.ring().one()
                    except AttributeError:
                        pass
                    except NotImplementedError:
                        pass
                return self.__class__(self._parent, tnum, tden,
                    coerce=False, reduce=False)
            except AttributeError:
                pass
            except NotImplementedError:
                pass
            except TypeError:
                pass

        rnum = self.__numerator
        rden = self.__denominator
        snum = (<FractionFieldElement> right).__numerator
        sden = (<FractionFieldElement> right).__denominator

        return self.__class__(self._parent, rnum * snum, rden * sden,
            coerce=False, reduce=False)

    cpdef RingElement _div_(self, RingElement right):
        """
        Computes the quotient of ``self`` and ``right``.

        INPUT:

        - ``right`` -- ``RingElement`` that is the divisor

        OUTPUT:

        Quotient of ``self`` and ``right``

        EXAMPLES::

            sage: K.<x,y,z> = Frac(ZZ['x,y,z'])
            sage: a = (x+1)*(x+y)/(z-3) # indirect doctest
            sage: b = (x+y)/(z-1)
            sage: a/b
            (x*z - x + z - 1)/(z - 3)
        """
        snum = (<FractionFieldElement> right).__numerator
        sden = (<FractionFieldElement> right).__denominator

        if snum.is_zero():
            raise ZeroDivisionError, "fraction field element division by zero"

        rightinv = self.__class__(self._parent, sden, snum,
            coerce=True, reduce=False)

        return self._mul_(rightinv)

    def __int__(self):
        """
        EXAMPLES::

            sage: K = Frac(ZZ['x'])
            sage: int(K(-3))
            -3
            sage: K.<x> = Frac(RR['x'])
            sage: x/x
            x/x
            sage: int(x/x)
            1
            sage: int(K(.5))
            0
        """
        if self.__denominator != 1:
            self.reduce()
        if self.__denominator == 1:
            return int(self.__numerator)
        else:
            raise TypeError, "denominator must equal 1"

    def _integer_(self, Z=ZZ):
        """
        EXAMPLES::

            sage: K = Frac(ZZ['x'])
            sage: K(5)._integer_()
            5
            sage: K.<x> = Frac(RR['x'])
            sage: ZZ(2*x/x)
            2
        """
        if self.__denominator != 1:
            self.reduce()
        if self.__denominator == 1:
            return Z(self.__numerator)
        raise TypeError, "no way to coerce to an integer."

    def _rational_(self, Q=QQ):
        """
        EXAMPLES::

            sage: K.<x> = Frac(QQ['x'])
            sage: K(1/2)._rational_()
            1/2
            sage: K(1/2 + x/x)._rational_()
            3/2
        """
        return Q(self.__numerator) / Q(self.__denominator)

    def __long__(self):
        """
        EXAMPLES::

            sage: K.<x> = Frac(QQ['x'])
            sage: long(K(3))
            3L
            sage: long(K(3/5))
            0L
        """
        return long(int(self))

    def __pow__(self, right, dummy):
        r"""
        Returns self raised to the `right^{th}` power.

        Note that we need to check whether or not right is negative so we
        don't set ``__numerator`` or ``__denominator`` to an element of the
        fraction field instead of the underlying ring.

        EXAMPLES::

            sage: R = QQ['x','y']
            sage: FR = R.fraction_field()
            sage: x,y = FR.gens()
            sage: a = x^2; a
            x^2
            sage: type(a.numerator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: type(a.denominator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: a = x^(-2); a
            1/x^2
            sage: type(a.numerator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: type(a.denominator())
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: x^0
            1
            sage: ((x+y)/(x-y))^2
            (x^2 + 2*x*y + y^2)/(x^2 - 2*x*y + y^2)
            sage: ((x+y)/(x-y))^-2
            (x^2 - 2*x*y + y^2)/(x^2 + 2*x*y + y^2)
            sage: ((x+y)/(x-y))^0
            1
        """
        snum = (<FractionFieldElement> self).__numerator
        sden = (<FractionFieldElement> self).__denominator
        if right == 0:
            R = self.parent().ring()
            return self.__class__(self.parent(),
                R.one(), R.one(),
                coerce=False, reduce=False)
        elif right > 0:
            return self.__class__(self.parent(),
                snum**right, sden**right,
                coerce=False, reduce=False)
        else:
            right = -right
            return self.__class__(self.parent(),
                sden**right, snum**right,
                coerce=False, reduce=False)

    def __neg__(self):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(5)['t'])
            sage: f = (t^2+t)/(t+2); f
            (t^2 + t)/(t + 2)
            sage: -f
            (4*t^2 + 4*t)/(t + 2)
        """
        return self.__class__(self._parent,
            -self.__numerator, self.__denominator,
            coerce=False, reduce=False)

    def __abs__(self):
        """
        EXAMPLES::

            sage: from sage.rings.fraction_field_element import FractionFieldElement
            sage: abs(FractionFieldElement(QQ, -2, 3, coerce=False))
            2/3
        """
        return abs(self.__numerator)/abs(self.__denominator)

    def __invert__(self):
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: f = (t^2+5)/(t-1)
            sage: ~f
            (t + 6)/(t^2 + 5)
        """
        if self.is_zero():
            raise ZeroDivisionError, "Cannot invert 0"
        return self.__class__(self._parent,
            self.__denominator, self.__numerator, coerce=False, reduce=False)

    def __float__(self):
        """
        EXAMPLES::

            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: float(x/x + y/y)
            2.0
        """
        return float(self.__numerator) / float(self.__denominator)

    cpdef int _cmp_(self, Element other) except -2:
        """
        EXAMPLES::

            sage: K.<t> = Frac(GF(7)['t'])
            sage: t/t == 1
            True
            sage: t+1/t == (t^2+1)/t
            True
            sage: t == t/5
            False

        ::

            sage: K.<x,y> = Frac(ZZ['x,y'])
            sage: x > y
            True
            sage: 1 > y
            False
        """
        return cmp(self.__numerator * \
                (<FractionFieldElement>other).__denominator,
                self.__denominator*(<FractionFieldElement>other).__numerator)

    def valuation(self, v=None):
        """
        Return the valuation of ``self``, assuming that the numerator and
        denominator have valuation functions defined on them.

        EXAMPLES::

            sage: x = PolynomialRing(RationalField(),'x').gen()
            sage: f = (x^3 + x)/(x^2 - 2*x^3)
            sage: f
            (x^2 + 1)/(-2*x^2 + x)
            sage: f.valuation()
            -1
            sage: f.valuation(x^2+1)
            1
        """
        return self.__numerator.valuation(v) - self.__denominator.valuation(v)

    def __nonzero__(self):
        """
        Return ``True`` if this element is nonzero.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: t = F(0)/x
            sage: t.__nonzero__()
            False

        ::

            sage: (1/x).__nonzero__()
            True
        """
        return not self.__numerator.is_zero()

    def is_zero(self):
        """
        Return ``True`` if this element is equal to zero.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: t = F(0)/x
            sage: t.is_zero()
            True
            sage: u = 1/x - 1/x
            sage: u.is_zero()
            True
            sage: u.parent() is F
            True
        """
        return self.__numerator.is_zero()

    def is_one(self):
        """
        Return ``True`` if this element is equal to one.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: (x/x).is_one()
            True
            sage: (x/y).is_one()
            False
        """
        return self.__numerator == self.__denominator

    def _symbolic_(self, ring):
        """
        Return ``self`` as a fraction in the ring ``ring``. Used for
        :func:`symbolic_expression` in creating a symbolic expression of
        ``self``.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: elt = (2*x + 2*y) / (3*x - 3*y); elt
            (2*x + 2*y)/(3*x - 3*y)
            sage: elt._symbolic_(SR)
            2/3*(x + y)/(x - y)
            sage: symbolic_expression(elt)
            2/3*(x + y)/(x - y)
        """
        return ring(self.__numerator)/ring(self.__denominator)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: F = ZZ['x,y'].fraction_field()
            sage: f = F.random_element()
            sage: loads(f.dumps()) == f
            True
        """
        return (make_element,
                (self._parent, self.__numerator, self.__denominator))

    def _evaluate_polynomial(self, pol):
        """
        Evaluate a univariate polynomial on this fraction.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: pol = x^3 + 1
            sage: pol(1/x)
            (x^3 + 1)/x^3

        TESTS::

            sage: R.<y,z> = ZZ[]
            sage: (~(y+z))._evaluate_polynomial(pol)
            (y^3 + 3*y^2*z + 3*y*z^2 + z^3 + 1)/(y^3 + 3*y^2*z + 3*y*z^2 + z^3)
            sage: rat = (y+z)/y
            sage: rat._evaluate_polynomial(pol)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: pol(rat)
            (2*y^3 + 3*y^2*z + 3*y*z^2 + z^3)/y^3
        """
        inverse = ~self
        if inverse.denominator().is_one():
            num = inverse.numerator()
            return pol.reverse()(num)/num**pol.degree()
        else:
            raise NotImplementedError

class FractionFieldElement_1poly_field(FractionFieldElement):
    """
    A fraction field element where the parent is the fraction field of a
    univariate polynomial ring.

    Many of the functions here are included for coherence with number fields.
    """
    def is_integral(self):
        """
        Returns whether this element is actually a polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: elt = (t^2 + t - 2) / (t + 2); elt # == (t + 2)*(t - 1)/(t + 2)
            t - 1
            sage: elt.is_integral()
            True
            sage: elt = (t^2 - t) / (t+2); elt # == t*(t - 1)/(t + 2)
            (t^2 - t)/(t + 2)
            sage: elt.is_integral()
            False
        """
        if self.denominator() != 1:
            self.reduce()
        return self.denominator() == 1

    def support(self):
        """
        Returns a sorted list of primes dividing either the numerator or
        denominator of this element.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: h = (t^14 + 2*t^12 - 4*t^11 - 8*t^9 + 6*t^8 + 12*t^6 - 4*t^5 - 8*t^3 + t^2 + 2)/(t^6 + 6*t^5 + 9*t^4 - 2*t^2 - 12*t - 18)
            sage: h.support()
            [t - 1, t + 3, t^2 + 2, t^2 + t + 1, t^4 - 2]
        """
        L = [fac[0] for fac in self.numerator().factor()] + [fac[0] for fac in self.denominator().factor()]
        L.sort()
        return L


def make_element(parent, numerator, denominator):
    """
    Used for unpickling :class:`FractionFieldElement` objects (and subclasses).

    EXAMPLES::

        sage: from sage.rings.fraction_field_element import make_element
        sage: R = ZZ['x,y']
        sage: x,y = R.gens()
        sage: F = R.fraction_field()
        sage: make_element(F, 1+x, 1+y)
        (x + 1)/(y + 1)
    """

    return parent._element_class(parent, numerator, denominator)

def make_element_old(parent, cdict):
    """
    Used for unpickling old :class:`FractionFieldElement` pickles.

    EXAMPLES::

        sage: from sage.rings.fraction_field_element import make_element_old
        sage: R.<x,y> = ZZ[]
        sage: F = R.fraction_field()
        sage: make_element_old(F, {'_FractionFieldElement__numerator':x+y,'_FractionFieldElement__denominator':x-y})
        (x + y)/(x - y)
    """
    return FractionFieldElement(parent,
            cdict['_FractionFieldElement__numerator'],
            cdict['_FractionFieldElement__denominator'],
            coerce=False, reduce=False)

