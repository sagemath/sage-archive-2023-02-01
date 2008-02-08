"""
Polynomial Interfaces to Singular

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-04-21)
     -- Robert Bradshaw: Re-factor to avoid multiple inheritance vs. Cython (2007-09)

TESTS:
    sage: R = MPolynomialRing(GF(2**8,'a'),10,'x', order='invlex')
    sage: R == loads(dumps(R))
    True
    sage: P.<a,b> = PolynomialRing(GF(7), 2)
    sage: f = (a^3 + 2*b^2*a)^7; f
    a^21 + 2*a^7*b^14

"""

#################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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
#
######################################################################

import sage.rings.finite_field
import sage.rings.number_field as number_field

from sage.interfaces.all import singular as singular_default, is_SingularElement
from sage.rings.complex_field import is_ComplexField
from sage.rings.real_mpfr import is_RealField
from sage.rings.complex_double import is_ComplexDoubleField
from sage.rings.real_double import is_RealDoubleField
from sage.rings.rational_field import is_RationalField
from sage.rings.integer_ring import ZZ
import sage.rings.arith
import sage.rings.ring


class PolynomialRing_singular_repr:
    """
    Implements methods to convert polynomial rings to Singular.

    This class is a base class for all univariate and multivariate
    polynomial rings which support conversion from and to Singular
    rings.
    """
    def _singular_(self, singular=singular_default, force=False):
        """
        Returns a singular ring for this polynomial ring over a field.
        Currently QQ, GF(p), and GF(p^n), CC, and RR are supported.

        INPUT:
            singular -- Singular instance
            force -- polynomials over ZZ may be coerced to Singular by
                     treating them as polynomials over RR. This is
                     inexact but works for some cases where the
                     coeffients are not considered (default: False).

        OUTPUT:
            singular ring matching this ring

        EXAMPLES:
            sage: r = MPolynomialRing(GF(2**8,'a'),10,'x', order='invlex')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 10
            //        block   1 : ordering rp
            //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
            //        block   2 : ordering C
            sage: r = MPolynomialRing(GF(127),2,'x', order='invlex')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r = MPolynomialRing(QQ,2,'x', order='invlex')
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r = PolynomialRing(QQ,'x')
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r = PolynomialRing(GF(127),'x')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r = PolynomialRing(GF(2**8,'a'),'y')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    y
            //        block   2 : ordering C
            sage: R.<x,y> = PolynomialRing(CC,'x',2)
            sage: R._singular_()
            //   characteristic : 0 (complex:15 digits, additional 0 digits)
            //   1 parameter    : I
            //   minpoly        : (I^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            sage: R.<x,y> = PolynomialRing(RealField(100),'x',2)
            sage: R._singular_()
            //   characteristic : 0 (real:29 digits, additional 0 digits)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: w = var('w')
            sage: R.<x,y> = PolynomialRing(NumberField(w^2+1,'s'))
            sage: R._singular_()
            //   characteristic : 0
            //   1 parameter    : s
            //   minpoly        : (s^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

        WARNING:
           If the base ring is a finite extension field or a number field
           the ring will not only be returned but also be set as the current
           ring in Singular.

        NOTE:
            Singular represents precision of floating point numbers base 10
            while SAGE represents floating point precision base 2.
        """
        try:
            R = self.__singular
            if not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring() is ZZ or self.base_ring().is_prime_field():
                return R
            if sage.rings.ring.is_FiniteField(self.base_ring()) or\
                    number_field.all.is_NumberField(self.base_ring()):
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                    singular.eval("minpoly=%s"%(self.__minpoly))
                    self.__minpoly = singular.eval('minpoly')[1:-1]

            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular, force)

    def _singular_init_(self, singular=singular_default, force=False):
        """
        Return a newly created Singular ring matching this ring.
        """
        if not self._can_convert_to_singular() and not force:
            raise TypeError, "no conversion of this ring to a Singular ring defined"

        if self.ngens()==1:
            _vars = str(self.gen())
            if "*" in _vars: # 1.000...000*x
                _vars = _vars.split("*")[1]
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        if is_RealField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(real,%d,0)"%digits, _vars, order=order, check=False)

        elif is_ComplexField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(complex,%d,0,I)"%digits, _vars,  order=order, check=False)

        elif is_RealDoubleField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(real,15,0)", _vars, order=order, check=False)

        elif is_ComplexDoubleField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(complex,15,0,I)", _vars,  order=order, check=False)

        elif self.base_ring().is_prime_field() or (self.base_ring() is ZZ and force):
            self.__singular = singular.ring(self.characteristic(), _vars, order=order, check=False)

        elif sage.rings.ring.is_FiniteField(self.base_ring()):
            # not the prime field!
            gen = str(self.base_ring().gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)
            self.__minpoly = (str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]

            self.__singular = r

        elif number_field.all.is_NumberField(self.base_ring()):
            # not the rationals!
            gen = str(self.base_ring().gen())
            poly=self.base_ring().polynomial()
            poly_gen=str(poly.parent().gen())
            poly_str=str(poly).replace(poly_gen,gen)
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)
            self.__minpoly = (poly_str).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]

            self.__singular = r

        else:
            raise TypeError, "no conversion to a Singular ring defined"

        return self.__singular

    def _can_convert_to_singular(self):
        """
        Returns True if this ring's base field or ring can be
        represented in Singular.  If this is True then this polynomial
        ring can be represented in Singular.

        The following base rings are supported: $GF(p)$, $GF(p^n)$,
        rationals, number fields, and real and complex fields.
        """
        base_ring = self.base_ring()
        return ( sage.rings.ring.is_FiniteField(base_ring)
                 or is_RationalField(base_ring)
                 or (base_ring.is_prime_field() and base_ring.characteristic() <= 2147483647)
                 or is_RealField(base_ring)
                 or is_ComplexField(base_ring)
                 or is_RealDoubleField(base_ring)
                 or is_ComplexDoubleField(base_ring)
                 or number_field.all.is_NumberField(base_ring)
                 or base_ring is ZZ )


class Polynomial_singular_repr:
    """
    Implements coercion of polynomials to Singular polynomials.

    This class is a base class for all (univariate and multivariate)
    polynomial classes which support conversion from and to
    Singular polynomials.

    Due to the incompatablity of Python extension classes and multiple inheritance,
    this just defers to module-level functions.
    """
    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        return _singular_func(self, singular, have_ring, force)
    def _singular_init_func(self, singular=singular_default, have_ring=False, force=False):
        return _singular_init_func(self, singular, have_ring, force)
    def lcm(self, singular=singular_default, have_ring=False):
        return lcm_func(self, singular, have_ring)
    def diff(self, variable, have_ring=False):
        return diff_func(self, variable, have_ring)
    def resultant(self, other, variable=None):
        return resultant_func(self, other, variable)


def _singular_func(self, singular=singular_default, have_ring=False, force=False):
    """
    Return Singular polynomial matching this polynomial.

    INPUT:
        singular -- Singular instance to use

        have_ring -- if True we will not attempt to set this
                     element's ring as the current Singular
                     ring. This is useful to speed up a batch of
                     f._singular_() calls. However, it's dangerous
                     as it might lead to wrong results if another
                     ring is singluar.current_ring().  (default:
                     False)

        force -- polynomials over ZZ may be coerced to Singular by
                 treating them as polynomials over QQ. This is
                 inexact but works for some cases where the
                 coeffients are not considered (default: False).


    EXAMPLES:
        sage: P.<a,b> = PolynomialRing(GF(7), 2)
        sage: f = (a^3 + 2*b^2*a)^7; f
        a^21 + 2*a^7*b^14
        sage: h = f._singular_(); h
        a^21+2*a^7*b^14
        sage: P(h)
        a^21 + 2*a^7*b^14
        sage: P(h^20) == f^20
        True

        sage: R.<x> = PolynomialRing(GF(7))
        sage: f = (x^3 + 2*x^2*x)^7
        sage: f
        3*x^21
        sage: h = f._singular_(); h
        3*x^21
        sage: R(h)
        3*x^21
        sage: R(h^20) == f^20
        True
    """
    if not have_ring:
        self.parent()._singular_(singular,force=force).set_ring() #this is expensive

    try:
        self.__singular._check_valid()
        if self.__singular.parent() is singular:
            return self.__singular
    except (AttributeError, ValueError):
        pass
#    return self._singular_init_(singular,have_ring=have_ring)
    return _singular_init_func(self, singular,have_ring=have_ring)

def _singular_init_func(self, singular=singular_default, have_ring=False, force=False):
    """
    Return corresponding Singular polynomial but enforce that a new
    instance is created in the Singular interpreter.

    Use self._singular_() instead.
    """
    if not have_ring:
        self.parent()._singular_(singular,force=force).set_ring() #this is expensive

    self.__singular = singular(str(self))

    return self.__singular

def lcm_func(self, right, have_ring=False):
    """
    Returns the least common multiple of this element and the right element.

    INPUT:
        right -- multivariate polynomial
        have_ring -- see self._singular_() (default:False)

    OUTPUT:
        multivariate polynomial representing the least common
        multiple of self and right

    ALGORITHM: Singular

    EXAMPLES:
        sage: r.<x,y> = MPolynomialRing(GF(2**8,'a'),2)
        sage: a = r.base_ring().0
        sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
        sage: f.lcm(x^4)
        (a^2 + a)*x^6*y + (a^4 + a^3 + a)*x^4*y + (a^5)*x^4

        sage: w = var('w')
        sage: r.<x,y> = MPolynomialRing(NumberField(w^4+1,'a'),2)
        sage: a = r.base_ring().0
        sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
        sage: f.lcm(x^4)
        (a^2 + a)*x^6*y + (a^3 + a - 1)*x^4*y + (-a)*x^4
    """
    lcm = self._singular_(have_ring=have_ring).lcm(right._singular_(have_ring=have_ring))
    return lcm.sage_poly(self.parent())

def diff_func(self, variable, have_ring=False):
    """
    Differentiates self with respect to the provided variable. This
    is completely symbolic so it is also defined over e.g. finite
    fields.

    INPUT:
        variable -- the derivative is taken with respect to variable
        have_ring -- see self._singular_() (default:False)

    EXAMPLES:
        sage: R.<x,y> = PolynomialRing(RR,2)
        sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
        sage: f.diff(x)
        9.00000000000000*x^2*y^2 + 3.00000000000000
        sage: f.diff(y)
        6.00000000000000*x^3*y + 10.0000000000000*y

        The derivate is also defined over finite fields:

        sage: R.<x,y> = PolynomialRing(GF(2**8, 'a'),2)
        sage: f = x^3*y^2 + y^2 + x + 2
        sage: f.diff(x)
        x^2*y^2 + 1

        The new coefficients are coerced to the base ring:

        sage: f.diff(y)
        0

        sage: w = var('w')
        sage: R.<x,y> = PolynomialRing(NumberField(w^3-2, 'a'),2)
        sage: a=R.base_ring().0
        sage: f = x^3*y^2 + y^2 + a*x + 2
        sage: f.diff(x)
        3*x^2*y^2 + a

    ALGORITHM: Singular

    """
    df = self._singular_(have_ring=have_ring).diff(variable._singular_(have_ring=have_ring))
    return df.sage_poly(self.parent())


def resultant_func(self, other, variable=None):
    """
    computes the resultant of self and the first argument with
    respect to the variable given as the second argument.

    If a second argument is not provide the first variable of
    self.parent() is chosen.

    INPUT:
        other -- polynomial in self.parent()
        variable -- optional variable (of type polynomial) in self.parent() (default: None)

    EXAMPLE:
        sage: P.<x,y> = PolynomialRing(QQ,2)
        sage: a = x+y
        sage: b = x^3-y^3
        sage: c = a.resultant(b); c
        -2*y^3
        sage: d = a.resultant(b,y); d
        2*x^3

    TESTS:
        sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
        sage: P.<x,y> = MPolynomialRing_polydict_domain(QQ,2,order='degrevlex')
        sage: a = x+y
        sage: b = x^3-y^3
        sage: c = a.resultant(b); c
        -2*y^3
        sage: d = a.resultant(b,y); d
        2*x^3

    """
    if variable is None:
        variable = self.parent().gen(0)
    rt = self._singular_().resultant(other._singular_(), variable._singular_())
    r = rt.sage_poly(self.parent())
    if self.parent().ngens() <= 1 and r.degree() <= 0:
        return self.parent().base_ring()(r[0])
    else:
        return r
