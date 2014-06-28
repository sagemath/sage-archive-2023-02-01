"""
Polynomial Interfaces to Singular

AUTHORS:

- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-04-21)
- Robert Bradshaw: Re-factor to avoid multiple inheritance vs. Cython (2007-09)
- Syed Ahmad Lavasani: Added function field to _singular_init_ (2011-12-16)
       Added non-prime finite fields to _singular_init_ (2012-1-22)

TESTS::

    sage: R = PolynomialRing(GF(2**8,'a'),10,'x', order='invlex')
    sage: R == loads(dumps(R))
    True
    sage: P.<a,b> = PolynomialRing(GF(7), 2)
    sage: f = (a^3 + 2*b^2*a)^7; f
    a^21 + 2*a^7*b^14

"""

#################################################################
#
#   Sage: System for Algebra and Geometry Experimentation
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

import sage.rings.fraction_field
import sage.rings.number_field as number_field

from sage.interfaces.all import singular as singular_default
from sage.rings.complex_field import is_ComplexField
from sage.rings.real_mpfr import is_RealField
from sage.rings.complex_double import is_ComplexDoubleField
from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
from sage.rings.real_double import is_RealDoubleField
from sage.rings.rational_field import is_RationalField
from sage.rings.function_field.function_field import is_RationalFunctionField
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.rings.integer_ring import ZZ

import sage.rings.arith
import sage.rings.finite_rings.constructor


class PolynomialRing_singular_repr:
    """
    Implements methods to convert polynomial rings to Singular.

    This class is a base class for all univariate and multivariate
    polynomial rings which support conversion from and to Singular
    rings.
    """
    def _singular_(self, singular=singular_default):
        r"""
        Returns a singular ring for this polynomial ring.
        Currently `\QQ`, `{\rm GF}(p), {\rm GF}(p^n)`, `\CC`, `\RR`, `\ZZ` and
        `\ZZ/n\ZZ` are supported.

        INPUT:

        - ``singular`` - Singular instance

        OUTPUT: Singular ring matching this ring

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(CC,'x',2)
            sage: singular(R)
            //   characteristic : 0 (complex:15 digits, additional 0 digits)
            //   1 parameter    : I
            //   minpoly        : (I^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            sage: R.<x,y> = PolynomialRing(RealField(100),'x',2)
            sage: singular(R)
            //   characteristic : 0 (real:29 digits, additional 0 digits)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: w = var('w')
            sage: R.<x> = PolynomialRing(NumberField(w^2+1,'s'))
            sage: singular(R)
            //   characteristic : 0
            //   1 parameter    : s
            //   minpoly        : (s^2+1)
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),1,'x')
            sage: singular(R)
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,1,'x')
            sage: singular(R)
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(QQ,'x')
            sage: singular(R)
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = PolynomialRing(GF(127),'x')
            sage: singular(R)
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

            sage: R = Frac(ZZ['a,b'])['x,y']
            sage: singular(R)
            //   characteristic : 0
            //   2 parameter    : a b
            //   minpoly        : 0
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C


            sage: R = IntegerModRing(1024)['x,y']
            sage: singular(R)
            //   coeff. ring is : Z/2^10
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = IntegerModRing(15)['x,y']
            sage: singular(R)
            //   coeff. ring is : Z/15
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: R = ZZ['x,y']
            sage: singular(R)
            //   coeff. ring is : Integers
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

            sage: k.<a> = FiniteField(25)
            sage: R = k['x']
            sage: K = R.fraction_field()
            sage: S = K['y']
            sage: singular(S)
            //   characteristic : 5
            //   1 parameter    : x
            //   minpoly        : 0
            //   number of vars : 2
            //        block   1 : ordering lp
            //                  : names    a y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=a2-a+2

        .. warning::

            - If the base ring is a finite extension field or a number field
              the ring will not only be returned but also be set as the current
              ring in Singular.
            - Singular represents precision of floating point numbers base 10
              while Sage represents floating point precision base 2.
        """
        try:
            R = self.__singular
            if not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring() is ZZ or self.base_ring().is_prime_field():
                return R
            if sage.rings.finite_rings.constructor.is_FiniteField(self.base_ring()) or\
                    (number_field.number_field_base.is_NumberField(self.base_ring()) and self.base_ring().is_absolute()):
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                    singular.eval("minpoly=%s"%(self.__minpoly))
                    self.__minpoly = singular.eval('minpoly')[1:-1]

            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
        Return a newly created Singular ring matching this ring.

        EXAMPLES::

            sage: PolynomialRing(QQ,'u_ba')._singular_init_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    u_ba
            //        block   2 : ordering C
        """
        if not can_convert_to_singular(self):
            raise TypeError("no conversion of this ring to a Singular ring defined")

        if self.ngens()==1:
            _vars = '(%s)'%self.gen()
            if "*" in _vars: # 1.000...000*x
                _vars = _vars.split("*")[1]
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        base_ring = self.base_ring()

        if is_RealField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = base_ring.precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(real,%d,0)"%digits, _vars, order=order, check=False)

        elif is_ComplexField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = base_ring.precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(complex,%d,0,I)"%digits, _vars,  order=order, check=False)

        elif is_RealDoubleField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(real,15,0)", _vars, order=order, check=False)

        elif is_ComplexDoubleField(base_ring):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(complex,15,0,I)", _vars,  order=order, check=False)

        elif base_ring.is_prime_field():
            self.__singular = singular.ring(self.characteristic(), _vars, order=order, check=False)

        elif sage.rings.finite_rings.constructor.is_FiniteField(base_ring):
            # not the prime field!
            gen = str(base_ring.gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)

            self.__minpoly = (str(base_ring.modulus()).replace("x",gen)).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]

            self.__singular = r

        elif number_field.number_field_base.is_NumberField(base_ring) and base_ring.is_absolute():
            # not the rationals!
            gen = str(base_ring.gen())
            poly=base_ring.polynomial()
            poly_gen=str(poly.parent().gen())
            poly_str=str(poly).replace(poly_gen,gen)
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)
            self.__minpoly = (poly_str).replace(" ","")
            if  singular.eval('minpoly') != "(" + self.__minpoly + ")":
                singular.eval("minpoly=%s"%(self.__minpoly) )
                self.__minpoly = singular.eval('minpoly')[1:-1]

            self.__singular = r

        elif sage.rings.fraction_field.is_FractionField(base_ring) and (base_ring.base_ring() is ZZ or base_ring.base_ring().is_prime_field() or is_FiniteField(base_ring.base_ring())):
            if base_ring.ngens()==1:
              gens = str(base_ring.gen())
            else:
              gens = str(base_ring.gens())

            if not (not base_ring.base_ring().is_prime_field() and is_FiniteField(base_ring.base_ring())) :
                self.__singular = singular.ring( "(%s,%s)"%(base_ring.characteristic(),gens), _vars, order=order, check=False)
            else:
                ext_gen = str(base_ring.base_ring().gen())
                _vars = '(' + ext_gen + ', ' + _vars[1:];

                R = self.__singular = singular.ring( "(%s,%s)"%(base_ring.characteristic(),gens), _vars, order=order, check=False)

                self.base_ring().__minpoly = (str(base_ring.base_ring().modulus()).replace("x",ext_gen)).replace(" ","")
                singular.eval('setring '+R._name);
                self.__singular = singular("std(ideal(%s))"%(self.base_ring().__minpoly),type='qring')

        elif sage.rings.function_field.function_field.is_RationalFunctionField(base_ring) and base_ring.constant_field().is_prime_field():
            gen = str(base_ring.gen())
            self.__singular = singular.ring( "(%s,%s)"%(base_ring.characteristic(),gen), _vars, order=order, check=False)

        elif is_IntegerModRing(base_ring):
            ch = base_ring.characteristic()
            if ch.is_power_of(2):
                exp = ch.nbits() -1
                self.__singular = singular.ring("(integer,2,%d)"%(exp,), _vars, order=order, check=False)
            else:
                self.__singular = singular.ring("(integer,%d)"%(ch,), _vars, order=order, check=False)

        elif base_ring is ZZ:
            self.__singular = singular.ring("(integer)", _vars, order=order, check=False)
        else:
            raise TypeError("no conversion to a Singular ring defined")

        return self.__singular

def can_convert_to_singular(R):
    """
    Returns True if this ring's base field or ring can be
    represented in Singular, and the polynomial ring has at
    least one generator.  If this is True then this polynomial
    ring can be represented in Singular.

    The following base rings are supported: finite fields, rationals, number
    fields, and real and complex fields.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
        sage: can_convert_to_singular(PolynomialRing(QQ, names=['x']))
        True

        sage: can_convert_to_singular(PolynomialRing(QQ, names=[]))
        False

    """
    if R.ngens() == 0:
        return False;

    base_ring = R.base_ring()
    return ( sage.rings.finite_rings.constructor.is_FiniteField(base_ring)
             or is_RationalField(base_ring)
             or (base_ring.is_prime_field() and base_ring.characteristic() <= 2147483647)
             or is_RealField(base_ring)
             or is_ComplexField(base_ring)
             or is_RealDoubleField(base_ring)
             or is_ComplexDoubleField(base_ring)
             or number_field.number_field_base.is_NumberField(base_ring)
             or ( sage.rings.fraction_field.is_FractionField(base_ring) and ( base_ring.base_ring().is_prime_field() or base_ring.base_ring() is ZZ or is_FiniteField(base_ring.base_ring()) ) )
             or base_ring is ZZ
             or is_IntegerModRing(base_ring)
             or (is_RationalFunctionField(base_ring) and base_ring.constant_field().is_prime_field()) )


class Polynomial_singular_repr:
    """
    Implements coercion of polynomials to Singular polynomials.

    This class is a base class for all (univariate and multivariate)
    polynomial classes which support conversion from and to
    Singular polynomials.

    Due to the incompatibility of Python extension classes and multiple inheritance,
    this just defers to module-level functions.
    """
    def _singular_(self, singular=singular_default, have_ring=False):
        return _singular_func(self, singular, have_ring)

    def _singular_init_func(self, singular=singular_default, have_ring=False):
        return _singular_init_func(self, singular, have_ring)

def _singular_func(self, singular=singular_default, have_ring=False):
    """
    Return Singular polynomial matching this polynomial.

    INPUT:

    - ``singular`` - Singular instance to use.
    - ``have_ring`` - if True we will not attempt to set this element's ring as
      the current Singular ring. This is useful to speed up a batch of
      ``f._singular_()`` calls. However, it's dangerous as it might lead to wrong
      results if another ring is ``singular.current_ring()``. (Default: False)

    EXAMPLES::

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
        self.parent()._singular_(singular).set_ring() #this is expensive

    try:
        self.__singular._check_valid()
        if self.__singular.parent() is singular:
            return self.__singular
    except (AttributeError, ValueError):
        pass
#    return self._singular_init_(singular,have_ring=have_ring)
    return _singular_init_func(self, singular,have_ring=have_ring)

def _singular_init_func(self, singular=singular_default, have_ring=False):
    """
    Return corresponding Singular polynomial but enforce that a new
    instance is created in the Singular interpreter.

    Use ``self._singular_()`` instead.
    """
    if not have_ring:
        self.parent()._singular_(singular).set_ring() #this is expensive

    self.__singular = singular(str(self))

    return self.__singular
