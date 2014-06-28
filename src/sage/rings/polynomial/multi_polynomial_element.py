"""
Generic Multivariate Polynomials

AUTHORS:

- David Joyner: first version

- William Stein: use dict's instead of lists

- Martin Albrecht malb@informatik.uni-bremen.de: some functions added

- William Stein (2006-02-11): added better __div__ behavior.

- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of some
  Singular features

- William Stein (2006-04-19): added e.g.,
  ``f[1,3]`` to get coeff of `xy^3`; added examples of the new
  ``R.x,y = PolynomialRing(QQ,2)`` notation.

- Martin Albrecht: improved singular coercions (restructured class
  hierarchy) and added ETuples

- Robert Bradshaw (2007-08-14): added support for coercion of
  polynomials in a subset of variables (including multi-level
  univariate rings)

- Joel B. Mohler (2008-03): Refactored interactions with ETuples.

EXAMPLES:

We verify Lagrange's four squares identity::

    sage: R.<a0,a1,a2,a3,b0,b1,b2,b3> = QQbar[]
    sage: (a0^2 + a1^2 + a2^2 + a3^2)*(b0^2 + b1^2 + b2^2 + b3^2) == (a0*b0 - a1*b1 - a2*b2 - a3*b3)^2 + (a0*b1 + a1*b0 + a2*b3 - a3*b2)^2 + (a0*b2 - a1*b3 + a2*b0 + a3*b1)^2 + (a0*b3 + a1*b2 - a2*b1 + a3*b0)^2
    True
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


from sage.structure.element import CommutativeRingElement, canonical_coercion, coerce_binop


import sage.misc.misc as misc
import sage.rings.integer

import polydict

from sage.structure.factorization import Factorization

from sage.rings.polynomial.polynomial_singular_interface import Polynomial_singular_repr

from sage.structure.sequence import Sequence


from multi_polynomial import MPolynomial

def is_MPolynomial(x):
    return isinstance(x, MPolynomial)

class MPolynomial_element(MPolynomial):
    def __init__(self, parent, x):
        """
        EXAMPLE::

            sage: K.<cuberoot2> = NumberField(x^3 - 2)
            sage: L.<cuberoot3> = K.extension(x^3 - 3)
            sage: S.<sqrt2> = L.extension(x^2 - 2)
            sage: S
            Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
            sage: P.<x,y,z> = PolynomialRing(S) # indirect doctest
        """
        CommutativeRingElement.__init__(self, parent)
        self.__element = x

    def _repr_(self):
        """
        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(QQbar)
            sage: x + QQbar.random_element() # indirect doctest
            x + 4
        """
        return "%s"%self.__element

    ####################

    def __call__(self, *x, **kwds):
        """
        Evaluate this multi-variate polynomial at `x`, where
        `x` is either the tuple of values to substitute in, or one
        can use functional notation `f(a_0,a_1,a_2, \ldots)` to
        evaluate `f` with the ith variable replaced by
        `a_i`.

        EXAMPLES::

            sage: R.<x,y> = CC[]
            sage: f = x^2 + y^2
            sage: f(1,2)
            5.00000000000000
            sage: f((1,2))
            5.00000000000000

        ::

            sage: x = PolynomialRing(CC,3,'x').gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            (-2.00000000000000)*x1*x2 + x0 + x1
            sage: f(1,2,0)
            3.00000000000000
            sage: f(1,2,5)
            -17.0000000000000

        AUTHORS:

        - David Kohel (2005-09-27)
        """
        if len(kwds) > 0:
            f = self.subs(**kwds)
            if len(x) > 0:
                return f(*x)
            else:
                return f
        if len(x) == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
        n = self.parent().ngens()
        if len(x) != n:
            raise TypeError("x must be of correct length")
        if n == 0:
            return self
        try:
            K = x[0].parent()
        except AttributeError:
            K = self.parent().base_ring()
        y = K(0)
        for (m,c) in self.element().dict().iteritems():
            y += c*misc.mul([ x[i]**m[i] for i in range(n) if m[i] != 0])
        return y

    def __cmp__(self, right):
        """
        Compares right to self with respect to the term order of
        self.parent().

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(QQbar,3,order='lex')
            sage: x^1*y^2 > y^3*z^4
            True
            sage: x^3*y^2*z^4 < x^3*y^2*z^1
            False

        ::

            sage: R.<x,y,z>=PolynomialRing(CC,3,order='deglex')
            sage: x^1*y^2*z^3 > x^3*y^2*z^0
            True
            sage: x^1*y^2*z^4 < x^1*y^1*z^5
            False

        ::

            sage: R.<x,y,z>=PolynomialRing(QQbar,3,order='degrevlex')
            sage: x^1*y^5*z^2 > x^4*y^1*z^3
            True
            sage: x^4*y^7*z^1 < x^4*y^2*z^3
            False
        """
        try:
            return self.__element.compare(right.__element,
                             self.parent().term_order().compare_tuples)
        except AttributeError:
            return self.__element.compare(right.__element)

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQbar, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            3*x^5 + y^2
        """
        n = self.parent().ngens()
        if n == 0:
            return codomain._coerce_(self)
        y = codomain(0)
        for (m,c) in self.element().dict().iteritems():
            y += codomain(c)*misc.mul([ im_gens[i]**m[i] for i in range(n) if m[i] ])
        return y


    def _add_(self, right):
        #return self.parent()(self.__element + right.__element)
        return self.__class__(self.parent(),self.__element + right.__element)

    def _sub_(self, right):
        # return self.parent()(self.__element - right.__element)
        return self.__class__(self.parent(),self.__element - right.__element)

    def _mul_(self, right):
        #return self.parent()(self.__element * right.__element)
        return self.__class__(self.parent(),self.__element * right.__element)

    def _lmul_(self, a):
        """
        Left Scalar Multiplication

        EXAMPLES:

        Note that it is not really possible to do a meaningful
        example since sage mpoly rings refuse to have non-commutative
        bases.

        ::

            sage: R.<x,y> = QQbar[]
            sage: f = (x + y)
            sage: 3*f
            3*x + 3*y
        """
        return self.__class__(self.parent(),self.__element.scalar_lmult(a))

    def _rmul_(self, a):
        """
        Right Scalar Multiplication

        EXAMPLES:

        Note that it is not really possible to do a meaningful
        example since sage mpoly rings refuse to have non-commutative
        bases.

        ::

            sage: R.<x,y> = QQbar[]
            sage: f = (x + y)
            sage: f*3
            3*x + 3*y
        """
        return self.__class__(self.parent(),self.__element.scalar_rmult(a))

    def _div_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x,y> = CC['x,y']
            sage: f = (x + y)/x; f
            (x + y)/x
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over
            Complex Field with 53 bits of precision

        If dividing by a scalar, there is no need to go to the fraction
        field of the polynomial ring::

            sage: f = (x + y)/2; f
            0.500000000000000*x + 0.500000000000000*y
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Complex Field with
            53 bits of precision

        TESTS:

        Ensure that :trac:`13704` is fixed.::

            sage: R.<t>=PolynomialRing(QQ)
            sage: S.<x,y>=PolynomialRing(R)
            sage: x/S(2)
            1/2*x
        """
        if right in self.base_ring():
            inv = self.base_ring().one()/self.base_ring()(right)
            return inv*self
        return self.parent().fraction_field()(self, right, coerce=False)

    def __rpow__(self, n):
        if not isinstance(n, (int, long, sage.rings.integer.Integer)):
            raise TypeError("The exponent must be an integer.")
        return self.parent()(self.__element**n)

    def element(self):
        return self.__element

    def change_ring(self, R):
        return self.parent().change_ring(R)(self)


#class MPolynomial_macaulay2_repr:
#    """
#    Multivariate polynomials that are representable in Macaulay2.
#    """
#    def _macaulay2_(self, macaulay2=macaulay2):
#        """
#        Return corresponding Macaulay2 polynomial.
#
#        EXAMPLES:
#            sage: R.<x,y> = ZZ[]
#            sage: f = (x^3 + 2*y^2*x)^7; f
#            x^21 + 14*x^19*y^2 + 84*x^17*y^4 + 280*x^15*y^6 + 560*x^13*y^8 + 672*x^11*y^10 + 448*x^9*y^12 + 128*x^7*y^14
#            sage: macaulay2(R)                      # optional - M2
#            ZZ [x, y, MonomialOrder => GRevLex, MonomialSize => 16]
#            sage: h = f._macaulay2_()               # optional
#            sage: R(h)                              # optional
#            x^21 + 14*x^19*y^2 + 84*x^17*y^4 + 280*x^15*y^6 + 560*x^13*y^8 + 672*x^11*y^10 + 448*x^9*y^12 + 128*x^7*y^14
#            sage: R(h^20) == f^20                   # optional
#            True
#        """
#        try:
#            if self.__macaulay2.parent() is macaulay2:
#                return self.__macaulay2
#        except AttributeError:
#            pass
#        self.parent()._macaulay2_set_ring(macaulay2)
#        self.__macaulay2 = macaulay2(repr(self))
#        return self.__macaulay2
#

class MPolynomial_polydict(Polynomial_singular_repr, MPolynomial_element):
    r"""
    Multivariate polynomials implemented in pure python using
    polydicts.
    """
    def __init__(self, parent, x):
        """
        EXAMPLES::

            sage: R, x = PolynomialRing(QQbar, 10, 'x').objgens()
            sage: x
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
            sage: loads(dumps(x)) == x
            True
        """
        if not isinstance(x, polydict.PolyDict):
            x = polydict.PolyDict(x, parent.base_ring()(0), remove_zero=True)
        MPolynomial_element.__init__(self, parent, x)

    def _new_constant_poly(self, x, P):
        """
        Quickly create a new constant polynomial with value x in parent P.

        ASSUMPTION:

        x must be an element of the base ring of P. That assumption is
        not verified.

        EXAMPLE::

            sage: R.<x,y> = QQ['t'][]
            sage: x._new_constant_poly(R.base_ring()(2),R)
            2

        """
        return MPolynomial_polydict(P, {P._zero_tuple:x})

    def __neg__(self):
        """
        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: -x
            -x
            sage: -(y-1)
            -y + 1
        """
        return self*(-1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: repr(-x^2-y+1)  # indirect doc-test
            '-x^2 - y + 1'
            sage: K.<I>=QuadraticField(-1)
            sage: R.<x,y>=K[]
            sage: repr(-I*y-x^2)  # indirect doc-test
            '-x^2 + (-I)*y'
        """
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self.element().poly_repr(self.parent().variable_names(),
                                        atomic_coefficients=atomic, cmpfn=cmpfn )

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: latex(-x^2-y+1)
            -x^{2} - y + 1
            sage: K.<I>=QuadraticField(-1)
            sage: R.<x,y>=K[]
            sage: latex(-I*y+I*x^2)
            \left(\sqrt{-1}\right) x^{2} + \left(-\sqrt{-1}\right) y
        """
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self.element().latex(self.parent().latex_variable_names(),
                                    atomic_coefficients=atomic, cmpfn=cmpfn)

    def _repr_with_changed_varnames(self, varnames):
        """
        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: f=-x^2-y+1
            sage: f._repr_with_changed_varnames(['jack','jill'])
            '-jack^2 - jill + 1'
        """
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self.element().poly_repr(varnames,
                                        atomic_coefficients=atomic, cmpfn=cmpfn)

    def degrees(self):
        r"""
        Returns a tuple (precisely - an ``ETuple``) with the
        degree of each variable in this polynomial. The list of degrees is,
        of course, ordered by the order of the generators.

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(QQbar)
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.degrees()
            (2, 2, 0)
            sage: f = x^2+z^2
            sage: f.degrees()
            (2, 0, 2)
            sage: f.total_degree()  # this simply illustrates that total degree is not the sum of the degrees
            2
            sage: R.<x,y,z,u>=PolynomialRing(QQbar)
            sage: f=(1-x)*(1+y+z+x^3)^5
            sage: f.degrees()
            (16, 5, 5, 0)
            sage: R(0).degrees()
            (0, 0, 0, 0)
        """
        if self.is_zero():
            return polydict.ETuple({},self.parent().ngens())
        else:
            return self._MPolynomial_element__element.max_exp()

    def degree(self, x=None):
        """
        Return the degree of self in x, where x must be one of the
        generators for the parent of self.

        INPUT:


        -  ``x`` - multivariate polynomial (a generator of the
           parent of self) If x is not specified (or is None), return the
           total degree, which is the maximum degree of any monomial.


        OUTPUT: integer

        EXAMPLE::

            sage: R.<x,y> = RR[]
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10
        """
        if x is None:
            return self.element().degree(None)
        if not (isinstance(x, MPolynomial) and x.parent() is self.parent() and x.is_generator()):
            raise TypeError("x must be one of the generators of the parent.")
        return self.element().degree(x.element())

    def total_degree(self):
        """
        Return the total degree of self, which is the maximum degree of any
        monomial in self.

        EXAMPLES::

            sage: R.<x,y,z> = QQbar[]
            sage: f=2*x*y^3*z^2
            sage: f.total_degree()
            6
            sage: f=4*x^2*y^2*z^3
            sage: f.total_degree()
            7
            sage: f=99*x^6*y^3*z^9
            sage: f.total_degree()
            18
            sage: f=x*y^3*z^6+3*x^2
            sage: f.total_degree()
            10
            sage: f=z^3+8*x^4*y^5*z
            sage: f.total_degree()
            10
            sage: f=z^9+10*x^4+y^8*x^2
            sage: f.total_degree()
            10
        """
        return self.degree()

    def monomial_coefficient(self, mon):
        """
        Return the coefficient in the base ring of the monomial mon in
        self, where mon must have the same parent as self.

        This function contrasts with the function
        ``coefficient`` which returns the coefficient of a
        monomial viewing this polynomial in a polynomial ring over a base
        ring having fewer variables.

        INPUT:

        -  ``mon`` - a monomial


        OUTPUT: coefficient in base ring

        .. seealso::

           For coefficients in a base ring of fewer variables, look
           at :meth:`coefficient`.

        EXAMPLES:

        The parent of the return is a member of the base ring.

        ::

            sage: R.<x,y>=QQbar[]

        The parent of the return is a member of the base ring.

        ::

            sage: f = 2 * x * y
            sage: c = f.monomial_coefficient(x*y); c
            2
            sage: c.parent()
            Algebraic Field

        ::

            sage: f = y^2 + y^2*x - x^9 - 7*x + 5*x*y
            sage: f.monomial_coefficient(y^2)
            1
            sage: f.monomial_coefficient(x*y)
            5
            sage: f.monomial_coefficient(x^9)
            -1
            sage: f.monomial_coefficient(x^10)
            0

        ::

            sage: var('a')
            a
            sage: K.<a> = NumberField(a^2+a+1)
            sage: P.<x,y> = K[]
            sage: f=(a*x-1)*((a+1)*y-1); f
            -x*y + (-a)*x + (-a - 1)*y + 1
            sage: f.monomial_coefficient(x)
            -a
        """
        if not (isinstance(mon, MPolynomial) and mon.parent() is self.parent() and mon.is_monomial()):
            raise TypeError("mon must be a monomial in the parent of self.")
        R = self.parent().base_ring()
        return R(self.element().monomial_coefficient(mon.element().dict()))

    def dict(self):
        """
        Return underlying dictionary with keys the exponents and values
        the coefficients of this polynomial.
        """
        return self.element().dict()

    #def __iter__(self):
    #    """
    #    Facilitates iterating over the monomials of self,
    #    returning tuples of the form (coeff, mon) for each
    #    non-zero monomial.
    #
    #    EXAMPLES:
    #        sage: R = ZZ['t']
    #        sage: P.<x,y,z> = PolynomialRing(R,3)
    #        sage: f = 3*x^3*y + 16*x + 7
    #        sage: [(c,m) for c,m in f]
    #        [(3, x^3*y), (16, x), (7, 1)]
    #        sage: f = P.random_element(10,10)
    #        sage: sum(c*m for c,m in f) == f
    #        True
    #    """
    #    exps = self.exponents()
    #    parent = self.parent()
    #    for exp in exps:
    #        yield self.element()[exp], MPolynomial_polydict(parent, {exp: 1})

    def __getitem__(self, x):
        """
        INPUT:


        -  ``x`` - a tuple or, in case of a single-variable
           MPolynomial ring x can also be an integer.


        EXAMPLES::

            sage: R.<x, y> = PolynomialRing(QQbar, 2)
            sage: f = -10*x^3*y + 17*x*y
            sage: f[3,1]
            -10
            sage: f[1,1]
            17
            sage: f[0,1]
            0

        ::

            sage: R.<x> = PolynomialRing(QQbar,1); R
            Multivariate Polynomial Ring in x over Algebraic Field
            sage: f = 5*x^2 + 3; f
            5*x^2 + 3
            sage: f[2]
            5
        """
        if isinstance(x, MPolynomial):
            return self.monomial_coefficient(x)
        if not isinstance(x, tuple):
            try:
                x = tuple(x)
            except TypeError:
                x = (x, )
        try:
            return self.element()[x]
        except KeyError:
            return self.parent().base_ring()(0)

    def coefficient(self, degrees):
        """
        Return the coefficient of the variables with the degrees specified
        in the python dictionary ``degrees``. Mathematically,
        this is the coefficient in the base ring adjoined by the variables
        of this ring not listed in ``degrees``. However, the
        result has the same parent as this polynomial.

        This function contrasts with the function
        ``monomial_coefficient`` which returns the coefficient
        in the base ring of a monomial.

        INPUT:


        -  ``degrees`` - Can be any of:

           -  a dictionary of degree restrictions

           -  a list of degree restrictions (with None in
              the unrestricted variables)

           -  a monomial (very fast, but not as flexible)


        OUTPUT: element of the parent of self

        .. seealso::

           For coefficients of specific monomials, look at
           :meth:`monomial_coefficient`.

        EXAMPLES::

            sage: R.<x, y> = QQbar[]
            sage: f = 2 * x * y
            sage: c = f.coefficient({x:1,y:1}); c
            2
            sage: c.parent()
            Multivariate Polynomial Ring in x, y over Algebraic Field
            sage: c in PolynomialRing(QQbar, 2, names = ['x','y'])
            True
            sage: f = y^2 - x^9 - 7*x + 5*x*y
            sage: f.coefficient({y:1})
            5*x
            sage: f.coefficient({y:0})
            -x^9 + (-7)*x
            sage: f.coefficient({x:0,y:0})
            0
            sage: f=(1+y+y^2)*(1+x+x^2)
            sage: f.coefficient({x:0})
            y^2 + y + 1
            sage: f.coefficient([0,None])
            y^2 + y + 1
            sage: f.coefficient(x)
            y^2 + y + 1
            sage: # Be aware that this may not be what you think!
            sage: # The physical appearance of the variable x is deceiving -- particularly if the exponent would be a variable.
            sage: f.coefficient(x^0) # outputs the full polynomial
            x^2*y^2 + x^2*y + x*y^2 + x^2 + x*y + y^2 + x + y + 1

        ::

            sage: R.<x,y> = RR[]
            sage: f=x*y+5
            sage: c=f.coefficient({x:0,y:0}); c
            5.00000000000000
            sage: parent(c)
            Multivariate Polynomial Ring in x, y over Real Field with 53 bits of precision

        AUTHORS:

        - Joel B. Mohler (2007-10-31)
        """
        looking_for = None
        if isinstance(degrees, MPolynomial) and degrees.parent() == self.parent() and degrees.is_monomial():
            looking_for = [e if e > 0 else None for e in degrees.exponents()[0]]
        elif isinstance(degrees, list):
            looking_for = degrees
        elif isinstance(degrees, dict):
            poly_vars = self.parent().gens()
            looking_for = [None] * len(poly_vars)
            for d, exp in degrees.items():
                for i in range(len(poly_vars)):
                    if d == poly_vars[i]:
                        looking_for[i] = exp
        if not looking_for:
            raise ValueError("You must pass a dictionary list or monomial.")
        return self.parent()(self.element().polynomial_coefficient(looking_for))

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES::

            sage: R.<a,b,c> = PolynomialRing(QQbar, 3)
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        try:
            return self.__exponents
        except AttributeError:
            self.__exponents = self.element().dict().keys()
            try:
                self.__exponents.sort(cmp=self.parent().term_order().compare_tuples, reverse=True)
            except AttributeError:
                pass
            return self.__exponents

    def is_unit(self):
        """
        Return True if self is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: (x+y).is_unit()
            False
            sage: R(0).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(-1 + x).is_unit()
            False
            sage: R(2).is_unit()
            True
        """
        d = self.element().dict()
        k = d.keys()
        if len(k) != 1:
            return False
        k = k[0]
        if k != polydict.ETuple([0]*self.parent().ngens()):
            return False
        return bool(d[k].is_unit())

    def inverse_of_unit(self):
        d = self.element().dict()
        k = d.keys()
        if len(k) != 1:
            raise ArithmeticError("is not a unit")
        k = k[0]
        if k != polydict.ETuple([0]*self.parent().ngens()):
            raise ArithmeticError("is not a unit")
        return ~d[k]

    def is_homogeneous(self):
        """
        Return True if self is a homogeneous polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: (x+y).is_homogeneous()
            True
            sage: (x.parent()(0)).is_homogeneous()
            True
            sage: (x+y^2).is_homogeneous()
            False
            sage: (x^2 + y^2).is_homogeneous()
            True
            sage: (x^2 + y^2*x).is_homogeneous()
            False
            sage: (x^2*y + y^2*x).is_homogeneous()
            True
        """
        return self.element().is_homogeneous()

    def _homogenize(self, var):
        r"""
        Return ``self`` if ``self`` is homogeneous.
        Otherwise return a homogenized polynomial constructed by modifying
        the degree of the variable with index ``var``.

        INPUT:


        -  ``var`` - an integer indicating which variable to
           use to homogenize (0 <= var < parent(self).ngens())


        OUTPUT: a multivariate polynomial

        EXAMPLES::

            sage: P.<x,y> = QQbar[]
            sage: f = x^2 + y + 1 + 5*x*y^1
            sage: g = f.homogenize('z'); g # indirect doctest
            x^2 + 5*x*y + y*z + z^2
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Algebraic Field

        SEE: ``self.homogenize``
        """
        if self.is_homogeneous():
            return self
        X = self.element().homogenize(var)
        R = self.parent()
        return R(X)

    def is_generator(self):
        """
        Returns True if self is a generator of it's parent.

        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: x.is_generator()
            True
            sage: (x+y-y).is_generator()
            True
            sage: (x*y).is_generator()
            False
        """
        d = self.element().dict()
        if len(d) == 1:
            e,c = d.items()[0]
            if c.is_one() and len(e.nonzero_positions()) == 1 and e.nonzero_values()[0] == 1:
                return True
        return False

    def is_monomial(self):
        """
        Returns True if self is a monomial, which we define to be a
        product of generators with coefficient 1.

        Use is_term to allow the coefficient to not be 1.

        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: x.is_monomial()
            True
            sage: (x+2*y).is_monomial()
            False
            sage: (2*x).is_monomial()
            False
            sage: (x*y).is_monomial()
            True

        To allow a non-1 leading coefficient, use is_term()::

            sage: (2*x*y).is_term()
            True
            sage: (2*x*y).is_monomial()
            False
        """
        term = (len(self.element().dict().keys()) == 1)
        if term:
            if self.coefficients()[0] == 1:
                return True
            else:
                return False
        else:
            return False

    def is_term(self):
        """
        Returns True if self is a term, which we define to be a
        product of generators times some coefficient, which need
        not be 1.

        Use is_monomial to require that the coefficent be 1.

        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: x.is_term()
            True
            sage: (x+2*y).is_term()
            False
            sage: (2*x).is_term()
            True
            sage: (7*x^5*y).is_term()
            True

        To require leading coefficient 1, use is_monomial()::

            sage: (2*x*y).is_monomial()
            False
            sage: (2*x*y).is_term()
            True
        """
        return len(self.element().dict().keys()) == 1

    def subs(self, fixed=None, **kw):
        """
        Fixes some given variables in a given multivariate polynomial and
        returns the changed multivariate polynomials. The polynomial itself
        is not affected. The variable,value pairs for fixing are to be
        provided as a dictionary of the form {variable:value}.

        This is a special case of evaluating the polynomial with some of
        the variables constants and the others the original variables.

        INPUT:


        -  ``fixed`` - (optional) dictionary of inputs

        -  ``**kw`` - named parameters


        OUTPUT: new MPolynomial

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
        """
        variables = list(self.parent().gens())
        for i in range(0,len(variables)):
            if str(variables[i]) in kw:
                variables[i]=kw[str(variables[i])]
            elif fixed and variables[i] in fixed:
                variables[i] = fixed[variables[i]]
        return self(tuple(variables))

    def monomials(self):
        """
        Returns the list of monomials in self. The returned list is
        decreasingly ordered by the term ordering of self.parent().

        OUTPUT: list of MPolynomials representing Monomials

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.monomials()
            [x^2*y^2, x^2, y, 1]

        ::

            sage: R.<fx,fy,gx,gy> = QQbar[]
            sage: F = ((fx*gy - fy*gx)^3)
            sage: F
            -fy^3*gx^3 + 3*fx*fy^2*gx^2*gy + (-3)*fx^2*fy*gx*gy^2 + fx^3*gy^3
            sage: F.monomials()
            [fy^3*gx^3, fx*fy^2*gx^2*gy, fx^2*fy*gx*gy^2, fx^3*gy^3]
            sage: F.coefficients()
            [-1, 3, -3, 1]
            sage: sum(map(mul,zip(F.coefficients(),F.monomials()))) == F
            True
        """
        ring = self.parent()
        one = ring.base_ring()(1)
        return [MPolynomial_polydict(ring, polydict.PolyDict({m:one}, force_int_exponents=False, force_etuples=False)) for m in self.exponents()]
        try:
            return self.__monomials
        except AttributeError:
            ring = self.parent()
            one = self.parent().base_ring()(1)
            self.__monomials = sorted([ MPolynomial_polydict(ring, polydict.PolyDict( {m:one}, force_int_exponents=False,  force_etuples=False ) ) \
                                for m in self._MPolynomial_element__element.dict().keys() ], reverse=True)
            return self.__monomials

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        #v = (0,)*int(self.parent().ngens())
        d = self.element().dict()
        try:
            return d[polydict.ETuple({},self.parent().ngens())]
        except KeyError:
            return self.parent().base_ring()(0)

    def is_univariate(self):
        """
        Returns True if this multivariate polynomial is univariate and
        False otherwise.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.is_univariate()
            False
            sage: g = f.subs({x:10}); g
            700*y^2 + (-2)*y + 305
            sage: g.is_univariate()
            True
            sage: f = x^0
            sage: f.is_univariate()
            True
        """
        mons = self.element().dict().keys()

        found = -1
        for mon in mons:
            for i in mon.nonzero_positions():
                if found != i:
                    if found != -1:
                        return False
                    else:
                        found = i
        return True

    def univariate_polynomial(self, R=None):
        """
        Returns a univariate polynomial associated to this multivariate
        polynomial.

        INPUT:


        -  ``R`` - (default: None) PolynomialRing


        If this polynomial is not in at most one variable, then a
        ValueError exception is raised. This is checked using the
        is_univariate() method. The new Polynomial is over the same base
        ring as the given MPolynomial.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f.subs({x:10}); g
            700*y^2 + (-2)*y + 305
            sage: g.univariate_polynomial ()
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305

        TESTS::

            sage: P = PolynomialRing(QQ, 0, '')
            sage: P(5).univariate_polynomial()
            5
        """
        if self.parent().ngens() == 0:
            if R is None:
                return self.base_ring()(self)
            else:
                return R(self)

        if not self.is_univariate():
            raise TypeError("polynomial must involve at most one variable")

        #construct ring if None
        if R is None:
            # constant, we just pick first variable from parent
            if self.is_constant():
                R = self.base_ring()[self.parent().variable_names()[0]]
            else:
                R = self.base_ring()[str(self.variables()[0])]

        monomial_coefficients = self._MPolynomial_element__element.dict()

        if( not self.is_constant() ):
            var_idx = self.degrees().nonzero_positions()[0] #variable
        else:
            var_idx = 0; #constant
            if( len(monomial_coefficients.keys())==0 ):
                return R(0)

        #construct list
        lookup = [int(0),]*len( monomial_coefficients.keys()[0] )
        coefficients = []
        for degree in range( 0 , max([ m[var_idx] for m in monomial_coefficients.keys() ])+1 ):
            lookup[var_idx]=int(degree);
            try:
                coefficients.append( monomial_coefficients[ polydict.ETuple(lookup) ] ) #if we find something, add the coefficient
            except KeyError:
                coefficients.append( 0 ) #else add zero

        #construct polynomial
        return R(coefficients)

    def variables(self):
        """
        Returns the tuple of variables occurring in this polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.variables()
            (x, y)
            sage: g = f.subs({x:10}); g
            700*y^2 + (-2)*y + 305
            sage: g.variables()
            (y,)

        TESTS:

        This shows that the issue at :trac:`7077` is fixed::

            sage: x,y,z=polygens(QQ,'x,y,z')
            sage: (x^2).variables()
            (x,)
        """
        return tuple([self.parent().gen(index) for index in self.degrees().nonzero_positions()])

    def variable(self,i):
        """
        Returns `i`-th variable occurring in this polynomial.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.variable(0)
            x
            sage: f.variable(1)
            y
        """
        return self.variables()[int(i)]

    def nvariables(self):
        """
        Number of variables in this polynomial

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.nvariables ()
            2
            sage: g = f.subs({x:10}); g
            700*y^2 + (-2)*y + 305
            sage: g.nvariables ()
            1
        """
        return len(self.degrees().nonzero_positions())

    def is_constant(self):
        """
        True if polynomial is constant, and False otherwise.

        EXAMPLES::

            sage: R.<x,y> = QQbar[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.is_constant()
            False
            sage: g = 10*x^0
            sage: g.is_constant()
            True
        """
        if len(self.dict()) <= 1 and len(self.degrees().nonzero_positions()) == 0:
            return True
        else:
            return False

    def lm(self):
        """
        Returns the lead monomial of self with respect to the term order of
        self.parent().

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: (x^1*y^2 + y^3*z^4).lm()
            x*y^2
            sage: (x^3*y^2*z^4 + x^3*y^2*z^1).lm()
            x^3*y^2*z^4

        ::

            sage: R.<x,y,z>=PolynomialRing(CC,3,order='deglex')
            sage: (x^1*y^2*z^3 + x^3*y^2*z^0).lm()
            x*y^2*z^3
            sage: (x^1*y^2*z^4 + x^1*y^1*z^5).lm()
            x*y^2*z^4

        ::

            sage: R.<x,y,z>=PolynomialRing(QQbar,3,order='degrevlex')
            sage: (x^1*y^5*z^2 + x^4*y^1*z^3).lm()
            x*y^5*z^2
            sage: (x^4*y^7*z^1 + x^4*y^2*z^3).lm()
            x^4*y^7*z

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
            sage: R.<x,y>=MPolynomialRing_polydict(GF(2),2,order='lex')
            sage: f=x+y
            sage: f.lm()
            x

        """
        try:
            return self.__lm
        except AttributeError:
            if self.is_zero():
                return self
            R = self.parent()
            f = self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )
            one = R.base_ring()(1)
            self.__lm = MPolynomial_polydict(R,polydict.PolyDict({f:one},zero=R.base_ring().zero_element(),force_int_exponents=False,  force_etuples=False))
            return self.__lm

    def lc(self):
        """
        Returns the leading coefficient of self i.e.,
        self.coefficient(self.lm())

        EXAMPLES::

            sage: R.<x,y,z>=QQbar[]
            sage: f=3*x^2-y^2-x*y
            sage: f.lc()
            3
        """
        try:
            return self.__lc
        except AttributeError:
            if self.is_zero():
                return self.base_ring()._zero_element
            R = self.parent()
            f = self._MPolynomial_element__element.dict()
            self.__lc = f[self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )]
            return self.__lc

    def lt(self):
        """
        Returns the leading term of self i.e., self.lc()\*self.lm(). The
        notion of "leading term" depends on the ordering defined in the
        parent ring.

        EXAMPLES::

            sage: R.<x,y,z>=PolynomialRing(QQbar)
            sage: f=3*x^2-y^2-x*y
            sage: f.lt()
            3*x^2
            sage: R.<x,y,z>=PolynomialRing(QQbar,order="invlex")
            sage: f=3*x^2-y^2-x*y
            sage: f.lt()
            -y^2

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
            sage: R.<x,y>=MPolynomialRing_polydict(GF(2),2,order='lex')
            sage: f=x+y
            sage: f.lt()
            x

        """
        try:
            return self.__lt
        except AttributeError:
            if self.is_zero():
                return self
            R = self.parent()
            f = self._MPolynomial_element__element.dict()
            res = self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )
            self.__lt = MPolynomial_polydict(R,polydict.PolyDict({res:f[res]},zero=R.base_ring().zero_element(),force_int_exponents=False, force_etuples=False))
            return self.__lt

    def __eq__(self,right):
        """

        """
        if not isinstance(right,MPolynomial_polydict):
            # we want comparison with zero to be fast
            if right == 0:
                if self._MPolynomial_element__element.dict()=={}:
                    return True
                else:
                    return False
            return self._richcmp_(right,2)
        return self._MPolynomial_element__element == right._MPolynomial_element__element

    def __ne__(self,right):
        """

        """
        if not isinstance(right,MPolynomial_polydict):
            # we want comparison with zero to be fast
            if right == 0:
                if self._MPolynomial_element__element.dict()=={}:
                    return False
                else:
                    return True
            # maybe add constant elements as well
            return self._richcmp_(right,3)
        return self._MPolynomial_element__element != right._MPolynomial_element__element

    def __nonzero__(self):
        """
        Returns True if self != 0

        .. note::

           This is much faster than actually writing ``self == 0``.
        """
        return self._MPolynomial_element__element.dict()!={}

    def __floordiv__(self,right):
        r"""
        Quotient of division of self by other. This is denoted //.

        .. note::

           It's not clear to me that this is well-defined if
           ``self`` is not exactly divisible by other.

        EXAMPLES::

            sage: R.<x,y>=QQbar[]
            sage: 2*x*y//y
            2*x
            sage: 2*x//y
            0
            sage: 2*x//4
            1/2*x
            sage: type(0//y)
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>
        """
        if not isinstance(self, type(right)) or self.parent() is not right.parent():
            self, right = canonical_coercion(self, right)
            return self // right  # this looks like recursion, but, in fact, it may be that self, right are a totally new composite type
        # handle division by monomials without using Singular
        if len(right.dict()) == 1:
            P = self.parent()
            ret = P(0)
            denC,denM = iter(right).next()
            for c,m in self:
                t = c*m
                if denC.divides(c) and P.monomial_divides(denM, m):
                    ret += P.monomial_quotient(t, right, coeff=True)
            return ret

        Q, _ = self.quo_rem(right)
        return Q

    def _derivative(self, var=None):
        r"""
        Differentiates ``self`` with respect to variable ``var``.

        If ``var`` is not one of the generators of this ring, _derivative(var)
        is called recursively on each coefficient of this polynomial.

        .. SEEALSO::

            :meth:`derivative`

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQbar)
            sage: S.<x, y> = PolynomialRing(R)
            sage: f = (t^2 + O(t^3))*x^2*y^3 + (37*t^4 + O(t^5))*x^3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Power Series Ring in t over Algebraic Field
            sage: f._derivative(x)   # with respect to x
            (2*t^2 + O(t^3))*x*y^3 + (111*t^4 + O(t^5))*x^2
            sage: f._derivative(x).parent()
            Multivariate Polynomial Ring in x, y over Power Series Ring in t over Algebraic Field
            sage: f._derivative(y)   # with respect to y
            (3*t^2 + O(t^3))*x^2*y^2
            sage: f._derivative(t)   # with respect to t (recurses into base ring)
            (2*t + O(t^2))*x^2*y^3 + (148*t^3 + O(t^4))*x^3
            sage: f._derivative(x)._derivative(y) # with respect to x and then y
            (6*t^2 + O(t^3))*x*y^2
            sage: f.derivative(y, 3) # with respect to y three times
            (6*t^2 + O(t^3))*x^2
            sage: f._derivative()    # can't figure out the variable
            Traceback (most recent call last):
            ...
            ValueError: must specify which variable to differentiate with respect to
        """
        if var is None:
            raise ValueError("must specify which variable to differentiate with respect to")

        gens = list(self.parent().gens())

        # check if var is one of the generators
        try:
            index = gens.index(var)
        except ValueError:
            # var is not a generator; do term-by-term differentiation recursively
            # var may be, for example, a generator of the base ring
            d = dict([(e, x._derivative(var)) for (e, x) in self.dict().iteritems()])
            d = polydict.PolyDict(d, self.parent().base_ring()(0), remove_zero=True)
            return MPolynomial_polydict(self.parent(), d)

        # differentiate w.r.t. indicated variable
        d = {}
        v = polydict.ETuple({index:1}, len(gens))
        for (exp, coeff) in self.dict().iteritems():
            if exp[index] > 0:
                d[exp.esub(v)] = coeff * exp[index]
        d = polydict.PolyDict(d, self.parent().base_ring()(0), remove_zero=True)
        return MPolynomial_polydict(self.parent(), d)

    def integral(self, var=None):
        r"""
        Integrates ``self`` with respect to variable ``var``.

        .. NOTE::

            The integral is always chosen so the constant term is 0.

        If ``var`` is not one of the generators of this ring, integral(var)
        is called recursively on each coefficient of this polynomial.

        EXAMPLES:

        On polynomials with rational coefficients::

            sage: x, y = PolynomialRing(QQ, 'x, y').gens()
            sage: ex = x*y + x - y
            sage: it = ex.integral(x); it
            1/2*x^2*y + 1/2*x^2 - x*y
            sage: it.parent() == x.parent()
            True

        On polynomials with coefficients in power series::

            sage: R.<t> = PowerSeriesRing(QQbar)
            sage: S.<x, y> = PolynomialRing(R)
            sage: f = (t^2 + O(t^3))*x^2*y^3 + (37*t^4 + O(t^5))*x^3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Power Series Ring in t over Algebraic Field
            sage: f.integral(x)   # with respect to x
            (1/3*t^2 + O(t^3))*x^3*y^3 + (37/4*t^4 + O(t^5))*x^4
            sage: f.integral(x).parent()
            Multivariate Polynomial Ring in x, y over Power Series Ring in t over Algebraic Field

            sage: f.integral(y)   # with respect to y
            (1/4*t^2 + O(t^3))*x^2*y^4 + (37*t^4 + O(t^5))*x^3*y
            sage: f.integral(t)   # with respect to t (recurses into base ring)
            (1/3*t^3 + O(t^4))*x^2*y^3 + (37/5*t^5 + O(t^6))*x^3

        TESTS::

            sage: f.integral()    # can't figure out the variable
            Traceback (most recent call last):
            ...
            ValueError: must specify which variable to integrate with respect to
        """
        if var is None:
            raise ValueError("must specify which variable to integrate "
                             "with respect to")

        gens = list(self.parent().gens())

        # check if var is one of the generators
        try:
            index = gens.index(var)
        except ValueError:
            # var is not a generator; do term-by-term integration recursively
            # var may be, for example, a generator of the base ring
            d = dict([(e, x.integral(var))
                      for (e, x) in self.dict().iteritems()])
            d = polydict.PolyDict(d, self.parent().base_ring()(0),
                                  remove_zero=True)
            return MPolynomial_polydict(self.parent(), d)

        # integrate w.r.t. indicated variable
        d = {}
        v = polydict.ETuple({index:1}, len(gens))
        for (exp, coeff) in self.dict().iteritems():
            d[exp.eadd(v)] = coeff / (1+exp[index])
        d = polydict.PolyDict(d, self.parent().base_ring()(0), remove_zero=True)
        return MPolynomial_polydict(self.parent(), d)

    def factor(self, proof=True):
        r"""
        Compute the irreducible factorization of this polynomial.

        INPUT:

        - ``proof'' - insist on provably correct results (ignored, always ``True``)

        ALGORITHM: Use univariate factorization code.

        If a polynomial is univariate, the appropriate univariate
        factorization code is called::

            sage: R.<z> = PolynomialRing(CC,1)
            sage: f = z^4 - 6*z + 3
            sage: f.factor()
            (z - 1.60443920904349) * (z - 0.511399619393097) * (z + 1.05791941421830 - 1.59281852704435*I) * (z + 1.05791941421830 + 1.59281852704435*I)

        TESTS:

        Check if we can handle polynomials with no variables, see :trac:`7950`::

            sage: P = PolynomialRing(ZZ,0,'')
            sage: res = P(10).factor(); res
            2 * 5
            sage: res[0][0].parent()
            Multivariate Polynomial Ring in no variables over Integer Ring
            sage: R = PolynomialRing(QQ,0,'')
            sage: res = R(10).factor(); res
            10
            sage: res.unit().parent()
            Rational Field
            sage: P(0).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: Prime factorization of 0 not defined.

        Check if we can factor a constant polynomial, see :trac:`8207`::

            sage: R.<x,y> = CC[]
            sage: R(1).factor()
            1.00000000000000

        Check that we prohibit too large moduli, :trac:`11829`::

            sage: R.<x,y> = GF(previous_prime(2^31))[]
            sage: factor(x+y+1,proof=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.
        """
        R = self.parent()

        # raise error if trying to factor zero
        if self == 0:
            raise ArithmeticError("Prime factorization of 0 not defined.")

        # if number of variables is zero ...
        if R.ngens() == 0:
            base_ring = self.base_ring()
            if base_ring.is_field():
                return Factorization([],unit=self.base_ring()(self))
            else:
                F = base_ring(self).factor()
                return Factorization([(R(f),m) for f,m in F], unit=F.unit())

        # try to use univariate factoring
        try:
            F = self.univariate_polynomial().factor()
            return Factorization([(R(f),m) for f,m in F], unit=F.unit())
        except TypeError:
            pass

        base_ring = self.base_ring()
        if base_ring.is_finite():
            if base_ring.characteristic() > 1<<29:
                raise NotImplementedError("Factorization of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.")
        if proof:
            raise NotImplementedError("proof = True factorization not implemented.  Call factor with proof=False.")

        R._singular_().set_ring()
        S = self._singular_().factorize()
        factors = S[1]
        exponents = S[2]
        v = sorted([(R(factors[i+1]), sage.rings.integer.Integer(exponents[i+1])) \
                        for i in range(len(factors))])
        unit = R(1)
        for i in range(len(v)):
            if v[i][0].is_unit():
                unit = unit * v[i][0]
                del v[i]
                break
        F = Factorization(v, unit=unit)
        F.sort()
        return F

    def lift(self,I):
        """
        given an ideal I = (f_1,...,f_r) and some g (== self) in I, find
        s_1,...,s_r such that g = s_1 f_1 + ... + s_r f_r

        ALGORITHM: Use Singular.

        EXAMPLE::

            sage: A.<x,y> = PolynomialRing(CC,2,order='degrevlex')
            sage: I = A.ideal([x^10 + x^9*y^2, y^8 - x^2*y^7 ])
            sage: f = x*y^13 + y^12
            sage: M = f.lift(I)
            sage: M
            [y^7, x^7*y^2 + x^8 + x^5*y^3 + x^6*y + x^3*y^4 + x^4*y^2 + x*y^5 + x^2*y^3 + y^4]
            sage: sum( map( mul , zip( M, I.gens() ) ) ) == f
            True
        """
        fs = self._singular_()
        Is = I._singular_()
        P = I.ring()
        try:
            M = Is.lift(fs)._sage_(P)
        except TypeError:
            raise ArithmeticError("f is not in I")
        return Sequence(M.list(), P, check=False, immutable=True)

    @coerce_binop
    def quo_rem(self, right):
        """
        Returns quotient and remainder of self and right.

        EXAMPLE::

            sage: R.<x,y> = CC[]
            sage: f = y*x^2 + x + 1
            sage: f.quo_rem(x)
            (x*y + 1.00000000000000, 1.00000000000000)

        ALGORITHM: Use Singular.
        """
        R = self.parent()
        R._singular_().set_ring()
        X = self._singular_().division(right._singular_())
        return R(X[1][1,1]), R(X[2][1])

    def resultant(self, other, variable=None):
        """
        Compute the resultant of ``self`` and ``other`` with respect
        to ``variable``.

        If a second argument is not provided, the first variable of
        ``self.parent()`` is chosen.

        INPUT:

        - ``other`` -- polynomial in ``self.parent()``

        - ``variable`` -- (optional) variable (of type polynomial) in
          ``self.parent()``

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ, 2)
            sage: a = x + y
            sage: b = x^3 - y^3
            sage: a.resultant(b)
            -2*y^3
            sage: a.resultant(b, y)
            2*x^3

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y> = MPolynomialRing_polydict_domain(QQ, 2, order='degrevlex')
            sage: a = x + y
            sage: b = x^3 - y^3
            sage: a.resultant(b)
            -2*y^3
            sage: a.resultant(b, y)
            2*x^3

        Check that :trac:`15061` is fixed::

            sage: R.<x, y> = AA[]
            sage: (x^2 + 1).resultant(x^2 - y)
            y^2 + 2*y + 1

        """
        R = self.parent()
        if variable is None:
            variable = R.gen(0)
        if R._has_singular:
            rt = self._singular_().resultant(other._singular_(), variable._singular_())
            r = rt.sage_poly(R)
        else:
            r = self.sylvester_matrix(other, variable).det()
        if R.ngens() <= 1 and r.degree() <= 0:
            return R.base_ring()(r[0])
        else:
            return r

    def reduce(self, I):
        """
        Reduce this polynomial by the the polynomials in I.

        INPUT:


        -  ``I`` - a list of polynomials or an ideal


        EXAMPLE::

            sage: P.<x,y,z> = QQbar[]
            sage: f1 = -2 * x^2 + x^3
            sage: f2 = -2 * y + x* y
            sage: f3 = -x^2 + y^2
            sage: F = Ideal([f1,f2,f3])
            sage: g = x*y - 3*x*y^2
            sage: g.reduce(F)
            (-6)*y^2 + 2*y
            sage: g.reduce(F.gens())
            (-6)*y^2 + 2*y

        ::

            sage: f = 3*x
            sage: f.reduce([2*x,y])
            0

        ::

            sage: k.<w> = CyclotomicField(3)
            sage: A.<y9,y12,y13,y15> = PolynomialRing(k)
            sage: J = [ y9 + y12]
            sage: f = y9 - y12; f.reduce(J)
            -2*y12
            sage: f = y13*y15; f.reduce(J)
            y13*y15
            sage: f = y13*y15 + y9 - y12; f.reduce(J)
            y13*y15 - 2*y12

        Make sure the remainder returns the correct type, fixing :trac:`13903`::

            sage: R.<y1,y2>=PolynomialRing(Qp(5),2, order='lex')
            sage: G=[y1^2 + y2^2, y1*y2 + y2^2, y2^3]
            sage: type((y2^3).reduce(G))
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>
        """
        from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal

        k = self.base_ring()
        P = self.parent()

        if isinstance(I, MPolynomialIdeal):
            I = I.gens()

        if not k.is_field():
            raise TypeError("Can only reduce polynomials over fields.")

        try:
            fs = self._singular_()
            Is = fs.parent().ideal(I)
            return P(fs.reduce(Is))
        except (NotImplementedError, TypeError):
            pass

        lI = len(I)
        I = list(I)
        r = P(0)
        p = self

        while p != 0:
            for i in xrange(lI):
                gi = I[i]
                plm = p.lm()
                gilm = gi.lm()
                if P.monomial_divides(gilm, plm):
                    quot = p.lc()/gi.lc() * P.monomial_quotient(plm, gilm)
                    p -= quot*I[i]
                    break
            else:
                plt = p.lt()
                r += plt
                p -= plt
        return r

###############################################################
# Useful for some geometry code.
###############################################################

def degree_lowest_rational_function(r,x):
    r"""
    INPUT:


    -  ``r`` - a multivariate rational function

    -  ``x`` - a multivariate polynomial ring generator x


    OUTPUT:


    -  ``integer`` - the degree of r in x and its "leading"
       (in the x-adic sense) coefficient.


    .. note::

       This function is dependent on the ordering of a python dict.
       Thus, it isn't really mathematically well-defined. I think that
       it should made a method of the FractionFieldElement class and
       rewritten.

    EXAMPLES::

        sage: R1 = PolynomialRing(FiniteField(5), 3, names = ["a","b","c"])
        sage: F = FractionField(R1)
        sage: a,b,c = R1.gens()
        sage: f = 3*a*b^2*c^3+4*a*b*c
        sage: g = a^2*b*c^2+2*a^2*b^4*c^7

    Consider the quotient
    `f/g = \frac{4 + 3 bc^{2}}{ac + 2 ab^{3}c^{6}}` (note the
    cancellation).

    ::

        sage: r = f/g; r
        (-2*b*c^2 - 1)/(2*a*b^3*c^6 + a*c)
        sage: degree_lowest_rational_function(r,a)
        (-1, 3)
        sage: degree_lowest_rational_function(r,b)
        (0, 4)
        sage: degree_lowest_rational_function(r,c)
        (-1, 4)
    """
    from sage.rings.fraction_field import FractionField
    R = r.parent()
    F = FractionField(R)
    r = F(r)
    if r == 0:
        return (0, F(0))
    L = x.dict().keys()[0]
    for ix in range(len(L)):
        if L[ix] != 0:
            break
    f = r.numerator()
    g = r.denominator()
    M = f.dict()
    numtermsf = len(M)
    degreesf = [M.keys()[j][ix] for j in range(numtermsf)]
    lowdegf = min(degreesf)
    cf = M[M.keys()[degreesf.index(lowdegf)]] ## constant coeff of lowest degree term
    M = g.dict()
    numtermsg = len(M)
    degreesg = [M.keys()[j][ix] for j in range(numtermsg)]
    lowdegg = min(degreesg)
    cg = M[M.keys()[degreesg.index(lowdegg)]] ## constant coeff of lowest degree term
    return (lowdegf-lowdegg,cf/cg)

