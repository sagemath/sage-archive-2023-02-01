"""
Multivariate Polynomials

AUTHORS:
    -- David Joyner: first version
    -- William Stein: use dict's instead of lists
    -- Martin Albrecht <malb@informatik.uni-bremen.de>: some functions added
    -- William Stein (2006-02-11): added better __div__ behavior.
    -- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of
              some Singular features
    -- William Stein (2006-04-19): added e.g., \code{f[1,3]} to get coeff of $xy^3$;
              added examples of the new \code{R.<x,y> = PolynomialRing(QQ,2) notation}.
    -- Martin Albrecht: improved singular coercions (restructed class hierarchy) and added
                        ETuples

EXAMPLES:
We verify Lagrange's four squares identity:
    sage: R.<a0,a1,a2,a3,b0,b1,b2,b3> = ZZ[]
    sage: (a0^2 + a1^2 + a2^2 + a3^2)*(b0^2 + b1^2 + b2^2 + b3^2) == (a0*b0 - a1*b1 - a2*b2 - a3*b3)^2 + (a0*b1 + a1*b0 + a2*b3 - a3*b2)^2 + (a0*b2 - a1*b3 + a2*b0 + a3*b1)^2 + (a0*b3 + a1*b2 - a2*b1 + a3*b0)^2
    True
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
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

import sage.rings.arith

from sage.structure.element import CommutativeRingElement, Element, is_Element

from sage.interfaces.all import singular, macaulay2

import sage.misc.misc as misc
import sage.rings.integer

import polydict

from sage.structure.factorization import Factorization

from sage.rings.polynomial.polynomial_singular_interface import Polynomial_singular_repr

from sage.structure.sequence import Sequence

import multi_polynomial_ring
import polynomial_ring

from sage.rings.integer_ring import ZZ

from multi_polynomial import MPolynomial

def is_MPolynomial(x):
    return isinstance(x, MPolynomial)

class MPolynomial_element(MPolynomial):
    def __init__(self, parent, x):
        CommutativeRingElement.__init__(self, parent)
        self.__element = x

    def _repr_(self):
        return "%s"%self.__element



    ####################

    def __call__(self, *x, **kwds):
        """
        Evaluate this multi-variate polynomial at $x$, where $x$ is
        either the tuple of values to substitute in, or one can use
        functional notation $f(a_0,a_1,a_2, \ldots)$ to evaluate $f$
        with the ith variable replaced by $a_i$.

        EXAMPLES:
            sage: R.<x, y> = MPolynomialRing(RationalField(),2)
            sage: f = x^2 + y^2
            sage: f(1,2)
            5
            sage: f((1,2))
            5

            sage: x = MPolynomialRing(RationalField(),'x',3).gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            -2*x1*x2 + x0 + x1
            sage: f(1,2,0)
            3
            sage: f(1,2,5)
            -17

        AUTHOR: David Kohel, 2005-09-27
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
            raise TypeError, "x must be of correct length"
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

        EXAMPLES:
             sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
             sage: x^1*y^2 > y^3*z^4
             True
             sage: x^3*y^2*z^4 < x^3*y^2*z^1
             False

             sage: R.<x,y,z>=PolynomialRing(QQ,3,order='deglex')
             sage: x^1*y^2*z^3 > x^3*y^2*z^0
             True
             sage: x^1*y^2*z^4 < x^1*y^1*z^5
             False

             sage: R.<x,y,z>=PolynomialRing(ZZ,3,order='degrevlex')
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
        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            3*x^5 + y^2
        """
        n = self.parent().ngens()
        if n == 0:
            return codomain._coerce_(self)
        y = codomain(0)
        for (m,c) in self.element().dict().iteritems():
            y += codomain(c)*misc.mul([ im_gens[i]**m[i] for i in range(n) ])
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

    def __div__(self, right):
        r"""
        EXAMPLES:
            sage: R.<x,y> = QQ['x,y']
            sage: f = (x + y)/3
            sage: f.parent()
            Polynomial Ring in x, y over Rational Field

        If we do the same over $\ZZ$ the result has to lie
        in the fraction field.

            sage: x,y = ZZ['x,y'].gens()
            sage: f = (x + y)/3
            sage: f.parent()
            Fraction Field of Polynomial Ring in x, y over Integer Ring

        """
        try:
            if not isinstance(right, Element) or right.parent() != self.parent():
                R = self.base_ring()
                x = R(right)
                return ~x * self
        except (TypeError, ValueError, ZeroDivisionError):
            pass
        return CommutativeRingElement.__div__(self, right)

    def _div_(self, right):
        return self.parent().fraction_field()(self.__element, right.__element)

    def __pow__(self, n):
        if not isinstance(n, (int, long, sage.rings.integer.Integer)):
            n = sage.rings.integer.Integer(n)
        if n < 0:
            return 1/(self**(-n))
        return self.parent()(self.__element**n)

    def __rpow__(self, n):
        if not isinstance(n, (int, long, sage.rings.integer.Integer)):
            raise TypeError, "The exponent must be an integer."
        return self.parent()(self.__element**n)

    def element(self):
        return self.__element

    def change_ring(self, R):
        return self.parent().change_ring(R)(self)


class MPolynomial_macaulay2_repr:
    """
    Multivariate polynomials that are representable in Macaulay2.
    """
    def _macaulay2_(self, macaulay2=macaulay2):
        """
        Return corresponding Macaulay2 polynomial.

        EXAMPLES:
            sage: R.<x,y> = GF(7)[]
            sage: f = (x^3 + 2*y^2*x)^7; f
            x^21 + 2*x^7*y^14
            sage: macaulay2(R)                      # optional
            ZZ/7 [x, y, MonomialOrder => GRevLex, MonomialSize => 16]
            sage: h = f._macaulay2_(); print h      # optional
             21     7 14
            x   + 2x y
            sage: R(h)                              # optional
            x^21 + 2*x^7*y^14
            sage: R(h^20) == f^20                   # optional
            True
        """
        try:
            if self.__macaulay2.parent() is macaulay2:
                return self.__macaulay2
        except AttributeError:
            pass
        self.parent()._macaulay2_set_ring(macaulay2)
        self.__macaulay2 = macaulay2(repr(self))
        return self.__macaulay2


class MPolynomial_polydict(Polynomial_singular_repr, MPolynomial_macaulay2_repr, MPolynomial_element):
    def __init__(self, parent, x):
        """
        EXAMPLES:
            sage: R, x = MPolynomialRing(QQ, 'x', 10).objgens()
            sage: x
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
            sage: loads(dumps(x)) == x
            True
        """
        if not isinstance(x, polydict.PolyDict):
            x = polydict.PolyDict(x, parent.base_ring()(0), remove_zero=True)
        MPolynomial_element.__init__(self, parent, x)

    def __neg__(self):
        return self*(-1)

    def _repr_(self):
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None

        return self.element().poly_repr(self.parent().variable_names(),
                                        atomic_coefficients=self.parent().base_ring().is_atomic_repr(),cmpfn=cmpfn )

    def _latex_(self):
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None

        return self.element().latex(self.parent().latex_variable_names(),
                                    atomic_coefficients=self.parent().base_ring().is_atomic_repr(), cmpfn=cmpfn)

    def _repr_with_changed_varnames(self, varnames):
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None

        return self.element().poly_repr(varnames,
                                        atomic_coefficients=self.parent().base_ring().is_atomic_repr(), cmpfn=cmpfn)


    def degree(self, x=None):
        """
        Return the degree of self in x, where x must be one of the
        generators for the parent of self.

        INPUT:
            x -- multivariate polynomial (a generator of the parent of self)
                 If x is not specified (or is None), return the total degree,
                 which is the maximum degree of any monomial.

        OUTPUT:
            integer

        EXAMPLE:
            sage: R.<x, y> = MPolynomialRing(QQ, 2)
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
        if not (isinstance(x, MPolynomial) and x.parent() == self.parent() and x.is_monomial()):
            raise TypeError, "x must be one of the generators of the parent."
        return self.element().degree(x.element())

    def newton_polytope(self):
        """
        Return the Newton polytope of this polynomial.

        You should have the optional polymake package installed.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: f = 1 + x*y + x^3 + y^3
            sage: P = f.newton_polytope()
            sage: P
            Convex hull of points [[1, 0, 0], [1, 0, 3], [1, 1, 1], [1, 3, 0]]
            sage: P.facets()
            [(0, 1, 0), (3, -1, -1), (0, 0, 1)]
            sage: P.is_simple()
            True
        """
        try:
            return self.__newton_polytope
        except AttributeError:
            from sage.geometry.all import polymake
            e = self.exponents()
            a = [[1] + list(v) for v in e]
            P = polymake.convex_hull(a)
            self.__newton_polytope = P
            return P

    def total_degree(self):
        """
        Return the total degree of self, which is the
        maximum degree of any monomial in self.

        EXAMPLES:
            sage: R.<x,y,z> = MPolynomialRing(QQ, 3)
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
        Return the coefficient of the monomial mon in self, where mon
        must have the same parent as self.

        INPUT:
            mon -- a monomial

        OUTPUT:
            ring element

        EXAMPLE:
            sage: x, y = MPolynomialRing(RationalField(), 2, names = ['x','y']).gens()

        The coefficient returned is an element of the base ring of self; in
        this case, QQ.
            sage: f = 2 * x * y
            sage: c = f.monomial_coefficient(x*y); c
            2
            sage: c in QQ
            True

            sage: f = y^2 - x^9 - 7*x + 5*x*y
            sage: f.monomial_coefficient(y^2)
            1
            sage: f.monomial_coefficient(x*y)
            5
            sage: f.monomial_coefficient(x^9)
            -1
            sage: f.monomial_coefficient(x^10)
            0
        """
        if not (isinstance(mon, MPolynomial) and mon.parent() == self.parent() and mon.is_monomial()):
            raise TypeError, "mon must be a monomial in the parent of self."
        R = self.parent().base_ring()
        return R(self.element().monomial_coefficient(mon.element().dict()))

    def dict(self):
        """
        Return underlying dictioniary with keys the exponents and
        values the coefficients of this polynomial.
        """
        return self.element().dict()


    def __getitem__(self, x):
        """
        INPUT:
            x -- a tuple or, in case of a single-variable MPolynomial
                 ring x can also be an integer.

        EXAMPLES:
            sage: R.<x, y> = PolynomialRing(QQ, 2)
            sage: f = -10*x^3*y + 17*x*y
            sage: f[3,1]
            -10
            sage: f[1,1]
            17
            sage: f[0,1]
            0

            sage: R.<x> = PolynomialRing(GF(7),1); R
            Polynomial Ring in x over Finite Field of size 7
            sage: f = 5*x^2 + 3; f
            -2*x^2 + 3
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


    def coefficient(self, mon):
        """
        Return the coefficient of mon in self, where mon must have the
        same parent as self.  The coefficient is defined as follows.
        If f is this polynomial, then the coefficient is the sum T/mon
        where the sum is over terms T in f that are exactly divisible
        by mon.

        INPUT:
            mon -- a monomial

        OUTPUT:
            element of the parent of self

        EXAMPLE:
            sage: x, y = MPolynomialRing(RationalField(), 2, names = ['x','y']).gens()

        The coefficient returned is an element of the parent of self; in
        this case, QQ[x, y].
            sage: f = 2 * x * y
            sage: c = f.coefficient(x*y); c
            2
            sage: c.parent()
            Polynomial Ring in x, y over Rational Field
            sage: c in MPolynomialRing(RationalField(), 2, names = ['x','y'])
            True

            sage: f = y^2 - x^9 - 7*x + 5*x*y
            sage: f.coefficient(y)
            5*x
            sage: f = y - x^9*y - 7*x + 5*x*y
            sage: f.coefficient(y)
            -x^9 + 5*x + 1

        The coefficient of 1 is also an element of the multivariate
        polynomial ring:
            sage: R.<x,y> = GF(389)[]
            sage: parent(R(x*y+5).coefficient(R(1)))
            Polynomial Ring in x, y over Finite Field of size 389
        """
        R = self.parent()
        if mon == 1:
            return R(self.constant_coefficient())
        if not (isinstance(mon, MPolynomial) and mon.parent() == self.parent() and mon.is_monomial()):
            raise TypeError, "mon must be a monomial in the parent of self."
        return R(self.element().coefficient(mon.element().dict()))

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES:
           sage: R.<a,b,c> = PolynomialRing(QQ, 3)
           sage: f = a^3 + b + 2*b^2
           sage: f.exponents()
           [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        return self.element().exponents()

    def is_unit(self):
        """
        Return True if self is a unit.

        EXAMPLES:
            sage: R = PolynomialRing(IntegerRing(), 2, ['x','y']); x,y = R.gens()
            sage: (x+y).is_unit()
            False
            sage: R(0).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(-1 + x).is_unit()
            False
            sage: R(2).is_unit()
            False
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
            raise ArithmeticError, "is not a unit"
        k = k[0]
        if k != polydict.ETuple([0]*self.parent().ngens()):
            raise ArithmeticError, "is not a unit"
        return ~d[k]

    def is_homogeneous(self):
        """
        Return True if self is a homogeneous polynomial.

        EXAMPLES:
            sage: x, y = MPolynomialRing(RationalField(), 2, names=['x', 'y']).gens()
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

    def homogenize(self, var="h"):
        """
        Return self is self is homogeneous.  Otherwise return a homogeneous
        polynomial in one more variable such that setting that variable
        equal to 1 yields self.

        INPUT:
            var -- string (default: "h"); a variable name for the new variable
                   to be added in when homogenizing.

        OUTPUT:
            a multivariate polynomial

        EXAMPLES:
            sage: x,y = MPolynomialRing(RationalField(),2,['x','y']).gens()
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: g = f.homogenize('z'); g
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Polynomial Ring in x, y, z over Rational Field
        """
        if self.is_homogeneous():
            return self
        X = self.element().homogenize()
        R = self.parent()
        S = multi_polynomial_ring.MPolynomialRing(
                        R.base_ring(),
                        R.ngens() + 1,
                        names=R.variable_names() + (var,),
                        order = R.term_order())
        return S(X)

    def is_monomial(self):
        return len(self.element().dict().keys()) == 1

    def subs(self, fixed=None, **kw):
        """
        Fixes some given variables in a given multivariate polynomial and
        returns the changed multivariate polynomials. The polynomial
        itself is not affected.  The variable,value pairs for fixing are
        to be provided as a dictionary of the form {variable:value}.

        This is a special case of evaluating the polynomial with some of
        the variables constants and the others the original variables.

        INPUT:
            fixed -- (optional) dictionary of inputs
            **kw  -- named parameters

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
        """
        variables = list(self.parent().gens())
        for i in range(0,len(variables)):
            if kw.has_key(str(variables[i])):
                variables[i]=kw[str(variables[i])]
            elif fixed and fixed.has_key(variables[i]):
                variables[i] = fixed[variables[i]]
        return self(tuple(variables))

    def monomials(self):
        """
        Returns list of all monomials which occure in this
        multivariate polynomial.

        OUTPUT:
            list of MPolynomials representing Monomials

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.monomials()
            [1, x^2*y^2, x^2, y]
        """
        try:
            return self.__monomials
        except AttributeError:
            ring = self.parent()
            one = self.parent().base_ring()(1)
            self.__monomials = [ MPolynomial_polydict(ring, polydict.PolyDict( {m:one}, force_int_exponents=False,  force_etuples=False ) ) \
                                for m in self._MPolynomial_element__element.dict().keys() ]
            return self.__monomials

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate polynomial.

        EXAMPLES:
            sage: x, y = ZZ['x,y'].gens()
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
        Returns True if this multivariate polynomial is univariate and False otherwise.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.is_univariate()
            False
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.is_univariate()
            True
            sage: f = x^0
            sage: f.is_univariate()
            True
        """
        mons = self.element().dict().keys()
        try:
            ngens = len(mons[0]) # number of generators
        except:
            return True # zero

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
        Returns a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:
            R -- (default: None) PolynomialRing

        If this polynomial is not in at most one variable, then a
        ValueError exception is raised.  This is checked using the
        is_univariate() method.  The new Polynomial is over the same
        base ring as the given MPolynomial and in the variable 'x' if
        no ring 'ring' is provided.

        EXAMPLES:
            sage: R.<x, y> = MPolynomialRing(ZZ,2,'xy')
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.univariate_polynomial ()
            700*x^2 - 2*x + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305
        """
        if not self.is_univariate():
            raise TypeError, "polynomial must involve at most one variable"

        #construct ring if none
        if R == None:
            R =  polynomial_ring.PolynomialRing(self.base_ring(),'x')

        monomial_coefficients = self._MPolynomial_element__element.dict()

        if( not self.is_constant() ):
            var_idx = self._variable_indices_()[0] #variable
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

    def _variable_indices_(self):

        ETuples = self._MPolynomial_element__element.dict().keys()

        idx = set()
        for e in ETuples:
            idx = idx.union(e.nonzero_positions())
        return sorted(idx)

    def variables(self):
        """
        Returns the list of variables occuring in this polynomial.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.variables()
            [x, y]
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.variables()
            [y]
        """
        return [self.parent().gen(index) for index in self._variable_indices_() ]


    def variable(self,i):
        """
        Returns $i$-th variable occuring in this polynomial.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
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

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ, 2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.nvariables ()
            2
            sage: g = f.subs({x:10}); g
            700*y^2 - 2*y + 305
            sage: g.nvariables ()
            1
        """
        return len(self._variable_indices_())

    def is_constant(self):
        """
        True if polynomial is constant, and False otherwise.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.is_constant()
            False
            sage: g = 10*x^0
            sage: g.is_constant()
            True
        """
        if( len(self._variable_indices_()) == 0 ):
            return True
        else:
            return False

    def __hash__(self):
        #requires base field elements are hashable!
        return hash(tuple(self._MPolynomial_element__element.dict().items()))

    def lm(self):
        """
        Returns the lead monomial of self with respect to the term order of
        self.parent(). Where 'lex', 'deglex', 'revlex', and 'degrevlex' are
        accepted.

        We say $a >_{lex} b$ if, in the vector difference $a-b \in Z^n$,
        the left-most nonzero entry is positive.

        We say $a >_{revlex} b$ if, in the vector difference $a-b \in Z^n$,
        the right-most nonzero entry is positive.

        We say $a >_{deglex} b$ if, $|a| > |b|$, or $|a| = |b|$ and the
        left-most nonzero entry of $a -b \in \ZZ^n$ is positive.

        We say $a >_{degrevlex} b$ if, $|a| > |b|$, or $|a| = |b|$ and the
        right-most nonzero entry of $a -b \in Z^n$ is negative.

        EXAMPLES:
             sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
             sage: (x^1*y^2 + y^3*z^4).lm()
             x*y^2
             sage: (x^3*y^2*z^4 + x^3*y^2*z^1).lm()
             x^3*y^2*z^4

             sage: R.<x,y,z>=PolynomialRing(QQ,3,order='deglex')
             sage: (x^1*y^2*z^3 + x^3*y^2*z^0).lm()
             x*y^2*z^3
             sage: (x^1*y^2*z^4 + x^1*y^1*z^5).lm()
             x*y^2*z^4

             sage: R.<x,y,z>=PolynomialRing(ZZ,3,order='degrevlex')
             sage: (x^1*y^5*z^2 + x^4*y^1*z^3).lm()
             x*y^5*z^2
             sage: (x^4*y^7*z^1 + x^4*y^2*z^3).lm()
             x^4*y^7*z

        """
        try:
            return self.__lm
        except AttributeError:
            if self.is_zero():
                return self
            R = self.parent()
            f = self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )
            one = R.base_ring()(1)
            self.__lm = MPolynomial_polydict(R,polydict.PolyDict({f:one},force_int_exponents=False,  force_etuples=False))
            return self.__lm

    def lc(self):
        """
        Returns the leading coefficent of self i.e.,
        self.coefficient(self.lm())
        """
        try:
            return self.__lc
        except AttributeError:
            if self.is_zero():
                return self
            R = self.parent()
            f = self._MPolynomial_element__element.dict()
            self.__lc = f[self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )]
            return self.__lc

    def lt(self):
        """
        Returns the leading term of self i.e., self.lc()*self.lm()
        """
        try:
            return self.__lt
        except AttributeError:
            if self.is_zero():
                return self
            R = self.parent()
            f = self._MPolynomial_element__element.dict()
            res = self._MPolynomial_element__element.lcmt( R.term_order().greater_tuple )
            self.__lt = MPolynomial_polydict(R,polydict.PolyDict({res:f[res]},force_int_exponents=False, force_etuples=False))
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

        \note{This is much faster than actually writing self == 0}
        """
        return self._MPolynomial_element__element.dict()!={}

    def __floordiv__(self,right):
        """
        Quotient of division of self by other.  This is denoted //.
        """
        Q, _ = self.quo_rem(right)
        return Q

    def factor(self):
        r"""
        Compute the irreducible factorization of this polynomial.

        ALGORITHM: Use Singular.

        EXAMPLES:
            sage: R.<x, y> = QQ[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 + 2*x^3*y^2 + x^4 + 2*x^2*y^2 + x^3 + 2*x*y^2
            sage: F = f.factor()
            sage: F
            x * (x^2 + x + 1) * (x^2 + 2*y^2)

        Next we factor the same polynomial, but over the finite field
        of order $3$.

            sage: R.<x, y> = GF(3)[]
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            x^5 - x^3*y^2 + x^4 - x^2*y^2 + x^3 - x*y^2
            sage: F = f.factor()
            sage: F # order is somewhat random
            (-1) * x * (-x + y) * (x + y) * (x - 1)^2

        Next we factor a polynomial over a number field.
            sage: p = var('p')
            sage: K.<s> = NumberField(p^3-2)
            sage: KXY.<x,y> = K[]
            sage: factor(x^3 - 2*y^3)
            (x + (-s)*y) * (x^2 + s*x*y + s^2*y^2)
            sage: k = (x^3-2*y^3)^5*(x+s*y)^2*(2/3 + s^2)
            sage: k.factor()
            (s^2 + 2/3) * (x + s*y)^2 * (x + (-s)*y)^5 * (x^2 + s*x*y + s^2*y^2)^5
        """
        # I do not think this applied anymore.  Or at least it's
        # more relevant to optimizing the NTL build.
        #\note{Singular multi-variate polynomial factorization is very
        #slow in \SAGE.  This \emph{not} a fault of Singular but of how
        #the \SAGE NTL is built.  If you download and install a
        #Singular binary from the Singular website it will not have
        #this problem (you can use it with \SAGE by putting it in
        #local/bin/).}
        R = self.parent()
        R._singular_().set_ring()
        S = self._singular_().factorize()
        factors = S[1]
        exponents = S[2]
        v = [(R(factors[i+1]), sage.rings.integer.Integer(exponents[i+1])) \
                        for i in range(len(factors))]
        v.sort()
        for i in range(len(v)):
            if str(v[i][0]) == '1':
                del v[i]
                break
        F = Factorization(v)
        F.sort()
        return F

    def lift(self,I):
        """
        given an ideal I = (f_1,...,f_r) and some g (== self) in I,
        find s_1,...,s_r such that g = s_1 f_1 + ... + s_r f_r

        ALGORITHM: Use Singular.

        EXAMPLE:
            sage: A.<x,y> = PolynomialRing(QQ,2,order='degrevlex')
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
            raise ArithmeticError, "f is not in I"
        return Sequence(M.list(), P, check=False, immutable=True)

    def gcd(self, f):
        """
        Compute the greatest common divisor of this polynomial and f.

        ALGORITHM: Use Singular.

        EXAMPLES:
            sage: x, y = QQ['x,y'].gens()
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2

        This also works correctly over ZZ:
            sage: R.<x,y> = ZZ[]
            sage: gcd(2*x,4*x)
            2*x
            sage: gcd(2*x,4*x)
            2*x
            sage: gcd(9*x*y*(x^2-y^2), 15*x*y^2*(x^2+y^2))
            3*x*y

        We compute a gcd over a finite field.
            sage: F.<u> = GF(31^2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
            sage: gcd(p,q)  # yes, twice -- tests that singular ring is properly set.
            x^3 + (u + 1)*y^3 + z^3

        We compute a gcd over a number field:
            sage: x = polygen(QQ)
            sage: F.<u> = NumberField(x^3 - 2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
        """
        if not isinstance(f, MPolynomial) and self.parent() is f.parent():
            raise TypeError, "self and f must have the same parent"


        # Singular ignores coefficents anyway, thus it is okay to work over Z here
        # PARI uses the coefficents btw.
        # TODO: This is slow

        P = self.parent()
        if P.base_ring() == ZZ:
            res = self.parent()(self._singular_(force=True).gcd(f._singular_(force=True)))
            coef = sage.rings.arith.gcd(self.element().dict().values() + f.element().dict().values(),True)
            return coef*res

        P._singular_().set_ring()
        return P(self._singular_().gcd(f._singular_()))

    def quo_rem(self, right):
        """
        Returns quotient and remainder of self and right.

        ALGORITHM: Use Singular.
        """
        if not isinstance(right, MPolynomial) or right.parent() != self.parent():
            right = self.parent()(right)
        R = self.parent()
        R._singular_().set_ring()
        X = self._singular_().division(right._singular_())
        return R(X[1][1,1]), R(X[2][1])

    def _magma_(self, magma=None):
        """
        Returns the MAGMA representation of self.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(GF(2),2)
            sage: f = y*x^2 + x +1
            sage: f._magma_() #optional
            x^2*y + x + 1

        """
        if magma == None:
            import sage.interfaces.magma
            magma = sage.interfaces.magma.magma

        try:
            m = self.__magma
            m._check_valid()
            if not m.parent() is magma:
                raise ValueError
            return m
        except (AttributeError,ValueError):
            magma_gens = [e.name() for e in self.parent()._magma_().gens()]
            f = self.element().poly_repr(magma_gens,atomic_coefficients=False)
            self.__magma = magma(f)
            return self.__magma



###############################################################
# Useful for some geometry code.
###############################################################

def degree_lowest_rational_function(r,x):
    r"""
    INPUT:
        r -- a multivariate rational function
        x -- a multivariate polynomial ring generator x

    OUTPUT:
        integer -- the degree of r in x and its "leading"
                   (in the x-adic sense) coefficient.

    EXAMPLES:
        sage: R1 = MPolynomialRing(FiniteField(5), 3, names = ["a","b","c"])
        sage: F = FractionField(R1)
        sage: a,b,c = R1.gens()
        sage: f = 3*a*b^2*c^3+4*a*b*c
        sage: g = a^2*b*c^2+2*a^2*b^4*c^7

    Consider the quotient $f/g = \frac{4 + 3 bc^{2}}{ac + 2 ab^{3}c^{6}}$ (note
    the cancellation).
        sage: r = f/g; r
        (-2*b*c^2 - 1)/(2*a*b^3*c^6 + a*c)
        sage: degree_lowest_rational_function(r,a)
              (-1, 4)
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

