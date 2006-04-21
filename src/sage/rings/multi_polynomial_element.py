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
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import arith

from sage.structure.element import CommutativeRingElement, Element_cmp_, Element
from coerce import bin_op, cmp as coerce_cmp

from sage.interfaces.all import singular, macaulay2

import sage.misc.misc as misc
import integer

import polydict

from sage.structure.factorization import Factorization

from sage.rings.polynomial_singular_interface import Polynomial_singular_repr

import multi_polynomial_ring
import polynomial_ring

def is_MPolynomialRingElement(x):
    return isinstance(x, MPolynomial)

class MPolynomial(Element_cmp_, CommutativeRingElement):
    def __init__(self, parent, x):
        CommutativeRingElement.__init__(self, parent)
        self.__element = x

    def _repr_(self):
        return "%s"%self.__element

    def __call__(self, *x):
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

            sage: x = MPolynomialRing(RationalField(),3).gens()
            sage: f = x[0] + x[1] - 2*x[1]*x[2]
            sage: f
            x1 - 2*x1*x2 + x0
            sage: f(1,2,0)
            3
            sage: f(1,2,5)
            -17

        AUTHOR: David Kohel, 2005-09-27
        """
        if len(x) == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
        n = self.parent().ngens()
        if len(x) != n:
            raise TypeError, "x (=%s) must be of length %s"%(x, n)
        if n == 0:
            return self
        try:
            K = x[0].parent()
        except AttributeError:
            K = self.parent().base_ring()
        y = K(0)
        for (m,c) in self.element().dict().iteritems():
            y += c*misc.mul([ x[i]**m[i] for i in range(n) ])
        return y

    def _cmp_(self, right):
        # MAJOR todo -- this should be relative to the term order, but isn't right now.
        return self.__element.__cmp__(right.__element)

    def _im_gens_(self, codomain, im_gens):
        """
        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = R.hom([y,x], R)
            sage: f(x^2 + 3*y^5)
            y^2 + 3*x^5
        """
        n = self.parent().ngens()
        if n == 0:
            return codomain._coerce_(self)
        y = codomain(0)
        for (m,c) in self.element().dict().iteritems():
            y += codomain(c)*misc.mul([ im_gens[i]**m[i] for i in range(n) ])
        return y


    def _add_(self, right):
        return self.parent()(self.__element + right.__element)

    def _sub_(self, right):
        return self.parent()(self.__element - right.__element)

    def _mul_(self, right):
        return self.parent()(self.__element * right.__element)

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

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent.
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: f = x^3 + y
            sage: g = R(3)
            sage: h = f/g; h
            1/3*y + 1/3*x^3
            sage: h.parent()
            Fraction Field of Polynomial Ring in x, y over Rational Field
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
        if not isinstance(n, (int, long, integer.Integer)):
            raise TypeError, "The exponent (n=%s) must be an integer."%n
        if n < 0:
            return 1/(self**(-n))
        return self.parent()(self.__element**n)

    def __rpow__(self, n):
        if not isinstance(n, (int, long, integer.Integer)):
            raise TypeError, "The exponent (n=%s) must be an integer."%n
        return self.parent()(self.__element**n)

    def element(self):
        return self.__element



class MPolynomial_polydict(Polynomial_singular_repr,MPolynomial):
    def __init__(self, parent, x):
        """
        EXAMPLES:
            sage: R, x = MPolynomialRing(QQ, 10).objgens()
            sage: x
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
            sage: loads(dumps(x)) == x
            True
        """
        if not isinstance(x, polydict.PolyDict):
            x = polydict.PolyDict(x, parent.base_ring()(0), remove_zero=True)
        MPolynomial.__init__(self, parent, x)

    def __neg__(self):
        return self*(-1)

    def _repr_(self):
        return self.element().poly_repr(self.parent().variable_names(),
                                        atomic_coefficients=self.parent().base_ring().is_atomic_repr())

    def _latex_(self):
        return self.element().latex(self.parent().latex_variable_names(),
                                    atomic_coefficients=self.parent().base_ring().is_atomic_repr())

    def _repr_with_changed_varnames(self, varnames):
        return self.element().poly_repr(varnames,
                                        atomic_coefficients=self.parent().base_ring().is_atomic_repr())


    def degree(self, x=None):
        """
        Return the degree of self in x, where x must be one of the
        generators for the parent of self.

        INPUT:
            x -- multivariate polynmial (a generator of the parent of self)
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
            raise TypeError, "x (=%s) must be one of the generators of the parent."%x
        return self.element().degree(x.element())

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
            raise TypeError, "mon (=%s) must be a monomial in the parent of self."%mon
        R = self.parent().base_ring()
        return R(self.element().monomial_coefficient(mon.element().dict()))

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

            sage: R.<x> = MPolynomialRing(GF(7)); R
            Polynomial Ring in x over Finite Field of size 7
            sage: f = 5*x^2 + 3; f
            3 + 5*x^2
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
            sage: f = y^2 - x^9 - 7*x + 5*x*y
            sage: f.coefficient(y)
            5*x
            sage: f = y - x^9*y - 7*x + 5*x*y
            sage: f.coefficient(y)
            1 + 5*x - x^9
        """
        if mon == 1:
            mon = self.parent().gen(0)**0
        if not (isinstance(mon, MPolynomial) and mon.parent() == self.parent() and mon.is_monomial()):
            raise TypeError, "mon (=%s) must be a monomial in the parent of self."%mon
        R = self.parent()
        return R(self.element().coefficient(mon.element().dict()))

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES:
           sage: R.<a,b,c> = PolynomialRing(QQ, 3)
           sage: f = a^3 + b + 2*b^2
           sage: f.exponents()
           [(0, 2, 0), (3, 0, 0), (0, 1, 0)]
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
        if k != tuple([0]*self.parent().ngens()):
            return False
        return bool(d[k].is_unit())

    def inverse_of_unit(self):
        d = self.element().dict()
        k = d.keys()
        if len(k) != 1:
            raise ArithmeticError, "is not a unit"
        k = k[0]
        if k != tuple([0]*self.parent().ngens()):
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
            z^11 + y*z^10 + 5*x*y^10 + x^2*z^9
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


    ############################################################################
    # Some functions added by Martin Albrecht <malb@informatik.uni-bremen.de>
    # (and documented by W. Stein)
    ############################################################################

    def fix(self, fixed):
        """
        Fixes some given variables in a given multivariate polynomial and
        returns the changed multivariate polynomials. The polynomial
        itself is not affected.  The variable,value pairs for fixing are
        to be provided as dictionary of the form {variable:value}.

        This is a special case of evaluating the polynomial with some of
        the variables constants and the others the original variables, but
        should be much faster.

        INPUT:
            fixed -- dict with variable:value pairs

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            30 + y + 25*y^2
            sage: f.fix({x:5})
            30 + y + 25*y^2
        """
        variables = list(self.parent().gens())
        for i in range(0,len(variables)):
            if fixed.has_key(variables[i]):
                variables[i] = fixed[variables[i]]
        return self(tuple(variables))


    def monomials(self):
        """
        Returns a list of all monomials which occure in this
        multivariate polynomial.

        OUTPUT:
            list of MPolynomials representing Monomials

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.monomials()
            [y, x^2, 1, x^2*y^2]
        """
        x = self.parent().gens()
        monomials = []
        for m in self.element().dict().keys():
            monomials.append(misc.mul([ x[i]**m[i] for i in range(self.parent().ngens()) ]))
        return monomials

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
        v = (0,)*int(self.parent().ngens())
        d = self.element().dict()
        try:
            return d[v]
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
            sage: g = f.fix({x:10}); g
            305 - 2*y + 700*y^2
            sage: g.is_univariate()
            True
            sage: f = x^0
            sage: f.is_univariate()
            True
        """
        degrees = self._MPolynomial__element.dict().keys()
        try:
            ngens = len(degrees[0]) # number of generators
        except:
            return True        # zero
        nmons = len(degrees) # number of monomials
        from sage.matrix.all import MatrixSpace
        from sage.rings.all import Z
        degrees = MatrixSpace(Z,nmons,ngens)(list(sum(degrees, ()))).transpose()
        all_one = MatrixSpace(Z,nmons,1)([1]*nmons)

        found = 0
        # use matrix multiplication to add all corresponding variables occurences
        #
        # keys of dictionary are exponent tuples so the keys() list is:
        #
        # [(2,0,0), |
        #  (1,1,0), | addition done by transposing and multiplication with (1,1,...)^T
        #  (0,1,0)]\/
        # ---------
        #   3,2,0  if more than one is not zero -> not univariate
        for elem in (degrees*all_one).list():
            if(elem!=0):
                if(found!=0):
                    return False
                else:
                    found = elem
        return True


    def univariate_polynomial(self, R=None):
        """
        Returns a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:
            R -- (defualt: None) PolynomialRing

        If this polynomial is not in at most one variable, then a
        ValueError exception is raised.  This is checked using the
        is_univariate() method.  The new Polynomial is over the same
        base ring as the given MPolynomial and in the variable 'x' if
        no ring 'ring' is provided.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: polynomial (=5 - 2*y + 3*x^2 + 7*x^2*y^2) must involve at most one variable
            sage: g = f.fix({x:10}); g
            305 - 2*y + 700*y^2
            sage: g.univariate_polynomial ()
            700*x^2 - 2*x + 305
            sage: g.univariate_polynomial(PolynomialRing(QQ,'z'))
            700*z^2 - 2*z + 305
            sage: R = PolynomialRing(QQ,'w')
            sage: R(g)
            700*w^2 - 2*w + 305
        """
        if not self.is_univariate():
            raise ValueError, "polynomial (=%s) must involve at most one variable"%self

        #construct ring if none
        if R == None:
            R =  polynomial_ring.PolynomialRing(self.base_ring(),'x')

        monomial_coefficients = self._MPolynomial__element.dict()

        if( not self.is_constant() ):
            var_idx = self._variable_indices_()[0] #variable
        else:
            var_idx = 0; #constant
            if( len(monomial_coefficients.keys())==0 ):
                return 0

        #construct list
        lookup = [int(0),]*len( monomial_coefficients.keys()[0] )
        coefficients = []
        for degree in range( 0 , max([ m[var_idx] for m in monomial_coefficients.keys() ])+1 ):
            lookup[var_idx]=int(degree);
            try:
                coefficients.append( monomial_coefficients[ tuple(lookup) ] ) #if we find something, add the coefficient
            except KeyError:
                coefficients.append( 0 ) #else add zero

        #construct polynomial
        return R(coefficients)

    def _variable_indices_(self):
        m_coefficients = self._MPolynomial__element.dict()
        variable_dict = dict()

        #get variables
        for exponents in m_coefficients.keys():
            for idx in range( 0 , len(exponents) ):
                if( exponents[idx] != 0 ):
                    variable_dict[idx] = 1

        return variable_dict.keys()

    def variables(self):
        """
        Returns the list of variables occuring in this polynomial.

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.variables()
            [x, y]
            sage: g = f.fix({x:10}); g
            305 - 2*y + 700*y^2
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
            sage: g = f.fix({x:10}); g
            305 - 2*y + 700*y^2
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


    ############################################################################
    # END: Some functions added by Martin Albrecht <malb@informatik.uni-bremen.de>
    ############################################################################

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
            sage: x, y = PolynomialRing(QQ, 2, ['x','y']).gens()
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            2*x*y^2 + 2*x^2*y^2 + x^3 + 2*x^3*y^2 + x^4 + x^5
            sage: F = f.factor()
            sage: F
            x * (2*y^2 + x^2) * (1 + x + x^2)

        Next we factor the same polynomial, but over the finite field
        of order $3$.

            sage: x, y = PolynomialRing(GF(3), 2, ['x','y']).gens()
            sage: f = (x^3 + 2*y^2*x) * (x^2 + x + 1); f
            2*x*y^2 + 2*x^2*y^2 + x^3 + 2*x^3*y^2 + x^4 + x^5
            sage: F = f.factor()
            sage: F
            2 * x * (2 + x)^2 * (y + x) * (y + 2*x)

        \note{Singular multi-variate polynomial factorization is very
        slow in \SAGE.  This \emph{not} a fault of Singular but of how
        the \SAGE NTL is built.  If you download and install a
        Singular binary from the Singular website it will not have
        this problem (you can use it with \SAGE by putting it in
        local/bin/).}
        """
        R = self.parent()
        S = self._singular_().factorize()
        factors = S[1]
        exponents = S[2]
        v = [(R(factors[i+1]), integer.Integer(exponents[i+1])) \
                        for i in range(len(factors))]
        v.sort()
        for i in range(len(v)):
            if str(v[i][0]) == '1':
                del v[i]
                break
        F = Factorization(v)
        F.sort()
        return F

    def gcd(self, f):
        """
        Compute the greatest common divisor of this polynomial and f.

        ALGORITHM: Use Singular.

        EXAMPLES:
            sage: x, y = PolynomialRing(RationalField(), 2, ['x','y']).gens()
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2
        """
        if not isinstance(f, MPolynomial) and self.parent() is f.parent():
            raise TypeError, "self (=%s) and f (=%s) must have the same parent"%(self, f)
        return self.parent()(self._singular_().gcd(f._singular_()))

    def quo_rem(self, right):
        """
        Returns quotient and remainder of self and right.

        ALGORITHM: Use Singular.
        """
        if not isinstance(right, MPolynomial) or right.parent() != self.parent():
            right = self.parent()(right)
        R = self.parent()
        X = self._singular_().division(right._singular_())
        return R(X[1][1,1]), R(X[2][1])


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
        (4 + 3*b*c^2)/(a*c + 2*a*b^3*c^6)
        sage: degree_lowest_rational_function(r,a)
              (-1, 4)
        sage: degree_lowest_rational_function(r,b)
              (0, 4)
        sage: degree_lowest_rational_function(r,c)
              (-1, 4)
    """
    from fraction_field import FractionField
    R = r.parent()
    F = FractionField(R)
    r = F(r)
    if r == 0:
        return (0, F(0))
    L = x.element().dict().keys()[0]
    for ix in range(len(L)):
        if L[ix] != 0:
            break
    f = r.numerator()
    g = r.denominator()
    M = f.element().dict()
    numtermsf = len(M)
    degreesf = [M.keys()[j][ix] for j in range(numtermsf)]
    lowdegf = min(degreesf)
    cf = M[M.keys()[degreesf.index(lowdegf)]] ## constant coeff of lowest degree term
    M = g.element().dict()
    numtermsg = len(M)
    degreesg = [M.keys()[j][ix] for j in range(numtermsg)]
    lowdegg = min(degreesg)
    cg = M[M.keys()[degreesg.index(lowdegg)]] ## constant coeff of lowest degree term
    return (lowdegf-lowdegg,cf/cg)



################################################
class MPolynomial_macaulay2_repr(MPolynomial_polydict):
    """
    Multivariate polynomials that are representable in Macaulay2.
    """
    def _macaulay2_(self, macaulay2=macaulay2):
        """
        Return corresponding Macaulay2 polynomial.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(GF(7), 2, macaulay2=True)   # optional
            sage: f = (x^3 + 2*y^2*x)^7; f          # optional
            2*x^7*y^14 + x^21
            sage: h = f._macaulay2_(); h            # optional
            x^21+2*x^7*y^14
            sage: R(h)                              # optional
            2*x^7*y^14 + x^21
            sage: R(h^20) == f^20                   # optional
            True
        """
        try:
            if self.__macaulay2.parent() is macaulay2:
                return self.__macaulay2
        except AttributeError:
            pass
        self.parent()._macaulay2_(macaulay2)
        self.__macaulay2 = macaulay2(str(self))
        return self.__macaulay2

