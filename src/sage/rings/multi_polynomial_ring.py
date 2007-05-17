r"""
Multivariate Polynomial Rings

AUTHORS:
    -- David Joyner and William Stein
    -- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of
              Singular features
    -- Martin Albrecht (2006-04-21): reorganize class hiearchy for singular rep

EXAMPLES:

We construct the Frobenius morphism on $\mbox{\rm F}_{5}[x,y,z]$ over $\F_5$:

    sage: R, (x,y,z) = PolynomialRing(GF(5), 3, 'xyz').objgens()
    sage: frob = R.hom([x^5, y^5, z^5])
    sage: frob(x^2 + 2*y - z^4)
    -z^20 + x^10 + 2*y^5
    sage: frob((x + 2*y)^3)
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15
    sage: (x^5 + 2*y^5)^3
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15

We make a polynomial ring in one variable over a polynomial ring in
two variables:
    sage: R.<x, y> = PolynomialRing(QQ, 2)
    sage: S.<t> = PowerSeriesRing(R)
    sage: t*(x+y)
    (x + y)*t
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

# Changed:
# Kiran Kedlaya (2006-02-12): added Macaulay2 names to TermOrder

import weakref
import re

import commutative_ring
import integral_domain

import fraction_field
import fraction_field_element

import finite_field

import multi_polynomial_element
import multi_polynomial_ideal
import polydict

import sage.misc.latex as latex

from sage.interfaces.all import singular as singular_default, is_SingularElement
from sage.interfaces.all import macaulay2 as macaulay2_default
from sage.interfaces.macaulay2 import is_Macaulay2Element

from sage.structure.sage_object import SageObject

from sage.rings.integer_ring import is_IntegerRing
from sage.rings.integer import Integer

from sage.rings.polynomial_singular_interface import PolynomialRing_singular_repr

import multi_polynomial_ideal

from sage.rings.polynomial_ring_constructor import PolynomialRing as MPolynomialRing

from sage.structure.parent_gens import ParentWithGens

from multi_polynomial_ring_generic import MPolynomialRing_generic, is_MPolynomialRing

from polydict import ETuple

class MPolynomialRing_macaulay2_repr:
    """
    """
    def _macaulay2_(self, macaulay2=None):
        if macaulay2 is None:
            macaulay2 = macaulay2_default
        try:
            R = self.__macaulay2
            if not (R.parent() is macaulay2):
                raise ValueError
            R._check_valid()
            return R
        except (AttributeError, ValueError):
            if self.base_ring().is_prime_field():
                if self.characteristic() == 0:
                    base_str = "QQ"
                else:
                    base_str = "ZZ/" + str(self.characteristic())
            elif is_IntegerRing(self.base_ring()):
                base_str = "ZZ"
            else:
                raise TypeError, "no conversion of to a Macaulay2 ring defined"
            self.__macaulay2 = macaulay2.ring(base_str, str(self.gens()), \
                                              self.term_order().macaulay2_str())
        return self.__macaulay2

    def is_exact(self):
        return self.base_ring().is_exact()

    def change_ring(self, R):
        from polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(R, self.variable_names(), order=self.term_order())

class MPolynomialRing_polydict( MPolynomialRing_macaulay2_repr, MPolynomialRing_generic):
    """
    Multivariable polynomial ring.

    EXAMPLES:
        sage: R = MPolynomialRing(Integers(12), 'x', 5); R
        Polynomial Ring in x0, x1, x2, x3, x4 over Ring of integers modulo 12
        sage.: loads(R.dumps()) == R     # TODO -- this currently hangs sometimes (??)
        True
    """
    def __init__(self, base_ring, n, names, order):
        MPolynomialRing_generic.__init__(self, base_ring, n, names, order)
        # Construct the generators
        v = [0 for _ in xrange(n)]
        one = base_ring(1);
        self._gens = []
        C = self._poly_class()
        for i in xrange(n):
            v[i] = 1  # int's!
            self._gens.append(C(self, {tuple(v):one}))
            v[i] = 0
        self._gens = tuple(self._gens)
        self._zero_tuple = tuple(v)

    def _monomial_order_function(self):
        return self.__monomial_order_function

    def _poly_class(self):
        return multi_polynomial_element.MPolynomial_polydict

    def __cmp__(left, right):
        if not is_MPolynomialRing(right):
            return cmp(type(left),type(right))
        else:
            return cmp((left.base_ring(), left.ngens(), left.variable_names(), left.term_order()),
                       (right.base_ring(), right.ngens(), right.variable_names(), right.term_order()))

    def __call__(self, x, check=True):
        """
        Coerce x into this multivariate polynomial ring.

        EXAMPLES:
        We create a Macaulay2 multivariate polynomial via ideal arithmetic,
        then coerce it into R.
            sage: R.<x,y> = PolynomialRing(QQ, 2)                        # optional
            sage: I = R.ideal([x^3 + y, y])                              # optional
            sage: S = I._macaulay2_()                                    # optional
            sage: T = S*S*S                                              # optional
            sage: U = T.gens().entries().flatten()                       # optional
            sage: f = U[2]; f                                            # optional
             3 2    3
            x y  + y
            sage: R(f)                                                   # optional
            y^3 + x^3*y^2

        Some other subtle coercions.  We create polynomial rings in 2 variables
        over the rationals, integers, and a finite field.
            sage: R.<x,y> = QQ[]
            sage: S.<x,y> = ZZ[]
            sage: T.<x,y> = GF(7)[]

        We coerce from the integer to the rationals, and back:
            sage: f = R(S.0^2 - 4*S.1^3); f
            -4*y^3 + x^2
            sage: parent(f)
            Polynomial Ring in x, y over Rational Field
            sage: parent(S(f))
            Polynomial Ring in x, y over Integer Ring

        We coerce from the finite field.
            sage: f = R(T.0^2 - 4*T.1^3); f
            3*y^3 + x^2
            sage: parent(f)
            Polynomial Ring in x, y over Rational Field

        We dump and load a the polynomial ring S:
            sage: S2 = loads(dumps(S))
            sage: S2 == S
            True

        Coerce works and gets the right parent.
            sage: parent(S2._coerce_(S.0)) is S2
            True

        Coercion to reduce modulo a prime between rings with different variable names:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = PolynomialRing(GF(7),2)
            sage: f = x^2 + 2/3*y^3
            sage: S(f)
            3*b^3 + a^2

        Coercion from symbolic variables:
            sage: x,y,z = var('x,y,z')
            sage: R = QQ[x,y,z]
            sage: type(x)
            <class 'sage.calculus.calculus.SymbolicVariable'>
            sage: type(R(x))
            <type 'sage.rings.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f = R(x^3 + y^3 - z^3); f
            x^3 + y^3 - z^3
            sage: type(f)
            <type 'sage.rings.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: parent(f)
            Polynomial Ring in x, y, z over Rational Field

        A more complicated symbolic and computational mix.  Behind the scenes
        Singular and Maxima are doing the real work.
            sage: R = QQ[x,y,z]
            sage: f = (x^3 + y^3 - z^3)^10; f
            (-z^3 + y^3 + x^3)^10
            sage: g = R(f); parent(g)
            Polynomial Ring in x, y, z over Rational Field
            sage: (f - g).expand()
            0

        """
        from sage.rings.multi_polynomial_libsingular import MPolynomial_libsingular

        if isinstance(x, multi_polynomial_element.MPolynomial_polydict):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return multi_polynomial_element.MPolynomial_polydict(self, x.element().dict())
            elif len(P.variable_names()) == len(self.variable_names()):
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                K = self.base_ring()
                D = x.element().dict()
                for i, a in D.iteritems():
                    D[i] = K(a)
                return multi_polynomial_element.MPolynomial_polydict(self, D)
            else:
                raise TypeError

        if isinstance(x, MPolynomial_libsingular):
            P = x.parent()
            if P == self:
                return multi_polynomial_element.MPolynomial_polydict(self, x.dict())
            elif len(P.variable_names()) == len(self.variable_names()):
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                K = self.base_ring()
                D = x.dict()
                for i, a in D.iteritems():
                    D[i] = K(a)
                return multi_polynomial_element.MPolynomial_polydict(self, D)
            else:
                raise TypeError

        elif isinstance(x, polydict.PolyDict):
            return multi_polynomial_element.MPolynomial_polydict(self, x)

        elif isinstance(x, fraction_field_element.FractionFieldElement) and x.parent().ring() == self:
            if x.denominator() == 1:
                return x.numerator()
            else:
                raise TypeError, "unable to coerce since the denominator is not 1"

        elif is_SingularElement(x) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except TypeError:
                raise TypeError, "unable to coerce singular object"

        elif hasattr(x, '_polynomial_'):
            return x._polynomial_(self)

        elif isinstance(x , str) and self._has_singular:
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except TypeError:
                raise TypeError,"unable to coerce string"

        elif is_Macaulay2Element(x):
            try:
                s = x.sage_polystring()
                if len(s) == 0:
                    raise TypeError
                # NOTE: It's CRUCIAL to use the eval command as follows,
                # i.e., with the gen dict as the third arg and the second
                # empty.  Otherwise pickling won't work after calls to this eval!!!
                # This took a while to figure out!
                return self(eval(s, {}, self.gens_dict()))
            except (AttributeError, TypeError, NameError):
                raise TypeError, "Unable to coerce macaulay2 object"
            return multi_polynomial_element.MPolynomial_polydict(self, x)

        if isinstance(x, dict):
            return multi_polynomial_element.MPolynomial_polydict(self, x)
        else:
            c = self.base_ring()(x)
            return multi_polynomial_element.MPolynomial_polydict(self, {self._zero_tuple:c})



class MPolynomialRing_polydict_domain(integral_domain.IntegralDomain,
                                      MPolynomialRing_polydict,
                                      PolynomialRing_singular_repr,
                                      MPolynomialRing_macaulay2_repr):
    def __init__(self, base_ring, n, names, order):
        MPolynomialRing_polydict.__init__(self, base_ring, n, names, order)
        self._has_singular = self._can_convert_to_singular()

    def is_integral_domain(self):
        return True

    def is_field(self):
        if self.ngens() == 0:
            return self.base_ring().is_field()
        return False

    def ideal(self, gens, coerce=True):
        """
        Create an ideal in this polynomial ring.
        """
        if not self._has_singular:
            # pass through
            MPolynomialRing_generic.ideal(self,gens,coerce)
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        if is_Macaulay2Element(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return multi_polynomial_ideal.MPolynomialIdeal(self, gens, coerce=False)


    def monomial_quotient(self,f, g, coeff=False):
        """
        Return f/g, where both f and g are treated as
        monomials. Coefficients are ignored by default.

        INPUT:
            f -- monomial
            g -- monomial
            coeff -- divide coefficents as well (default: False)

        EXAMPLE:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ, 3, order=TermOrder('degrevlex'))
            sage: P.monomial_quotient(3/2*x*y,x)
            y

            sage: P.monomial_quotient(3/2*x*y,2*x,coeff=True)
            3/4*y

        TESTS:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: R.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_quotient(x*y,x)
            y

            sage: P.monomial_quotient(x*y,R.gen())
            y

            sage: P.monomial_quotient(P(0),P(1))
            0

            sage: P.monomial_quotient(P(1),P(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

            sage: P.monomial_quotient(P(3/2),P(2/3), coeff=True)
            9/4

            sage: P.monomial_quotient(x,y) # Note the wrong result
            x*y^-1

            sage: P.monomial_quotient(x,P(1))
            x

        NOTE: Assumes that the head term of f is a multiple of the
        head term of g and return the multiplicant m. If this rule is
        violated, funny things may happen.

        """

        if not f:
          return f
        if not g:
          raise ZeroDivisionError

        if not coeff:
          coeff= self.base_ring()(1)
        else:
          coeff = f.dict().values()[0] /  g.dict().values()[0]

        f = f.dict().keys()[0]
        g = g.dict().keys()[0]

        res = f.esub(g)

        return multi_polynomial_element.MPolynomial_polydict(self,polydict.PolyDict({res:coeff},\
                                                                                    force_int_exponents=False, \
                                                                                    force_etuples=False))

    def monomial_lcm(self, f, g):
        """
        LCM for monomials. Coefficients are ignored.

        INPUT:
            f -- monomial
            g -- monomial

        EXAMPLE:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_lcm(3/2*x*y,x)
            x*y

        TESTS:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: R.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_lcm(x*y,R.gen())
            x*y

            sage: P.monomial_lcm(P(3/2),P(2/3))
            1

            sage: P.monomial_lcm(x,P(1))
            x

        """
        one = self.base_ring()(1)

        f=f.dict().keys()[0]
        g=g.dict().keys()[0]


        length = len(f)

        res = {}

        nonzero = []

        for i in f.common_nonzero_positions(g):
            res[i] = max([f[i],g[i]])

        res =  self(polydict.PolyDict({ETuple(res,length):one},\
                                   force_int_exponents=False,force_etuples=False))
        return res

    def monomial_reduce(self, f, G):
        """
        Try to find a g in G where g.lm() divides f. If found (g,flt)
        is returned, (0,0) otherwise, where flt is f/g.lm().

        It is assumed that G is iterable and contains ONLY elements in
        self.

        INPUT:
            f -- monomial
            G -- list/set of mpolynomials

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, P(1/2)  ]
            sage: P.monomial_reduce(f,G)
            (1/4*x*y + 2/7, y)

        TESTS:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: f = x*y^2
            sage: G = [ 3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, P(1/2)  ]

            sage: P.monomial_reduce(P(0),G)
            (0, 0)

            sage: P.monomial_reduce(f,[P(0)])
            (0, 0)

        """
        if not f:
            return 0,0
        for g in G:
            t = g.lm()
            if self.monomial_is_divisible_by(f,t):
                return g,self.monomial_quotient(f,t)
        return 0,0


    def monomial_is_divisible_by(self, a, b):
        """
        Return False if b does not divide a and True otherwise.

        INPUT:
            a -- monomial
            b -- monomial

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_is_divisible_by(x^3*y^2*z^4, x*y*z)
            True
            sage: P.monomial_is_divisible_by(x*y*z, x^3*y^2*z^4)
            False

        TESTS:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_is_divisible_by(P(0),P(1))
            True
            sage: P.monomial_is_divisible_by(x,P(1))
            True

        """

        if not a:
            return True
        if not b:
            return False

        one = self.base_ring()(1)

        a=a.dict().keys()[0]
        b=b.dict().keys()[0]

        for i in a.common_nonzero_positions(b):
          if a[i] - b[i] < 0:
            return False
        return True

    def monomial_pairwise_prime(self, h, g):
        """
        Return True if h and g are pairwise prime. Both are treated as monomials.

        INPUT:
            h -- monomial
            g -- monomial

        EXAMPLES:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: Q.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_pairwise_prime(x^2*z^3, Q('y^4'))
            True

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, Q(0))
            True

            sage: P.monomial_pairwise_prime(P(1/2),x)
            False
        """
        if not g:
            if not h:
                return False #GCD(0,0) = 0
            else:
                return True #GCD(x,0) = 1

        elif not h:
            return True # GCD(0,x) = 1

        return self.monomial_lcm(g,h) == g*h

    def monomial_all_divisors(self,t):
        """
        Return a list of all monomials that divide t, coefficients are
        ignored.

        INPUT:
            t -- a monomial

        OUTPUT:
            a list of monomials


        EXAMPLE:
            sage: from sage.rings.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order=TermOrder('degrevlex'))
            sage: P.monomial_all_divisors(x^2*z^3)
            [x, x^2, z, x*z, x^2*z, z^2, x*z^2, x^2*z^2, z^3, x*z^3, x^2*z^3]

        ALGORITHM: addwithcarry idea by Toon Segers
        """

        def addwithcarry(tempvector, maxvector, pos):
            if tempvector[pos] < maxvector[pos]:
              tempvector[pos] += 1
            else:
              tempvector[pos] = 0
              tempvector = addwithcarry(tempvector, maxvector, pos + 1)
            return tempvector

        if not t.is_monomial():
          raise TypeError, "Only monomials are supported"

        R = self
        one = self.base_ring()(1)
        M = list()

        maxvector = list(t.dict().keys()[0])

        tempvector =[0,]*len(maxvector)

        pos = 0

        while tempvector != maxvector:
          tempvector = addwithcarry(list(tempvector) , maxvector, pos)
          M.append(R(polydict.PolyDict({ETuple(tempvector):one}, \
                                       force_int_exponents=False,force_etuples=False)))
        return M



#######################

singular_name_mapping = {'lex':'lp', \
                'revlex':'rp', \
                'degrevlex':'dp', \
                'deglex':'Dp'}

m2_name_mapping = {'lex':'Lex', \
                   'revlex':'RevLex', \
                   'degrevlex':'GRevLex', \
                   'deglex':'GLex'}

magma_name_mapping = {'lex': '"lex"', \
                      'revlex' : '"revlex"', \
                      'deglex' : '"glex"', \
                      'degrevlex' : '"grevlex"'}

class TermOrder(SageObject):
    """
    EXAMPLES:
        sage: t = TermOrder('lex')
        sage: t
        Lexicographic term order
        sage: loads(dumps(t)) == t
        True
    """
    def __init__(self, name='lex'):
        if isinstance(name, TermOrder):
            name = name.__name
        name = name.lower()
        self.__name = name

        if singular_name_mapping.has_key(name):
            singular_name = singular_name_mapping[name]
            self.__singular_str = singular_name
        else:
            self.__singular_str = name

        if m2_name_mapping.has_key(name):
            macaulay2_name = m2_name_mapping[name]
            self.__macaulay2_str = macaulay2_name
        else:
            self.__macaulay2_str = name

        if magma_name_mapping.has_key(name):
            magma_name = magma_name_mapping[name]
            self.__magma_str = magma_name
        else:
            self.__magma_str = name


    def __getattr__(self,name):
        if name=='compare_tuples':
            return getattr(self,'compare_tuples_'+self.__singular_str)
        elif name=='greater_tuple':
            return getattr(self,'greater_tuple_'+self.__singular_str)
        else:
            raise AttributeError,name

    def compare_tuples_lp(self,f,g):
        """
        Compares two exponent tuples with respect to the
        lexicographical term order.
        """

        if f>g:
            return 1
        elif f<g:
            return -1
        else:
            return 0

    def compare_tuples_rp(self,f,g):
        """
        Compares two exponent tuples with respect to the reversed
        lexicographical term order.
        """
        return (-1)*self.compare_tuples_lp(f.reversed(),g.reversed())

    def compare_tuples_Dp(self,f,g):
        """
        Compares two exponent tuples with respect to the
        degree lexicographical term order.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return 1
        elif sf<sg:
            return -1
        elif sf == sg:
            return self.compare_tuples_lp(f,g)

    def compare_tuples_dp(self,f,g):
        """
        Compares two exponent tuples with respect to the degree
        reversed lexicographical term order.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return 1
        elif sf<sg:
            return -1
        elif sf == sg:
            return (-1)*self.compare_tuples_lp(f.reversed(),g.reversed())

    def greater_tuple_lp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        lexicographical term order.
        """
        return f > g and f or g

    def greater_tuple_rp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        reversed lexicographical term order.
        """
        return f.reversed() < g.reversed()   and f or g

    def greater_tuple_Dp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the total
        degree lexicographical term order.
        """
        return (sum(f.nonzero_values(sort=False))>sum(g.nonzero_values(sort=False))
                or (sum(f.nonzero_values(sort=False))==sum(g.nonzero_values(sort=False)) and f  > g )) and f or g

    def greater_tuple_dp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the total
        degree reversed lexicographical term order.
        """
        return (sum(f.nonzero_values(sort=False))>sum(g.nonzero_values(sort=False))
                or (sum(f.nonzero_values(sort=False))==sum(g.nonzero_values(sort=False)) and f.reversed() < g.reversed())) and f or g

    def _repr_(self):
        if self.__name == 'lex':
            s = 'Lexicographic'
        elif self.__name == 'degrevlex':
            s = 'Degree reverse lexicographic'
        else:
            s = self.__name
        return '%s term order'%s

    def singular_str(self):
        return self.__singular_str

    def macaulay2_str(self):
        return self.__macaulay2_str

    def magma_str(self):
        return self.__magma_str

    def __cmp__(self, other):
        if not isinstance(other, TermOrder):
            if isinstance(other, str):
                other = TermOrder(other)
            else:
                return cmp(type(self), type(other))
        return cmp(self.__singular_str, other.__singular_str)



