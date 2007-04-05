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
    4*z^20 + 2*y^5 + x^10
    sage: frob((x + 2*y)^3)
    3*y^15 + 2*x^5*y^10 + x^10*y^5 + x^15
    sage: (x^5 + 2*y^5)^3
    3*y^15 + 2*x^5*y^10 + x^10*y^5 + x^15

We make a polynomial ring in one variable over a polynomial ring in
two variables:
    sage: R.<x, y> = PolynomialRing(QQ, 2)
    sage: S.<t> = PowerSeriesRing(R)
    sage: t*(x+y)
    (y + x)*t
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

        We create an equal but not identical copy of the integer ring
        by dumping and loading:
            sage: S2 = loads(dumps(S))
            sage: S2 is S
            False
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
        """
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
            except:
                raise TypeError, "Unable to coerce singular object"
        elif isinstance(x , str) and self._has_singular:
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except:
                raise TypeError,"Unable to coerce string"
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


#class MPolynomialRing_macaulay2_repr_domain(MPolynomialRing_macaulay2_repr, integral_domain.IntegralDomain):
#    pass


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



