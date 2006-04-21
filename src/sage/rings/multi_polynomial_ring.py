r"""
Multivariate Polynomial Rings

AUTHORS:
    -- David Joyner and William Stein
    -- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of
              Singular features

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

from sage.ext.sage_object import SageObject

from sage.rings.integer_ring import IntegerRing

import multi_polynomial_ideal

#_cache = {}

def MPolynomialRing(base_ring, n=1, names=None,
                    order='lex', macaulay2=False):
    r"""
    Create a Multivariate polynomial ring over a commutative base ring.

    INPUT:
        base_ring -- CommutativeRing
        n -- int, number of variables  (default: 1)
        names -- tuple or string:
                   - tuple of n variable names
                   - if string, names the variables the characters in the string.
                 default: names variables x0, x1, etc.

        order -- string; the term order, or an object of type TermOrder:
                 'lex' (default) -- lexicographic
                 'revlex' -- reverse lexicographic
                 'degrevlex' -- degree reverse lexicographic
                 'deglex' -- degree lexicographic
                 'wp(w1,...,wn)' -- weight reverse lexicographic
                 'Wp(w1,...,wn)' -- weight lexicographic

        macaulay2 -- boolean (default: False); Use Macaulay2 for
                     internal representations; provides some additional
                     functionality. (Currently only supported when the
                     base ring is ZZ or a prime field.)

    EXAMPLES:
        sage: R = MPolynomialRing(RationalField(), 3)
        sage: R
        Polynomial Ring in x0, x1, x2 over Rational Field
        sage: x0,x1,x2 = R.gens()
        sage: x0.element()
        PolyDict with representation {(1, 0, 0): 1}
        sage: x0 + x1 + x2
        x2 + x1 + x0
        sage: (x0 + x1 + x2)**2
        x2^2 + 2*x1*x2 + x1^2 + 2*x0*x2 + 2*x0*x1 + x0^2

    This example illustrates the quick shorthand for naming several
    variables one-letter names.
        sage: MPolynomialRing(ZZ, 4, 'xyzw')
        Polynomial Ring in x, y, z, w over Integer Ring

    To obtain both the ring and its generators, use the \code{objgens} function.
        sage: R, (x,y,z,w) = MPolynomialRing(ZZ, 4, 'xyzw').objgens()
        sage: (x+y+z+w)^2
        w^2 + 2*z*w + z^2 + 2*y*w + 2*y*z + y^2 + 2*x*w + 2*x*z + 2*x*y + x^2

    We can construct multi-variate polynomials rings over completely
    arbitrary SAGE rings.  In this example, we construct a polynomial
    ring S in 3 variables over a polynomial ring in 2 variables over
    GF(9).  Then we construct a polynomial ring in 20 variables over S!

        sage: R, (n1,n2) = MPolynomialRing(GF(9),2, names=['n1','n2']).objgens()
        sage: n1^2 + 2*n2
        2*n2 + n1^2
        sage: S = MPolynomialRing(R,3, names='a'); a0,a1,a2=S.gens()
        sage: S
        Polynomial Ring in a0, a1, a2 over Polynomial Ring in n1, n2 over Finite Field in a of size 3^2
        sage: x = (n1+n2)*a0 + 2*a1**2
        sage: x
        2*a1^2 + (n2 + n1)*a0
        sage: x**3
        2*a1^6 + (n2^3 + n1^3)*a0^3
        sage: T = MPolynomialRing(S, 20)
        sage: T
        Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19 over Polynomial Ring in a0, a1, a2 over Polynomial Ring in n1, n2 over Finite Field in a of size 3^2

    We create a polynomial ring that uses the Macaualy2 interface.
        sage: R, (x,y,z) = MPolynomialRing(ZZ, 3, 'xyz', macaulay2=True).objgens()     # optional
        sage: type(R)                                                                  # optional
        <class 'sage.rings.multi_polynomial_ring.MPolynomialRing_macaulay2_repr_domain'>
    """
    global _cache
    T = TermOrder(order)
    if isinstance(names, list):
        names = tuple(names)

    #elif isinstance(names, str):
    #    if len(names) > 1:
    #        names = tuple(names)
    #key = (base_ring, n, names, T, macaulay2)
    #if _cache.has_key(key):
    #    R = _cache[key]()
    #    if not (R is None):
    #        return R

    if not isinstance(base_ring, commutative_ring.CommutativeRing):
        raise TypeError, "Base ring (=%s) must be a commutative ring."%base_ring

    if macaulay2:
        if integral_domain.is_IntegralDomain(base_ring):
            R = MPolynomialRing_macaulay2_repr_domain(base_ring, n, names, T)
        else:
            R = MPolynomialRing_macaulay2_repr(base_ring, n, names, T)
    elif finite_field.is_FiniteField(base_ring) or base_ring.is_prime_field():
        R = MPolynomialRing_singular_repr_domain(base_ring, n, names, T)
    else:
        if integral_domain.is_IntegralDomain(base_ring):
            R = MPolynomialRing_polydict_domain(base_ring, n, names, T)
        else:
            R = MPolynomialRing_polydict(base_ring, n, names, T)

    #_cache[key] = weakref.ref(R)

    return R

def is_MPolynomialRing(x):
    return isinstance(x, MPolynomialRing_generic)

class MPolynomialRing_generic(commutative_ring.CommutativeRing):
    def __init__(self, base_ring, n, names, order):
        if not isinstance(base_ring, commutative_ring.CommutativeRing):
            raise TypeError, "Base ring (=%s) must be a commutative ring."%base_ring
        self.__base_ring = base_ring
        n = int(n)
        if n < 0:
            raise ValueError, "Multivariate Polynomial Rings must " + \
                  "have more than 0 variables."
        self.__ngens = n
        self.assign_names(names)
        self.__term_order = order

    def __cmp__(self, right):
        if not is_MPolynomialRing(right):
            return -1
        return cmp((self.__base_ring, self.__ngens, self.variable_names(), self.__term_order),
                   (right.__base_ring, right.__ngens, right.variable_names(), right.__term_order))

    def __contains__(self, x):
        """
        This definition of containment does not involve a natural
        inclusion from rings with less variables into rings with more.
        """
        try:
            return x.parent() == self
        except AttributeError:
            return False

    def _repr_(self):
        return "Polynomial Ring in %s over %s"%(", ".join(self.variable_names()), self.base_ring())

    def _latex_(self):
        vars = str(self.latex_variable_names()).replace('\n','').replace("'",'')
        return "%s[%s]"%(latex.latex(self.base_ring()), vars[1:-1])


    def _ideal_class_(self):
        return multi_polynomial_ideal.MPolynomialIdeal

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # all that is needed is that elements of the base ring
            # of the polynomial ring canonically coerce into codomain.
            # Since poly rings are free, any image of the gen
            # determines a homomorphism
            codomain._coerce_(self.base_ring()(1))
        except TypeError:
            return False
        return True

    def var_dict(self):
        """
        Return dictionary of paris varname:var of the variables
        of this multivariate polynomial ring.
        """
        return dict([(str(g),g) for g in self.gens()])


    def is_finite(self):
        if self.ngens() == 0:
            return self.base_ring().is_finite()
        return False

    def is_field(self):
        """
        Return True if this multivariate polynomial ring is a field, i.e.,
        it is a ring in 0 generators over a field.
        """
        if self.ngens() == 0:
            return self.base_ring().is_field()
        return False

    def term_order(self):
        return self.__term_order

    def base_ring(self):
        return self.__base_ring

    def characteristic(self):
        """
        Return the characteristic of this polynomial ring.

        EXAMPLES:
            sage: R = MPolynomialRing(RationalField(), 3)
            sage: R.characteristic()
            0
            sage: R = MPolynomialRing(GF(7),20)
            sage: R.characteristic()
            7
        """
        return self.__base_ring.characteristic()

    def gen(self, n=0):
        if n < 0 or n >= self.__ngens:
            raise ValueError, "Generator %s not defined."%n
        return self._gens[int(n)]

    def gens(self):
        return self._gens

    def krull_dimension(self):
        return self.base_ring().krull_dimension() + self.ngens()

    def ngens(self):
        return self.__ngens

    def _monomial_order_function(self):
        raise NotImplementedError

    def latex_variable_names(self):
        """
        Returns the list of variable names suitable for latex output.

        All '_SOMETHING' substrings are replaced by '_{SOMETHING}' recursively
        so that subscripts of subscripts work.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ,12).objgens();
            sage: x
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
            sage: print R.latex_variable_names ()
            ['x_{0}', 'x_{1}', 'x_{2}', 'x_{3}', 'x_{4}', 'x_{5}', 'x_{6}', 'x_{7}', 'x_{8}', 'x_{9}', 'x_{10}', 'x_{11}']
            sage: f = x[0]^3 + 15/3 * x[1]^10
            sage: print latex(f)
            5 x_{1}^{10} + x_{0}^{3}
        """
        try:
            return self.__latex_variable_names
        except AttributeError:
            pass
        names = []
        for g in self.variable_names():
            i = len(g)-1
            while i >= 0 and g[i].isdigit():
                i -= 1
            if i < len(g)-1:
                g = '%s_{%s}'%(g[:i+1], g[i+1:])
            names.append(g)
        self.__latex_variable_names = names
        return names

    def assign_names(self, names=None):
        try:
            del self.__latex_variable_names
        except AttributeError:
            pass
        commutative_ring.CommutativeRing.assign_names(self, names)

class MPolynomialRing_polydict(MPolynomialRing_generic):
    """
    Multivariable polynomial ring.

    EXAMPLES:
        sage: R = MPolynomialRing(Integers(12),5); R
        Polynomial Ring in x0, x1, x2, x3, x4 over Ring of integers modulo 12
        sage: loads(R.dumps()) == R
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

    def __call__(self, x, check=True):
        """
        Coerce x into this multivariate polynomial ring.
        """
        if isinstance(x, multi_polynomial_element.MPolynomial_polydict) and x.parent() == self:
            return x
        elif isinstance(x, polydict.PolyDict):
            return multi_polynomial_element.MPolynomial_polydict(self, x)
        elif isinstance(x, fraction_field_element.FractionFieldElement) and x.parent().ring() == self:
            if x.denominator() == 1:
                return x.numerator()
            else:
                raise TypeError, "unable to coerce in %s since the denominator is not 1"%x
        c = self.base_ring()(x)
        return multi_polynomial_element.MPolynomial_polydict(self, {self._zero_tuple:c})

class MPolynomialRing_polydict_domain(integral_domain.IntegralDomain, MPolynomialRing_polydict):
    pass


class MPolynomialRing_singular_repr(MPolynomialRing_polydict):
    def _singular_(self, singular=singular_default):
        """
        Returns a singular ring for a given MPolynomialRing
        over a finite field.

        INPUT:
           singular -- Singular instance

        OUTPUT:
           singular ring matching this ring

        EXAMPLES:
           sage: r=MPolynomialRing(GF(2**8),10,'x')
           sage: r._singular_()
           //   characteristic : 2
           //   1 parameter    : a
           //   minpoly        : (a^8+a^4+a^3+a^2+1)
           //   number of vars : 10
           //        block   1 : ordering lp
           //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
           //        block   2 : ordering C
           sage: r=MPolynomialRing(GF(127),2,'x')
           sage: r._singular_()
           //   characteristic : 127
           //   number of vars : 2
           //        block   1 : ordering lp
           //                  : names    x0 x1
           //        block   2 : ordering C
           sage: r=MPolynomialRing(QQ,2,'x')
           sage: r._singular_()
           //   characteristic : 0
           //   number of vars : 2
           //        block   1 : ordering lp
           //                  : names    x0 x1
           //        block   2 : ordering C

        WARNING:
           If the base ring is a finite extension field the ring will not only be
           returned but also be set as the current ring in Singular.
        """
        try:
            R = self.__singular
            if not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite():
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != self.__minpoly:
                    singular.eval("minpoly=%s"%(self.__minpoly))
            return R
        except (AttributeError, ValueError):
            if not self.base_ring().is_finite() and not self.base_ring().is_prime_field():
                raise TypeError, "no conversion of %s to a Singular ring defined"%self
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
        Return a newly created singular ring matching this ring.
        """
        if self.base_ring().is_prime_field():
            self.__singular = singular.ring(self.characteristic(), \
                                            str(self.gens()), \
                                            self.term_order().singular_str())
            return self.__singular
        if self.base_ring().is_finite(): #must be extension field
            gen = str(self.base_ring().gen())
            if self.ngens()==1:
              _vars = str(self.gen())
            else:
              _vars = str(self.gens())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen),
                               _vars,
                               self.term_order().singular_str() )
            self.__minpoly = "("+(str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")+")"
            singular.eval("minpoly=%s"%(self.__minpoly) )
            self.__singular = r
            return self.__singular
        if isinstance(self.base_ring(),IntegerRing): #char=0
            self.__singular = singular.ring(0, str(self.gens()), \
                                                self.term_order().singular_str())
            return self.__singular
        raise TypeError, "no conversion of %s to a Singular ring defined"%self

    def __call__(self, x, check=True):
        """
        Coerce x into this multivariate polynomial ring.

        EXAMPLES:
        We create a Singular multivariate polynomial via ideal arithmetic,
        then coerce it into R.
            sage: R, (x,y) = PolynomialRing(QQ, 2, ['x','y']).objgens()
            sage: I = R.ideal([x^3 + y, y])
            sage: S = singular(I)
            sage: f = (S*S*S)[2]
            sage: f
            x^3*y^2+y^3
            sage: R(f)
            y^3 + x^3*y^2
        """
        if isinstance(x, multi_polynomial_element.MPolynomial_singular_repr) and x.parent() == self:
            return x
        elif isinstance(x, polydict.PolyDict):
            return multi_polynomial_element.MPolynomial_singular_repr(self, x)
        elif is_SingularElement(x):
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except:
                raise TypeError, "Unable to coerce singular object %s to %s (string='%s')"%(x, self, str(x))
        elif isinstance(x, fraction_field_element.FractionFieldElement) and x.parent().ring() == self:
            if x.denominator() == 1:
                return x.numerator()
            else:
                raise TypeError, "unable to coerce in %s since the denominator is not 1"%x
        elif isinstance(x , str):
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except:
                raise TypeError,"Unable to coerce string %s to %s"%(x,self)
        c = self.base_ring()(x)
        return multi_polynomial_element.MPolynomial_singular_repr(self, {self._zero_tuple:c})

    def _poly_class(self):
        return multi_polynomial_element.MPolynomial_singular_repr

    def ideal(self, gens, coerce=True):
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return multi_polynomial_ideal.MPolynomialIdeal_singular_repr(self, gens, coerce=False)

class MPolynomialRing_singular_repr_domain(MPolynomialRing_singular_repr, integral_domain.IntegralDomain):
    pass

class MPolynomialRing_macaulay2_repr(MPolynomialRing_polydict):
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
            elif isinstance(self.base_ring(), IntegerRing):
                base_str = "ZZ"
            else:
                raise TypeError, "no conversion of %s to a Macaulay2 ring defined"%self
            self.__macaulay2 = macaulay2.ring(base_str, str(self.gens()), \
                                              self.term_order().macaulay2_str())
        return self.__macaulay2

    def __call__(self, x, check=True):
        """
        Coerce x into this multivariate polynomial ring.

        EXAMPLES:
        We create a Macaulay2 multivariate polynomial via ideal arithmetic,
        then coerce it into R.
            sage: R, (x,y) = PolynomialRing(QQ, 2, ['x','y'], macaulay2=True).objgens()       # optional
            sage: I = R.ideal([x^3 + y, y])                                                   # optional
            sage: S = I._macaulay2_()                                                         # optional
            sage: T = S*S*S                                                                   # optional
            sage: U = T.gens().entries().flatten()                                            # optional
            sage: f = U[2]; f                                                                 # optional
            x^3*y^2+y^3
            sage: R(f)                                                                        # optional
            y^3 + x^3*y^2
        """
        if isinstance(x, multi_polynomial_element.MPolynomial_macaulay2_repr) and x.parent() == self:
            return x
        elif isinstance(x, polydict.PolyDict):
            return multi_polynomial_element.MPolynomial_macaulay2_repr(self, x)
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
                raise TypeError, "Unable to coerce macaulay2 object %s to %s (string='%s')"%(x, self, s)
            return multi_polynomial_element.MPolynomial_macaulay2_repr(self, x)
        elif isinstance(x, fraction_field_element.FractionFieldElement) and x.parent().ring() == self:
            if x.denominator() == 1:
                return x.numerator()
            else:
                raise TypeError, "unable to coerce in %s since the denominator is not 1"%x
        c = self.base_ring()(x)
        return multi_polynomial_element.MPolynomial_macaulay2_repr(self, {self._zero_tuple:c})

    def _poly_class(self):
        return multi_polynomial_element.MPolynomial_macaulay2_repr

    def ideal(self, gens, coerce=True):
        if is_Macaulay2Element(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # will this coerce from macaulay2 ideals correctly?
        return multi_polynomial_ideal.MPolynomialIdeal_macaulay2_repr(self, gens, coerce=False)


class MPolynomialRing_macaulay2_repr_domain(MPolynomialRing_macaulay2_repr, integral_domain.IntegralDomain):
    pass


#######################

name_mapping = {'lex':'lp', \
                'revlex':'rp', \
                'degrevlex':'dp', \
                'deglex':'Dp'}

m2_name_mapping = {'lex':'Lex', \
                   'revlex':'RevLex', \
                   'degrevlex':'GRevLex', \
                   'deglex':'GLex'}

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
        if name_mapping.has_key(name):
            singular_name = name_mapping[name]
            self.__singular_str = singular_name
        else:
            self.__singular_str = name
        if m2_name_mapping.has_key(name):
            macaulay2_name = m2_name_mapping[name]
            self.__macaulay2_str = macaulay2_name
        else:
            self.__macaulay2_str = name

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

    def __cmp__(self, other):
        if not isinstance(other, TermOrder):
            if isinstance(other, str):
                other = TermOrder(other)
            else:
                return -1
        return cmp(self.__name, other.__name)

