"""
Univariate Polynomial Rings

AUTHOR:
   -- William Stein
   -- Kiran Kedlaya (2006-02-13): added macaulay2 option
"""


#################################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random
import weakref

import commutative_ring
import ring
import ring_element
import field
import integral_domain
import multi_polynomial_ring
import multi_polynomial_element
import principal_ideal_domain
import polynomial_element as polynomial
import rational_field
import integer_ring
import integer
import integer_mod_ring
from sage.libs.all import pari
import sage.misc.defaults
import sage.misc.latex as latex

from sage.libs.ntl.all import ZZ as ntl_ZZ, set_modulus

from sage.interfaces.all import singular as singular_default, is_SingularElement

from sage.rings.polynomial_singular_interface import PolynomialRing_singular_repr

#_objsPolynomialRing = {}

def PolynomialRing(base_ring, name=None, sparse=False, names=None, order=None, macaulay2=False):
    """
    Return a univariate or multivariate polynomial ring.

    INPUT:
        base_ring -- the base ring
        name -- (str) the name of the generator
        sparse -- (bool; default: False) whether or not elements are represented using
                  sparse methods; note that multivariate polynomials are always sparse
        names -- names of the generators (for multivariate poly)
        order -- term order of ring
        macaulay2 (bool; default: False) -- whether or not to use Macaulay2 (multivariate only)

    EXAMPLES:
        sage: PolynomialRing(ZZ)
        Univariate Polynomial Ring in x over Integer Ring
        sage: PolynomialRing(ZZ, 'y')
        Univariate Polynomial Ring in y over Integer Ring
        sage: PolynomialRing(PolynomialRing(QQ,'z'), 'y')
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in z over Rational Field
        sage: PolynomialRing(QQ, name='abc')
        Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, name='abc', sparse=True)
        Sparse Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, 3, sparse=True)
        Polynomial Ring in x0, x1, x2 over Rational Field
        sage: PolynomialRing(QQ, 3, macaulay2=True)
        Polynomial Ring in x0, x1, x2 over Rational Field
    """
    if isinstance(name, (int,long,integer.Integer)):
        if isinstance(sparse, (list, str)):
            if order is None:
                order = names
            names = sparse
        if order is None:
            order = 'lex'
        return multi_polynomial_ring.MPolynomialRing(base_ring, n=name, names=names, order=order, macaulay2=macaulay2)

    #global _objsPolynomialRing
    #key = (base_ring, name, sparse)
    #if _objsPolynomialRing.has_key(key):
    #    R = _objsPolynomialRing[key]()
    #    if R != None:
    #        return R

    if integer_mod_ring.is_IntegerModRing(base_ring) and not sparse:
        n = base_ring.order()
        if n.is_prime():
            R = PolynomialRing_dense_mod_p(base_ring, name)
        else:
            R = PolynomialRing_dense_mod_n(base_ring, name)
    elif isinstance(base_ring, field.Field):
        R = PolynomialRing_field(base_ring, name, sparse)
    elif isinstance(base_ring, integral_domain.IntegralDomain):
        R = PolynomialRing_integral_domain(base_ring, name, sparse)
    else:
        R = PolynomialRing_generic(base_ring, name, sparse)

    #_objsPolynomialRing[key] = weakref.ref(R)
    return R

def is_PolynomialRing(x):
    return isinstance(x, PolynomialRing_generic)

class PolynomialRing_generic(commutative_ring.CommutativeRing):
    """
    Univariate polynomial ring over a commutative ring.
    """
    def __init__(self, base_ring, name=None, sparse=False):
        """
        EXAMPLES:
            sage: R, x = Q['x'].objgen()
            sage: R(-1) + R(1)
            0
            sage: (x - Q('2/3'))*(x**2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3
        """
        if not isinstance(base_ring, commutative_ring.CommutativeRing):
            raise TypeError, "Base ring (=%s) must be a commutative ring."%base_ring
        self.__base_ring = base_ring
        self.assign_names(name)
        self.__is_sparse = sparse
        ring.Ring.__init__(self)
        self.__set_polynomial_class()
        self.__generator = self([0,1], is_gen=True)
        self.__cyclopoly_cache = {}
        self._has_singular = False

    def __call__(self, x=None, check=True, is_gen = False, construct=False):
        C = self.__polynomial_class
        if isinstance(x, C) and x.parent() is self:
            return x
        elif is_SingularElement(x) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except:
                raise TypeError, "Unable to coerce singular object %s to %s (string='%s')"%(x, self, str(x))
        elif isinstance(x , str) and self._has_singular:
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except:
                raise TypeError,"Unable to coerce string %s to %s"%(x,self)
        if isinstance(x, multi_polynomial_element.MPolynomial_polydict):
            return x.univariate_polynomial(self)
        return C(self, x, check, is_gen, construct=construct)

    def __reduce__(self):
        return sage.rings.polynomial_ring.PolynomialRing, \
               (self.__base_ring, self.variable_name(), self.__is_sparse)

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

    def _coerce_(self, x):
        if isinstance(x, polynomial.Polynomial) and x.parent() is self:
            return x
        if isinstance(x, ring_element.RingElement) and x.parent() == self.base_ring():
            return self([x])
        return self([self.base_ring()._coerce_(x)])

    def __cmp__(self, other):
        if not isinstance(other, PolynomialRing_generic):
            return -1
        if self.variable_name() == other.variable_name() and \
               self.__base_ring == other.__base_ring:
            return 0
        elif self.variable_name() < other.variable_name():
            return -1
        return +1

    def __repr__(self):
        s = "Univariate Polynomial Ring in %s over %s"%(
                self.variable_name(), self.__base_ring)
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        return "%s[%s]"%(latex.latex(self.base_ring()), latex.latex(self.variable_name()))

    def __set_polynomial_class(self, cls=None):
        if not (cls is None):
            self.__polynomial_class = cls
            return
        if isinstance(self.base_ring(), rational_field.RationalField) and not self.is_sparse():
            self.__polynomial_class = polynomial.Polynomial_rational_dense
        elif isinstance(self.base_ring(), integer_ring.IntegerRing) and not self.is_sparse():
            self.__polynomial_class = polynomial.Polynomial_integer_dense
        elif isinstance(self.base_ring(), field.Field):
            if self.__is_sparse:
                self.__polynomial_class = polynomial.Polynomial_generic_sparse_field
            else:
                self.__polynomial_class = polynomial.Polynomial_generic_dense_field
        elif self.__is_sparse:
            self.__polynomial_class = polynomial.Polynomial_generic_sparse
        else:
            self.__polynomial_class = polynomial.Polynomial_generic_dense

    def base_ring(self):
        return self.__base_ring

    def base_extend(self, R):
        return PolynomialRing(R, name=self.variable_name(), sparse=self.is_sparse())

    def characteristic(self):
        return self.__base_ring.characteristic()

    def cyclotomic_polynomial(self, n):
        """
        The nth cyclotomic polynomial.

        EXAMPLES:
            sage: R = Q['x']
            sage: R.cyclotomic_polynomial(8)
            x^4 + 1
            sage: R.cyclotomic_polynomial(12)
            x^4 - x^2 + 1
            sage: S = PolynomialRing(FiniteField(7))
            sage: S.cyclotomic_polynomial(12)
            x^4 + 6*x^2 + 1
        """
        if n <= 0:
            raise ArithmeticError, "n=%s must be positive"%n
        f = pari.polcyclo(n)
        C = self.__polynomial_class
        if C == polynomial.Polynomial_rational_dense:
            return self(f, construct=True)
        coeffs = str(f.Vec())
        if C == polynomial.Polynomial_integer_dense:
            return self(coeffs, construct=True)

        coeffs = eval(coeffs)
        return self(coeffs, check=True)

    def gen(self, n=0):
        """
        If this is R[x], return x.
        """
        if n != 0:
            raise IndexError, "Generator %s not defined."%n
        return self.__generator

    def parameter(self):
        return self.gen()

    def is_field(self):
        return False

    def is_sparse(self):
        return self.__is_sparse

    def krull_dimension(self):
        return self.base_ring().krull_dimension() + 1

    def ngens(self):
        return 1

    def quotient_by_principal_ideal(self, f, names=None):
        """
        Return the quotient of this polynomial ring by the principal
        ideal generated by $f$.
        """
        import polynomial_quotient_ring
        return polynomial_quotient_ring.PolynomialQuotientRing(self, f, names)

    #def quotient(self, I, name=None):
    #    Q = commutative_ring.CommutativeRing.quotient(self, I)
    #Q.assign_names([name])
    #    return Q

    def random_element(self, degree, bound=0):
        """
        Return a random polynomial.

        INPUT:
            degree -- an int
            bound -- an int (default: 0, which tries to spread choice across ring, if implemented)

        OUTPUT:
            Polynomial -- A polynomial such that the coefficient of x^i,
            for i up to degree, are coercisions to the base ring of
            random integers between -bound and bound.

        """
        R = self.base_ring()
        return self([R.random_element(bound) for _ in xrange(degree+1)])


class PolynomialRing_integral_domain(PolynomialRing_generic, integral_domain.IntegralDomain):
    def __init__(self, base_ring, name="x", sparse=False):
        PolynomialRing_generic.__init__(self, base_ring, name, sparse)

class PolynomialRing_field(PolynomialRing_integral_domain,
                           PolynomialRing_singular_repr,
                           principal_ideal_domain.PrincipalIdealDomain,
                           ):
    def __init__(self, base_ring, name="x", sparse=False):
        PolynomialRing_generic.__init__(self, base_ring, name, sparse)
        self._has_singular = self._can_convert_to_singular()

class PolynomialRing_dense_mod_n(PolynomialRing_generic):
    def __init__(self, base_ring, name="x"):
        self.__modulus = ntl_ZZ(base_ring.order())
        PolynomialRing_generic.__init__(self, base_ring, name)

    def _ntl_set_modulus(self):
        set_modulus(self.__modulus)

    def __call__(self, x=None, check=True, is_gen = False, construct=False):
        set_modulus(self.__modulus)
        return polynomial.Polynomial_dense_mod_n(self, x, check, is_gen, construct=construct)

class PolynomialRing_dense_mod_p(PolynomialRing_dense_mod_n,
                                 PolynomialRing_singular_repr,
                                 principal_ideal_domain.PrincipalIdealDomain):
    def __init__(self, base_ring, name="x"):
        self.__modulus = ntl_ZZ(base_ring.order())
        PolynomialRing_dense_mod_n.__init__(self, base_ring, name)
        self._has_singular = self._can_convert_to_singular()

    def __call__(self, x=None, check=True, is_gen = False, construct=False):
        set_modulus(self.__modulus)
        if is_SingularElement(x) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except:
                raise TypeError, "Unable to coerce singular object %s to %s (string='%s')"%(x, self, str(x))
        elif isinstance(x , str) and self._has_singular:
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except:
                raise TypeError,"Unable to coerce string %s to %s"%(x,self)
        return polynomial.Polynomial_dense_mod_p(self, x, check, is_gen,construct=construct)


def polygen(base_ring, name="x"):
    if isinstance(name, list):
        return polygens(base_ring, name)
    return PolynomialRing(base_ring, name).gen()
