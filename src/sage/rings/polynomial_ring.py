"""
Univariate Polynomial Rings

AUTHOR:
   -- William Stein
   -- Kiran Kedlaya (2006-02-13): added macaulay2 option
   -- Martin Albrecht (2006-08-25): removed it again as it isn't needed anymore
"""


#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random, weakref

import commutative_ring
import commutative_algebra
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

from sage.interfaces.all import singular as singular_default, is_SingularElement, is_MagmaElement

from sage.rings.polynomial_singular_interface import PolynomialRing_singular_repr

def is_PolynomialRing(x):
    """
    sage: ???
    """
    return isinstance(x, PolynomialRing_generic)

#########################################################################################
# Factory function for making polynomial rings
#########################################################################################

_cache = {}

def PolynomialRing(base_ring, arg1=None, arg2=None,
                   sparse=False, order='degrevlex'):
    r"""
    Return the globally unique univariate or multivariate polynomial
    ring with given properties and variable name or names.

    There are three ways to call the polynomial ring constructor:
          1. PolynomialRing(base_ring, name,    sparse=False)
          2. PolynomialRing(base_ring, names,   order='degrevlex')
          3. PolynomialRing(base_ring, name, n, order='degrevlex')
          4. PolynomialRing(base_ring, n, name, order='degrevlex')

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:
         base_ring -- a commutative ring
         name -- a string
         names -- a list or tuple of names, or a comma separated string
         n -- an integer
         sparse -- bool (default: False), whether or not elements are sparse
         order -- string or TermOrder, e.g.,
                 'degrevlex' (default) -- degree reverse lexicographic
                 'revlex' -- reverse lexicographic
                 'lex'  -- lexicographic
                 'deglex' -- degree lexicographic
                 'wp(w1,...,wn)' -- weight reverse lexicographic
                 'Wp(w1,...,wn)' -- weight lexicographic

    OUTPUT:
        PolynomialRing(base_ring, name, sparse=False) returns a univariate
        polynomial ring; all other input formats return a multivariate
        polynomial ring.

    UNIQUENESS and IMMUTABILITY: In SAGE there is exactly one
    single-variate polynomial ring over each base ring in each choice
    of variable and sparsenes.  There is also exactly one multivariate
    polynomial ring over each base ring for each choice of names of
    variables and term order.  The names of the generators cannot be
    changed after the ring has been created (immutability).

    SQUARE BRACKETS NOTATION: You can alternatively create a single or
    multivariate polynomial ring over a ring $R$ by writing
    \code{R['varname']} or \code{R['var1,var2,var3,...']}.  This
    square brackets notation doesn't allow for setting any of the
    optional arguments.  Also, it must *never* be used in library
    code, since it injects variables into the global scope.

    EXAMPLES:
        sage: PolynomialRing(ZZ, 'y')
        Univariate Polynomial Ring in y over Integer Ring
        sage: PolynomialRing(QQ['z'], 'y')
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in z over Rational Field
        sage: PolynomialRing(QQ, 'abc')
        Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, 'abc', sparse=True)
        Sparse Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, 'y', 3, sparse=True)
        Polynomial Ring in y0, y1, y2 over Rational Field
    """
    import polynomial_ring as m

    if isinstance(arg1, (int, long, m.integer.Integer)):
        arg1, arg2 = arg2, arg1

    if not m.ring.is_Ring(base_ring):
        raise TypeError, 'base_ring must be a ring'

    R = None
    if isinstance(arg2, (int, long, m.integer.Integer)):
        # 3. PolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError, "You *must* specify the names of the variables (as a list of strings or a comma-seperated string)."
        n = int(arg2)
        names = arg1
        R = _multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str):
        if not ',' in arg1:
            # 1. PolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError, "if second arguments is a string with no commas, then there must be no other non-optional arguments"
            name = arg1
            R = _single_variate(base_ring, name, sparse)
        else:
            # 2-4. PolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
            names = arg1.split(',')
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)
    elif isinstance(arg1, (list, tuple)):
            # PolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
            names = arg1
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)

    if R is None:
        raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"

    return R


def _get_from_cache(key):
    try:
        if _cache.has_key(key):
            return _cache[key]()
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)
    return None

def _save_in_cache(key, R):
    try:
        _cache[key] = weakref.ref(R)
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)


def _single_variate(base_ring, name, sparse):
    import polynomial_ring as m
    key = (base_ring, name, sparse)
    R = _get_from_cache(key)
    if not R is None: return R

    if m.integer_mod_ring.is_IntegerModRing(base_ring) and not sparse:
        n = base_ring.order()
        if n.is_prime():
            R = m.PolynomialRing_dense_mod_p(base_ring, name)
        else:
            R = m.PolynomialRing_dense_mod_n(base_ring, name)

    elif base_ring.is_field():
        R = m.PolynomialRing_field(base_ring, name, sparse)

    elif base_ring.is_integral_domain():
        R = m.PolynomialRing_integral_domain(base_ring, name, sparse)
    else:
        R = m.PolynomialRing_generic(base_ring, name, sparse)

    _save_in_cache(key, R)
    return R

def _multi_variate(base_ring, names, n, sparse, order):
    import multi_polynomial_ring as m

    order = m.TermOrder(order)

    if isinstance(names, list):
        names = tuple(names)

    elif isinstance(names, str):
        if len(names) > 1:
            names = tuple(names)

    key = (base_ring, names, n, sparse, order)
    R = _get_from_cache(key)
    if not R is None:
        return R

    if m.integral_domain.is_IntegralDomain(base_ring):
        R = m.MPolynomialRing_polydict_domain(base_ring, n, names, order)
    else:
        R = m.MPolynomialRing_polydict(base_ring, n, names, order)

    _save_in_cache(key, R)
    return R

#########################################################################################
# END (Factory function for making polynomial rings)
#########################################################################################


#########################################################################################

class PolynomialRing_generic(commutative_algebra.CommutativeAlgebra):
    """
    Univariate polynomial ring over a commutative ring.
    """
    def __init__(self, base_ring, name=None, sparse=False):
        """
        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: R(-1) + R(1)
            0
            sage: (x - 2/3)*(x^2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3
        """
        if not isinstance(base_ring, commutative_ring.CommutativeRing):
            raise TypeError, "Base ring must be a commutative ring."
        commutative_algebra.CommutativeAlgebra.__init__(self, base_ring)
        self.__is_sparse = sparse
        ring.Ring.__init__(self)
        self.__set_polynomial_class()
        self.__generator = self([0,1], is_gen=True)
        self.__cyclopoly_cache = {}
        self._has_singular = False
        self._assign_names(name)

    def __call__(self, x=None, check=True, is_gen = False, construct=False):
        C = self.__polynomial_class
        if isinstance(x, C) and x.parent() is self:
            return x
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
        elif isinstance(x, multi_polynomial_element.MPolynomial_polydict):
            return x.univariate_polynomial(self)
        elif is_MagmaElement(x):
            x = list(x.Eltseq())
        return C(self, x, check, is_gen, construct=construct)

    def __reduce__(self):
        return sage.rings.polynomial_ring.PolynomialRing, \
               (self.base_ring(), self.variable_name(), self.__is_sparse)

    def _magma_(self, G=None):
        """
        Used in converting this ring to the corresponding ring in MAGMA.

        EXAMPLES:
            sage: R.<y> = PolynomialRing(QQ)
            sage: S = magma(R) #optional
            sage: print S #optional
            Univariate Polynomial Ring in y over Rational Field
            sage: S.1 #optional
            y

            sage: magma(PolynomialRing(GF(7))) #optional
            Univariate Polynomial Ring in x over GF(7)

            sage: magma(PolynomialRing(GF(49))) #optional
            Univariate Polynomial Ring in x over GF(7^2)

            sage: magma(PolynomialRing(PolynomialRing(ZZ,'w'))) #optional
            Univariate Polynomial Ring in x over Univariate Polynomial
            Ring in w over Integer Ring
        """
        if G is None:
            import sage.interfaces.magma
            G = sage.interfaces.magma.magma
        B = G(self.base_ring())
        R = G('PolynomialRing(%s)'%B.name())
        R._assign_names(self.variable_names())
        return R

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
        return self([self.base_ring()._coerce_(x)])

    def __cmp__(self, other):
        if not isinstance(other, PolynomialRing_generic):
            return -1
        if self.variable_name() == other.variable_name() and \
               self.base_ring() == other.base_ring():
            return 0
        elif self.variable_name() < other.variable_name():
            return -1
        return +1

    def __repr__(self):
        s = "Univariate Polynomial Ring in %s over %s"%(
                self.variable_name(), self.base_ring())
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

    def base_extend(self, R):
        return PolynomialRing(R, name=self.variable_name(), sparse=self.is_sparse())

    def characteristic(self):
        return self.base_ring().characteristic()

    def cyclotomic_polynomial(self, n):
        """
        The nth cyclotomic polynomial.

        EXAMPLES:
            sage: R = QQ['x']
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
            raise IndexError, "generator n not defined"
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
    #Q._assign_names([name])
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

    def _monics_degree( self, of_degree ):
        base = self.base_ring()
        x = self.gen()
        for lt1 in sage.misc.mrange.xmrange_iter([[base(1)]]+[base]*of_degree):
            yield sum([x**i*lt1[of_degree-i] for i in range(len(lt1))])

    def _monics_max( self, max_degree ):
        for degree in xrange(max_degree + 1):
            for m in self._monics_degree( degree ):
                yield m

    def monics( self, of_degree = None, max_degree = None ):
        """
        Return an iterator over the monic polynomials of specified degree.

        INPUT:
            Pass exactly one of:
            max_degree -- an int; the iterator will generate all monic polynomials which have degree less than or equal to max_degree
            of_degree -- an int; the iterator will generate all monic polynomials which have degree of_degree

        OUTPUT:
            an iterator

        EXAMPLES:
            sage: P = PolynomialRing(GF(4),'y')
            sage: for p in P.monics( of_degree = 2 ): print p
            y^2
            y^2 + 1
            y^2 + a
            y^2 + a + 1
            y^2 + y
            y^2 + y + 1
            y^2 + y + a
            y^2 + y + a + 1
            y^2 + a*y
            y^2 + a*y + 1
            y^2 + a*y + a
            y^2 + a*y + a + 1
            y^2 + (a + 1)*y
            y^2 + (a + 1)*y + 1
            y^2 + (a + 1)*y + a
            y^2 + (a + 1)*y + a + 1
            sage: for p in P.monics( max_degree = 1 ): print p
            1
            y
            y + 1
            y + a
            y + a + 1
            sage: for p in P.monics( max_degree = 1, of_degree = 3 ): print p
            Traceback (most recent call last):
            ...
            ValueError
        """

        if self.base_ring().order() is sage.rings.infinity.Infinity:
            raise NotImplementedError
        if of_degree is not None and max_degree is None:
            return self._monics_degree( of_degree )
        if max_degree is not None and of_degree is None:
            return self._monics_max( max_degree )
        raise ValueError # You should pass exactly one of of_degree and max_degree

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

    def lagrange_polynomial(self, points):
        """
        Return the Lagrange interpolation polynomial in self
        associated to the given list of points.

        Given a list of points, i.e. tuples of elements of self's base
        ring, this function returns the interpolation polynomial in
        the Lagrange form.

        INPUT:
            points -- a list of tuples representing points through which
                      the polynomial returned by this function must pass.

        EXAMPLE:
            sage: R = PolynomialRing(QQ)
            sage: f = R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)]);f
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1
            sage: f(0)
            1
            sage: f(2)
            2
            sage: f(3)
            -2
            sage: f(-4)
            9
            sage: R = PolynomialRing(GF(2**3))
            sage: a = R.base_ring().gen()
            sage: f = R.lagrange_polynomial([(a^2+a,a),(a,1),(a^2,a^2+a+1)]); f
            a^2*x^2 + a^2*x + a^2
            sage: f(a^2+a)
            a
            sage: f(a)
            1
            sage: f(a^2)
            a^2 + a + 1

        NOTE:  This is a straight forward Lagrange construction and no
        measures were taken to optimize it.   (If you need something
        that is highly optimized, consider implementing it and including
        it with SAGE.)
        """
        var = self.gen()

        def Pj(j):
            denom = 1
            divis = 1
            for i in range(len(points)):
                if i!=j:
                    denom *= (var          - points[i][0])
                    divis *= (points[j][0] - points[i][0])
            return denom/divis

        P = 0
        for j in range(len(points)):
            P += Pj(j)*points[j][1]
        return P


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
                raise TypeError, "Unable to coerce singular object"
        elif isinstance(x , str) and self._has_singular:
            self._singular_().set_ring()
            try:
                return self._singular_().parent(x).sage_poly(self)
            except:
                raise TypeError,"Unable to coerce string"
        return polynomial.Polynomial_dense_mod_p(self, x, check, is_gen,construct=construct)


def polygen(base_ring, name="x"):
    if isinstance(name, list):
        return polygens(base_ring, name)
    return PolynomialRing(base_ring, name).gen()
