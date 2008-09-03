r"""
Univariate Power Series Rings

EXAMPLES:
Power series rings are constructed in the standard SAGE fashion.
    sage: R.<t> = PowerSeriesRing(QQ)
    sage: R.random_element(6)
    -4 - 1/2*t^2 - 1/95*t^3 + 1/2*t^4 - 12*t^5 + O(t^6)

The default precision is specified at construction, but does not bound
the precision of created elements.
    sage: R.<t> = PowerSeriesRing(QQ, default_prec=5)
    sage: R.random_element(6)
    1/2 - 1/4*t + 2/3*t^2 - 5/2*t^3 + 2/3*t^5 + O(t^6)

    sage: S = R([1, 3, 5, 7]); S  # XXX + O(t^5)
    1 + 3*t + 5*t^2 + 7*t^3

    sage: S.truncate(3)
    5*t^2 + 3*t + 1

    sage: S.<w> = PowerSeriesRing(QQ)
    sage: S.base_ring()
    Rational Field

An iterated example:
    sage: R.<t> = PowerSeriesRing(ZZ)
    sage: S.<t2> = PowerSeriesRing(R)
    sage: S
    Power Series Ring in t2 over Power Series Ring in t over Integer Ring
    sage: S.base_ring()
    Power Series Ring in t over Integer Ring

We compute with power series over the symbolic ring.
    sage: K.<t> = PowerSeriesRing(SR, 5)
    sage: a, b, c = var('a,b,c')
    sage: f = a + b*t + c*t^2 + O(t^3)
    sage: f*f
    a^2 + ((b + a)^2 - b^2 - a^2)*t + ((c + b + a)^2 - (c + b)^2 - (b + a)^2 + 2*b^2)*t^2 + O(t^3)
    sage: f = sqrt(2) + sqrt(3)*t + O(t^3)
    sage: f^2
    2 + ((sqrt(3) + sqrt(2))^2 - 5)*t + 3*t^2 + O(t^3)

Elements are first coerced to constants in base_ring, then coerced into the
PowerSeriesRing:
    sage: R.<t> = PowerSeriesRing(ZZ)
    sage: f = Mod(2, 3) * t; (f, f.parent())
    (2*t, Power Series Ring in t over Ring of integers modulo 3)

We make a sparse power series.
    sage: R.<x> = PowerSeriesRing(QQ, sparse=True); R
    Sparse Power Series Ring in x over Rational Field
    sage: f = 1 + x^1000000
    sage: g = f*f
    sage: g.degree()
    2000000

We make a sparse Laurent series from a power series generator:
    sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
    sage: latex(-2/3*(1/t^3) + 1/t + 3/5*t^2 + O(t^5))
    \frac{-\frac{2}{3}}{t^{3}} + \frac{1}{t} + \frac{3}{5}t^{2} + O(t^{5})
    sage: S = parent(1/t); S
    Sparse Laurent Series Ring in t over Rational Field

AUTHOR:
    -- William Stein: the code
    -- Jeremy Cho (2006-05-17): some examples (above)

TESTS:
    sage: R.<t> = PowerSeriesRing(QQ)
    sage: R == loads(dumps(R))
    True

    sage: R.<x> = PowerSeriesRing(QQ, sparse=True)
    sage: R == loads(dumps(R))
    True

"""

import weakref
import power_series_poly
import power_series_mpoly
import power_series_ring_element

from polynomial.all import is_MPolynomialRing, is_PolynomialRing
from polynomial.polynomial_ring_constructor import PolynomialRing
import laurent_series_ring
import commutative_ring
import integral_domain
import field
import integer
import sage.structure.parent_gens as gens
from infinity import infinity
import sage.misc.latex as latex
from sage.structure.nonexact import Nonexact

from sage.interfaces.magma import MagmaElement
from sage.misc.sage_eval import sage_eval

from sage.structure.parent_gens import ParentWithGens


_cache = {}

def PowerSeriesRing(base_ring, name=None, default_prec=20, names=None,
                    sparse=False):
    """
    Create a power series ring.

    INPUT:
        base_ring -- a commutative ring
        name -- name of the indeterminate
        default_prec -- (efault: 20) the default precision used if an exact object
            must be changed to an approximate object in order to do an
            arithmetic operation.
        sparse -- (default: False) whether power series are represented as sparse objects.

    There is a unique power series ring over each base ring with given
    variable name.  Two power series over the same base ring with
    different variable names are not equal or isomorphic.

    EXAMPLES:
        sage: R = PowerSeriesRing(QQ, 'x'); R
        Power Series Ring in x over Rational Field

        sage: S = PowerSeriesRing(QQ, 'y'); S
        Power Series Ring in y over Rational Field

        sage: R = PowerSeriesRing(QQ, 10)
        Traceback (most recent call last):
        ...
        ValueError: first letter of variable name must be a letter

        sage: S = PowerSeriesRing(QQ, 'x', default_prec = 15); S
        Power Series Ring in x over Rational Field
        sage: S.default_prec()
        15
    """
    if isinstance(name, (int,long,integer.Integer)):
        default_prec = name
    if not names is None:
        name = names
    try:
        name = gens.normalize_names(1, name)
    except TypeError:
        raise TypeError, "illegal variable name"

    if name is None:
        raise TypeError, "You must specify the name of the indeterminate of the Power series ring."

    key = (base_ring, name, default_prec, sparse)
    if _cache.has_key(key):
        R = _cache[key]()
        if not R is None:
            return R

    if isinstance(name, (tuple, list)):
        assert len(name) == 1
        name = name[0]

    if not (name is None or isinstance(name, str)):
        raise TypeError, "variable name must be a string or None"


    if isinstance(base_ring, field.Field):
        R = PowerSeriesRing_over_field(base_ring, name, default_prec, sparse=sparse)
    elif isinstance(base_ring, integral_domain.IntegralDomain):
        R = PowerSeriesRing_domain(base_ring, name, default_prec, sparse=sparse)
    elif isinstance(base_ring, commutative_ring.CommutativeRing):
        R = PowerSeriesRing_generic(base_ring, name, default_prec, sparse=sparse)
    else:
        raise TypeError, "base_ring must be a commutative ring"
    _cache[key] = weakref.ref(R)
    return R

def is_PowerSeriesRing(R):
    """
    Return True if R is a power series ring.

    EXAMPLES:
        sage: is_PowerSeriesRing(10)
        False
        sage: is_PowerSeriesRing(QQ[['x']])
        True
    """
    return isinstance(R, PowerSeriesRing_generic)

class PowerSeriesRing_generic(commutative_ring.CommutativeRing, Nonexact):
    """
    A power series ring.
    """
    def __init__(self, base_ring, name=None, default_prec=20, sparse=False,
                 use_lazy_mpoly_ring=False):
        """
        Initializes a power series ring.

        INPUT:
            base_ring -- a commutative ring
            name -- name of the indeterminate
            default_prec -- the default precision
            sparse -- whether or not power series are sparse

            use_lazy_mpoly_ring -- if base ring is a poly ring compute
                    with multivariate polynomials instead of a
                    univariate poly over the base ring.  Only use this
                    for dense power series where you won't do too much
                    arithmetic, but the arithmetic you do must be
                    fast.  You must explicitly call
                    \code{f.do_truncation()} on an element for it to
                    truncate away higher order terms (this is called
                    automatically before printing).
        """
        ParentWithGens.__init__(self, base_ring, name)
        Nonexact.__init__(self, default_prec)
        R = PolynomialRing(base_ring, name, sparse=sparse)
        self.__poly_ring = R
        self.__power_series_class = power_series_poly.PowerSeries_poly
        self.__generator = self.__power_series_class(self, R.gen(), check=True, is_gen=True)
        self.__is_sparse = sparse
        self.__params = (base_ring, name, default_prec, sparse)

        if use_lazy_mpoly_ring and (is_MPolynomialRing(base_ring) or \
                                    is_PolynomialRing(base_ring)):
            K = base_ring
            names = K.variable_names() + (name,)
            self.__mpoly_ring = PolynomialRing(K.base_ring(), names=names)
            assert is_MPolynomialRing(self.__mpoly_ring)
            self.__power_series_class = power_series_mpoly.PowerSeries_mpoly

    def __reduce__(self):
        """
        TESTS:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: S = loads(dumps(R)); S
            Power Series Ring in t over Integer Ring
            sage: type(S)
            <class 'sage.rings.power_series_ring.PowerSeriesRing_domain'>
            sage: R.<t> = PowerSeriesRing(QQ, default_prec=10, sparse=True); R
            Sparse Power Series Ring in t over Rational Field
            sage: S = loads(dumps(R)); S
            Sparse Power Series Ring in t over Rational Field
            sage: type(S)
            <class 'sage.rings.power_series_ring.PowerSeriesRing_over_field'>
        """
        return unpickle_power_series_ring_v0, self.__params

    def _repr_(self):
        """
        Prints out a power series ring.

        EXAMPLES:
            sage: R = GF(17)[['y']]
            sage: R
            Power Series Ring in y over Finite Field of size 17
            sage: R.__repr__()
            'Power Series Ring in y over Finite Field of size 17'
            sage: R.rename('my power series ring')
            sage: R
            my power series ring
        """
        s = "Power Series Ring in %s over %s"%(self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def is_sparse(self):
        return self.__is_sparse

    def is_dense(self):
        return not self.__is_sparse

    def _latex_(self):
        r"""
        Display latex representation of this power series ring.

        EXAMPLES:
            sage: R = GF(17)[['y']]
            sage: latex(R)
            \mathbf{F}_{17}[[y]]
            sage: R = GF(17)[['y12']]
            sage: latex(R)
            \mathbf{F}_{17}[[y_{12}]]
        """
        return "%s[[%s]]"%(latex.latex(self.base_ring()), self.latex_variable_names()[0])

    def __call__(self, f, prec=infinity, check=True):
        """
        Coerce object to this power series ring.

        Returns a new instance unless the parent of f is self, in
        which case f is returned (since f is immutable).

        INPUT:
             f -- object, e.g., a power series ring element
             prec -- (default: infinity); truncation precision for coercion
             check -- bool (default: True), whether to verify that the coefficients,
                      etc., coerce in correctly.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R(t+O(t^5))
            t + O(t^5)
            sage: R(13)
            13
            sage: R(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: R([1,2,3])
            1 + 2*t + 3*t^2
            sage: S.<w> = PowerSeriesRing(QQ)
            sage: R(w + 3*w^2 + O(w^3))
            t + 3*t^2 + O(t^3)
            sage: x = polygen(QQ,'x')
            sage: R(x + x^2 + x^3 + x^5, 3)
            t + t^2 + O(t^3)
        """
        if isinstance(f, power_series_ring_element.PowerSeries) and f.parent() is self:
            if prec >= f.prec():
                return f
            f = f.truncate(prec)
        elif isinstance(f, MagmaElement) and str(f.Type()) == 'RngSerPowElt':
            v = sage_eval(f.Eltseq())
            return self(v) * (self.gen(0)**f.Valuation())
        return self.__power_series_class(self, f, prec, check=check)

    def construction(self):
        """
        Returns the functorial construction of self, namely, completion of
        the univariate polynomial ring with respect to the indeterminate
        (to a given precision).

        EXAMPLE:
            sage: R = PowerSeriesRing(ZZ, 'x')
            sage: c, S = R.construction(); S
            Univariate Polynomial Ring in x over Integer Ring
            sage: R == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return CompletionFunctor(self._names[0], self.default_prec()),  self._poly_ring()

    def _coerce_impl(self, x):
        """
        Return canonical coercion of x into self.

        Rings that canonically coerce to this power series ring R:

           * R itself
           * Any power series ring in the same variable whose base ring canonically coerces to
             the base ring of R.
           * Any ring that canonically coerces to the polynomial ring over the base ring of R.
           * Any ring that canonically coerces to the base ring of R

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R._coerce_(t + t^2)
            t + t^2
            sage: R._coerce_(1/t)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
            sage: R._coerce_(5)
            5
            sage: tt = PolynomialRing(ZZ,'t').gen()
            sage: R._coerce_(tt^2 + tt - 1)
            -1 + t + t^2
            sage: R._coerce_(1/2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
            sage: S.<s> = PowerSeriesRing(ZZ)
            sage: R._coerce_(s)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self

        We illustrate canonical coercion between power series rings with compatible
        base rings:
            sage: R.<t> = PowerSeriesRing(GF(7)['w'])
            sage: S = PowerSeriesRing(ZZ, 't')
            sage: f = S([1,2,3,4]); f
            1 + 2*t + 3*t^2 + 4*t^3
            sage: g = R._coerce_(f); g
            1 + 2*t + 3*t^2 + 4*t^3
            sage: parent(g)
            Power Series Ring in t over Univariate Polynomial Ring in w over Finite Field of size 7
            sage: S._coerce_(g)
            Traceback (most recent call last):
            ...
            TypeError: no natural map between bases of power series rings
        """
        try:
            P = x.parent()
            if is_PowerSeriesRing(P):
                if P.variable_name() == self.variable_name():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)
                    else:
                        raise TypeError, "no natural map between bases of power series rings"

        except AttributeError:
            pass
        return self._coerce_try(x, [self.base_ring(), self.__poly_ring])



    def _is_valid_homomorphism_(self, codomain, im_gens):
        r"""
        This gets called implicitly when one constructs a ring
        homomorphism from a power series ring.

        EXAMPLE:
            sage: S = RationalField(); R.<t>=PowerSeriesRing(S)
            sage: f = R.hom([0])
            sage: f(3)
            3
            sage: g = R.hom([t^2])
            sage: g(-1 + 3/5 * t)
            -1 + 3/5*t^2

        NOTE: There are no ring homomorphisms from the ring of all
        formal power series to most rings, e.g, the p-adic field,
        since you can always (mathematically!) construct some power
        series that doesn't converge.  Note that 0 is not a
        \emph{ring} homomorphism.
        """
        if im_gens[0] == 0:
            return True   # this is allowed.
        from laurent_series_ring import is_LaurentSeriesRing
        if is_PowerSeriesRing(codomain) or is_LaurentSeriesRing(codomain):
            return im_gens[0].valuation() > 0
        return False

    def _poly_ring(self):
        """
        Return the underlying polynomial ring used to represent
        elements of this power series ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R._poly_ring()
            Univariate Polynomial Ring in t over Integer Ring
        """
        return self.__poly_ring

    def _mpoly_ring(self):
        return self.__mpoly_ring

    def base_extend(self, R):
        """
        Returns the power series ring over R in the same variable as
        self, assuming there is a canonical coerce map from the base
        ring of self to R.

        EXAMPLES:
            sage: R.<T> = GF(7)[[]]; R
            Power Series Ring in T over Finite Field of size 7
            sage: R.change_ring(ZZ)
            Power Series Ring in T over Integer Ring
            sage: R.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError, "no base extension defined"

    def change_ring(self, R):
        """
        Returns the power series ring over R in the same variable as
        self.

        EXAMPLES:
            sage: R.<T> = QQ[[]]; R
            Power Series Ring in T over Rational Field
            sage: R.change_ring(GF(7))
            Power Series Ring in T over Finite Field of size 7
            sage: R.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: R.base_extend(QuadraticField(3,'a'))
            Power Series Ring in T over Number Field in a with defining polynomial x^2 - 3
        """
        return PowerSeriesRing(R, name = self.variable_name(), default_prec = self.default_prec())

    def is_exact(self):
        return False

    def gen(self, n=0):
        """
        Return the generator of this power series ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.gen()
            t
            sage: R.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: generator n>0 not defined
        """
        if n != 0:
            raise IndexError, "generator n>0 not defined"
        return self.__generator

    def ngens(self):
        """
        Return the number of generators of this power series ring.

        This is always 1.

        EXAMPLES:
            sage: R.<t> = ZZ[[]]
            sage: R.ngens()
            1
        """
        return 1

    def random_element(self, prec, bound=None):
        r"""
        Return a random power series.

        INPUT:
            prec -- an integer
            bound -- an integer (default: None, which tries to spread choice across
                         ring, if implemented)

        OUTPUT:
            power series -- a power series such that the coefficient
            of $x^i$, for $i$ up to \var{degree}, are coercions to the base
            ring of random integers between -\var{bound} and \var{bound}.

        IMPLEMENTATION: Call the random_element method on the underlying polynomial ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ)
            sage: R.random_element(5)
            -4 - 1/2*t^2 - 1/95*t^3 + 1/2*t^4 + O(t^5)
            sage: R.random_element(5,20)
            1/15 + 19/17*t + 10/3*t^2 + 5/2*t^3 + 1/2*t^4 + O(t^5)
        """
        return self(self.__poly_ring.random_element(prec, bound), prec)

    def __cmp__(self, other):
        """
        Compare this power series ring to something else.

        Power series rings are considered equal if the base ring,
        variable names, and default truncation precision are the same.

        First the base rings are compared, then the variable names,
        then the default precision.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: S.<t> = PowerSeriesRing(ZZ)
            sage: R is S
            True
            sage: R == S
            True
            sage: S.<t> = PowerSeriesRing(ZZ, default_prec=10)
            sage: R == S
            False
            sage: PowerSeriesRing(QQ,'t') == PowerSeriesRing(ZZ,'t')
            False
            sage: PowerSeriesRing(QQ,'t') == 5
            False
        """
        if not isinstance(other, PowerSeriesRing_generic):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
        c = cmp(self.default_prec(), other.default_prec())
        if c: return c
        return 0

    def __contains__(self, x):
        """
        Returns true if x is an element of this power series ring or canonically
        coerces to this ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t + t^2 in R
            True
            sage: 1/t in R
            False
            sage: 5 in R
            True
            sage: 1/3 in R
            False
            sage: S.<s> = PowerSeriesRing(ZZ)
            sage: s in R
            False
        """
        if x.parent() == self:
            return True
        try:
            self._coerce_(x)
        except TypeError:
            return False
        return True

    def is_atomic_repr(self):
        """
        Return False since power objects do not appear atomically, i.e., they have plus and spaces.
        """
        return False

    def is_field(self):
        """
        Return False since the ring of power series over any ring is never a field.
        """
        return False

    def is_finite(self):
        """
        Return False since the ring of power series over any ring is never finite.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Return the characteristic of this power series ring, which is
        the same as the characteristic of the base ring of the power
        series ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.characteristic()
            0
            sage: R.<w> = Integers(2^50)[[]]; R
            Power Series Ring in w over Ring of integers modulo 1125899906842624
            sage: R.characteristic()
            1125899906842624
        """
        return self.base_ring().characteristic()

    def laurent_series_ring(self):
        """
        If this is the power series ring $R[[t]]$, this function returns the Laurent
        series ring $R((t))$.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.laurent_series_ring()
            Laurent Series Ring in t over Integer Ring
        """
        try:
            return self.__laurent_series_ring
        except AttributeError:
            self.__laurent_series_ring = laurent_series_ring.LaurentSeriesRing(
                                                 self.base_ring(), self.variable_name(), sparse=self.is_sparse())
            return self.__laurent_series_ring

class PowerSeriesRing_domain(PowerSeriesRing_generic, integral_domain.IntegralDomain):
      pass

class PowerSeriesRing_over_field(PowerSeriesRing_domain):
    def fraction_field(self):
        """
        Return the fraction field of this power series ring, which is defined since
        this is over a field.

        This fraction field is just the Laurent series ring over the base field.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(GF(7))
            sage: R.fraction_field()
            Laurent Series Ring in t over Finite Field of size 7
            sage: Frac(R)
            Laurent Series Ring in t over Finite Field of size 7
        """
        return self.laurent_series_ring()

def unpickle_power_series_ring_v0(base_ring, name, default_prec, sparse):
    return PowerSeriesRing(base_ring, name=name, default_prec = default_prec, sparse=sparse)

