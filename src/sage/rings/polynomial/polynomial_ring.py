
"""
Univariate Polynomial Rings

Sage implements sparse and dense polynomials over commutative and
non-commutative rings.  In the non-commutative case, the polynomial
variable commutes with the elements of the base ring.

AUTHOR:

- William Stein

- Kiran Kedlaya (2006-02-13): added macaulay2 option

- Martin Albrecht (2006-08-25): removed it again as it isn't needed anymore

- Simon King (2011-05): Dense and sparse polynomial rings must not be equal.

- Simon King (2011-10): Choice of categories for polynomial rings.

EXAMPLES:

Creating a polynomial ring injects the variable into the interpreter namespace::

    sage: z = QQ['z'].0
    sage: (z^3 + z - 1)^3
    z^9 + 3*z^7 - 3*z^6 + 3*z^5 - 6*z^4 + 4*z^3 - 3*z^2 + 3*z - 1

Saving and loading of polynomial rings works::

    sage: loads(dumps(QQ['x'])) == QQ['x']
    True
    sage: k = PolynomialRing(QQ['x'],'y'); loads(dumps(k))==k
    True
    sage: k = PolynomialRing(ZZ,'y'); loads(dumps(k)) == k
    True
    sage: k = PolynomialRing(ZZ,'y', sparse=True); loads(dumps(k))
    Sparse Univariate Polynomial Ring in y over Integer Ring

Rings with different variable names are not equal; in fact,
by trac ticket #9944, poynomial rings are equal if and only
if they are identic (which should be the  case for all parent
structures in Sage)::

    sage: QQ['y'] != QQ['x']
    True
    sage: QQ['y'] != QQ['z']
    True

We create a polynomial ring over a quaternion algebra::

    sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
    sage: R.<w> = PolynomialRing(A,sparse=True)
    sage: f = w^3 + (i+j)*w + 1
    sage: f
    w^3 + (i + j)*w + 1
    sage: f^2
    w^6 + (2*i + 2*j)*w^4 + 2*w^3 - 2*w^2 + (2*i + 2*j)*w + 1
    sage: f = w + i ; g = w + j
    sage: f * g
    w^2 + (i + j)*w + k
    sage: g * f
    w^2 + (i + j)*w - k

Trac ticket #9944 introduced some changes related with
coercion. Previously, a dense and a sparse polynomial ring with the
same variable name over the same base ring evaluated equal, but of
course they were not identical.Coercion maps are cached - but if a
coercion to a dense ring is requested and a coercion to a sparse ring
is returned instead (since the cache keys are equal!), all hell breaks
loose.

Therefore, the coercion between rings of sparse and dense polynomials
works as follows::

    sage: R.<x> = PolynomialRing(QQ, sparse=True)
    sage: S.<x> = QQ[]
    sage: S == R
    False
    sage: S.has_coerce_map_from(R)
    True
    sage: R.has_coerce_map_from(S)
    False
    sage: (R.0+S.0).parent()
    Univariate Polynomial Ring in x over Rational Field
    sage: (S.0+R.0).parent()
    Univariate Polynomial Ring in x over Rational Field

It may be that one has rings of dense or sparse polynomials over
different base rings. In that situation, coercion works by means of
the :func:`~sage.categories.pushout.pushout` formalism::

    sage: R.<x> = PolynomialRing(GF(5), sparse=True)
    sage: S.<x> = PolynomialRing(ZZ)
    sage: R.has_coerce_map_from(S)
    False
    sage: S.has_coerce_map_from(R)
    False
    sage: S.0 + R.0
    2*x
    sage: (S.0 + R.0).parent()
    Univariate Polynomial Ring in x over Finite Field of size 5
    sage: (S.0 + R.0).parent().is_sparse()
    False

Similarly, there is a coercion from the (non-default) NTL
implementation for univariate polynomials over the integers
to the default FLINT implementation, but not vice versa::

    sage: R.<x> = PolynomialRing(ZZ, implementation = 'NTL')
    sage: S.<x> = PolynomialRing(ZZ, implementation = 'FLINT')
    sage: (S.0+R.0).parent() is S
    True
    sage: (R.0+S.0).parent() is S
    True

TESTS::

    sage: K.<x>=FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x

Check that :trac:`5562` has been fixed::

    sage: R.<u> = PolynomialRing(RDF, 1, 'u')
    sage: v1 = vector([u])
    sage: v2 = vector([CDF(2)])
    sage: v1 * v2
    2.0*u

These may change over time::

    sage: type(ZZ['x'].0)
    <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    sage: type(QQ['x'].0)
    <type 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>
    sage: type(RR['x'].0)
    <type 'sage.rings.polynomial.polynomial_real_mpfr_dense.PolynomialRealDense'>
    sage: type(Integers(4)['x'].0)
    <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>
    sage: type(Integers(5*2^100)['x'].0)
    <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>
    sage: type(CC['x'].0)
    <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>
    sage: type(CC['t']['x'].0)
    <type 'sage.rings.polynomial.polynomial_element.Polynomial_generic_dense'>
"""


#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element
from sage.structure.category_object import check_default_category
import sage.algebras.algebra
import sage.categories.basic as categories
import sage.rings.commutative_ring as commutative_ring
import sage.rings.commutative_algebra as commutative_algebra
import sage.rings.ring as ring
import sage.rings.ring_element as ring_element
import sage.rings.integral_domain as integral_domain
import sage.rings.principal_ideal_domain as principal_ideal_domain
import sage.rings.polynomial.polynomial_element_generic as polynomial_element_generic
import sage.rings.rational_field as rational_field
from sage.rings.integer_ring import is_IntegerRing, IntegerRing
from sage.rings.integer import Integer
from sage.libs.pari.all import pari_gen
from sage.rings.polynomial.polynomial_ring_constructor import polynomial_default_category

import sage.misc.latex as latex
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method

from sage.rings.real_mpfr import is_RealField
from polynomial_real_mpfr_dense import PolynomialRealDense
from sage.rings.polynomial.polynomial_singular_interface import PolynomialRing_singular_repr
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.finite_rings.element_base import FiniteRingElement

from polynomial_element import PolynomialBaseringInjection

from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

import cyclotomic

ZZ_sage = IntegerRing()

from sage.interfaces.singular import SingularElement


def is_PolynomialRing(x):
    """
    Return True if x is a *univariate* polynomial ring (and not a
    sparse multivariate polynomial ring in one variable).

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        sage: from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
        sage: is_PolynomialRing(2)
        False

    This polynomial ring is not univariate.

    ::

        sage: is_PolynomialRing(ZZ['x,y,z'])
        False
        sage: is_MPolynomialRing(ZZ['x,y,z'])
        True

    ::

        sage: is_PolynomialRing(ZZ['w'])
        True

    Univariate means not only in one variable, but is a specific data
    type. There is a multivariate (sparse) polynomial ring data type,
    which supports a single variable as a special case.

    ::

        sage: is_PolynomialRing(PolynomialRing(ZZ,1,'w'))
        False
        sage: R = PolynomialRing(ZZ,1,'w'); R
        Multivariate Polynomial Ring in w over Integer Ring
        sage: is_PolynomialRing(R)
        False
        sage: type(R)
        <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
    """
    return isinstance(x, PolynomialRing_general)


#########################################################################################

class PolynomialRing_general(sage.algebras.algebra.Algebra):
    """
    Univariate polynomial ring over a ring.
    """
    _no_generic_basering_coercion = True
    def __init__(self, base_ring, name=None, sparse=False, element_class=None, category=None):
        """
        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: R(-1) + R(1)
            0
            sage: (x - 2/3)*(x^2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3

            sage: category(ZZ['x'])
            Join of Category of unique factorization domains and
                    Category of commutative algebras over euclidean domains
            sage: category(GF(7)['x'])
            Join of Category of euclidean domains and
                    Category of commutative algebras over (finite fields and
                        subquotients of monoids and quotients of semigroups)
        """
        # We trust that, if category is given, it is useful and does not need to be joined
        # with the default category
        if category is None:
            category = polynomial_default_category(base_ring.category(),False)
        self.__is_sparse = sparse
        if element_class:
            self._polynomial_class = element_class
        else:
            if sparse:
                self._polynomial_class = polynomial_element_generic.Polynomial_generic_sparse
            else:
                from sage.rings.polynomial import polynomial_element
                self._polynomial_class = polynomial_element.Polynomial_generic_dense
        self.__cyclopoly_cache = {}
        self._has_singular = False
        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which
        # tries to provide a conversion from the base ring, if it does not exist.
        # This is for algebras that only do the generic stuff in their initialisation.
        # But the attribute _no_generic_basering_coercion prevents that from happening,
        # since we want to use PolynomialBaseringInjection.
        sage.algebras.algebra.Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._polynomial_class(self, [0,1], is_gen=True)
        self._populate_coercion_lists_(
                #coerce_list = [base_inject],
                #convert_list = [list, base_inject],
                convert_method_name = '_polynomial_')
        if is_PolynomialRing(base_ring):
            self._Karatsuba_threshold = 0
        else:
            from sage.matrix.matrix_space import is_MatrixSpace
            if is_MatrixSpace(base_ring):
                self._Karatsuba_threshold = 0
            else:
                self._Karatsuba_threshold = 8

    def __reduce__(self):
        import sage.rings.polynomial.polynomial_ring_constructor
        return (sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing,
                (self.base_ring(), self.variable_name(), None, self.is_sparse()))


    def _element_constructor_(self, x=None, check=True, is_gen = False, construct=False, **kwds):
        r"""
        Convert ``x`` into this univariate polynomial ring,
        possibly non-canonically.

        Conversion from power series::

            sage: R.<x> = QQ[]
            sage: R(1 + x + x^2 + O(x^3))
            x^2 + x + 1

        Stacked polynomial rings coerce into constants if possible. First,
        the univariate case::

            sage: R.<x> = QQ[]
            sage: S.<u> = R[]
            sage: S(u + 2)
            u + 2
            sage: S(x + 3)
            x + 3
            sage: S(x + 3).degree()
            0

        Second, the multivariate case::

            sage: R.<x,y> = QQ[]
            sage: S.<u> = R[]
            sage: S(x + 2*y)
            x + 2*y
            sage: S(x + 2*y).degree()
            0
            sage: S(u + 2*x)
            u + 2*x
            sage: S(u + 2*x).degree()
            1

        Foreign polynomial rings coerce into the highest ring; the point
        here is that an element of T could coerce to an element of R or an
        element of S; it is anticipated that an element of T is more likely
        to be "the right thing" and is historically consistent.

        ::

            sage: R.<x> = QQ[]
            sage: S.<u> = R[]
            sage: T.<a> = QQ[]
            sage: S(a)
            u

        Coercing in pari elements::

            sage: QQ['x'](pari('[1,2,3/5]'))
            3/5*x^2 + 2*x + 1
            sage: QQ['x'](pari('(-1/3)*x^10 + (2/3)*x - 1/5'))
            -1/3*x^10 + 2/3*x - 1/5

        Coercing strings::

            sage: QQ['y']('-y')
            -y

        TESTS:

        This shows that the issue at trac #4106 is fixed::

            sage: x = var('x')
            sage: R = IntegerModRing(4)
            sage: S = PolynomialRing(R, x)
            sage: S(x)
            x

        Throw a TypeError if any of the coefficients cannot be coerced
        into the base ring (trac #6777)::

            sage: RealField(300)['x']( [ 1, ComplexField(300).gen(), 0 ])
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='1.00...00*I') to real number.

        Check that the bug in :trac:`11239` is fixed::

            sage: K.<a> = GF(5^2, conway=True, prefix='z')
            sage: L.<b> = GF(5^4, conway=True, prefix='z')
            sage: f = K['x'].gen() + a
            sage: L['x'](f)
            x + b^3 + b^2 + b + 3

        """
        C = self._polynomial_class
        if isinstance(x, list):
            return C(self, x, check=check, is_gen=False,construct=construct)
        if isinstance(x, Element):
            P = x.parent()
            if P is self:
                return x
            elif P is self.base_ring():
                # It *is* the base ring, hence, we should not need to check.
                # Moreover, if x is equal to zero then we usually need to
                # provide [] to the polynomial class, not [x], if we don't want
                # to check (normally, polynomials like to strip trailing zeroes).
                # However, in the padic case, we WANT that trailing
                # zeroes are not stripped, because O(5)==0, but still it must
                # not be forgotten. It should be the job of the __init__ method
                # to decide whether to strip or not to strip.
                return C(self, [x], check=False, is_gen=False,
                         construct=construct)
            elif P == self.base_ring():
                return C(self, [x], check=True, is_gen=False,
                         construct=construct)
        if isinstance(x, SingularElement) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except Exception:
                raise TypeError("Unable to coerce singular object")
        elif isinstance(x , str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(x))
            except NameError:
                raise TypeError("Unable to coerce string")
        elif isinstance(x, FractionFieldElement):
            if x.denominator().is_unit():
                x = x.numerator() * x.denominator().inverse_of_unit()
            else:
                raise TypeError("denominator must be a unit")
        elif isinstance(x, pari_gen):
            if x.type() == 't_RFRAC':
                raise TypeError("denominator must be a unit")
            if x.type() != 't_POL':
                x = x.Polrev()
        elif isinstance(x, FiniteRingElement):
            try:
                return self(x.polynomial())
            except AttributeError:
                pass
        elif isinstance(x, sage.rings.power_series_ring_element.PowerSeries):
            x = x.truncate()
        return C(self, x, check, is_gen, construct=construct, **kwds)

    def is_integral_domain(self, proof = True):
        """
        EXAMPLES::

            sage: ZZ['x'].is_integral_domain()
            True
            sage: Integers(8)['x'].is_integral_domain()
            False
        """
        return self.base_ring().is_integral_domain(proof)

    def is_noetherian(self):
        return self.base_ring().is_noetherian()

    def construction(self):
        from sage.categories.pushout import PolynomialFunctor
        return PolynomialFunctor(self.variable_name(), sparse=self.__is_sparse), self.base_ring()

    def completion(self, p, prec=20, extras=None):
        """
        Return the completion of self with respect to the irreducible
        polynomial p. Currently only implemented for p=self.gen(), i.e. you
        can only complete R[x] with respect to x, the result being a ring
        of power series in x. The prec variable controls the precision used
        in the power series ring.

        EXAMPLES::

            sage: P.<x>=PolynomialRing(QQ)
            sage: P
            Univariate Polynomial Ring in x over Rational Field
            sage: PP=P.completion(x)
            sage: PP
            Power Series Ring in x over Rational Field
            sage: f=1-x
            sage: PP(f)
            1 - x
            sage: 1/f
            1/(-x + 1)
            sage: 1/PP(f)
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11 + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + O(x^20)
        """
        if str(p) == self._names[0]:
            from sage.rings.power_series_ring import PowerSeriesRing
            return PowerSeriesRing(self.base_ring(), name=self._names[0], default_prec=prec)
        else:
            raise TypeError("Cannot complete %s with respect to %s" % (self, p))

    def _coerce_map_from_(self, P):
        """
        The rings that canonically coerce to this polynomial ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring.

        - polynomial rings in the same variable over any base ring that
          canonically coerces to the base ring of this ring.

        - a multivariate polynomial ring P such that self's variable name
          is among the variable names of P, and the ring obtained by
          removing that variable is different from the base ring of self,
          but coerces into it. (see trac ticket #813 for a discussion of this)

        Caveat: There is no coercion from a dense into a sparse
        polynomial ring. So, when adding a dense and a sparse
        polynomial, the result will be dense. See trac ticket #9944.

        EXAMPLES::

            sage: R = QQ['x']
            sage: R.has_coerce_map_from(QQ)
            True
            sage: R.has_coerce_map_from(ZZ)
            True
            sage: R.has_coerce_map_from(GF(7))
            False
            sage: R.has_coerce_map_from(ZZ['x'])
            True
            sage: R.has_coerce_map_from(ZZ['y'])
            False

            sage: R.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Univariate Polynomial Ring in x over Rational Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Polynomial base injection morphism:
                      From: Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field

        Here we test against the change in the coercions introduced
        in trac ticket #9944::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: S.<x> = QQ[]
            sage: (R.0+S.0).parent()
            Univariate Polynomial Ring in x over Rational Field
            sage: (S.0+R.0).parent()
            Univariate Polynomial Ring in x over Rational Field

        Here we test a feature that was implemented in trac ticket #813::

            sage: P = QQ['x','y']
            sage: Q = Frac(QQ['x'])['y']
            sage: Q.has_coerce_map_from(P)
            True
            sage: P.0+Q.0
            y + x

        In order to avoid bidirectional coercions (which are generally
        problematic), we only have a coercion from P to Q if the base
        ring of Q is more complicated than "P minus one variable"::

            sage: Q = QQ['x']['y']
            sage: P.has_coerce_map_from(Q)
            True
            sage: Q.has_coerce_map_from(P)
            False
            sage: Q.base_ring() is P.remove_var(Q.variable_name())
            True
        """
        # In the first place, handle the base ring
        base_ring = self.base_ring()
        if P is base_ring:
            return PolynomialBaseringInjection(base_ring, self)
        # handle constants that canonically coerce into self.base_ring()
        # first, if possible
        try:
            connecting = base_ring.coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(base_ring) * connecting
        except TypeError:
            pass

        # polynomial rings in the same variable over a base that canonically
        # coerces into self.base_ring()
        try:
            if is_PolynomialRing(P):
                if self.__is_sparse and not P.is_sparse():
                    return False
                if P.variable_name() == self.variable_name():
                    if P.base_ring() is base_ring and \
                            base_ring is ZZ_sage:
                        # We're trying to coerce from FLINT->NTL
                        # or vice versa.  Only allow coercions from
                        # NTL->FLINT, not vice versa.
                        # Unfortunately this doesn't work, because
                        # the parents for ZZ[x]-with-NTL and
                        # ZZ[x]-with-FLINT are equal, and the coercion model
                        # believes this means that both coercions are valid;
                        # but we'll probably change that in the
                        # coercion model, at which point this code will
                        # become useful.
                        if self._implementation_names == ('NTL',):
                            return False
                    f = base_ring.coerce_map_from(P.base_ring())
                    if f is not None:
                        from sage.rings.homset import RingHomset
                        from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
                        return PolynomialRingHomomorphism_from_base(RingHomset(P, self), f)
        except AttributeError:
            pass

        # Last, we consider multivariate polynomial rings:
        from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
        if is_MPolynomialRing(P) and self.variable_name() in P.variable_names():
            P_ = P.remove_var(self.variable_name())
            return self.base_ring()!=P_ and self.base_ring().has_coerce_map_from(P_)

    def _magma_init_(self, magma):
        """
        Used in converting this ring to the corresponding ring in MAGMA.

        EXAMPLES::

            sage: R = QQ['y']
            sage: R._magma_init_(magma)                     # optional - magma
            'SageCreateWithNames(PolynomialRing(_sage_ref...),["y"])'
            sage: S = magma(R)                              # optional - magma
            sage: print S                                   # optional - magma
            Univariate Polynomial Ring in y over Rational Field
            sage: S.1                                       # optional - magma
            y
            sage: magma(PolynomialRing(GF(7), 'x'))         # optional - magma
            Univariate Polynomial Ring in x over GF(7)
            sage: magma(PolynomialRing(GF(49,'a'), 'x'))    # optional - magma
            Univariate Polynomial Ring in x over GF(7^2)
            sage: magma(PolynomialRing(PolynomialRing(ZZ,'w'), 'x')) # optional - magma
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in w over Integer Ring

        Watch out, Magma has different semantics from Sage, i.e., in Magma
        there is a unique univariate polynomial ring, and the variable name
        has no intrinsic meaning (it only impacts printing), so can't be
        reliably set because of caching.

        ::

            sage: m = Magma()            # new magma session; optional - magma
            sage: m(QQ['w'])                                # optional - magma
            Univariate Polynomial Ring in w over Rational Field
            sage: m(QQ['x'])                                # optional - magma
            Univariate Polynomial Ring in x over Rational Field
            sage: m(QQ['w'])   # same magma object, now prints as x; optional - magma
            Univariate Polynomial Ring in x over Rational Field

        A nested example over a Givaro finite field::

            sage: k.<a> = GF(9)
            sage: R.<x> = k[]
            sage: magma(a^2*x^3 + (a+1)*x + a)              # optional - magma
            a^2*x^3 + a^2*x + a
        """
        B = magma(self.base_ring())
        Bref = B._ref()
        s = 'PolynomialRing(%s)'%(Bref)
        return magma._with_names(s, self.variable_names())

    def _gap_(self, G=None):
        """
        Used in converting this ring to the corresponding ring in GAP.

        EXAMPLES::

            sage: R.<z> = ZZ[]
            sage: gap(R)
            PolynomialRing( Integers, ["z"] )
            sage: gap(z^2 + z)
            z^2+z
        """
        if G is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
        R = G(self._gap_init_())
        v = self.variable_name()
        G.eval('%s := IndeterminatesOfPolynomialRing(%s)[1]'%(v, R.name()))
        return R

    def _gap_init_(self):
        return 'PolynomialRing(%s, ["%s"])'%(self.base_ring()._gap_init_(), self.variable_name())

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: sage_input(GF(5)['x']['y'], verify=True)
            # Verified
            GF(5)['x']['y']
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ZZ['z']._sage_input_(SageInputBuilder(), False)
            {constr_parent: {subscr: {atomic:ZZ}[{atomic:'z'}]} with gens: ('z',)}
        """
        base = sib(self.base_ring())
        sie = base[self.variable_name()]
        gens_syntax = sib.empty_subscript(base)
        return sib.parent_with_gens(self, sie, self.variable_names(), 'R',
                                    gens_syntax=gens_syntax)

    def _macaulay2_(self, m2=None):
        """
        EXAMPLES::

            sage: R = QQ['x']
            sage: macaulay2(R) # optional - macaulay2
            QQ[x, Degrees => {1}, Heft => {1}, MonomialOrder => {MonomialSize => 32}, DegreeRank => 1]
                                                                {GRevLex => {1}    }
                                                                {Position => Up    }
        """
        if m2 is None:
            import sage.interfaces.macaulay2
            m2 = sage.interfaces.macaulay2.macaulay2
        base_ring = m2( self.base_ring() )
        var = self.gen()
        return m2("%s[symbol %s]"%(base_ring.name(), var))


    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # all that is needed is that elements of the base ring
            # of the polynomial ring canonically coerce into codomain.
            # Since poly rings are free, any image of the gen
            # determines a homomorphism
            codomain.coerce(self.base_ring().one_element())
        except TypeError:
            return False
        return True

#    Polynomial rings should be unique parents. Hence,
#    no need for __cmp__. Or actually, having a __cmp__
#    method that identifies a dense with a sparse ring
#    is a bad bad idea!
#    def __cmp__(left, right):
#        c = cmp(type(left),type(right))
#        if c: return c
#        return cmp((left.base_ring(), left.variable_name()), (right.base_ring(), right.variable_name()))

    def __hash__(self):
        # should be faster than just relying on the string representation
        try:
            return self._cached_hash
        except AttributeError:
            pass
        h = self._cached_hash = hash((self.base_ring(),self.variable_name()))
        return h

    def _repr_(self):
        try:
            return self._cached_repr
        except AttributeError:
            pass
        s = "Univariate Polynomial Ring in %s over %s"%(
                self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = "Sparse " + s
        self._cached_repr = s
        return s

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: S.<alpha12>=ZZ[]
            sage: latex(S)
            \Bold{Z}[\alpha_{12}]
        """
        return "%s[%s]"%(latex.latex(self.base_ring()), self.latex_variable_names()[0])

    def base_extend(self, R):
        """
        Return the base extension of this polynomial ring to R.

        EXAMPLES::

            sage: R.<x> = RR[]; R
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: R.base_extend(CC)
            Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: R.base_extend(QQ)
            Traceback (most recent call last):
            ...
            TypeError: no such base extension
            sage: R.change_ring(QQ)
            Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if R.has_coerce_map_from(self.base_ring()):
            return PolynomialRing(R, names=self.variable_name(), sparse=self.is_sparse())
        else:
            raise TypeError("no such base extension")

    def change_ring(self, R):
        """
        Return the polynomial ring in the same variable as self over R.

        EXAMPLES::

            sage: R.<ZZZ> = RealIntervalField() []; R
            Univariate Polynomial Ring in ZZZ over Real Interval Field with 53 bits of precision
            sage: R.change_ring(GF(19^2,'b'))
            Univariate Polynomial Ring in ZZZ over Finite Field in b of size 19^2
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        return PolynomialRing(R, names=self.variable_name(), sparse=self.is_sparse())

    def change_var(self, var):
        r"""
        Return the polynomial ring in variable var over the same base
        ring.

        EXAMPLES::

            sage: R.<x> = ZZ[]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: R.change_var('y')
            Univariate Polynomial Ring in y over Integer Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        return PolynomialRing(self.base_ring(), names = var, sparse=self.is_sparse())

    def extend_variables(self, added_names, order = 'degrevlex'):
        r"""
        Returns a multivariate polynomial ring with the same base ring but
        with added_names as additional variables.

        EXAMPLES::

            sage: R.<x> = ZZ[]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: R.extend_variables('y, z')
            Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: R.extend_variables(('y', 'z'))
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if isinstance(added_names, str):
            added_names = added_names.split(',')
        return PolynomialRing(self.base_ring(), names = self.variable_names() + tuple(added_names), order = order)

    def variable_names_recursive(self, depth=sage.rings.infinity.infinity):
        r"""
        Returns the list of variable names of this and its base rings, as if
        it were a single multi-variate polynomial.

        EXAMPLES::

            sage: R = QQ['x']['y']['z']
            sage: R.variable_names_recursive()
            ('x', 'y', 'z')
            sage: R.variable_names_recursive(2)
            ('y', 'z')
        """
        if depth <= 0:
            return ()
        elif depth == 1:
            return self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
               return self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                return my_vars

    def _mpoly_base_ring(self, variables=None):
        r"""
        Returns the base ring if this is viewed as a polynomial ring over
        ``variables``. See also
        Polynomial._mpoly_dict_recursive
        """
        if variables is None:
            variables = self.variable_names_recursive()
        variables = list(variables)
        var = self.variable_name()
        if not var in variables:
            return self
        else:
            try:
                return self.base_ring()._mpoly_base_ring(variables[:variables.index(var)])
            except AttributeError:
                return self.base_ring()

    def characteristic(self):
        """
        Return the characteristic of this polynomial ring, which is the
        same as that of its base ring.

        EXAMPLES::

            sage: R.<ZZZ> = RealIntervalField() []; R
            Univariate Polynomial Ring in ZZZ over Real Interval Field with 53 bits of precision
            sage: R.characteristic()
            0
            sage: S = R.change_ring(GF(19^2,'b')); S
            Univariate Polynomial Ring in ZZZ over Finite Field in b of size 19^2
            sage: S.characteristic()
            19
        """
        return self.base_ring().characteristic()

    def cyclotomic_polynomial(self, n):
        """
        Return the nth cyclotomic polynomial as a polynomial in this
        polynomial ring. For details of the implementation, see the
        documentation for
        :func:`sage.rings.polynomial.cyclotomic.cyclotomic_coeffs`.

        EXAMPLES::

            sage: R = ZZ['x']
            sage: R.cyclotomic_polynomial(8)
            x^4 + 1
            sage: R.cyclotomic_polynomial(12)
            x^4 - x^2 + 1
            sage: S = PolynomialRing(FiniteField(7), 'x')
            sage: S.cyclotomic_polynomial(12)
            x^4 + 6*x^2 + 1
            sage: S.cyclotomic_polynomial(1)
            x + 6

        TESTS:

        Make sure it agrees with other systems for the trivial case::

            sage: ZZ['x'].cyclotomic_polynomial(1)
            x - 1
            sage: gp('polcyclo(1)')
            x - 1
        """
        if n <= 0:
            raise ArithmeticError("n=%s must be positive"%n)
        elif n == 1:
            return self.gen() - 1
        else:
            return self(cyclotomic.cyclotomic_coeffs(n), check=True)

    def gen(self, n=0):
        """
        Return the indeterminate generator of this polynomial ring.

        EXAMPLES::

            sage: R.<abc> = Integers(8)[]; R
            Univariate Polynomial Ring in abc over Ring of integers modulo 8
            sage: t = R.gen(); t
            abc
            sage: t.is_gen()
            True

        An identical generator is always returned.

        ::

            sage: t is R.gen()
            True
        """
        if n != 0:
            raise IndexError("generator n not defined")
        return self.__generator

    def gens_dict(self):
        """
        Returns a dictionary whose keys are the variable names of this
        ring as strings and whose values are the corresponding
        generators.

        EXAMPLES::

            sage: R.<x> = RR[]
            sage: R.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def parameter(self):
        """
        Return the generator of this polynomial ring.

        This is the same as ``self.gen()``.
        """
        return self.gen()

    def is_finite(self):
        """
        Return False since polynomial rings are not finite (unless the base
        ring is 0.)

        EXAMPLES::

            sage: R = Integers(1)['x']
            sage: R.is_finite()
            True
            sage: R = GF(7)['x']
            sage: R.is_finite()
            False
            sage: R['x']['y'].is_finite()
            False
        """
        R = self.base_ring()
        if R.is_finite() and R.order() == 1:
            return True
        return False

    def is_exact(self):
        return self.base_ring().is_exact()

    def is_field(self, proof = True):
        """
        Return False, since polynomial rings are never fields.

        EXAMPLES::

            sage: R.<z> = Integers(2)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 2 (using NTL)
            sage: R.is_field()
            False
        """
        return False

    def is_sparse(self):
        """
        Return true if elements of this polynomial ring have a sparse
        representation.

        EXAMPLES::

            sage: R.<z> = Integers(8)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 8
            sage: R.is_sparse()
            False
            sage: R.<W> = PolynomialRing(QQ, sparse=True); R
            Sparse Univariate Polynomial Ring in W over Rational Field
            sage: R.is_sparse()
            True
        """
        return self.__is_sparse

    def krull_dimension(self):
        """
        Return the Krull dimension of this polynomial ring, which is one
        more than the Krull dimension of the base ring.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.krull_dimension()
            1
            sage: R.<z> = GF(9,'a')[]; R
            Univariate Polynomial Ring in z over Finite Field in a of size 3^2
            sage: R.krull_dimension()
            1
            sage: S.<t> = R[]
            sage: S.krull_dimension()
            2
            sage: for n in range(10):
            ...    S = PolynomialRing(S,'w')
            sage: S.krull_dimension()
            12
        """
        return self.base_ring().krull_dimension() + 1

    def ngens(self):
        """
        Return the number of generators of this polynomial ring, which is 1
        since it is a univariate polynomial ring.

        EXAMPLES::

            sage: R.<z> = Integers(8)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 8
            sage: R.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, *args, **kwds):
        r"""
        Return a random polynomial.

        INPUT:

        -  ``degree`` - Integer with degree (default: 2)
           or a tuple of integers with minimum and maximum degrees

        -  ``*args, **kwds`` - Passed on to the ``random_element`` method for
           the base ring

        OUTPUT:

        -  Polynomial such that the coefficients of `x^i`, for `i` up to
           ``degree``, are random elements from the base ring, randomized
           subject to the arguments ``*args`` and ``**kwds``

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R.random_element(10, 5,10)
            9*x^10 + 8*x^9 + 6*x^8 + 8*x^7 + 8*x^6 + 9*x^5 + 8*x^4 + 8*x^3 + 6*x^2 + 8*x + 8
            sage: R.random_element(6)
            x^6 - 3*x^5 - x^4 + x^3 - x^2 + x + 1
            sage: R.random_element(6)
            -2*x^5 + 2*x^4 - 3*x^3 + 1
            sage: R.random_element(6)
            x^4 - x^3 + x - 2

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree::

            sage: R.random_element(degree=(0,8))
            2*x^7 - x^5 + 4*x^4 - 5*x^3 + x^2 + 14*x - 1
            sage: R.random_element(degree=(0,8))
            -2*x^3 + x^2 + x + 4

        TESTS::

            sage: R.random_element(degree=[5])
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)

            sage: R.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: minimum degree must be less or equal than maximum degree
        """
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            degree = randint(*degree)
        R = self.base_ring()
        return self([R.random_element(*args, **kwds) for _ in xrange(degree+1)])

    def _monics_degree( self, of_degree ):
        """
        Refer to monics() for full documentation.
        """
        base = self.base_ring()
        for coeffs in sage.misc.mrange.xmrange_iter([[base.one_element()]]+[base]*of_degree):
            # Each iteration returns a *new* list!
            # safe to mutate the return
            coeffs.reverse()
            yield self(coeffs)

    def _monics_max( self, max_degree ):
        """
        Refer to monics() for full documentation.
        """
        for degree in xrange(max_degree + 1):
            for m in self._monics_degree( degree ):
                yield m

    def _polys_degree( self, of_degree ):
        """
        Refer to polynomials() for full documentation.
        """
        base = self.base_ring()
        base0 = base.zero_element()
        for leading_coeff in base:
            if leading_coeff != base0:
                for lt1 in sage.misc.mrange.xmrange_iter([base]*(of_degree)):
                    # Each iteration returns a *new* list!
                    # safe to mutate the return
                    coeffs = [leading_coeff] + lt1
                    coeffs.reverse()
                    yield self(coeffs)

    def _polys_max( self, max_degree ):
        """
        Refer to polynomials() for full documentation.
        """
        base = self.base_ring()
        for coeffs in sage.misc.mrange.xmrange_iter([base]*(max_degree+1)):
            # Each iteration returns a *new* list!
            # safe to mutate the return
            coeffs.reverse()
            yield self(coeffs)

    def karatsuba_threshold(self):
        """
        Return the Karatsuba threshold used for this ring by the method
        _mul_karatsuba to fall back to the schoolbook algorithm.

        EXAMPLES::

            sage: K = QQ['x']
            sage: K.karatsuba_threshold()
            8
            sage: K = QQ['x']['y']
            sage: K.karatsuba_threshold()
            0
        """
        return self._Karatsuba_threshold

    def set_karatsuba_threshold(self, Karatsuba_threshold):
        """
        Changes the default threshold for this ring in the method _mul_karatsuba
        to fall back to the schoolbook algorithm.

        .. warning::

           This method may have a negative performance impact in polynomial
           arithmetic. So use it at your own risk.

        EXAMPLES::

            sage: K = QQ['x']
            sage: K.karatsuba_threshold()
            8
            sage: K.set_karatsuba_threshold(0)
            sage: K.karatsuba_threshold()
            0
        """
        self._Karatsuba_threshold = ZZ_sage(Karatsuba_threshold)


    def polynomials( self, of_degree = None, max_degree = None ):
        """
        Return an iterator over the polynomials of specified degree.

        INPUT: Pass exactly one of:


        -  ``max_degree`` - an int; the iterator will generate
           all polynomials which have degree less than or equal to
           max_degree

        -  ``of_degree`` - an int; the iterator will generate
           all polynomials which have degree of_degree


        OUTPUT: an iterator

        EXAMPLES::

            sage: P = PolynomialRing(GF(3),'y')
            sage: for p in P.polynomials( of_degree = 2 ): print p
            y^2
            y^2 + 1
            y^2 + 2
            y^2 + y
            y^2 + y + 1
            y^2 + y + 2
            y^2 + 2*y
            y^2 + 2*y + 1
            y^2 + 2*y + 2
            2*y^2
            2*y^2 + 1
            2*y^2 + 2
            2*y^2 + y
            2*y^2 + y + 1
            2*y^2 + y + 2
            2*y^2 + 2*y
            2*y^2 + 2*y + 1
            2*y^2 + 2*y + 2
            sage: for p in P.polynomials( max_degree = 1 ): print p
            0
            1
            2
            y
            y + 1
            y + 2
            2*y
            2*y + 1
            2*y + 2
            sage: for p in P.polynomials( max_degree = 1, of_degree = 3 ): print p
            Traceback (most recent call last):
            ...
            ValueError: you should pass exactly one of of_degree and max_degree

        AUTHORS:

        - Joel B. Mohler
        """

        if self.base_ring().order() is sage.rings.infinity.infinity:
            raise NotImplementedError
        if of_degree is not None and max_degree is None:
            return self._polys_degree( of_degree )
        if max_degree is not None and of_degree is None:
            return self._polys_max( max_degree )
        raise ValueError("you should pass exactly one of of_degree and max_degree")

    def monics( self, of_degree = None, max_degree = None ):
        """
        Return an iterator over the monic polynomials of specified degree.

        INPUT: Pass exactly one of:


        -  ``max_degree`` - an int; the iterator will generate
           all monic polynomials which have degree less than or equal to
           max_degree

        -  ``of_degree`` - an int; the iterator will generate
           all monic polynomials which have degree of_degree


        OUTPUT: an iterator

        EXAMPLES::

            sage: P = PolynomialRing(GF(4,'a'),'y')
            sage: for p in P.monics( of_degree = 2 ): print p
            y^2
            y^2 + a
            y^2 + a + 1
            y^2 + 1
            y^2 + a*y
            y^2 + a*y + a
            y^2 + a*y + a + 1
            y^2 + a*y + 1
            y^2 + (a + 1)*y
            y^2 + (a + 1)*y + a
            y^2 + (a + 1)*y + a + 1
            y^2 + (a + 1)*y + 1
            y^2 + y
            y^2 + y + a
            y^2 + y + a + 1
            y^2 + y + 1
            sage: for p in P.monics( max_degree = 1 ): print p
            1
            y
            y + a
            y + a + 1
            y + 1
            sage: for p in P.monics( max_degree = 1, of_degree = 3 ): print p
            Traceback (most recent call last):
            ...
            ValueError: you should pass exactly one of of_degree and max_degree

        AUTHORS:

        - Joel B. Mohler
        """

        if self.base_ring().order() is sage.rings.infinity.infinity:
            raise NotImplementedError
        if of_degree is not None and max_degree is None:
            return self._monics_degree( of_degree )
        if max_degree is not None and of_degree is None:
            return self._monics_max( max_degree )
        raise ValueError("you should pass exactly one of of_degree and max_degree")

class PolynomialRing_commutative(PolynomialRing_general, commutative_algebra.CommutativeAlgebra):
    """
    Univariate polynomial ring over a commutative ring.
    """
    def __init__(self, base_ring, name=None, sparse=False, element_class=None, category=None):
        if base_ring not in _CommutativeRings:
            raise TypeError("Base ring %s must be a commutative ring."%repr(base_ring))
        # We trust that, if a category is given, that it is useful.
        if category is None:
            category = polynomial_default_category(base_ring.category(),False)
        PolynomialRing_general.__init__(self, base_ring, name=name,
                sparse=sparse, element_class=element_class, category=category)

    def quotient_by_principal_ideal(self, f, names=None):
        """
        Return the quotient of this polynomial ring by the principal
        ideal (generated by) `f`.

        INPUT:

        - ``f`` - either a polynomial in ``self``, or a principal
          ideal of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: I = (x^2-1)*R
            sage: R.quotient_by_principal_ideal(I)
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 - 1

        The same example, using the polynomial instead of the ideal,
        and customizing the variable name::

            sage: R.<x> = QQ[]
            sage: R.quotient_by_principal_ideal(x^2-1, names=('foo',))
            Univariate Quotient Polynomial Ring in foo over Rational Field with modulus x^2 - 1

        TESTS:

        Quotienting by the zero ideal returns ``self`` (:trac:`5978`)::

            sage: R = QQ['x']
            sage: R.quotient_by_principal_ideal(R.zero_ideal()) is R
            True
            sage: R.quotient_by_principal_ideal(0) is R
            True
        """
        from sage.rings.ideal import Ideal
        I = Ideal(f)
        if I.is_zero():
            return self
        f = I.gen()
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing
        return PolynomialQuotientRing(self, f, names)



class PolynomialRing_integral_domain(PolynomialRing_commutative, integral_domain.IntegralDomain):
    def __init__(self, base_ring, name="x", sparse=False, implementation=None,
            element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x'); R
            Univariate Polynomial Ring in x over Integer Ring
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>

            sage: R = PRing(ZZ, 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>
        """
        self._implementation_repr = ''
        if not element_class:
            if is_IntegerRing(base_ring) and not sparse:
                if implementation == 'NTL':
                    from sage.rings.polynomial.polynomial_integer_dense_ntl \
                            import Polynomial_integer_dense_ntl
                    element_class = Polynomial_integer_dense_ntl
                    self._implementation_names = ('NTL',)
                    self._implementation_repr = ' (using NTL)'
                elif implementation == 'FLINT' or implementation is None:
                    from sage.rings.polynomial.polynomial_integer_dense_flint \
                            import Polynomial_integer_dense_flint
                    element_class = Polynomial_integer_dense_flint
                    self._implementation_names = (None, 'FLINT')
                else:
                    raise ValueError("Unknown implementation %s for ZZ[x]"%implementation)
        PolynomialRing_commutative.__init__(self, base_ring, name=name,
                sparse=sparse, element_class=element_class)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
        """
        s = PolynomialRing_commutative._repr_(self)
        return s + self._implementation_repr


class PolynomialRing_field(PolynomialRing_integral_domain,
                           PolynomialRing_singular_repr,
                           principal_ideal_domain.PrincipalIdealDomain,
                           ):
    def __init__(self, base_ring, name="x", sparse=False, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_field as PRing
            sage: R = PRing(QQ, 'x'); R
            Univariate Polynomial Ring in x over Rational Field
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>
            sage: R = PRing(QQ, 'x', sparse=True); R
            Sparse Univariate Polynomial Ring in x over Rational Field
            sage: type(R.gen())
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_sparse_field'>
            sage: R = PRing(CC, 'x'); R
            Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: type(R.gen())
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>

            #Demonstrate that Trac #8762 is fixed
            sage: R.<x> = PolynomialRing(GF(next_prime(10^20)), sparse=True)
            sage: x^(10^20) # this should be fast
            x^100000000000000000000
        """
        from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular

        if not element_class:
            if sparse:
                element_class = polynomial_element_generic.Polynomial_generic_sparse_field
            elif isinstance(base_ring, rational_field.RationalField):
                from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
                element_class = Polynomial_rational_flint
            elif is_RealField(base_ring):
                element_class = PolynomialRealDense
            else:
                element_class = polynomial_element_generic.Polynomial_generic_dense_field

        PolynomialRing_integral_domain.__init__(self, base_ring, name=name, sparse=sparse, element_class=element_class)

        self._has_singular = can_convert_to_singular(self)

    def _ideal_class_(self, n=0):
        """
        Returns the class representing ideals in univariate polynomial rings over fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: R._ideal_class_()
            <class 'sage.rings.polynomial.ideal.Ideal_1poly_field'>
        """
        from sage.rings.polynomial.ideal import Ideal_1poly_field
        return Ideal_1poly_field

    def divided_difference(self, points, full_table=False):
        """
        Return the Newton divided-difference coefficients of the `n`-th
        Lagrange interpolation polynomial of ``points``.

        If ``points`` are `n+1` distinct points
        `(x_0, f(x_0)), (x_1, f(x_1)), \dots, (x_n, f(x_n))`, then `P_n(x)`
        is the `n`-th Lagrange interpolation polynomial of `f(x)` that
        passes through the points `(x_i, f(x_i))`. This method returns
        the coefficients `F_{i,i}` such that

        .. MATH::

            P_n(x) = \sum_{i=0}^n F_{i,i} \prod_{j=0}^{i-1} (x - x_j)


        INPUT:

        - ``points`` -- a list of tuples
          `(x_0, f(x_0)), (x_1, f(x_1)), \dots, (x_n, f(x_n))` where each
          `x_i \\neq x_j` for `i \\neq j`

        - ``full_table`` -- (default: ``False``) if ``True`` then return the
          full divided-difference table; if ``False`` then only return
          entries along the main diagonal. The entries along the main
          diagonal are the Newton divided-difference coefficients `F_{i,i}`.


        OUTPUT:

        - The Newton divided-difference coefficients of the `n`-th Lagrange
          interpolation polynomial that passes through the points in
          ``points``.


        EXAMPLES:

        Only return the divided-difference coefficients `F_{i,i}`. This
        example is taken from Example 1, p.121 of [BF05]_::

            sage: points = [(1.0, 0.7651977), (1.3, 0.6200860), (1.6, 0.4554022), (1.9, 0.2818186), (2.2, 0.1103623)]
            sage: R = PolynomialRing(QQ, "x")
            sage: R.divided_difference(points)
            <BLANKLINE>
            [0.765197700000000,
            -0.483705666666666,
            -0.108733888888889,
            0.0658783950617283,
            0.00182510288066044]

        Now return the full divided-difference table::

            sage: points = [(1.0, 0.7651977), (1.3, 0.6200860), (1.6, 0.4554022), (1.9, 0.2818186), (2.2, 0.1103623)]
            sage: R = PolynomialRing(QQ, "x")
            sage: R.divided_difference(points, full_table=True)
            <BLANKLINE>
            [[0.765197700000000],
            [0.620086000000000, -0.483705666666666],
            [0.455402200000000, -0.548946000000000, -0.108733888888889],
            [0.281818600000000,
            -0.578612000000000,
            -0.0494433333333339,
            0.0658783950617283],
            [0.110362300000000,
            -0.571520999999999,
            0.0118183333333349,
            0.0680685185185209,
            0.00182510288066044]]

        The following example is taken from Example 4.12, p.225 of [MF99]_::

            sage: points = [(1, -3), (2, 0), (3, 15), (4, 48), (5, 105), (6, 192)]
            sage: R = PolynomialRing(RR, "x")
            sage: R.divided_difference(points)
            [-3, 3, 6, 1, 0, 0]
            sage: R.divided_difference(points, full_table=True)
            <BLANKLINE>
            [[-3],
            [0, 3],
            [15, 15, 6],
            [48, 33, 9, 1],
            [105, 57, 12, 1, 0],
            [192, 87, 15, 1, 0, 0]]


        REFERENCES:

        .. [MF99] J.H. Mathews and K.D. Fink. *Numerical Methods Using MATLAB*.
          3rd edition, Prentice-Hall, 1999.
        """
        n = len(points)
        F = [[points[i][1]] for i in xrange(n)]
        for i in xrange(1, n):
            for j in xrange(1, i+1):
                numer = F[i][j-1] - F[i-1][j-1]
                denom = points[i][0] - points[i-j][0]
                F[i].append(numer / denom)
        if full_table:
            return F
        else:
            return [F[i][i] for i in xrange(n)]

    def lagrange_polynomial(self, points, algorithm="divided_difference", previous_row=None):
        """
        Return the Lagrange interpolation polynomial in ``self`` associated to
        the given list of points.

        Given a list of points, i.e. tuples of elements of ``self``'s base
        ring, this function returns the interpolation polynomial in the
        Lagrange form.


        INPUT:

        - ``points`` -- a list of tuples representing points  through which
          the polynomial returned by this function must pass.

        - ``algorithm`` -- (default: ``'divided_difference'``) the available
          values for this option are ``'divided_difference'`` and ``neville``.

          - If ``algorithm='divided_difference'`` then use the method of
            divided difference.

          - If ``algorithm='neville'`` then adapt Neville's method as described
            on page 144 of [BF05]_ to recursively generate the Lagrange
            interpolation polynomial. Neville's method generates
            a table of approximating polynomials, where the last row of that
            table contains the `n`-th Lagrange interpolation polynomial. The
            adaptation implemented by this method is to only generate the
            last row of this table, instead of the full table itself.
            Generating the full table can be memory inefficient.

        - ``previous_row`` -- (default: ``None``) This option is only relevant
          if used together with ``algorithm='neville'``. If provided, this
          should be the last row of the table resulting from a previous use of
          Neville's method. If such a row is passed in, then ``points`` should
          consist of both previous and new interpolating points. Neville's
          method will then use that last row and the interpolating points to
          generate a new row which contains a better Lagrange interpolation
          polynomial.


        EXAMPLES:

        By default, we use the method of divided-difference::

            sage: R = PolynomialRing(QQ, 'x')
            sage: f = R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)]); f
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1
            sage: f(0)
            1
            sage: f(2)
            2
            sage: f(3)
            -2
            sage: f(-4)
            9
            sage: R = PolynomialRing(GF(2**3,'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: f = R.lagrange_polynomial([(a^2+a,a),(a,1),(a^2,a^2+a+1)]); f
            a^2*x^2 + a^2*x + a^2
            sage: f(a^2+a)
            a
            sage: f(a)
            1
            sage: f(a^2)
            a^2 + a + 1

        Now use a memory efficient version of Neville's method::

            sage: R = PolynomialRing(QQ, 'x')
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm="neville")
            <BLANKLINE>
            [9,
            -11/7*x + 19/7,
            -17/42*x^2 - 83/42*x + 53/7,
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1]
            sage: R = PolynomialRing(GF(2**3,'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: R.lagrange_polynomial([(a^2+a,a),(a,1),(a^2,a^2+a+1)], algorithm="neville")
            [a^2 + a + 1, x + a + 1, a^2*x^2 + a^2*x + a^2]

        Repeated use of Neville's method to get better Lagrange interpolation
        polynomials::

            sage: R = PolynomialRing(QQ, 'x')
            sage: p = R.lagrange_polynomial([(0,1),(2,2)], algorithm="neville")
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm="neville", previous_row=p)[-1]
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1
            sage: R = PolynomialRing(GF(2**3,'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: p = R.lagrange_polynomial([(a^2+a,a),(a,1)], algorithm="neville")
            sage: R.lagrange_polynomial([(a^2+a,a),(a,1),(a^2,a^2+a+1)], algorithm="neville", previous_row=p)[-1]
            a^2*x^2 + a^2*x + a^2


        TESTS:

        The value for ``algorithm`` must be either ``'divided_difference'`` (by
        default it is), or ``'neville'``::

            sage: R = PolynomialRing(QQ, "x")
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm="abc")
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm="divided difference")
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm="")
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'

        Make sure that ticket #10304 is fixed. The return value should always
        be an element of ``self`` in the case of ``divided_difference``, or
        a list of elements of ``self`` in the case of ``neville``. ::

            sage: R = PolynomialRing(QQ, "x")
            sage: R.lagrange_polynomial([]).parent() == R
            True
            sage: R.lagrange_polynomial([(2, 3)]).parent() == R
            True
            sage: row = R.lagrange_polynomial([], algorithm='neville')
            sage: all(poly.parent() == R for poly in row)
            True
            sage: row = R.lagrange_polynomial([(2, 3)], algorithm='neville')
            sage: all(poly.parent() == R for poly in row)
            True


        REFERENCES:

        .. [BF05] R.L. Burden and J.D. Faires. *Numerical Analysis*.
          Thomson Brooks/Cole, 8th edition, 2005.
        """
        var = self.gen()

        # use the method of divided-difference
        if algorithm == "divided_difference":
            # Evaluate in nested form, similar to Horner's method. This is
            # more efficient than evaluation using the definition of
            # Lagrange interpolation polynomial by means of divided
            # difference.
            n = len(points)
            if n == 0:
                return self.zero()

            F = self.divided_difference(points)
            P = self(F[n-1])
            for i in xrange(n-2, -1, -1):
                P *= (var - points[i][0])
                P += F[i]
            return P

            # Evaluate using the definition of Lagrange interpolation
            # polynomial by means of divided difference. This is slow
            # compared to that above, which is in nested form.
#             P = 0
#             for i in xrange(n):
#                 prod = 1
#                 for j in xrange(i):
#                     prod *= (var - points[j][0])
#                 P += (F[i] * prod)
#             return P

        # using Neville's method for recursively generating the
        # Lagrange interpolation polynomial
        elif algorithm == "neville":
            if previous_row is None:
                previous_row = []
            N = len(points)
            M = len(previous_row)
            # During the computation, P keeps track of the previous row,
            # and Q keeps track of the current row
            P = previous_row + [None] * (N - M) # use results of previous computation if available
            Q = [None] * N
            for i in xrange(M, N):
                Q[0] = self(points[i][1]) # start populating the current row
                for j in xrange(1, 1 + i):
                    numer = (var - points[i - j][0]) * Q[j - 1] - (var - points[i][0]) * P[j - 1]
                    denom = points[i][0] - points[i - j][0]
                    Q[j] = numer / denom
                P, Q = Q, P # the current row is complete, reuse the old P to hold the next row
            return P # return the last row in the Neville table

#        # use the definition of Lagrange interpolation polynomial
#        elif algorithm == "definition":
#            def Pj(j):
#                denom = 1
#                divis = 1
#                for i in range(len(points)):
#                    if i!=j:
#                        denom *= (var          - points[i][0])
#                        divis *= (points[j][0] - points[i][0])
#            return denom/divis
#
#            P = 0
#            for j in range(len(points)):
#                P += Pj(j)*points[j][1]
#            return P

        else:
            raise ValueError("algorithm must be one of 'divided_difference' or 'neville'")

    def fraction_field(self):
        """
        Returns the fraction field of self.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: R.fraction_field()
            Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
        """
        try:
            return self._fraction_field
        except AttributeError:
            R = self.base_ring()
            p = R.characteristic()
            if p != 0 and R.is_prime_field() and 2 < p and p < 2**16:
                from sage.rings.fraction_field_FpT import FpT
                self._fraction_field = FpT(self)
            else:
                from sage.rings.fraction_field import FractionField_1poly_field
                self._fraction_field = FractionField_1poly_field(self)
            return self._fraction_field

class PolynomialRing_dense_finite_field(PolynomialRing_field):
    """
    Univariate polynomial ring over a finite field.

    EXAMPLE::

        sage: R = PolynomialRing(GF(27, 'a'), 'x')
        sage: type(R)
        <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_finite_field_with_category'>
    """
    def __init__(self, base_ring, name="x", element_class=None, implementation=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_finite_field
            sage: R = PolynomialRing_dense_finite_field(GF(5), implementation='generic')
            sage: type(R(0))
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>

            sage: S = PolynomialRing_dense_finite_field(GF(25, 'a'), implementation='NTL')
            sage: type(S(0))
            <type 'sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX'>
        """
        if implementation is None:
            implementation = "NTL"

        if implementation == "NTL":
            from sage.libs.ntl.ntl_ZZ_pEContext import ntl_ZZ_pEContext
            from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX
            from sage.rings.polynomial.polynomial_zz_pex import Polynomial_ZZ_pEX

            p=base_ring.characteristic()
            self._modulus = ntl_ZZ_pEContext(ntl_ZZ_pX(list(base_ring.polynomial()), p))
            element_class = Polynomial_ZZ_pEX

        PolynomialRing_field.__init__(self, base_ring, sparse=False, name=name,
                                      element_class=element_class)

    def irreducible_element(self, n, algorithm=None):
        """
        Construct an irreducible polynomial of degree `n`.

        INPUT:

        - ``n`` -- integer: degree of the polynomial to construct

        - ``algorithm`` -- string: algorithm to use, or ``None``

          - ``'random'``: try random polynomials until an irreducible
            one is found.

          - ``'first_lexicographic'``: try polynomials in
            lexicographic order until an irreducible one is found.

        OUTPUT:

        A monic irreducible polynomial of degree `n` in ``self``.

        EXAMPLES::

            sage: GF(5^3, 'a')['x'].irreducible_element(2)
            x^2 + (4*a^2 + a + 4)*x + 2*a^2 + 2
            sage: GF(19)['x'].irreducible_element(21, algorithm="first_lexicographic")
            x^21 + x + 5
            sage: GF(5**2, 'a')['x'].irreducible_element(17, algorithm="first_lexicographic")
            x^17 + a*x + 4*a + 3

        AUTHORS:

        - Peter Bruin (June 2013)
        - Jean-Pierre Flori (May 2014)
        """
        if n < 1:
            raise ValueError("degree must be at least 1")

        if algorithm is None or algorithm == "random":
            while True:
                f = self.gen()**n + self.random_element(n - 1)
                if f.is_irreducible():
                    return f
        elif algorithm == "first_lexicographic":
            for g in self.polynomials(max_degree=n-1):
                f = self.gen()**n + g
                if f.is_irreducible():
                    return f
        else:
            raise ValueError("no such algorithm for finding an irreducible polynomial: %s" % algorithm)

class PolynomialRing_dense_padic_ring_generic(PolynomialRing_integral_domain):
    pass

class PolynomialRing_dense_padic_field_generic(PolynomialRing_field):
    pass

class PolynomialRing_dense_padic_ring_capped_relative(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_capped_relative as PRing
            sage: R = PRing(Zp(13), name='t'); R
            Univariate Polynomial Ring in t over 13-adic Ring with capped relative precision 20
            sage: type(R.gen())
            <class 'sage.rings.polynomial.padics.polynomial_padic_capped_relative_dense.Polynomial_padic_capped_relative_dense'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.\
                    polynomial_padic_capped_relative_dense import \
                    Polynomial_padic_capped_relative_dense
            element_class = Polynomial_padic_capped_relative_dense
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_padic_ring_capped_absolute(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_capped_absolute as PRing
            sage: R = PRing(Zp(13, type='capped-abs'), name='t'); R
            Univariate Polynomial Ring in t over 13-adic Ring with capped absolute precision 20
            sage: type(R.gen())
            <class 'sage.rings.polynomial.padics.polynomial_padic_flat.Polynomial_padic_flat'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.polynomial_padic_flat import \
                    Polynomial_padic_flat
            element_class = Polynomial_padic_flat
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_padic_ring_fixed_mod(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_fixed_mod as PRing
            sage: R = PRing(Zp(13, type='fixed-mod'), name='t'); R
            Univariate Polynomial Ring in t over 13-adic Ring of fixed modulus 13^20

            sage: type(R.gen())
            <class 'sage.rings.polynomial.padics.polynomial_padic_flat.Polynomial_padic_flat'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.polynomial_padic_flat import \
                    Polynomial_padic_flat
            element_class = Polynomial_padic_flat
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_padic_ring_lazy(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_lazy as PRing
            sage: R = PRing(Zp(13, type='lazy'), name='t')
            Traceback (most recent call last):
            ...
            NotImplementedError: lazy p-adics need more work.  Sorry.

            #sage: type(R.gen())

        """
        if element_class is None:
            element_class = polynomial_element_generic.Polynomial_generic_dense
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_padic_field_capped_relative(PolynomialRing_dense_padic_field_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_field_capped_relative as PRing
            sage: R = PRing(Qp(13), name='t'); R
            Univariate Polynomial Ring in t over 13-adic Field with capped relative precision 20
            sage: type(R.gen())
            <class 'sage.rings.polynomial.padics.polynomial_padic_capped_relative_dense.Polynomial_padic_capped_relative_dense'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.\
                    polynomial_padic_capped_relative_dense import \
                    Polynomial_padic_capped_relative_dense
            element_class = Polynomial_padic_capped_relative_dense
        PolynomialRing_dense_padic_field_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_padic_field_lazy(PolynomialRing_dense_padic_field_generic):
    def __init__(self, base_ring, name=None, element_class=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_field_lazy as PRing
            sage: R = PRing(Qp(13, type='lazy'), name='t')
            Traceback (most recent call last):
            ...
            NotImplementedError: lazy p-adics need more work.  Sorry.

            #sage: type(R.gen())
        """
        if element_class is None:
            element_class = polynomial_element_generic.Polynomial_generic_dense
        PolynomialRing_dense_padic_field_generic.__init__(self, base_ring,
                name=name, element_class=element_class)

class PolynomialRing_dense_mod_n(PolynomialRing_commutative):
    def __init__(self, base_ring, name=None, element_class=None,
            implementation=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_mod_n as PRing
            sage: R = PRing(Zmod(15), 'x'); R
            Univariate Polynomial Ring in x over Ring of integers modulo 15
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

            sage: R = PRing(Zmod(15), 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Ring of integers modulo 15 (using NTL)
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_zz'>

            sage: R = PRing(Zmod(2**63*3), 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Ring of integers modulo 27670116110564327424 (using NTL)
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>

            sage: R = PRing(Zmod(2**63*3), 'x', implementation='FLINT')
            Traceback (most recent call last):
            ...
            ValueError: FLINT does not support modulus 27670116110564327424

            sage: R = PRing(Zmod(2**63*3), 'x'); R
            Univariate Polynomial Ring in x over Ring of integers modulo 27670116110564327424 (using NTL)
            sage: type(R.gen())
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>
        """
        from sage.rings.polynomial.polynomial_zmod_flint import \
                Polynomial_zmod_flint
        import sage.rings.polynomial.polynomial_modn_dense_ntl as \
                modn_dense_ntl
        self.__modulus = base_ring.order()
        if not element_class:
            if implementation is None or implementation == 'FLINT':
                import sys
                if self.__modulus < sys.maxsize:
                    element_class = Polynomial_zmod_flint
                    self._implementation_names = (None, 'FLINT')
                    self._implementation_repr = ''
                elif implementation == 'FLINT':
                    raise ValueError("FLINT does not support modulus %s"%(self.__modulus))
            if not element_class:
                self._implementation_names = ('NTL',)
                self._implementation_repr = ' (using NTL)'
                if self.__modulus < ZZ_sage(modn_dense_ntl.zz_p_max):
                    element_class = modn_dense_ntl.Polynomial_dense_modn_ntl_zz
                else:
                    element_class = modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ
        PolynomialRing_commutative.__init__(self, base_ring, name=name,
                element_class=element_class)

    @cached_method
    def modulus(self):
        """
        EXAMPLES::

            sage: R.<x> = Zmod(15)[]
            sage: R.modulus()
            15
        """
        return self.base_ring().characteristic()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
        """
        s = PolynomialRing_commutative._repr_(self)
        return s + self._implementation_repr


class PolynomialRing_dense_mod_p(PolynomialRing_dense_finite_field,
                                 PolynomialRing_dense_mod_n,
                                 PolynomialRing_singular_repr):
    def __init__(self, base_ring, name="x", implementation=None):
        """
        TESTS::

            sage: P = GF(2)['x']; P
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
            sage: type(P.gen())
            <type 'sage.rings.polynomial.polynomial_gf2x.Polynomial_GF2X'>

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_mod_p
            sage: P = PolynomialRing_dense_mod_p(GF(5), 'x'); P
            Univariate Polynomial Ring in x over Finite Field of size 5
            sage: type(P.gen())
            <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

            sage: P = PolynomialRing_dense_mod_p(GF(5), 'x', implementation='NTL'); P
            Univariate Polynomial Ring in x over Finite Field of size 5 (using NTL)
            sage: type(P.gen())
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>

            sage: P = PolynomialRing_dense_mod_p(GF(9223372036854775837), 'x')
            sage: P
            Univariate Polynomial Ring in x over Finite Field of size 9223372036854775837 (using NTL)
            sage: type(P.gen())
            <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>


        """
        from sage.rings.polynomial.polynomial_zmod_flint import \
                Polynomial_zmod_flint
        __modulus = base_ring.characteristic()
        element_class = None
        if __modulus == 2:
            import sage.rings.polynomial.polynomial_gf2x as polynomial_gf2x
            element_class = polynomial_gf2x.Polynomial_GF2X
            self._implementation_repr = ' (using NTL)'
        elif implementation is None or implementation == 'FLINT':
            import sys
            if __modulus < sys.maxsize:
                self._implementation_names = (None, 'FLINT')
                self._implementation_repr = ''
                element_class = Polynomial_zmod_flint
            elif implementation == 'FLINT':
                raise ValueError("FLINT does not support modulus %s"%(__modulus))
        if not element_class:
            from sage.rings.polynomial.polynomial_modn_dense_ntl import \
                    Polynomial_dense_mod_p
            element_class = Polynomial_dense_mod_p
            self._implementation_names = ('NTL',)
            self._implementation_repr = ' (using NTL)'
        PolynomialRing_dense_mod_n.__init__(self, base_ring, name=name,
                element_class=element_class)

        from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
        self._has_singular = can_convert_to_singular(self)

    def irreducible_element(self, n, algorithm=None):
        """
        Construct an irreducible polynomial of degree `n`.

        INPUT:

        - ``n`` -- integer: the degree of the polynomial to construct

        - ``algorithm`` -- string: algorithm to use, or ``None``.
          Currently available options are:

          - ``'adleman-lenstra'``: a variant of the Adleman--Lenstra
              algorithm as implemented in PARI.

          - ``'conway'``: look up the Conway polynomial of degree `n`
            over the field of `p` elements in the database; raise a
            ``RuntimeError`` if it is not found.

          - ``'first_lexicographic'``: return the lexicographically
            smallest irreducible polynomial of degree `n`.

          - ``'minimal_weight'``: return an irreducible polynomial of
            degree `n` with minimal number of non-zero coefficients.
            Only implemented for `p = 2`.

          - ``'random'``: try random polynomials until an irreducible
            one is found.

          If ``algorithm`` is ``None``, the Conway polynomial is used
          if it is found in the database.  If no Conway polynomial is
          found, the algorithm ``minimal_weight`` is used if `p = 2`,
          and the algorithm ``adleman-lenstra`` if `p > 2`.

        OUTPUT:

        A monic irreducible polynomial of degree `n` in ``self``.

        EXAMPLES::

            sage: GF(5)['x'].irreducible_element(2)
            x^2 + 4*x + 2
            sage: GF(5)['x'].irreducible_element(2, algorithm="adleman-lenstra")
            x^2 + x + 1

            sage: GF(2)['x'].irreducible_element(33)
            x^33 + x^13 + x^12 + x^11 + x^10 + x^8 + x^6 + x^3 + 1
            sage: GF(2)['x'].irreducible_element(33, algorithm="minimal_weight")
            x^33 + x^10 + 1

        AUTHORS:

        - Peter Bruin (June 2013)
        """
        from sage.libs.pari.all import pari
        from sage.rings.finite_rings.conway_polynomials import (conway_polynomial,
                                                                exists_conway_polynomial)
        from polynomial_gf2x import (GF2X_BuildIrred_list,
                                     GF2X_BuildSparseIrred_list,
                                     GF2X_BuildRandomIrred_list)

        p = self.characteristic()
        n = int(n)
        if n < 1:
            raise ValueError("degree must be at least 1")
        if algorithm is None:
            if exists_conway_polynomial(p, n):
                algorithm = "conway"
            elif p == 2:
                algorithm = "minimal_weight"
            else:
                algorithm = "adleman-lenstra"

        if algorithm == "adleman-lenstra":
            return self(pari(p).ffinit(n))
        elif algorithm == "conway":
            return self(conway_polynomial(p, n))
        elif algorithm == "first_lexicographic":
            if p == 2:
                return self(GF2X_BuildIrred_list(n))
            else:
                # Fallback to PolynomialRing_dense_finite_field.irreducible_element
                pass
        elif algorithm == "minimal_weight":
            if p == 2:
                return self(GF2X_BuildSparseIrred_list(n))
            else:
                raise NotImplementedError("'minimal_weight' option only implemented for p = 2")
        elif algorithm == "random":
            if p == 2:
                return self(GF2X_BuildRandomIrred_list(n))
            else:
                pass

        # No suitable algorithm found, try algorithms from the base class.
        return PolynomialRing_dense_finite_field.irreducible_element(self, n, algorithm)

def polygen(ring_or_element, name="x"):
    """
    Return a polynomial indeterminate.

    INPUT:

    - polygen(base_ring, name="x")

    - polygen(ring_element, name="x")

    If the first input is a ring, return a polynomial generator over
    that ring. If it is a ring element, return a polynomial generator
    over the parent of the element.

    EXAMPLES::

        sage: z = polygen(QQ,'z')
        sage: z^3 + z +1
        z^3 + z + 1
        sage: parent(z)
        Univariate Polynomial Ring in z over Rational Field

    .. note::

       If you give a list or comma separated string to polygen, you'll
       get a tuple of indeterminates, exactly as if you called
       polygens.
    """
    if ring_element.is_RingElement(ring_or_element):
        base_ring = ring_or_element.parent()
    elif ring.is_Ring(ring_or_element):
        base_ring = ring_or_element
    else:
        raise TypeError("input must be a ring or ring element")
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    t = PolynomialRing(base_ring, name)
    if t.ngens() > 1:
        return t.gens()
    return t.gen()

def polygens(base_ring, names="x"):
    """
    Return indeterminates over the given base ring with the given
    names.

    EXAMPLES::

        sage: x,y,z = polygens(QQ,'x,y,z')
        sage: (x+y+z)^2
        x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2
        sage: parent(x)
        Multivariate Polynomial Ring in x, y, z over Rational Field
        sage: t = polygens(QQ,['x','yz','abc'])
        sage: t
        (x, yz, abc)
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(base_ring, names).gens()
