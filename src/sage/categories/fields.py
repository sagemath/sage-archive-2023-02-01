r"""
Fields
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                2012      Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton, Category_contains_method_by_parent_class
from sage.categories.euclidean_domains import EuclideanDomains
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.division_rings import DivisionRings

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
import sage.rings.ring
from sage.structure.element import coerce_binop

class Fields(Category_singleton):
    """
    The category of (commutative) fields, i.e. commutative rings where
    all non-zero elements have multiplicative inverses

    EXAMPLES::

        sage: K = Fields()
        sage: K
        Category of fields
        sage: Fields().super_categories()
        [Category of euclidean domains, Category of unique factorization domains, Category of division rings]

        sage: K(IntegerRing())
        Rational Field
        sage: K(PolynomialRing(GF(3), 'x'))
        Fraction Field of Univariate Polynomial Ring in x over
        Finite Field of size 3
        sage: K(RealField())
        Real Field with 53 bits of precision

    TESTS::

        sage: TestSuite(Fields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Fields().super_categories()
            [Category of euclidean domains, Category of unique factorization domains, Category of division rings]

        """
        return [EuclideanDomains(), UniqueFactorizationDomains(), DivisionRings()]

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in Fields()
            True
            sage: QQ in Fields()
            True
            sage: ZZ in Fields()
            False
            sage: IntegerModRing(4) in Fields()
            False
            sage: InfinityRing in Fields()
            False

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`Fields`().

        Caveat: this should eventually be fixed::

            sage: gap.Rationals in Fields()
            False

        typically by implementing the method :meth:`category`
        appropriately for Gap objects::

            sage: GR = gap.Rationals
            sage: GR.category = lambda : Fields()
            sage: GR in Fields()
            True

        The following tests against a memory leak fixed in :trac:`13370`::

            sage: import gc
            sage: _ = gc.collect()
            sage: n = len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)])
            sage: for i in prime_range(100):
            ...     R = ZZ.quotient(i)
            ...     t = R in Fields()
            sage: _ = gc.collect()
            sage: len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)]) - n
            1

        """
        try:
            return self._contains_helper(x) or sage.rings.ring._is_Field(x)
        except StandardError:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of fields.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        fields. There are, however, rings that are initialised
        as plain commutative rings and found out to be fields
        only afterwards. Hence, this helper alone is not enough
        for a proper containment test.

        TESTS::

            sage: P.<x> = QQ[]
            sage: Q = P.quotient(x^2+2)
            sage: Q.category()
            Join of Category of commutative algebras over Rational Field and Category of subquotients of monoids and Category of quotients of semigroups
            sage: F = Fields()
            sage: F._contains_helper(Q)
            False
            sage: Q in F  # This changes the category!
            True
            sage: F._contains_helper(Q)
            True

        """
        return Category_contains_method_by_parent_class(cls())

    def _call_(self, x):
        """
        Construct a field from the data in ``x``

        EXAMPLES::

            sage: K = Fields()
            sage: K
            Category of fields
            sage: Fields().super_categories()
            [Category of euclidean domains, Category of unique factorization domains, Category of division rings]

            sage: K(IntegerRing()) # indirect doctest
            Rational Field
            sage: K(PolynomialRing(GF(3), 'x')) # indirect doctest
            Fraction Field of Univariate Polynomial Ring in x over
            Finite Field of size 3
            sage: K(RealField())
            Real Field with 53 bits of precision
        """
        try:
            return x.fraction_field()
        except AttributeError:
            raise TypeError, "unable to associate a field to %s"%x

    class ParentMethods:
        def is_field(self):
            """
            Return True, since this in an object of the category of fields.

            EXAMPLES::

                sage: Parent(QQ,category=Fields()).is_field()
                True

            """
            return True

        def is_integrally_closed(self):
            r"""

            Return ``True``, as per :meth:`IntegralDomain.is_integraly_closed`:
            for every field `F`, `F` is its own field of fractions,
            hence every element of `F` is integral over `F`.

            EXAMPLES::

                sage: QQ.is_integrally_closed()
                True
                sage: QQbar.is_integrally_closed()
                True
                sage: Z5 = GF(5); Z5
                Finite Field of size 5
                sage: Z5.is_integrally_closed()
                True
            """
            return True

        def _test_characteristic_fields(self, **options):
            """
            Run generic tests on the method :meth:`.characteristic`.

            EXAMPLES::

                sage: QQ._test_characteristic_fields()

            .. NOTE::

                We cannot call this method ``_test_characteristic`` since that
                would overwrite the method in the super category, and for
                cython classes just calling
                ``super(sage.categories.fields.Fields().parent_class,
                self)._test_characteristic`` doesn't have the desired effect.

            .. SEEALSO::

                :meth:`sage.categories.rings.Rings.ParentMethods._test_characteristic`
            """
            tester = self._tester(**options)
            try:
                char = self.characteristic()
                tester.assertTrue(char.is_zero() or char.is_prime())
            except AttributeError:
                return
                # raised when self.one() does not have a additive_order() [or when char is an int and not an Integer which is already checked by _test_characteristic for rings]
            except NotImplementedError:
                return

        def is_integral_domain(self):
            r"""

            Returns ``True``, as fields are integral domains.

            EXAMPLES::

                sage: QQ.is_integral_domain()
                True
            """
            return True

        def is_field( self, proof=True ):
            r"""
            Returns True as ``self`` is a field.

            EXAMPLES::

                sage: QQ.is_field()
                True
            """
            return True

        def fraction_field(self):
            r"""
            Returns the *fraction field* of ``self``, which is ``self``.

            EXAMPLES::

                sage: QQ.fraction_field() is QQ
                True
            """
            return self

        def __pow__(self, n):
            r"""
            Returns the vector space of dimension `n` over ``self``.

            EXAMPLES::

                sage: QQ^4
                Vector space of dimension 4 over Rational Field
            """
            from sage.modules.all import FreeModule
            return FreeModule(self, n)

        def _xgcd_univariate_polynomial(self, other):
            r"""
            Extended gcd of ``self`` and ``other``.

            INPUT:

                - ``other`` -- a polynomial in the same ring as ``self``

            OUTPUT:

            Polynomials ``g``, ``u``, and ``v`` such that ``g = u*self + v*other``

            .. NOTE::

                This is a helper method for
                :meth:`sage.rings.polynomial.polynomial_element.xgcd`

            EXAMPLES::

                sage: P.<x> = QQ[]
                sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
                sage: g, u, v = QQ._xgcd_univariate_polynomial(F,G)
                sage: g, u, v
                (x^2 + 2, 1/27, -1/27*x^2 - 1/9*x - 1/3)
                sage: u*F + v*G
                x^2 + 2

            ::

                sage: g, u, v = QQ._xgcd_univariate_polynomial(x,P(0)); g, u, v
                (x, 1, 0)
                sage: g == u*x + v*P(0)
                True
                sage: g, u, v = QQ._xgcd_univariate_polynomial(P(0),x); g, u, v
                (x, 0, 1)
                sage: g == u*P(0) + v*x
                True

            """
            R = self.parent()
            if other.is_zero():
                return self, R.one_element(), R.zero_element()
            # Algorithm 3.2.2 of Cohen, GTM 138
            A = self
            B = other
            U = R.one_element()
            G = A
            V1 = R.zero_element()
            V3 = B
            while not V3.is_zero():
                Q, R = G.quo_rem(V3)
                T = U - V1*Q
                U = V1
                G = V3
                V1 = T
                V3 = R
            V = (G-A*U)//B
            lc = G.leading_coefficient()
            return G/lc, U/lc, V/lc


    class ElementMethods:

        def is_unit( self ):
            r"""
            Returns True if ``self`` has a multiplicative inverse.

            EXAMPLES::

                sage: QQ(2).is_unit()
                True
                sage: QQ(0).is_unit()
                False
            """
            return not self.is_zero()

        # Fields are unique factorization domains, so, there is gcd and lcm
        # Of course, in general gcd and lcm in a field are not very interesting.
        # However, they should be implemented!
        def gcd(self,other):
            """
            Greatest common divisor.

            NOTE:

            Since we are in a field and the greatest common divisor is
            only determined up to a unit, it is correct to either return
            zero or one. Note that fraction fields of unique factorization
            domains provide a more sophisticated gcd.

            EXAMPLES::

                sage: GF(5)(1).gcd(GF(5)(1))
                1
                sage: GF(5)(1).gcd(GF(5)(0))
                1
                sage: GF(5)(0).gcd(GF(5)(0))
                0

            For fields of characteristic zero (i.e., containing the
            integers as a sub-ring), evaluation in the integer ring is
            attempted. This is for backwards compatibility::

                sage: gcd(6.0,8); gcd(6.0,8).parent()
                2
                Integer Ring

            If this fails, we resort to the default we see above::

                sage: gcd(6.0*CC.0,8*CC.0); gcd(6.0*CC.0,8*CC.0).parent()
                1.00000000000000
                Complex Field with 53 bits of precision

            AUTHOR:

            - Simon King (2011-02): Trac ticket #10771

            """
            P = self.parent()
            try:
                other = P(other)
            except (TypeError, ValueError):
                raise ArithmeticError, "The second argument can not be interpreted in the parent of the first argument. Can't compute the gcd"
            from sage.rings.integer_ring import ZZ
            if ZZ.is_subring(P):
                try:
                    return ZZ(self).gcd(ZZ(other))
                except TypeError:
                    pass
            # there is no custom gcd, so, we resort to something that always exists
            # (that's new behaviour)
            if self==0 and other==0:
                return P.zero()
            return P.one()

        def lcm(self,other):
            """
            Least common multiple.

            NOTE:

            Since we are in a field and the least common multiple is
            only determined up to a unit, it is correct to either return
            zero or one. Note that fraction fields of unique factorization
            domains provide a more sophisticated lcm.

            EXAMPLES::

                sage: GF(2)(1).lcm(GF(2)(0))
                0
                sage: GF(2)(1).lcm(GF(2)(1))
                1

            If the field contains the integer ring, it is first
            attempted to compute the gcd there::

                sage: lcm(15.0,12.0); lcm(15.0,12.0).parent()
                60
                Integer Ring

            If this fails, we resort to the default we see above::

                sage: lcm(6.0*CC.0,8*CC.0); lcm(6.0*CC.0,8*CC.0).parent()
                1.00000000000000
                Complex Field with 53 bits of precision
                sage: lcm(15.2,12.0)
                1.00000000000000

            AUTHOR:

            - Simon King (2011-02): Trac ticket #10771

            """
            P = self.parent()
            try:
                other = P(other)
            except (TypeError, ValueError):
                raise ArithmeticError, "The second argument can not be interpreted in the parent of the first argument. Can't compute the lcm"
            from sage.rings.integer_ring import ZZ
            if ZZ.is_subring(P):
                try:
                    return ZZ(self).lcm(ZZ(other))
                except TypeError:
                    pass
            # there is no custom lcm, so, we resort to something that always exists
            if self==0 or other==0:
                return P.zero()
            return P.one()

        @coerce_binop
        def xgcd(self, other):
            """
            Compute the extended gcd of ``self`` and ``other``.

            INPUT:

                - ``other`` -- an element with the same parent as ``self``

            OUTPUT:

                A tuple ``r,s,t`` of elements in the parent of ``self`` such
                that ``r = s*self + t*other``. Since the computations are done
                over a field, ``r`` is zero if ``self`` and ``other`` are zero,
                and one otherwise.

            AUTHORS:

            - Julian Rueth (2012-10-19): moved here from
              :class:`sage.structure.element.FieldElement`

            EXAMPLES::

                sage: (1/2).xgcd(2)
                (1, 2, 0)
                sage: (0/2).xgcd(2)
                (1, 0, 1/2)
                sage: (0/2).xgcd(0)
                (0, 0, 0)

            """
            R = self.parent()
            if not self.is_zero():
                return R.one(), ~self, R.zero()
            elif not other.is_zero():
                return R.one(), R.zero(), ~other
            else: # both are 0
                return R.zero(), R.zero(), R.zero()
