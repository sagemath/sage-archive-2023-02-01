r"""
Fields
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
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
from sage.rings.field import is_Field

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
        """
        try:
            return self._contains_helper(x) or is_Field(x)
        except:
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

    class ElementMethods:
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
