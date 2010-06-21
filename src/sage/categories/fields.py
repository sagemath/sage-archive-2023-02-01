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
from sage.misc.cachefunc import cached_method

class Fields(Category):
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

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Fields().super_categories()
            [Category of euclidean domains, Category of unique factorization domains, Category of division rings]
        """
        from sage.categories.basic import EuclideanDomains, UniqueFactorizationDomains, DivisionRings
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
        from sage.rings.field import is_Field
        try:
            return super(Fields, self).__contains__(x) or is_Field(x)
        except:
            return False

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

    class ElementMethods:
        pass
