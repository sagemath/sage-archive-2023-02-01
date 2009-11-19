r"""
Schemes
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category, HomCategory
from sage.categories.category_types import Category_over_base
from sage.misc.cachefunc import cached_method

def Schemes(X=None):
    """
    Construct a category of schemes.

    EXAMPLES::

        sage: Schemes()
        Category of Schemes

        sage: Schemes(Spec(ZZ))
        Category of schemes over Spectrum of Integer Ring

        sage: Schemes(ZZ)
        Category of schemes over Spectrum of Integer Ring
    """
    if X is None:
        return Schemes_abstract()
    from sage.schemes.all import is_Scheme
    if not is_Scheme(X):
        X = Schemes()(X)
    return Schemes_over_base(X)

# TODO: rename into AbstractSchemes ???
class Schemes_abstract(Category):
    """
    The category of all abstract schemes.

    EXAMPLES::

        sage: Schemes()
        Category of Schemes
    """
    def __init__(self):
        """
        TESTS::

            sage: C = Schemes()
            sage: C
            Category of Schemes
            sage: TestSuite(C).run()
        """
        Category.__init__(self, "Schemes")

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Schemes().super_categories()
            [Category of sets]
        """
        from sets_cat import Sets
        return [Sets()]

    def _call_(self, x):
        """
        Construct a scheme from the data in ``x``

        EXAMPLES:

        Let us first construct the category of schemes::

            sage: S = Schemes(); S
            Category of Schemes

        We create a scheme from a ring::

            sage: X = S(ZZ); X                  # indirect doctest
            Spectrum of Integer Ring

        We create a scheme from a scheme (do nothing)::

            sage: S(X)
            Spectrum of Integer Ring

        We create a scheme morphism from a ring homomorphism.x::
\
            sage: phi = ZZ.hom(QQ); phi
            Ring Coercion morphism:
              From: Integer Ring
              To:   Rational Field
            sage: f = S(phi); f                 # indirect doctest
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

            sage: f.domain()
            Spectrum of Rational Field
            sage: f.codomain()
            Spectrum of Integer Ring
            sage: S(f)                          # indirect doctest
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        """
        from sage.rings.all import is_CommutativeRing, is_RingHomomorphism
        from sage.schemes.all import is_Scheme, Spec, is_SchemeMorphism
        if is_Scheme(x) or is_SchemeMorphism(x):
            return x
        elif is_CommutativeRing(x):
            return Spec(x)
        elif is_RingHomomorphism(x):
            A = Spec(x.codomain())
            return A.hom(x)
        else:
            raise TypeError, "No way to create an object or morphism in %s from %s"%(self, x)


    class HomCategory(HomCategory):
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Schemes().hom_category().extra_super_categories()
                []
                sage: Schemes().hom_category().super_categories()
                [Category of hom sets in Category of sets]

            FIXME: what category structure is there on Homsets of schemes?
            The result above is wrong, and should be fixed during the next
            homsets overhaul.
            """
            return []

        class ParentMethods:

            def __new__(cls, R, S, category):
                """
                TESTS::

                    sage: E = EllipticCurve('37a1')
                    sage: Hom(E, E).__class__
                    <class 'sage.schemes.generic.homset.SchemeHomset_generic_with_category'>

                If both schemes R and S are actually specs, we want
                the parent for Hom(R, S) to be in a different class::

                    sage: Hom(Spec(ZZ), Spec(ZZ)).__class__
                    <class 'sage.schemes.generic.homset.SchemeHomset_spec_with_category'>

                Currently, and to minimize the changes, this is done
                by delegating the job to SchemeHomset. This is not
                very robust: for example, only one category can do
                this hack.

                FIXME: this might be better handled by an extra Spec category
                """
                from sage.schemes.generic.homset import SchemeHomset
                return SchemeHomset(R, S, category=category)


#############################################################
# Schemes over a given base scheme.
#############################################################
class Schemes_over_base(Category_over_base):
    """
    The category of schemes over a given base scheme.

    EXAMPLES::

        sage: Schemes(Spec(ZZ))
        Category of schemes over Spectrum of Integer Ring

    TESTS::

        sage: C = Schemes(ZZ)
        sage: TestSuite(C).run()
    """

    def base_scheme(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)).base_scheme()
            Spectrum of Integer Ring
        """
        return self.base()

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)).super_categories()
            [Category of Schemes]
        """
        return [Schemes_abstract()]

    def _repr_(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)) # indirect doctest
            Category of schemes over Spectrum of Integer Ring
        """
        return "Category of schemes over %s"%self.base_scheme()
