r"""
Rings
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
from category import HomCategory
from sage.misc.cachefunc import cached_method

class Rings(Category):
    """
    The category of rings

    Associative rings with unit, not necessarily commutative

    EXAMPLES::

      sage: Rings()
      Category of rings
      sage: Rings().super_categories()
      [Category of rngs, Category of monoids]

    TESTS::

        sage: TestSuite(Rings()).run()

    TODO (see: http://trac.sagemath.org/sage_trac/wiki/CategoriesRoadMap)

     - Make Rings() into a subcategory or alias of Algebras(ZZ);

     - A parent P in the category ``Rings()`` should automatically be
       in the category ``Algebras(P)``.
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Rings().super_categories()
            [Category of rngs, Category of monoids]
        """
        from sage.categories.rngs import Rngs
        from sage.categories.monoids import Monoids
        return [Rngs(), Monoids()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass


    class HomCategory(HomCategory):
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Rings().hom_category().extra_super_categories()
                [Category of sets]
            """
            from sage.categories.sets_cat import Sets
            return [Sets()]

#         def get_Parent(self, X, Y):
#             """
#             Given two objects X and Y in this category, returns the parent
#             class to be used for the collection of the morphisms of this
#             category between X and Y.

#             Returns self.ParentMethods by default.

#             Rationale: some categories, like Rings or Schemes, currently
#             use different classes for their homset, depending on some
#             specific properties of X or Y which do not fit in the category
#             hierarchy. For example, if X is a quotient field, morphisms
#             can be defined by the image of the generators, even if Y
#             itself is not a quotient field.

#             Design question: should this really concern the parent for the
#             homset, or just the possible classes for the elements?
#             """
#             category = self.base_category
#             assert(X in category and Y in category)
#             # return self.hom_category()(X, Y)?
#             #print self.hom_category(), X, Y, self.hom_category().parent_class, self.hom_category().parent_class.mro()
#             return self.ParentMethods

        class ParentMethods:
            # Design issue: when X is a quotient field, we can build
            # morphisms from X to Y by specifying the images of the
            # generators. This is not something about the category,
            # because Y need not be a quotient field.

            # Currently, and to minimize the changes, this is done by
            # delegating the job to RingHomset. This is not very robust:
            # for example, only one category can do this hack.

            # This should be cleaned up upon the next homset overhaul

            def __new__(cls, X, Y, category):
                """
                    sage: Hom(QQ, QQ, category = Rings()).__class__                  # indirect doctest
                    <class 'sage.rings.homset.RingHomset_generic_with_category'>

                    sage: Hom(CyclotomicField(3), QQ, category = Rings()).__class__  # indirect doctest
                    <class 'sage.rings.number_field.morphism.CyclotomicFieldHomset_with_category'>
                """
                from sage.rings.homset import RingHomset
                return RingHomset(X, Y, category = category)

            def __getnewargs__(self):
                """
                Note: without this method, :meth:`.__new__` gets called with no
                argument upon unpickling. Maybe it would be preferable to
                have :meth:`.__new__` accept to be called without arguments.

                TESTS::

                    sage: Hom(QQ, QQ, category = Rings()).__getnewargs__()
                    (Rational Field, Rational Field, Category of hom sets in Category of rings)
                    sage: TestSuite(Hom(QQ, QQ, category = Rings())).run() # indirect doctest
                """
                return (self.domain(), self.codomain(), self.category())
