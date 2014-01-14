"""
Isomorphic Objects Functorial Construction

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

# This is Category.IsomorphicObjects
def IsomorphicObjects(self):
    """
    INPUT:
     - ``self`` -- a concrete category

    Given a concrete category ``As()`` (i.e. a subcategory of
    ``Sets()``), ``As().IsomorphicObjects()`` returns the category of
    objects of ``As()`` endowed with a distinguished description as
    the image of some other object of ``As()`` by an isomorphism.

    See :func:`~sage.categories.subquotients.Subquotients` for background.

    EXAMPLES::

        sage: C = Sets().IsomorphicObjects(); C
        Category of isomorphic objects of sets

        sage: C.super_categories()
        [Category of subobjects of sets, Category of quotients of sets]

        sage: C.all_super_categories()
        [Category of isomorphic objects of sets,
         Category of subobjects of sets,
         Category of quotients of sets,
         Category of subquotients of sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    Unless something specific about isomorphic objects is implemented
    for this category, one actually get an optimized super category::

        sage: C = Semigroups().IsomorphicObjects(); C
        Join of Category of quotients of semigroups and Category of isomorphic objects of sets

    TESTS::

        sage: TestSuite(Sets().IsomorphicObjects()).run()
    """
    return IsomorphicObjectsCategory.category_of(self)

Category.IsomorphicObjects = IsomorphicObjects

class IsomorphicObjectsCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "IsomorphicObjects"

    @classmethod
    def default_super_categories(cls, category):
        """
        Returns the default super categories of ``category.IsomorphicObjects()``

        Mathematical meaning: if `A` is the image of `B` by an
        isomorphism in the category `C`, then `A` is both a subobject
        of `B` and a quotient of `B` in the category `C`.

        INPUT:

         - ``cls`` -- the class ``IsomorphicObjectsCategory``
         - ``category`` -- a category `Cat`

        OUTPUT: a (join) category

        In practice, this returns ``category.Subobjects()`` and
        ``category.Quotients()``, joined together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and
        ``cat.IsomorphicObjects()`` for each ``cat`` in the super
        categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``, which has ``cat=Monoids()`` as
        super category. Then, the image of a group `G'` by a group
        isomorphism is simultaneously a subgroup of `G`, a subquotient
        of `G`, a group by itself, and the image of `G` by a monoid
        isomorphism::

            sage: Groups().IsomorphicObjects().super_categories()
            [Category of groups,
             Category of subquotients of monoids,
             Category of quotients of semigroups,
             Category of isomorphic objects of sets]

        Mind the last item above: there is indeed currently nothing
        implemented about isomorphic objects of monoids.

        This resulted from the following call::

            sage: sage.categories.isomorphic_objects.IsomorphicObjectsCategory.default_super_categories(Groups())
            Join of Category of groups and
            Category of subquotients of monoids and
            Category of quotients of semigroups and
            Category of isomorphic objects of sets
        """
        return Category.join([category.Subobjects(), category.Quotients(),
                              super(IsomorphicObjectsCategory, cls).default_super_categories(category)])

