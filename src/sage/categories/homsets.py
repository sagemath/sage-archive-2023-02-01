## -*- encoding: utf-8 -*-
r"""
Homset categories
"""
#*****************************************************************************
#  Copyright (C) 2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category import Category, JoinCategory
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.covariant_functorial_construction import FunctorialConstructionCategory

import sage.categories.category_with_axiom

class HomsetsCategory(FunctorialConstructionCategory):

    _functor_category = "Homsets"

    @classmethod
    @cached_function
    def category_of(cls, category, *args):
        functor_category = getattr(category.__class__, cls._functor_category)
        if isinstance(functor_category, type) and issubclass(functor_category, Category):
            return functor_category(category, *args)
        elif category.full_super_categories():
            return cls.default_super_categories(category, *args)
        else:
            return HomsetsOf(Category.join(category.super_structure_categories()))

    @classmethod
    def default_super_categories(cls, category):
        """
        Return the default super categories of ``category.Homsets()``

        INPUT:

         - ``cls`` -- the category class for the functor `F`
         - ``category`` -- a category `Cat`

        OUTPUT: a category

        .. TODO:: adapt everything below

        The default implementation is to return the join of the
        categories of `F(A,B,...)` for `A,B,...` in turn in each of
        the super categories of ``category``.

        This is implemented as a class method, in order to be able to
        reconstruct the functorial category associated to each of the
        super categories of ``category``.

        EXAMPLES::

            sage: AdditiveMagmas().Homsets().super_categories()
            [Category of additive magmas]
            sage: AdditiveMagmas().AdditiveUnital().Homsets().super_categories()
            [Category of additive unital additive magmas]

        For now nothing specific is implemented for homsets of
        additive groups compared to homsets of monoids::

            sage: from sage.categories.additive_groups import AdditiveGroups
            sage: AdditiveGroups().Homsets()
            Category of homsets of additive monoids

        Similarly for rings; so a ring homset is a homset of unital
        magmas and additive magmas::

            sage: Rings().Homsets()
            Category of homsets of unital magmas and additive unital additive magmas
        """
        return Category.join([getattr(cat, cls._functor_category)()
                              for cat in category.full_super_categories()])


    def extra_super_categories(self):
        """
        Return the extra super categories of ``self``.

        EXAMPLES::

            sage: Objects().Homsets().super_categories()
            [Category of homsets]
            sage: Sets().Homsets().super_categories()
            [Category of homsets]
            sage: (Magmas() & Posets()).Homsets().super_categories()
            [Category of homsets]
        """
        return [Homsets()]


    def _test_homsets_category(self, **options):
        r"""
        Run generic tests on this homsets category

        .. SEEALSO:: :class:`TestSuite`.

        EXAMPLES::

            sage: Sets().Homsets()._test_homsets_category()
        """
        # TODO: remove if unneeded
        #from sage.categories.objects    import Objects
        #from sage.categories.sets_cat import Sets
        tester = self._tester(**options)
        tester.assert_(self.is_subcategory(Category.join(self.base_category().super_structure_categories()).Homsets()))

    @cached_method
    def base(self):
        """
        If this homsets category is subcategory of a category with a base, return that base.

        .. TODO:: Is this really useful?

        EXAMPLES::

            sage: ModulesWithBasis(ZZ).Homsets().base()
            Integer Ring

        """
        from sage.categories.category_types import Category_over_base
        for C in self._all_super_categories_proper:
            if isinstance(C,Category_over_base):
                return C.base()
        raise AttributeError("This hom category has no base")

class HomsetsOf(HomsetsCategory):
    """
    Default class for homsets of a category.

    This is used for structure categories that do not define a homsets
    category of their own.

    This is needed because, unlike for covariant functorial
    constructions, we cannot represent the homset category of such a
    category by just the join of the homset categories of its super
    categories.

    EXAMPLES::

        sage: C = (Magmas() & Posets()).Homsets(); C
        Category of homsets of magmas and posets
        sage: type(C)
        <class 'sage.categories.homsets.HomsetsOf_with_category'>

    TESTS::

        sage: TestSuite(C).run()
        sage: C = Rings().Homsets()
        sage: TestSuite(C).run(skip=['_test_category_graph'])
        sage: TestSuite(C).run()
    """

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Semigroups().Homsets()
            Category of homsets of magmas
            sage: (Magmas() & AdditiveMagmas() & Posets()).Homsets()
            Category of homsets of magmas and additive magmas and posets
            sage: Rings().Homsets()
            Category of homsets of unital magmas and additive unital additive magmas
        """
        base_category = self.base_category()
        try:
            object_names = base_category._repr_object_names()
        except ValueError:
            assert isinstance(base_category, JoinCategory)
            object_names = ' and '.join(cat._repr_object_names() for cat in base_category.super_categories())
        return "homsets of %s"%(object_names)

class Homsets(Category_singleton):
    """
    The category of all homsets.

    EXAMPLES::

        sage: from sage.categories.homsets import Homsets
        sage: Homsets()
        Category of homsets

    This is a subcategory of ``Sets()``::

        sage: Homsets().super_categories()
        [Category of sets]

    By this, we assume that all homsets implemented in Sage are sets,
    or equivalently that we only implement locally small categories.
    See :wikipedia:`Category_(mathematics)`.

    .. TODO::

        We would want a more general mechanism. See also
        :meth:`Monoids.Homsets.extra_super_categories`.
    """
    def super_categories(self):
        """
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.categories.homsets import Homsets
            sage: Homsets()
            Category of homsets
        """
        from sets_cat import Sets
        return [Sets()]

    class Endset(CategoryWithAxiom):
        """
        The category of all end sets

        At this point, the main purpose of this class is to keep the
        information that a category is a Endset category, even if the
        base category implements nothing specific about those.
        """
        def super_categories(self):
            from monoids import Monoids
            return [Monoids()]
