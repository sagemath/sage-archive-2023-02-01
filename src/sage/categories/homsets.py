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
        """
        Return the homsets category of ``category``.

        This classmethod centralizes the construction of all homset
        categories.

        The ``cls`` and ``args`` arguments below are essentially
        unused. Their purpose is solely to let the code deviate as
        little as possible from the generic implementation of this
        method: :meth:`FunctorialConstructionCategory.category_of`.
        The reason for this deviation is that, unlike in the other
        functorial constructions which are covariant, we recurse only
        on *full* supercategories; then, we need a special treatment
        for the base case were a category neither defines the
        ``Homsets`` construction, nor inherits it from its full
        supercategories.

        INPUT:

         - ``cls`` -- :class:`HomsetsCategory` or a subclass thereof
         - ``category`` -- a category `Cat`
         - ``*args`` -- (unused)

        EXAMPLES:

        If ``category`` implements a ``Homsets`` class, then this
        class is used to build the homset category::

            sage: from sage.categories.homsets import HomsetsCategory
            sage: H = HomsetsCategory.category_of(Modules(ZZ)); H
            Category of homsets of modules over Integer Ring
            sage: type(H)
            <class 'sage.categories.modules.Modules.Homsets_with_category'>

        Otherwise, if ``category`` has one or more full super
        categories, then the join of their respective homsets category
        is returned. In this example, the join consists of a single
        category::

            sage: C = Modules(ZZ).WithBasis().FiniteDimensional()
            sage: C.full_super_categories()
            [Category of modules with basis over Integer Ring,
             Category of finite dimensional modules over Integer Ring]
            sage: H = HomsetsCategory.category_of(C); H
            Category of homsets of modules with basis over Integer Ring
            sage: type(H)
            <class 'sage.categories.modules_with_basis.ModulesWithBasis.Homsets_with_category'>

        As a last resort, a :class:`HomsetsOf` of the categories
        forming the structure of ``self`` is constructed::

            sage: H = HomsetsCategory.category_of(Magmas()); H
            Category of homsets of magmas
            sage: type(H)
            <class 'sage.categories.homsets.HomsetsOf_with_category'>

            sage: HomsetsCategory.category_of(Rings())
            Category of homsets of unital magmas and additive unital additive magmas
        """
        functor_category = getattr(category.__class__, cls._functor_category)
        if isinstance(functor_category, type) and issubclass(functor_category, Category):
            return functor_category(category, *args)
        elif category.full_super_categories():
            return cls.default_super_categories(category, *args)
        else:
            return HomsetsOf(Category.join(category.structure()))

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
        tester.assert_(self.is_subcategory(Category.join(self.base_category().structure()).Homsets()))

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

    This is used when a category `C` defines some additional structure
    but not a homset category of its own. Indeed, unlike for covariant
    functorial constructions, we cannot represent the homset category
    of `C` by just the join of the homset categories of its super
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

    class SubcategoryMethods:

        def Endset(self):
            """
            Return the subcategory of the homsets of ''self'' that are endomorphism sets.

            EXAMPLES::

                sage: Sets().Homsets().Endset()
                Category of endsets of sets

                sage: Posets().Homsets().Endset()
                Category of endsets of posets
            """
            return self._with_axiom("Endset")

    class Endset(CategoryWithAxiom):
        """
        The category of all endomorphism sets.

        This category serves too purposes: making sure that the
        ``Endset`` axiom is implemented in the category where it's
        defined, namely ``Homsets``, and specifying that ``Endsets``
        are monoids.

        EXAMPLES::

            sage: from sage.categories.homsets import Homsets
            sage: Homsets().Endset()
            Category of endsets
        """
        def extra_super_categories(self):
            """
            Implement the fact that endsets are monoids.

            .. SEEALSO:: :meth:`CategoryWithAxiom.extra_super_categories`

            EXAMPLES::

                sage: from sage.categories.homsets import Homsets
                sage: Homsets().Endset().extra_super_categories()
                [Category of monoids]
            """
            from monoids import Monoids
            return [Monoids()]
