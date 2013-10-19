r"""
Covariant Functorial Constructions

A *functorial construction* is a collection of functors
`(F_{Cat})_{Cat}` which associate to parents `A, B, ...` in some
category `Cat` a parent `F_{Cat}(A, B, ...)`. Examples of functorial
constructions are `F` = cartesian product, tensor product,
`\QQ`-algebras, ... The category of `F_{Cat}(A, B, ...)`, which only
depends on `Cat`, is called the (functorial) construction category.

A functorial construction is *(category)-covariant* if for every
categories `Cat` and `SuperCat`, the category of `F_{Cat}(A, B, ...)`
is a subcategory of the category of `F_{SuperCat}(A, B, ...)` whenever
`Cat` is a subcategory of `SuperCat`. A functorial construction is
*(category)-regressive* if the category of `F_{Cat}(A, B, ...)` is a
subcategory of `Cat`.

The goal of this module is to provide generic support for covariant
functorial constructions. In particular, given some parent `A`, `B`,
`C` of respective categories `Cat_A`, `Cat_B`, `Cat_C`, it provide
tools for calculating the best known category for the parent
`F(A,B,C)`. For examples, knowing that cartesian products of
semigroups (resp. monoids, groups), have a structure of semigroups
(resp. monoids, groups), and given a group `B` and two monoids `A` and
`C` it can calculate that `A \times B \times C` is naturally endowed
with a monoid structure.

See :class:`CovariantFunctorialConstruction`,
:class:`CovariantConstructionCategory` and
:class:`RegressiveCovariantConstructionCategory` for more details.

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.cachefunc import cached_function, cached_method
from sage.categories.category import Category
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

class CovariantFunctorialConstruction(UniqueRepresentation, SageObject):
    r"""
    An abstract class for construction functors `F` (eg `F` = cartesian
    product, tensor product, `\QQ`-algebras, ...) such that:

     - Each category `Cat` (eg `Cat=` ``Groups()``) can provide a category
       `F_{Cat}` for parents constructed via this functor (e.g.
       `F_{Cat} =` ``CartesianProductsOf(Groups())``).

     - For every category `Cat`, `F_{Cat}` is a subcategory of
       `F_{SuperCat}` for every super category `SuperCat` of
       `Cat` (the functorial construction is (category)-covariant).

     - For parents `A`, `B`, `C` respectively in the categories `Cat_A`,
       `Cat_B`, `Cat_C`, the category of `F(A,B,C)` is `F_{Cat}`
       where `Cat` is the meet of the categories `Cat_A`, `Cat_B`, `Cat_C`.

    This covers two slightly different use cases:

     - In the first use case, one uses directly the construction
       functor to create new parents::

         sage: tensor()  # todo: not implemented (add an example)

       or even new elements, which indirectly constructs the
       corresponding parent::

         sage: tensor(...) # todo: not implemented

     - In the second use case, one implements a parent, and then put
       it in the category `F_{Cat}` to specify supplementary
       mathematical information about that parent.

    The main purpose of this class is to handle automatically the
    trivial part of the category hiearchy. For example,
    ``CartesianProductsOf(Groups())`` is set automatically as a
    subcategory of ``CartesianProductsOf(Monoids())``.

    In practice, each subclass of this class should provide the
    following attributes:

     - ``_functor_category`` - a string which should match the name of
       the nested category class to be used in each category to
       specify information and generic operations for elements of this
       category.

     - ``_functor_name`` - an string which specifies the name of the
       functor, and also (when relevant) of the method on parents and
       elements used for calling the construction.


    TODO: What syntax do we want for `F_{Cat}`? For example, for the
    tensor product construction, which one of the followings do we want (see
    chat on IRC, on 07/12/2009):

     - ``tensor(Cat)``
     - ``tensor((Cat, Cat))``
     - ``tensor.of((Cat, Cat))``
     - ``tensor.category_from_categories((Cat, Cat, Cat))``
     - ``Cat.TensorProducts()``

    The syntax ``Cat.TensorProducts()`` does not supports well multivariate
    constructions like ``tensor.of([Algebras(), HopfAlgebras(), ...])``.
    Also it forces every category to be (somehow) aware of all the
    tensorial construction that could apply to it, even those which
    are only induced from super categories.

    Note: for each functorial construction, there probably is one (or several)
    largest categories on which it applies. For example, the
    :func:`~sage.categories.cartesian_product.CartesianProducts` construction makes
    only sense for concrete categories, that is subcategories of
    ``Sets()``. Maybe we want to model this one way or the other.
    """

    def category_from_parents(self, parents):
        """
        Returns the category of `F(A,B,C)` for `A,B,C` parents.

        INPUT:

         - self: a functor F
         - parents: a list (or iterable) of parents.

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"])
            sage: tensor.category_from_parents((E, E, E)) # todo: not implemented (see upcoming category patch #5985)
        """
        from sage.structure.parent import Parent
        assert(all(isinstance(parent, Parent) for parent in parents))
        # Should we pass a set of categories to reduce the cache size?
        # But then this would impose that, for any constructor, the
        # category of the result does not depend on the order/repetition
        # of the categories of the parents
        return self.category_from_categories(tuple(set(parent.category() for parent in parents)))


    @cached_method
    def category_from_categories(self, categories):
        """
        Returns the category of `F(A,B,C)` for `A,B,C` parents in the given categories

        INPUT:

         - ``self``: a functor `F`
         - ``categories``: a non empty tuple of categories

        EXAMPLES::

            sage: Cat1 = Rings()
            sage: Cat2 = Groups()
            sage: cartesian_product.category_from_categories((Cat1, Cat1, Cat1))
            Category of Cartesian products of monoids

            sage: cartesian_product.category_from_categories((Cat1, Cat2))
            Category of Cartesian products of monoids
        """
        assert(len(categories) > 0)
        return self.category_from_category(Category.meet(categories))

    def category_from_category(self, category):
        """
        Returns the category of `F(A,B,C)` for `A,B,C` parents in ``category``

        INPUT:

         - ``self``: a functor `F`
         - ``category``: a category

        EXAMPLES::

            sage: tensor.category_from_category(ModulesWithBasis(QQ))
            Category of tensor products of modules with basis over Rational Field

        # TODO: add support for parametrized functors
        """
        return getattr(category, self._functor_category)()

    def _repr_(self):
        """
        EXAMPLES::

            sage: tensor                       # indirect doctest
            The tensor functorial construction
        """
        return "The %s functorial construction"%self._functor_name

    def __call__(self, args):
        """
        Functorial construction application

        INPUT:
         - ``self``: a covariant functorial construction `F`
         - ``args``: a tuple (or iterable) of parents or elements

        Returns `F(args)`

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"]); E.rename("E")
            sage: tensor((E, E, E))
            E # E # E
        """
        args = tuple(args) # a bit brute force; let's see if this becomes a bottleneck later
        assert(all( hasattr(arg, self._functor_name) for arg in args))
        assert(len(args) > 0)
        return getattr(args[0], self._functor_name)(*args[1:])

class CovariantConstructionCategory(Category): # Should this be CategoryWithBase?
    """
    Abstract class for categories `F_{Cat}` obtained through a
    covariant functorial construction
    """

    @staticmethod
    def __classget__(cls, category, owner):
        r"""
        Special binding for covariant constructions

        This implements a hack allowing e.g. ``category.Subquotients``
        to recover the default ``Subquotients`` method defined in
        ``Category``, even if it has been overriden by a
        ``Subquotients`` class.

        TESTS::

            sage: Sets.Subquotients
            <class 'sage.categories.sets_cat.Sets.Subquotients'>
            sage: Sets().Subquotients
            Cached version of <function Subquotients at ...>
        """
        if category is None:
            return cls
        return getattr(super(category.__class__.__base__, category), cls._functor_category)

    @classmethod
    @cached_function
    def category_of(cls, category, *args):
        """
        This is the main entry point for constructing the category `F_{Cat}`
        of parents `F(A,B,C)` constructed from parents `A,B,C` in `Cat`.

        INPUT:

         - ``cls`` -- the category class for the functorial construction `F`
         - ``category`` -- a category `Cat`
         - ``*args`` -- further arguments for the functor

        OUTPUT: the image category of the functor `F_{Cat}`

        EXAMPLES::

            sage: sage.categories.tensor.TensorProductsCategory.category_of(ModulesWithBasis(QQ))
            Category of tensor products of modules with basis over Rational Field

            sage: sage.categories.algebra_functor.AlgebrasCategory.category_of(FiniteMonoids(), QQ)
            Join of Category of finite dimensional algebras with basis over Rational Field and Category of monoid algebras over Rational Field and Category of finite set algebras over Rational Field
        """
        # TODO: find a better test
        # the purpose is to test whether ``category`` implements that functor
        functor_category = getattr(category.__class__, cls._functor_category)
        if isinstance(functor_category, type) and issubclass(functor_category, Category):
            return functor_category(category, *args)
        else:
            return cls.default_super_categories(category, *args)

    @classmethod
    def default_super_categories(cls, category, *args):
        """
        Returns the default super categories of `F_{Cat}(A,B,C)` for `A,B,C`
        parents in `Cat`.

        INPUT:

         - ``cls`` -- the category class for the functor `F`
         - ``category`` -- a category `Cat`
         - ``*args`` -- further arguments for the functor

        OUTPUT: a (join) category

        The default implementation is to return the join of the
        categories of `F(A,B,C)` for `A,B,C` in turn in each of the
        super categories of ``category``.

        This is implemented as a class method, in order to be able to
        reconstruct the functorial category associated to each of the
        super categories of ``category``.

        EXAMPLES:

        Bialgebras are both algebras and coalgebras::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field, Category of coalgebras over Rational Field]

        Hence tensor products of bialgebras are tensor products of
        algebras and tensor products of coalgebras::

            sage: Bialgebras(QQ).TensorProducts().super_categories()
            [Category of tensor products of algebras over Rational Field, Category of tensor products of coalgebras over Rational Field]

        Here is how :meth:`default_super_categories` was called internally::

            sage: sage.categories.tensor.TensorProductsCategory.default_super_categories(Bialgebras(QQ))
            Join of Category of tensor products of algebras over Rational Field and Category of tensor products of coalgebras over Rational Field

        We now show a similar example, with the ``Algebra`` functor
        which takes a parameter `\QQ`::

            sage: FiniteMonoids().super_categories()
            [Category of monoids, Category of finite semigroups]
            sage: sorted(FiniteMonoids().Algebras(QQ).super_categories(), key=str)
            [Category of finite dimensional algebras with basis over Rational Field,
             Category of finite set algebras over Rational Field,
             Category of monoid algebras over Rational Field]

        Note that neither the category of *finite* semigroup algebras
        nor that of monoid algebras appear in the result; this is
        because there is currently nothing specific implemented about them.

        Here is how :meth:`default_super_categories` was called internally::

            sage: sage.categories.algebra_functor.AlgebrasCategory.default_super_categories(FiniteMonoids(), QQ)
            Join of Category of finite dimensional algebras with basis over Rational Field and Category of monoid algebras over Rational Field and Category of finite set algebras over Rational Field
        """
        return Category.join([getattr(cat, cls._functor_category)(*args) for cat in category._super_categories if hasattr(cat, cls._functor_category)])

    def __init__(self, category, *args):
        """
        TESTS::

            sage: from sage.categories.covariant_functorial_construction import CovariantConstructionCategory
            sage: class FooBars(CovariantConstructionCategory):
            ...       _functor_category = "FooBars"
            sage: Category.FooBars = lambda self: FooBars.category_of(self)
            sage: C = FooBars(ModulesWithBasis(ZZ))
            sage: C
            Category of foo bars of modules with basis over Integer Ring
            sage: C.base_category()
            Category of modules with basis over Integer Ring
            sage: latex(C)
            \mathbf{FooBars}(\mathbf{ModulesWithBasis}_{\Bold{Z}})
            sage: import __main__; __main__.FooBars = FooBars # Fake FooBars being defined in a python module
            sage: TestSuite(C).run()
        """
        assert isinstance(category, Category)
        self._base_category = category
        self._args = args
        super(CovariantConstructionCategory, self).__init__(*args)

    def base_category(self):
        """
        Returns the base category of the category ``self``

        For any category ``B`` = `F_{Cat}` obtained through a functorial
        construction `F`, the call ``B.base_category()`` returns the
        category `Cat`

        EXAMPLES::

            sage: Semigroups().Quotients().base_category()
            Category of semigroups
        """
        return self._base_category

    def extra_super_categories(self):
        """
        Returns the extra super categories of a construction category

        Default implementation which returns ``[]``

        EXAMPLES::

            sage: Sets().Subquotients().extra_super_categories()
            []
            sage: Semigroups().Quotients().extra_super_categories()
            []
        """
        return []

    def super_categories(self):
        """
        Returns the super categories of a construction category

        EXAMPLES::

            sage: Sets().Subquotients().super_categories()
            [Category of sets]
            sage: Semigroups().Quotients().super_categories()
            [Category of subquotients of semigroups, Category of quotients of sets]
        """
        return Category.join([self.__class__.default_super_categories(self.base_category(), *self._args)] +
                             self.extra_super_categories(),
                             as_list = True)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Semigroups().Subquotients()  # indirect doctest
            Category of subquotients of semigroups
        """
        return "%s of %s"%(Category._repr_object_names(self), self.base_category()._repr_object_names())

    def _latex_(self):
        """
        EXAMPLES::

            sage: latex(Semigroups().Subquotients())   # indirect doctest
            \mathbf{Subquotients}(\mathbf{Semigroups})
            sage: latex(ModulesWithBasis(QQ).TensorProducts())
            \mathbf{TensorProducts}(\mathbf{ModulesWithBasis}_{\Bold{Q}})
            sage: latex(Semigroups().Algebras(QQ))
            \mathbf{Algebras}(\mathbf{Semigroups})
        """
        from sage.misc.latex import latex
        return "\\mathbf{%s}(%s)"%(self._short_name(), latex(self.base_category()))

class RegressiveCovariantConstructionCategory(CovariantConstructionCategory):
    """
    Abstract class for categories `F_{Cat}` obtained through a
    regressive covariant functorial construction
    """

    @classmethod
    def default_super_categories(cls, category, *args):
        """
        Returns the default super categories of ``...``

        INPUT:

         - ``category`` -- a category

        OUTPUT: a join category

        This implements the property that an induced subcategory is a
        subcategory.

        EXAMPLES:

        A subquotient of a monoid is a monoid, and a subquotient of
        semigroup::

            sage: Monoids().Subquotients().super_categories()
            [Category of monoids, Category of subquotients of semigroups]

        TESTS::

            sage: C = Monoids().Subquotients()
            sage: C.__class__.default_super_categories(C.base_category(), *C._args)
            Category of unital subquotients of semigroups
        """
        return Category.join([category, super(RegressiveCovariantConstructionCategory, cls).default_super_categories(category, *args)])
