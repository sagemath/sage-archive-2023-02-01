r"""
Categories

AUTHORS:

- David Kohel, William Stein and Nicolas M. Thiery

Every Sage object lies in a category. Categories in Sage are
modeled on the mathematical idea of category, and are distinct from
Python classes, which are a programming construct.

In most cases, typing ``x.category()`` returns the category to which ``x``
belongs. If ``C`` is a category and ``x`` is any object, ``C(x)`` tries to
make an object in ``C`` from ``x``. Checking if ``x`` belongs to ``C`` is done
as usually by ``x in C``.

See :class:`Category` and :mod:`sage.categories.primer` for more details.

EXAMPLES:

We create a couple of categories::

    sage: Sets()
    Category of sets
    sage: GSets(AbelianGroup([2,4,9]))
    Category of G-sets for Multiplicative Abelian group isomorphic to C2 x C4 x C9
    sage: Semigroups()
    Category of semigroups
    sage: VectorSpaces(FiniteField(11))
    Category of vector spaces over Finite Field of size 11
    sage: Ideals(IntegerRing())
    Category of ring ideals in Integer Ring

Let's request the category of some objects::

    sage: V = VectorSpace(RationalField(), 3)
    sage: V.category()
    Category of vector spaces over Rational Field
    sage: G = SymmetricGroup(9)
    sage: G.category()
    Join of Category of finite permutation groups and Category of finite weyl groups
    sage: P = PerfectMatchings(3)
    sage: P.category()
    Category of finite enumerated sets

Let's check some memberships::

    sage: V in VectorSpaces(QQ)
    True
    sage: V in VectorSpaces(FiniteField(11))
    False
    sage: G in Monoids()
    True
    sage: P in Rings()
    False

For parametrized categories one can use the following shorthand::

    sage: V in VectorSpaces
    True
    sage: G in VectorSpaces
    False

A parent ``P`` is in a category ``C`` if ``P.category()`` is a subcategory of
``C``.

.. note::

    Any object of a category should be an instance of
    :class:`~sage.structure.category_object.CategoryObject`.

    For backward compatibilty this is not yet enforced::

        sage: class A:
        ...     def category(self):
        ...         return Fields()
        sage: A() in Rings()
        True

    By default, the category of an element `x` of a parent `P` is the category
    of all objects of `P` (this is dubious an may be deprecated)::

        sage: V = VectorSpace(RationalField(), 3)
        sage: v = V.gen(1)
        sage: v.category()
        Category of elements of Vector space of dimension 3 over Rational Field
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import inspect
from sage.misc.abstract_method import abstract_method, abstract_methods_of_class
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.c3_controlled import C3_sorted_merge, category_sort_key, _cmp_key, _cmp_key_named
from sage.misc.unknown import Unknown

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.dynamic_class import DynamicMetaclass, dynamic_class

from weakref import WeakValueDictionary
_join_cache = WeakValueDictionary()

def _join(categories, as_list):
    """
    This is an auxiliary function for :meth:`Category.join`

    INPUT:

    - ``categories``: A tuple (no list) of categories.
    - ``as_list`` (boolean): Whether or not the result should be represented as a list.

    EXAMPLES::

        sage: Category.join((Groups(), CommutativeAdditiveMonoids()))  # indirect doctest
        Join of Category of groups and Category of commutative additive monoids
        sage: Category.join((Modules(ZZ), FiniteFields()), as_list=True)
        [Category of finite fields, Category of modules over Integer Ring]

    """
    # Since Objects() is the top category, it is the neutral element of join
    if len(categories) == 0:
        from objects import Objects
        return Objects()

    if not as_list:
        try:
            return _join_cache[categories]
        except KeyError:
            pass

    # Ensure associativity by flattening JoinCategory's
    # Invariant: the super categories of a JoinCategory are not JoinCategories themselves
    categories = sum( (tuple(category._super_categories) if isinstance(category, JoinCategory) else (category,)
                       for category in categories), ())

    # canonicalize, by removing redundant categories which are super
    # categories of others, and by sorting
    result = ()
    for category in categories:
        if any(cat.is_subcategory(category) for cat in result):
            continue
        result = tuple( cat for cat in result if not category.is_subcategory(cat) ) + (category,)
    result = tuple(sorted(result, key = category_sort_key, reverse=True))
    if as_list:
        return list(result)
    if len(result) == 1:
        out = _join_cache[categories] = result[0]
    else:
        out = _join_cache[categories] = JoinCategory(result)
    return out


class Category(UniqueRepresentation, SageObject):
    r"""
    The base class for modeling mathematical categories, like for example:

    - Groups(): the category of groups
    - EuclideanRings(): the category of euclidean rings
    - VectorSpaces(QQ): the category of vector spaces over the field of rational

    See :mod:`sage.categories.primer` for an introduction to
    categories in Sage, their relevance, purpose and usage. The
    documentation below focus on their implementation.

    Technically, a category is an instance of the class
    :class:`Category` or some of its subclasses. Some categories, like
    VectorSpaces, are parametrized: ``VectorSpaces(QQ)`` is one of
    many instances of the class :class:`VectorSpaces`. On the other
    hand, ``EuclideanRings()`` is the single instance of the class
    :class:`EuclideanRings`.

    Recall that an algebraic structure (say the ring QQ[x]) is
    modelled in Sage by an object which is called a parent. This
    object belongs to certain categories (here EuclideanRings() and
    Algebras()). The elements of the ring are themselves objects.

    The class of a category (say EuclideanRings) can define simultaneously:

    - Operations on the category itself (what are its super categories, its category of
      morphisms?, its dual category)
    - Generic operations on parents in this category, like the ring QQ[x]
    - Generic operations on elements of this ring (Euclide algorithm for computing gcds)

    This is achieved as follows::

        sage: from sage.categories.all import Category
        sage: class EuclideanRings(Category):
        ...       # operations on the category itself
        ...       def super_categories(self):
        ...           [Rings()]
        ...
        ...       def dummy(self): # TODO: find some good examples
        ...            pass
        ...
        ...       class ParentMethods: # holds the generic operations on parents
        ...            # find a good example of operation
        ...            pass
        ...
        ...       class ElementMethods:# holds the generic operations on elements
        ...            def gcd(x,y):
        ...                # Euclid algorithms
        ...                pass

    Note that the EuclideanRings.ParentMethods and .Element class above do
    not inherit from anything. They are merely containers of
    operations. The hierarchy between the different categories is
    defined once at the level of the categories. Behind the scene, a
    parallel hierarchy of classes is built automatically from all the
    .ParentMethods classes. Then, a parent in a category receives the
    appropriate operations from all the super categories by usual
    class inheritance. Similarly, a third hierarchy of classes is
    built for elements from the .Elements.

    EXAMPLES:

    We define a hierarchy of four categories As(), Bs(), Cs(), Ds()
    with a diamond inheritance. Think for example:

    - As(): the category of sets
    - Bs(): the category of additive groups
    - Cs(): the category of multiplicative monoids
    - Ds(): the category of rings

    ::

        sage: from sage.categories.all import Category
        sage: from sage.misc.lazy_attribute import lazy_attribute
        sage: class As (Category):
        ...       def super_categories(self):
        ...           return []
        ...
        ...       class ParentMethods:
        ...           def fA(self):
        ...               return "A"
        ...           f = fA
        ...
        sage: class Bs (Category):
        ...       def super_categories(self):
        ...           return [As()]
        ...
        ...       class ParentMethods:
        ...           def fB(self):
        ...               return "B"
        ...
        sage: class Cs (Category):
        ...       def super_categories(self):
        ...           return [As()]
        ...
        ...       class ParentMethods:
        ...           def fC(self):
        ...               return "C"
        ...           f = fC
        ...
        sage: class Ds (Category):
        ...       def super_categories(self):
        ...           return [Bs(),Cs()]
        ...
        ...       class ParentMethods:
        ...           def fD(self):
        ...               return "D"
        ...

    Categories should always have unique representation; by trac ticket
    #12215, this means that it will be kept in cache, but only if there
    is still some strong reference to it.

    We check this before proceeding::

        sage: import gc
        sage: idAs = id(As())
        sage: _ = gc.collect()
        sage: n == id(As())
        False
        sage: a = As()
        sage: id(As()) == id(As())
        True
        sage: As().parent_class == As().parent_class
        True

    We construct a parent in the category Ds() (that is an instance of
    Ds().parent_class), and check that it has access to all the
    methods provided by all the categories, with the appropriate
    inheritance order.
    ::

        sage: D = Ds().parent_class()
        sage: [ D.fA(), D.fB(), D.fC(), D.fD() ]
        ['A', 'B', 'C', 'D']
        sage: D.f()
        'C'

    ::

        sage: C = Cs().parent_class()
        sage: [ C.fA(), C.fC() ]
        ['A', 'C']
        sage: C.f()
        'C'

    Here is the parallel hierarchy of classes which has been built
    automatically, together with the method resolution order (``.mro()``)::

        sage: As().parent_class
        <class '__main__.As.parent_class'>
        sage: As().parent_class.__bases__
        (<type 'object'>,)
        sage: As().parent_class.mro()
        [<class '__main__.As.parent_class'>, <type 'object'>]

    ::

        sage: Bs().parent_class
        <class '__main__.Bs.parent_class'>
        sage: Bs().parent_class.__bases__
        (<class '__main__.As.parent_class'>,)
        sage: Bs().parent_class.mro()
        [<class '__main__.Bs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>]

    ::

        sage: Cs().parent_class
        <class '__main__.Cs.parent_class'>
        sage: Cs().parent_class.__bases__
        (<class '__main__.As.parent_class'>,)
        sage: Cs().parent_class.__mro__
        (<class '__main__.Cs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>)

    ::

        sage: Ds().parent_class
        <class '__main__.Ds.parent_class'>
        sage: Ds().parent_class.__bases__
        (<class '__main__.Cs.parent_class'>, <class '__main__.Bs.parent_class'>)
        sage: Ds().parent_class.mro()
        [<class '__main__.Ds.parent_class'>, <class '__main__.Cs.parent_class'>, <class '__main__.Bs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>]

    Note that that two categories in the same class need not have the
    same super_categories. For example, Algebras(QQ) has
    VectorSpaces(QQ) as super category, whereas Algebras(ZZ) only has
    Modules(ZZ) as super category. In particular, the constructed
    parent class and element class will differ (inheriting, or not,
    methods specific for vector spaces)::

        sage: Algebras(QQ).parent_class is Algebras(ZZ).parent_class
        False
        sage: issubclass(Algebras(QQ).parent_class, VectorSpaces(QQ).parent_class)
        True

    On the other hand, identical hierarchies of classes are,
    preferably, built only once (e.g. for categories over a base ring)::

        sage: Algebras(GF(5)).parent_class is Algebras(GF(7)).parent_class
        True
        sage: Coalgebras(QQ).parent_class is Coalgebras(FractionField(QQ[x])).parent_class
        True

    We now construct a parent in the usual way::

        sage: class myparent(Parent):
        ...       def __init__(self):
        ...           Parent.__init__(self, category=Ds())
        ...       def g(self):
        ...           return "myparent"
        ...       class Element:
        ...           pass
        sage: D = myparent()
        sage: D.__class__
        <class '__main__.myparent_with_category'>
        sage: D.__class__.__bases__
        (<class '__main__.myparent'>, <class '__main__.Ds.parent_class'>)
        sage: D.__class__.mro()
        [<class '__main__.myparent_with_category'>,
        <class '__main__.myparent'>,
        <type 'sage.structure.parent.Parent'>,
        <type 'sage.structure.category_object.CategoryObject'>,
        <type 'sage.structure.sage_object.SageObject'>,
        <class '__main__.Ds.parent_class'>,
        <class '__main__.Cs.parent_class'>,
        <class '__main__.Bs.parent_class'>,
        <class '__main__.As.parent_class'>,
        <type 'object'>]
        sage: D.fA()
        'A'
        sage: D.fB()
        'B'
        sage: D.fC()
        'C'
        sage: D.fD()
        'D'
        sage: D.f()
        'C'
        sage: D.g()
        'myparent'

    ::

        sage: D.element_class
        <class '__main__.myparent_with_category.element_class'>
        sage: D.element_class.mro()
        [<class '__main__.myparent_with_category.element_class'>,
        <class __main__.Element at ...>,
        <class 'sage.categories.category.Ds.element_class'>,
        <class 'sage.categories.category.Cs.element_class'>,
        <class 'sage.categories.category.Bs.element_class'>,
        <class 'sage.categories.category.As.element_class'>,
        <type 'object'>]


    TESTS::

        sage: import __main__
        sage: __main__.myparent = myparent
        sage: __main__.As = As
        sage: __main__.Bs = Bs
        sage: __main__.Cs = Cs
        sage: __main__.Ds = Ds
        sage: loads(dumps(Ds)) is Ds
        True
        sage: loads(dumps(Ds())) is Ds()
        True
        sage: loads(dumps(Ds().element_class)) is Ds().element_class
        True

    """
    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Input mangling for unique representation.

        Let ``C = Cs(...)`` be a category. Since :trac:`12895`, the
        class of ``C`` is a dynamic subclass ``Cs_with_category`` of
        ``Cs`` in order for ``C`` to inherit code from the
        ``SubcategoryMethods`` nested classes of its super categories.

        The purpose of this ``__classcall__`` method is to ensure that
        reconstructing ``C`` from its class with
        ``Cs_with_category(...)`` actually calls properly ``Cs(...)``
        and gives back ``C``.

        .. SEEALSO:: :meth:`subcategory_class`

        EXAMPLES::

            sage: A = Algebras(QQ)
            sage: A.__class__
            <class 'sage.categories.algebras.Algebras_with_category'>
            sage: A is Algebras(QQ)
            True
            sage: A is A.__class__(QQ)
            True
        """
        if isinstance(cls, DynamicMetaclass):
            cls = cls.__base__
        return super(Category, cls).__classcall__(cls, *args, **options)

    def __init__(self, s=None):
        """
        Initializes this category.

        INPUT:

        - ``s`` -- (Default: ``None``) A string giving the name of this
          category. If ``None``, the name is determined from the name of
          the class.

        EXAMPLES::

            sage: class SemiprimitiveRings(Category):
            ...       def super_categories(self):
            ...           return [Rings()]
            ...
            ...       class ParentMethods:
            ...           def jacobson_radical(self):
            ...               return self.ideal(0)
            ...
            sage: C = SemiprimitiveRings("SPR")
            sage: C
            Category of SPR
            sage: C.__class__
            <class '__main__.SemiprimitiveRings_with_category'>
        """
        if s is not None:
            if isinstance(s, str):
                self._label = s
                self.__repr_object_names = s
            else:
                raise TypeError, "Argument string must be a string."
        self.__class__ = dynamic_class("%s_with_category"%self.__class__.__name__,
                                       (self.__class__, self.subcategory_class, ),
                                       cache = False, reduction = None,
                                       doccls=self.__class__)

    @lazy_attribute
    def _label(self):
        """
        A short name of self, obtained from its type.

        EXAMPLES::

            sage: Rings()._label
            'Rings'

        """
        t = str(self.__class__.__base__)
        t = t[t.rfind('.')+1:]
        return t[:t.rfind("'")]

    @lazy_attribute
    def __repr_object_names(self):
        """
        Determine the name of the objects of this category
        from its type, if it has not been explicitly given
        at initialisation.

        EXAMPLES::

            sage: Rings()._Category__repr_object_names
            'rings'
            sage: PrincipalIdealDomains()._Category__repr_object_names
            'principal ideal domains'

            sage: Rings()
            Category of rings
        """
        i = -1
        s = self._label
        while i < len(s)-1:
            for i in range(len(s)):
                if s[i].isupper():
                    s = s[:i] + " " + s[i].lower() + s[i+1:]
                    break
        return s.lstrip()

    def _repr_object_names(self):
        """
        Returns the name of the objects of this category

        EXAMPLES::

            sage: FiniteGroups()._repr_object_names()
            'finite groups'
            sage: AlgebrasWithBasis(QQ)._repr_object_names()
            'algebras with basis over Rational Field'
        """
        return self.__repr_object_names

    def _short_name(self):
        """
        Returns a CamelCase name for this category

        EXAMPLES::

            sage: FiniteGroups()._short_name()
            'FiniteGroups'

            sage: AlgebrasWithBasis(QQ)._short_name()
            'AlgebrasWithBasis'

        Conventions for short names should be discussed at the level
        of Sage, and only then applied accordingly here.
        """
        return self._label

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: Rings.an_instance()
            Category of rings

        Parametrized categories should overload this default
        implementation to provide appropriate arguments::

            sage: Algebras.an_instance()
            Category of algebras over Rational Field
            sage: Bimodules.an_instance()
            Category of bimodules over Rational Field on the left and Real Field with 53 bits of precision on the right
            sage: AlgebraIdeals.an_instance()
            Category of algebra ideals in Univariate Polynomial Ring in x over Rational Field
        """
        return cls()

    def __call__(self, x, *args, **opts):
        """
        Constructs an object in this category from the data in ``x``,
        or throws ``TypeError`` or ``NotImplementedError``.

        If ``x`` is readily in ``self`` it is returned unchanged.
        Categories wishing to extend this minimal behavior should
        implement :meth:`._call_`.

        EXAMPLES::

            sage: Rings()(ZZ)
            Integer Ring
        """
        if x in self:
            return x
        return self._call_(x, *args, **opts)

    def _call_(self, x):
        """
        Constructs an object in this category from the data in ``x``,
        or throws NotImplementedError.

        EXAMPLES::

            sage: Semigroups()._call_(3)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_(self):
        """
        Returns the print representation of this category.

        EXAMPLES::

            sage: Sets() #indirect doctest
            Category of sets
        """
        return "Category of %s"%self._repr_object_names()

    def _latex_(self):
        r"""
        Returns the latex representation of this category.

        EXAMPLES::

            sage: latex(Sets()) #indirect doctest
            \mathbf{Sets}
        """
        return "\\mathbf{%s}"%self._short_name()

#   The convention for which hash function to use should be decided at the level of UniqueRepresentation
#   The implementation below is bad (hash independent of the base ring)
#     def __hash__(self):
#         """
#         Returns a hash for this category.
#
#         Currently this is just the hash of the string representing the category.
#
#         EXAMPLES::
#
#             sage: hash(Algebras(QQ)) #indirect doctest
#             699942203
#             sage: hash(Algebras(ZZ))
#             699942203
#         """
#         return hash(self.__category) # Any reason not to use id?

    def _subcategory_hook_(self, category):
        """
        Quick subcategory check.

        INPUT:

        - ``category`` -- a category

        OUTPUT:

        - ``True``, if ``category`` is a subcategory of ``self``.
        - ``False``, if ``category`` is not a subcategory of ``self``.
        - ``Unknown``, if a quick check was not enough to determine
          whether ``category`` is a subcategory of ``self`` or not.

        The aim of this method is to offer a framework to add cheap
        tests for subcategories. When doing
        ``category.is_subcategory(self)`` (note the reverse order of
        ``self`` and ``category``), this method is usually called
        first.  Only if it returns ``Unknown``, :meth:`is_subcategory`
        will build the list of super categories of ``category``.

        This method need not to handle the case where ``category`` is
        ``self``; this is the first test that is done in
        :meth:`is_subcategory`.

        This default implementation tests whether the parent class of
        ``category`` is a subclass of the parent class of ``self``.
        This is most of the time a complete subcategory test.

        .. WARNING::

            This test is incomplete for categories in
            :class:`CategoryWithParameters`, as introduced by
            :trac:`11935`. This method is therefore overwritten by
            :meth:`~sage.categories.category.CategoryWithParameters._subcategory_hook_`.

        EXAMPLES::

            sage: Rings()._subcategory_hook_(Rings())
            True
        """
        return issubclass(category.parent_class, self.parent_class)

    def __contains__(self, x):
        """
        Membership testing

        Returns whether ``x`` is an object in this category, that is
        if the category of ``x`` is a subcategory of ``self``.

        EXAMPLES::

            sage: ZZ in Sets()
            True
        """
        try:
            c = x.category()
        except AttributeError:
            return False
        return c.is_subcategory(self)

    @staticmethod
    def __classcontains__(cls, x):
        """
        Membership testing, without arguments

        INPUT:

        - ``cls`` -- a category class
        - ``x`` -- any object

        Returns whether ``x`` is an object of a category which is an instance
        of ``cls``.

        EXAMPLES:

        This method makes it easy to test if an object is, say, a
        vector space, without having to specify the base ring::

            sage: F = FreeModule(QQ,3)
            sage: F in VectorSpaces
            True

            sage: F = FreeModule(ZZ,3)
            sage: F in VectorSpaces
            False

            sage: F in Algebras
            False

        TESTS:

        Non category objects shall be handled properly::

            sage: [1,2] in Algebras
            False
        """
        try:
            c = x.categories()
        except AttributeError:
            return False
        return any(isinstance(cat, cls) for cat in c)

    def is_abelian(self):
        """
        Returns whether this category is abelian.

        An abelian category is a category satisfying:

        - It has a zero object;
        - It has all pullbacks and pushouts;
        - All monomorphisms and epimorphisms are normal.

        Equivalently, one can define an increasing sequence of conditions:

        - A category is pre-additive if it is enriched over abelian groups (all homsets are abelian groups and composition is bilinear);
        - A pre-additive category is additive if every finite set of objects has a biproduct (we can form direct sums and direct products);
        - An additive category is pre-abelian if every morphism has both a kernel and a cokernel;
        - A pre-abelian category is abelian if every monomorphism is the kernel of some morphism and every epimorphism is the cokernel of some morphism.

        EXAMPLES::

            sage: Modules(ZZ).is_abelian()
            True
            sage: FreeModules(ZZ).is_abelian()
            False
            sage: FreeModules(QQ).is_abelian()
            True
            sage: CommutativeAdditiveGroups().is_abelian()
            True
            sage: Semigroups().is_abelian()
            Traceback (most recent call last):
            NotImplementedError: is_abelian
        """
        raise NotImplementedError("is_abelian")

    def category_graph(self):
         r"""
         Returns the graph of all super categories of this category

         EXAMPLES::

             sage: C = Algebras(QQ)
             sage: G = C.category_graph()
             sage: G.is_directed_acyclic()
             True
             sage: G.girth()
             4
         """
         return category_graph([self])

    @abstract_method
    def super_categories(self):
        """
        Returns the *immediate* super categories of ``self``

        Every category should implement this method.

        EXAMPLES::

            sage: Groups().super_categories()
            [Category of monoids]
            sage: Objects().super_categories()
            []

        .. note::

            Mathematically speaking, the order of the super categories
            should be irrelevant. However, in practice, this order
            influences the result of :meth:`all_super_categories`, and
            accordingly of the method resolution order for parent and
            element classes. Namely, since ticket 11943, Sage uses the
            same `C3` algorithm for determining the order on the list
            of *all* super categories as Python is using for the
            method resolution order of new style classes.

        .. note::

            Whenever speed matters, developers are advised to use the
            lazy attribute :meth:`_super_categories` instead of
            calling this method.
        """

    @lazy_attribute
    def _all_super_categories(self):
        r"""
        All the super categories of this category, including this category.

        Since :trac:`11943`, the order of super categories is
        determined by Python's method resolution order C3 algorithm.

        .. seealso:: :meth:`all_super_categories`

        .. note:: this attribute is likely to eventually become a tuple.

        .. note:: this sets :meth:`_super_categories_for_classes` as a side effect

        EXAMPLES::

            sage: C = Rings(); C
            Category of rings
            sage: C._all_super_categories
            [Category of rings,
             Category of rngs,
             Category of semirings,
             Category of monoids,
             Category of semigroups,
             Category of magmas,
             Category of commutative additive groups,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """
        (result, bases) = C3_sorted_merge([cat._all_super_categories
                                           for cat in self._super_categories] +
                                          [self._super_categories],
                                          category_sort_key)
        self._super_categories_for_classes = bases
        return [self] + result

    @lazy_attribute
    def _all_super_categories_proper(self):
        r"""
        All the proper super categories of this category.

        Since :trac:`11943`, the order of super categories is
        determined by Python's method resolution order C3 algorithm.

        .. seealso:: :meth:`all_super_categories`

        .. note:: this attribute is likely to eventually become a tuple.

        EXAMPLES::

            sage: C = Rings(); C
            Category of rings
            sage: C._all_super_categories_proper
            [Category of rngs,
             Category of semirings,
             Category of monoids,
             Category of semigroups,
             Category of magmas,
             Category of commutative additive groups,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """
        return self._all_super_categories[1:]

    @lazy_attribute
    def _set_of_super_categories(self):
        """
        The frozen set of all proper super categories of this category.

        .. note:: this is used for speeding up category containment tests.

        .. seealso:: :meth:`all_super_categories`

        EXAMPLES::

            sage: Groups()._set_of_super_categories
            frozenset([...])
            sage: sorted(Groups()._set_of_super_categories, key=repr)
            [Category of magmas, Category of monoids, Category of objects, Category of semigroups,
             Category of sets, Category of sets with partial maps]

        TESTS::

            sage: C = HopfAlgebrasWithBasis(GF(7))
            sage: C._set_of_super_categories == frozenset(C._all_super_categories_proper)
            True
        """
        return frozenset(self._all_super_categories_proper)

    def all_super_categories(self, proper=False):
        """
        Returns the list of all super categories of this category.

        INPUT:

         - ``proper`` -- a boolean (default: ``False``); whether to exclude this category.

        Since :trac:`11943`, the order of super categories is
        determined by Python's method resolution order C3 algorithm.

        .. note::

            Whenever speed matters, the developers are advised to use
            instead the lazy attributes :meth:`_all_super_categories`,
            :meth:`_all_super_categories_proper`, or
            :meth:`_set_of_all_super_categories`, as
            appropriate. Simply because lazy attributes are much
            faster than any method.

        EXAMPLES::

            sage: C = Rings(); C
            Category of rings
            sage: C.all_super_categories()
            [Category of rings,
             Category of rngs,
             Category of semirings,
             Category of monoids,
             Category of semigroups,
             Category of magmas,
             Category of commutative additive groups,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]

            sage: C.all_super_categories(proper = True)
            [Category of rngs,
             Category of semirings,
             Category of monoids,
             Category of semigroups,
             Category of magmas,
             Category of commutative additive groups,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]

            sage: Sets().all_super_categories()
            [Category of sets, Category of sets with partial maps, Category of objects]
            sage: Sets().all_super_categories(proper=True)
            [Category of sets with partial maps, Category of objects]
            sage: Sets().all_super_categories() is Sets()._all_super_categories
            True
            sage: Sets().all_super_categories(proper=True) is Sets()._all_super_categories_proper
            True

        """
        if proper:
            return self._all_super_categories_proper
        return self._all_super_categories

    @lazy_attribute
    def _super_categories(self):
        """
        The immediate super categories of this category.

        This lazy attributes caches the result of the mandatory method
        :meth:`super_categories` for speed. It also does some mangling
        (flattening join categories, sorting, ...).

        Whenever speed matters, developers are advised to use this
        lazy attribute rather than calling :meth:`super_categories`.

        .. NOTE:: this attribute is likely to eventually become a tuple.

        EXAMPLES::

            sage: Rings()._super_categories
            [Category of rngs, Category of semirings]
        """
        return sorted(Category._flatten_categories(self.super_categories()), key = category_sort_key, reverse=True)

    @lazy_attribute
    def _super_categories_for_classes(self):
        """
        The super categories of this category used for building classes.

        This is a close variant of :meth:`_super_categories` used for
        constructing the list of the bases for :meth:`parent_class`,
        :meth:`element_class`, and friends. The purpose is ensure that
        Python will find a proper Method Resolution Order for those
        classes. For background, see :mod:`sage.misc.c3_controlled`.

        .. SEEALSO:: :meth:`_cmp_key`.

        .. NOTE::

            This attribute is calculated as a by-product of computing
            :meth:`_all_super_categories`.

        EXAMPLES::

            sage: Rings()._super_categories_for_classes
            [Category of rngs, Category of semirings]
        """
        self._all_super_categories
        return self._super_categories_for_classes

    def _test_category_graph(self, **options):
        """
        Check that the category graph matches with Python's method resolution order

        .. note::

            By :trac:`11943`, the list of categories returned by
            :meth:`all_super_categories` is supposed to match with the
            method resolution order of the parent and element
            classes. This method checks this.

        .. todo:: currently, this won't work for hom categories.

        EXAMPLES::

            sage: C = HopfAlgebrasWithBasis(QQ)
            sage: C.parent_class.mro() == [X.parent_class for X in C._all_super_categories] + [object]
            True
            sage: C.element_class.mro() == [X.element_class for X in C._all_super_categories] + [object]
            True
            sage: TestSuite(C).run()    # indirect doctest

        """
        tester = self._tester(**options)
        tester.assert_(self.parent_class.mro() == [C.parent_class for C in self._all_super_categories] + [object])
        tester.assert_(self.element_class.mro() == [C.element_class for C in self._all_super_categories] + [object])

#    def construction(self):
#        return (self.__class__,)

#    def __reduce__(self):
#        construction = self.construction()
#        return (construction[0], construction[1:])

    class ParentMethods:
        """
        Put methods for parents here.
        """
        pass

    class ElementMethods:
        """
        Put methods for elements here.
        """
        pass

    _cmp_key = _cmp_key

    def _make_named_class(self, name, method_provider, cache=False, picklable=True):
        """
        Construction of the parent/element/... class of ``self``.

        INPUT:

        - ``name`` -- a string; the name of the class as an attribute of
          ``self``. E.g. "parent_class"
        - ``method_provider`` -- a string; the name of an attribute of
          ``self`` that provides methods for the new class (in
          addition to those coming from the super categories).
          E.g. "ParentMethods"
        - ``cache`` -- a boolean or ``ignore_reduction`` (default: ``False``)
          (passed down to dynamic_class; for internal use only)
        - ``picklable`` -- a boolean (default: ``True``)

        ASSUMPTION:

        It is assumed that this method is only called from a lazy
        attribute whose name coincides with the given ``name``.

        OUTPUT:

        A dynamic class with bases given by the corresponding named
        classes of ``self``'s super_categories, and methods taken from
        the class ``getattr(self,method_provider)``.

        .. NOTE::

            - In this default implementation, the reduction data of
              the named class makes it depend on ``self``. Since the
              result is going to be stored in a lazy attribute of
              ``self`` anyway, we may as well disable the caching in
              ``dynamic_class`` (hence the default value
              ``cache=False``).

            - :class:`CategoryWithParameters` overrides this method so
              that the same parent/element/... classes can be shared
              between closely related categories.

            - The bases of the named class may also contain the named
              classes of some indirect super categories, according to
              :meth:`_super_categories_for_classes`. This is to
              guarantee that Python will build consistent method
              resolution orders. For background, see
              :mod:`sage.misc.c3_controlled`.

        .. SEEALSO:: :meth:`CategoryWithParameters._make_named_class`

        EXAMPLES::

            sage: PC = Rings()._make_named_class("parent_class", "ParentMethods"); PC
            <class 'sage.categories.rings.Rings.parent_class'>
            sage: type(PC)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>
            sage: PC.__bases__
            (<class 'sage.categories.rngs.Rngs.parent_class'>,
             <class 'sage.categories.semirings.Semirings.parent_class'>)

        Note that, by default, the result is not cached::

            sage: PC is Rings()._make_named_class("parent_class", "ParentMethods")
            False

        Indeed this method is only meant to construct lazy attributes
        like ``parent_class`` which already handle this caching::

            sage: Rings().parent_class
            <class 'sage.categories.rings.Rings.parent_class'>

        Reduction for pickling also assumes the existence of this lazy
        attribute::

            sage: PC._reduction
            (<built-in function getattr>, (Category of rings, 'parent_class'))
            sage: loads(dumps(PC)) is Rings().parent_class
            True

        TESTS::

            sage: class A: pass
            sage: class BrokenCategory(Category):
            ....:     def super_categories(self): return []
            ....:     ParentMethods = 1
            ....:     class ElementMethods(A):
            ....:         pass
            ....:     class MorphismMethods(object):
            ....:         pass
            sage: C = BrokenCategory()
            sage: C._make_named_class("parent_class",   "ParentMethods")
            Traceback (most recent call last):
            ...
            AssertionError: BrokenCategory.ParentMethods should be a class
            sage: C._make_named_class("element_class",  "ElementMethods")
            doctest:...: UserWarning: BrokenCategory.ElementMethods should not have a super class
            <class '__main__.BrokenCategory.element_class'>
            sage: C._make_named_class("morphism_class", "MorphismMethods")
            <class '__main__.BrokenCategory.morphism_class'>
        """
        cls = self.__class__
        if isinstance(cls, DynamicMetaclass):
            cls = cls.__base__
        class_name = "%s.%s"%(cls.__name__, name)
        method_provider_cls = getattr(self, method_provider, None)
        if method_provider_cls is None:
            # If the category provides no XXXMethods class,
            # point to the documentation of the category itself
            doccls = cls
        else:
            # Otherwise, check XXXMethods
            assert inspect.isclass(method_provider_cls),\
                "%s.%s should be a class"%(cls.__name__, method_provider)
            mro = inspect.getmro(method_provider_cls)
            if len(mro) > 2 or (len(mro) == 2 and mro[1] is not object):
                from warnings import warn
                warn("%s.%s should not have a super class"%(cls.__name__, method_provider))
            # and point the documentation to it
            doccls = method_provider_cls
        if picklable:
            reduction = (getattr, (self, name))
        else:
            reduction = None
        return dynamic_class(class_name,
                             tuple(getattr(cat,name) for cat in self._super_categories_for_classes),
                             method_provider_cls, prepend_cls_bases = False, doccls = doccls,
                             reduction = reduction, cache = cache)


    @lazy_attribute
    def subcategory_class(self):
        """
        A common superclass for all subcategories of this category (including this one).

        This class derives from ``D.subcategory_class`` for each super
        category `D` of ``self``, and includes all the methods from
        the nested class ``self.SubcategoryMethods``, if it exists.

        .. SEEALSO::

            - :trac:`12895`
            - :meth:`parent_class`
            - :meth:`element_class`
            - :meth:`_make_named_class`

        EXAMPLES::

            sage: cls = Rings().subcategory_class; cls
            <class 'sage.categories.rings.Rings.subcategory_class'>
            sage: type(cls)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>

        ``Rings()`` is an instance of this class, as well as all its subcategories::

            sage: isinstance(Rings(), cls)
            True
            sage: isinstance(AlgebrasWithBasis(QQ), cls)
            True

        TESTS::

            sage: cls = Algebras(QQ).subcategory_class; cls
            <class 'sage.categories.algebras.Algebras.subcategory_class'>
            sage: type(cls)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>

        """
        return self._make_named_class('subcategory_class', 'SubcategoryMethods',
                                      cache=False, picklable=False)

    @lazy_attribute
    def parent_class(self):
        """
        A common super class for all parents in this category.

        EXAMPLES::

            sage: C = Algebras(QQ).parent_class; C
            <class 'sage.categories.algebras.Algebras.parent_class'>
            sage: type(C)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>

        By :trac:`11935`, some categories share their parent
        classes. For example, the parent class of an algebra only
        depends on the category of the base ring. A typical example is
        the category of algebras over a finite field versus algebras
        over a non-field::

            sage: Algebras(GF(7)).parent_class is Algebras(GF(5)).parent_class
            True
            sage: Algebras(QQ).parent_class is Algebras(ZZ).parent_class
            False
            sage: Algebras(ZZ['t']).parent_class is Algebras(ZZ['t','x']).parent_class
            True

        See :class:`CategoryWithParameters` for an abstract base class for
        categories that depend on parameters, even though the parent
        and element classes only depend on the parent or element
        classes of its super categories. It is used in
        :class:`~sage.categories.bimodules.Bimodules`,
        :class:`~sage.categories.category_types.Category_over_base` and
        :class:`sage.categories.category.JoinCategory`.
        """
        return self._make_named_class('parent_class', 'ParentMethods')

    @lazy_attribute
    def element_class(self):
        """
        A common super class for all elements of parents in this category.

        EXAMPLES::

            sage: C = Algebras(QQ).element_class; C
            <class 'sage.categories.algebras.Algebras.element_class'>
            sage: type(C)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>

        By :trac:`11935`, some categories share their element
        classes. For example, the element class of an algebra only
        depends on the category of the base. A typical example is the
        category of algebras over a field versus algebras over a
        non-field::

            sage: Algebras(GF(5)).element_class is Algebras(GF(3)).element_class
            True
            sage: Algebras(QQ).element_class is Algebras(ZZ).element_class
            False
            sage: Algebras(ZZ['t']).element_class is Algebras(ZZ['t','x']).element_class
            True

        .. SEEALSO:: :meth:`parent_class`
        """
        return self._make_named_class('element_class', 'ElementMethods')

    def required_methods(self):
        """
        Returns the methods that are required and optional for parents
        in this category and their elements.

        EXAMPLES::

            sage: Algebras(QQ).required_methods()
            {'parent': {'required': ['__contains__'], 'optional': []}, 'element': {'required': ['__nonzero__'], 'optional': ['_add_', '_mul_']}}
        """
        return { "parent"  : abstract_methods_of_class(self.parent_class),
                 "element" : abstract_methods_of_class(self.element_class) }


    # Operations on the lattice of categories
    def is_subcategory(self, c):
        """
        Returns True if self is naturally embedded as a subcategory of c.

        EXAMPLES::

            sage: AbGrps = CommutativeAdditiveGroups()
            sage: Rings().is_subcategory(AbGrps)
            True
            sage: AbGrps.is_subcategory(Rings())
            False

        The ``is_subcategory`` function takes into account the
        base.

        ::

            sage: M3 = VectorSpaces(FiniteField(3))
            sage: M9 = VectorSpaces(FiniteField(9, 'a'))
            sage: M3.is_subcategory(M9)
            False

        Join categories are properly handled::

            sage: CatJ = Category.join((CommutativeAdditiveGroups(), Semigroups()))
            sage: Rings().is_subcategory(CatJ)
            True

        ::

            sage: V3 = VectorSpaces(FiniteField(3))
            sage: POSet = PartiallyOrderedSets()
            sage: PoV3 = Category.join((V3, POSet))
            sage: A3 = AlgebrasWithBasis(FiniteField(3))
            sage: PoA3 = Category.join((A3, POSet))
            sage: PoA3.is_subcategory(PoV3)
            True
            sage: PoV3.is_subcategory(PoV3)
            True
            sage: PoV3.is_subcategory(PoA3)
            False
        """
        if c is self:
            return True
        subcat_hook = c._subcategory_hook_(self)
        if subcat_hook is Unknown:
            return c in self._set_of_super_categories
        return subcat_hook

    def or_subcategory(self, category = None):
        """
        INPUT:

        - ``category`` - a sub category of ``self``, tuple/list thereof, or ``None``

        OUTPUT:

        - a category

        Returns ``category`` or ``self`` if ``category`` is None.

        EXAMPLES::

            sage: Monoids().or_subcategory(Groups())
            Category of groups
            sage: Monoids().or_subcategory(None)
            Category of monoids

        If category is a list/tuple, then a join category is returned::

            sage: Monoids().or_subcategory((FiniteEnumeratedSets(), Groups()))
            Join of Category of groups and Category of finite enumerated sets

        An error if raised if category is not a subcategory of ``self``.
        ::

            sage: Monoids().or_subcategory(EnumeratedSets())
            Traceback (most recent call last):
            ...
            AssertionError: Subcategory of `Category of enumerated sets` required; got `Category of monoids`
        """
        if category is None:
            return self
        if isinstance(category, (tuple, list)):
            category = Category.join(category)
        assert isinstance(category, Category)
        assert category.is_subcategory(self), "Subcategory of `%s` required; got `%s`"%(category, self)
        return category

    def _is_subclass(self, c):
        """
        Same as is_subcategory, but c may also be the class of a
        category instead of a category.

        EXAMPLES::

            sage: Fields()._is_subclass(Rings)
            True
            sage: Algebras(QQ)._is_subclass(Modules)
            True
            sage: Algebras(QQ)._is_subclass(ModulesWithBasis)
            False
        """
        assert( isinstance(c, Category) or (issubclass(c.__class__, type) and issubclass(c, Category)) )
        if isinstance(c, Category):
            return self.is_subcategory(c)
        else:
            return any(isinstance(cat, c) for cat in self._all_super_categories)

    @cached_method
    def _meet_(self, other):
        """
        Returns the largest common subcategory of self and other:

        EXAMPLES::

            sage: Monoids()._meet_(Monoids())
            Category of monoids
            sage: Rings()._meet_(Rings())
            Category of rings
            sage: Rings()._meet_(Monoids())
            Category of monoids
            sage: Monoids()._meet_(Rings())
            Category of monoids

            sage: VectorSpaces(QQ)._meet_(Modules(ZZ))
            Category of commutative additive groups
            sage: Algebras(ZZ)._meet_(Algebras(QQ))
            Category of rings
            sage: Groups()._meet_(Rings())
            Category of monoids
            sage: Algebras(QQ)._meet_(Category.join([Fields(), ModulesWithBasis(QQ)]))
            Join of Category of rings and Category of vector spaces over Rational Field

        Note: abstractly, the category poset is a distributive
        lattice, so this is well defined; however, the subset of those
        categories actually implemented is not: we need to also
        include their join-categories.

        For example, the category of rings is *not* the join of the
        category of abelian groups and that of semi groups, just a
        subcategory of their join, since rings further require
        distributivity.

        For the meet computation, there may be several lowest common
        sub categories of self and other, in which case, we need to
        take the join of them all.

        FIXME:

        - If A is a subcategory of B, A has *more* structure than B,
          but then *less* objects in there. We should choose an
          appropriate convention for A<B. Using subcategory calls
          for A<B, but the current meet and join call for A>B.
        """
        if self is other: # useful? fast pathway
            return self
        elif self.is_subcategory(other):
            return other
        elif other.is_subcategory(self):
            # Useful fast pathway; try:
            # %time L = EllipticCurve('960d1').prove_BSD()
            return self
        else:
            return Category.join(self._meet_(sup) for sup in other._super_categories)

    @staticmethod
    def meet(categories):
        """
        Returns the meet of a list of categories

        INPUT:

        - ``categories`` - a non empty list (or iterable) of categories

        EXAMPLES::

            sage: Category.meet([Algebras(ZZ), Algebras(QQ), Groups()])
            Category of monoids

        That meet of an empty list should be a category which is a
        subcategory of all categories, which does not make practical sense::

            sage: Category.meet([])
            Traceback (most recent call last):
            ...
            ValueError: The meet of an empty list of categories is not implemented
        """
        categories = tuple(categories)
        if len(categories) == 0:
            raise ValueError, "The meet of an empty list of categories is not implemented"
        result = categories[0]
        for category in categories[1:]:
            result = result._meet_(category)
        return result

    @staticmethod
    def _flatten_categories(categories):
        """
        INPUT:

        - ``categories`` -- a list (or iterable) of categories

        Returns the tuple of categories in ``categories``, while
        flattening join categories

        EXAMPLES::

            sage: Category._flatten_categories([Algebras(QQ), Category.join([Monoids(), Coalgebras(QQ)]), Sets()])
            (Category of algebras over Rational Field, Category of monoids, Category of coalgebras over Rational Field, Category of sets)
        """
        # Invariant: the super categories of a category are not JoinCategories
        return tuple(cat
                     for category in categories
                     for cat in (category._super_categories if isinstance(category, JoinCategory) else (category,)))

    @staticmethod
    def join(categories, as_list = False):
        """
        Returns the join of the input categories in the lattice of categories

        INPUT:

        - a sequence of categories (FIXME: should this be a list or iterable?)
        - as_list: a boolean, False by default (keyword only)

        EXAMPLES::

            sage: J = Category.join((Groups(), CommutativeAdditiveMonoids())); J
            Join of Category of groups and Category of commutative additive monoids
            sage: J.super_categories()
            [Category of groups, Category of commutative additive monoids]
            sage: J.all_super_categories(proper=True)
            [Category of groups,
             Category of monoids,
             Category of semigroups,
             Category of magmas,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]

        This is an associative operation::

            sage: Category.join((Objects(), Sets(), Category.join((Monoids(), Sets(), Monoids())), Category.join((Objects(), CommutativeAdditiveGroups()))))
            Join of Category of monoids and Category of commutative additive groups

        The join of a single category is the category itself::

            sage: Category.join((Monoids(),))
            Category of monoids

        Similarly, the join of several mutually comparable categories is the smallest one::

            sage: Category.join((Sets(), Rings(), Monoids()))
            Category of rings

        If the optional parameter ``as_list`` is ``True``, this just
        returns the super categories of the join as a list, without
        constructing the join category itself::

            sage: Category.join((Groups(), CommutativeAdditiveMonoids()), as_list=True)
            [Category of groups, Category of commutative additive monoids]
            sage: Category.join((Sets(), Rings(), Monoids()), as_list=True)
            [Category of rings]

        """
        return _join(tuple(categories), as_list)

    def category(self):
        """
        Returns the category of this category. So far all categories
        are in the category of objects.

        EXAMPLES::

            sage: Sets().category()
            Category of objects
            sage: VectorSpaces(QQ).category()
            Category of objects
        """
        from objects import Objects
        return Objects()

    # For better code locality and to avoid import loops, those
    # assignements are achieved by the respective modules:
    #
    # TensorProducts    = tensor.TensorProducts
    # CartesianProducts = cartesian_products.CartesianProducts
    #
    # Subquotients      = subquotients.Subquotients
    # Subobjects        = subobjects.Subobjects
    # Quotients         = quotients.Quotients
    # IsomorphicObjects = isomorphic_objects.IsomorphicObjects
    #
    # DualObjects       = dual.DualObjects
    # Algebras          = algebra_functor.Algebras
    @cached_method
    def hom_category(self):
        """
        Returns the category for homsets between objects this category.

        A category which needs to give specific information about this
        category should provide a HomCategory class.

        To avoid generating billions of categories, if there is
        nothing specific for homsets of this category, then this just
        returns the join of the categories of homsets of the super
        categories.

        EXAMPLES::

            sage: Sets().hom_category()
            Category of hom sets in Category of sets

        """
        try: #if hasattr(self, "HomCategory"):
            return self.HomCategory(self)
        except AttributeError:
            return Category.join((category.hom_category() for category in self._super_categories))

    def example(self, *args, **keywords):
        """
        Returns an object in this category. Most of the time, this is a parent.

        This serves three purposes:

         - Give a typical example to better explain what the category is all about.
           (and by the way prove that the category is non empty :-) )
         - Provide a minimal template for implementing other objects in this category
         - Provide an object on which to test generic code implemented by the category

        For all those applications, the implementation of the object
        shall be kept to a strict minimum. The object is therefore not
        meant to be used for other applications; most of the time a
        full featured version is available elsewhere in Sage, and
        should be used insted.

        Technical note: by default ``FooBar(...).example()`` is
        constructed by looking up
        ``sage.categories.examples.foo_bar.Example`` and calling it as
        ``Example()``. Extra positional or named parameters are also
        passed down. For a category over base ring, the base ring is
        further passed down as an optional argument.

        Categories are welcome to override this default implementation.

        EXAMPLES::

            sage: Semigroups().example()
            An example of a semigroup: the left zero semigroup

            sage: Monoids().Subquotients().example()
            NotImplemented
        """
        if '.' in self.__class__.__name__:
            # this magic should not apply to nested categories like Monoids.Subquotients
            return NotImplemented
        module_name = self.__module__.replace("sage.categories", "sage.categories.examples")
        import sys
        try:
            __import__(module_name)
            module = sys.modules[module_name]
        except ImportError:
            return NotImplemented
        try:
            cls = module.Example
        except AttributeError:
            return NotImplemented
        # Add the base ring as optional argument if this is a category over base ring
        # This really should be in Category_over_base_ring.example,
        # but that would mean duplicating the documentation above.
        from category_types import Category_over_base_ring
        if isinstance(self, Category_over_base_ring): # Huh, smelly Run Time Type Checking, isn't it?
            if "base_ring" not in keywords:
                keywords["base_ring"]=self.base_ring()
        return cls(*args, **keywords)


def is_Category(x):
    """
    Returns True if x is a category.

    EXAMPLES::

        sage: sage.categories.category.is_Category(CommutativeAdditiveSemigroups())
        True
        sage: sage.categories.category.is_Category(ZZ)
        False
    """
    return isinstance(x, Category)

@cached_function
def category_sample():
    r"""
    Return a sample of categories.

    It is constructed by looking for all concrete category classes declared in
    ``sage.categories.all``, calling :meth:`Category.an_instance` on those and
    taking all their super categories.

    EXAMPLES::

        sage: from sage.categories.category import category_sample
        sage: sorted(category_sample(), key=str)
        [Category of G-sets for Symmetric group of order 8! as a permutation group,
         Category of Hecke modules over Rational Field,
         Category of additive magmas, ...,
         Category of fields, ...,
         Category of graded hopf algebras with basis over Rational Field, ...,
         Category of modular abelian varieties over Rational Field, ...,
         Category of simplicial complexes, ...,
         Category of vector spaces over Rational Field, ...,
         Category of weyl groups,...
    """
    import sage.categories.all
    abstract_classes_for_categories = [Category, HomCategory]
    return tuple(cls.an_instance()
                 for cls in sage.categories.all.__dict__.values()
                 if isinstance(cls, type) and issubclass(cls, Category) and cls not in abstract_classes_for_categories)

def category_graph(categories = None):
    """
    Returns the graph of the categories in Sage

    INPUT:

    - ``categories`` -- a list (or iterable) of categories

    If ``categories`` is specified, then the graph will contain the
    mentionned categories together with all their super
    categories. Otherwise the graph will contain (an instance of) each
    category in :mod:`sage.categories.all` (e.g. ``Algebras(QQ)`` for
    algebras).

    For readability, the names of the category are shortened, and in
    particular do not contain base rings.

    EXAMPLES::

        sage: G = sage.categories.category.category_graph(categories = [Rings()])
        sage: G.vertices()
        ['additive magmas',
         'commutative additive groups',
         'commutative additive monoids',
         'commutative additive semigroups',
         'magmas',
         'monoids',
         'objects',
         'rings',
         'rngs',
         'semigroups',
         'semirings',
         'sets',
         'sets with partial maps']
        sage: G.plot()

        sage: sage.categories.category.category_graph().plot()
    """
    from sage import graphs
    if categories is None:
        categories = category_sample()
    cats = set()
    for category in categories:
        for cat in category.all_super_categories():
            cats.add(cat)

    categories = cats
    g = graphs.digraph.DiGraph()
    for cat in categories:
        g.add_vertex(cat._repr_object_names())
        for source in categories:
            for target in source.super_categories():
                g.add_edge([source._repr_object_names(), target._repr_object_names()])
    return g

#############################################################
# Homsets categories
#############################################################

class HomCategory(Category):
    """
    An abstract base class for all categories of homsets

    .. todo::

        Get a consistent hierarchy of homset categories. Currently, it
        is built in parallel to that of their base categories (which
        is plain wrong!!!)

    """
    def __init__(self, category, name=None):
        """
        Initializes this HomCategory

        INPUT:
         - ``category`` -- the category whose Homsets are the objects of this category.
         - ``name`` -- An optional name for this category.

        EXAMPLES:

        We need to skip one test, since the hierarchy of hom categories isn't
        consistent yet::

            sage: C = sage.categories.category.HomCategory(Rings()); C
            Category of hom sets in Category of rings
            sage: TestSuite(C).run(skip=['_test_category_graph'])
        """
        self.base_category = category
        Category.__init__(self, name)

    def _repr_object_names(self): # improve?
        """
        Print representation.

        EXAMPLES::

            sage: Sets().hom_category() #indirect doctest
            Category of hom sets in Category of sets
        """
        return "hom sets in %s"%self.base_category

    @cached_method
    def base(self):
        """
        If this hom-category is subcategory of a category with a base, return that base.

        EXAMPLES::

            sage: ModulesWithBasis(ZZ).hom_category().base()
            Integer Ring

        """
        from sage.categories.category_types import Category_over_base
        for C in self._all_super_categories_proper:
            if isinstance(C,Category_over_base):
                return C.base()
        raise AttributeError, "This hom category has no base"

    def super_categories(self):
        """
        Returns the immediate super categories, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: HomCategory(Sets()).super_categories()
            [Category of hom sets in Category of sets with partial maps]
        """
        return Category.join(self.extra_super_categories() +
                             [category.hom_category()
                              for category in self.base_category._super_categories],
                             as_list=True)
    @cached_method
    def extra_super_categories(self):
        """
        The super categories of self that are not derived from the
        inheritance diagram of the base category, as a list.

        EXAMPLES::

            sage: HomCategory(Sets()).extra_super_categories()
            []
        """
        return []


##############################################################################
# Parametrized categories whose parent/element class depend only on
# the super categories
##############################################################################

class CategoryWithParameters(Category):
    """
    A parametrized category whose parent/element classes depend only on
    its super categories.

    Many categories in Sage are parametrized, like ``C = Algebras(K)``
    which takes a base ring as parameter. In many cases, however, the
    operations provided by ``C`` in the parent class and element class
    depend only on the super categories of ``C``. For example, the
    vector space operations are provided if and only if ``K`` is a
    field, since ``VectorSpaces(K)`` is a super category of ``C`` only
    in that case. In such cases, and as an optimization (see :trac:`11935`),
    we want to use the same parent and element class for all fields.
    This is the purpose of this abstract class.

    Currently, :class:`~sage.categories.category.JoinCategory`,
    :class:`~sage.categories.category_types.Category_over_base` and
    :class:`~sage.categories.bimodules.Bimodules` inherit from this
    class.

    EXAMPLES::

        sage: C1 = Algebras(GF(5))
        sage: C2 = Algebras(GF(3))
        sage: C3 = Algebras(ZZ)
        sage: from sage.categories.category import CategoryWithParameters
        sage: isinstance(C1, CategoryWithParameters)
        True
        sage: C1.parent_class is C2.parent_class
        True
        sage: C1.parent_class is C3.parent_class
        False
    """
    def _make_named_class(self, name, method_provider, cache = False, **options):
        """
        Return the parent/element/... class of ``self``.

        INPUT:

        - ``name`` -- a string; the name of the class as an attribute
          of ``self``
        - ``method_provider`` -- a string; the name of an attribute of
          ``self`` that provides methods for the new class (in
          addition to what comes from the super categories)

        ASSUMPTION:

        It is assumed that this method is only called from a lazy
        attribute whose name coincides with the given ``name``.

        OUTPUT:

        A dynamic class that has the corresponding named classes of
        the super categories of ``self`` as bases and contains the
        methods provided by ``getattr(self, method_provider)``.

        .. NOTE::

            This method overrides :meth:`Category._make_named_class`
            so that the returned class *only* depends on the
            corresponding named classes of the super categories and on
            the provided methods. This allows for sharing the named
            classes across closely related categories providing the
            same code to their parents, elements and so on.

        EXAMPLES:

        The categories of bimodules over the fields ``CC`` or ``RR``
        provide the same methods to their parents and elements::

            sage: Bimodules(ZZ,RR).parent_class is Bimodules(ZZ,RDF).parent_class #indirect doctest
            True
            sage: Bimodules(CC,ZZ).element_class is Bimodules(RR,ZZ).element_class
            True

        On the other hand, modules over a field have more methods than
        modules over a ring::

            sage: Modules(GF(3)).parent_class is Modules(ZZ).parent_class
            False
            sage: Modules(GF(3)).element_class is Modules(ZZ).element_class
            False

        For a more subtle example, one could possibly share the classes for
        ``GF(3)`` and ``GF(2^3, 'x')``, but this is not currently the case::

            sage: Modules(GF(3)).parent_class is Modules(GF(2^3,'x')).parent_class
            False

        This is because those two fields do not have the exact same category::

            sage: GF(3).category()
            Join of Category of finite fields and Category of subquotients of monoids and Category of quotients of semigroups
            sage: GF(2^3,'x').category()
            Category of finite fields

        Similarly for ``QQ`` and ``RR``::

            sage: QQ.category()
            Category of quotient fields
            sage: RR.category()
            Category of fields
            sage: Modules(QQ).parent_class is Modules(RR).parent_class
            False

        Some other cases where one could potentially share those classes::

            sage: Modules(GF(3),dispatch=False).parent_class  is Modules(ZZ).parent_class
            False
            sage: Modules(GF(3),dispatch=False).element_class is Modules(ZZ).element_class
            False

        TESTS::

            sage: PC = Algebras(QQ).parent_class; PC   # indirect doctest
            <class 'sage.categories.algebras.Algebras.parent_class'>
            sage: type(PC)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>
            sage: PC.__bases__
            (<class 'sage.categories.rings.Rings.parent_class'>,
             <class 'sage.categories.vector_spaces.VectorSpaces.parent_class'>)
            sage: loads(dumps(PC)) is PC
            True
        """
        cls = self.__class__
        if isinstance(cls, DynamicMetaclass):
            cls = cls.__base__
        key = (cls, name, self._make_named_class_key(name))
        try:
            return self._make_named_class_cache[key]
        except KeyError:
            pass
        result = Category._make_named_class(self, name, method_provider,
                                            cache=cache, **options)
        self._make_named_class_cache[key] = result
        return result


    @abstract_method
    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... class depend on.

        INPUT:

        - ``name`` -- a string; the name of the class as an attribute
          of ``self``

        .. SEEALSO::

            - :meth:`_make_named_class`
            - :meth:`sage.categories.category_types.Category_over_base._make_named_class_key`
            - :meth:`sage.categories.bimodules.Bimodules._make_named_class_key`
            - :meth:`JoinCategory._make_named_class_key`

        EXAMPLES:

        The parent class of an algebra depends only on the category of the base ring::

            sage: Algebras(ZZ)._make_named_class_key("parent_class")
            Category of euclidean domains

        The morphism class of a bimodule depends only on the category
        of the left and right base rings::

            sage: Bimodules(QQ, ZZ)._make_named_class_key("morphism_class")
            (Category of quotient fields, Category of euclidean domains)

        The element class of a join category depends only on the
        element class of its super categories::

            sage: Category.join([Groups(), Posets()])._make_named_class_key("element_class")
            (<class 'sage.categories.groups.Groups.element_class'>,
             <class 'sage.categories.posets.Posets.element_class'>)
        """

    _make_named_class_cache = dict()

    _cmp_key = _cmp_key_named

    def _subcategory_hook_(self, C):
        """
        A quick but partial test whether ``C`` is a subcategory of ``self``.

        INPUT:

        - ``C`` -- a category

        OUTPUT:

        ``False``, if ``C.parent_class`` is not a subclass of
        ``self.parent_class``, and :obj:`~sage.misc.unknown.Unknown`
        otherwise.

        EXAMPLES::

            sage: Bimodules(QQ,QQ)._subcategory_hook_(Modules(QQ))
            Unknown
            sage: Bimodules(QQ,QQ)._subcategory_hook_(Rings())
            False
        """
        if not issubclass(C.parent_class, self.parent_class):
            return False
        return Unknown


#############################################################
# Join of several categories
#############################################################

class JoinCategory(CategoryWithParameters):
    """
    A class for joins of several categories. Do not use directly;
    see Category.join instead.

    EXAMPLES::

        sage: from sage.categories.category import JoinCategory
        sage: J = JoinCategory((Groups(), CommutativeAdditiveMonoids())); J
        Join of Category of groups and Category of commutative additive monoids
        sage: J.super_categories()
        [Category of groups, Category of commutative additive monoids]
        sage: J.all_super_categories(proper=True)
        [Category of groups, Category of monoids, Category of semigroups, Category of magmas, Category of commutative additive monoids, Category of commutative additive semigroups, Category of additive magmas, Category of sets, Category of sets with partial maps, Category of objects]

    By :trac:`11935`, join categories and categories over base
    rings inherit from :class:`CategoryWithParameters`. This allows
    for sharing parent and element classes between similar
    categories. For example, since polynomial rings belong to a join
    category and since the underlying implementation is the same for
    all finite fields, we have::

        sage: GF(3)['x'].category()
        Join of Category of euclidean domains and Category of commutative algebras over Finite Field of size 3
        sage: type(GF(3)['x']) is type(GF(5)['z'])
        True
    """

    def __init__(self, super_categories, **kwds):
        """
        Initializes this JoinCategory

        INPUT:

        - super_categories -- Categories to join.  This category will
          consist of objects and morphisms that lie in all of these
          categories.

        - name -- An optional name for this category.

        TESTS::

            sage: from sage.categories.category import JoinCategory
            sage: C = JoinCategory((Groups(), CommutativeAdditiveMonoids())); C
            Join of Category of groups and Category of commutative additive monoids
            sage: TestSuite(C).run()

        """
        assert(len(super_categories) >= 2)
        assert(all(not isinstance(category, JoinCategory) for category in super_categories))
        # Use __super_categories to not overwrite the lazy attribute Category._super_categories
        # Maybe this would not be needed if the flattening/sorting is does consistently?
        self.__super_categories = list(super_categories)
        if kwds.has_key('name'):
            Category.__init__(self, kwds['name'])
        else:
            Category.__init__(self)

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        Since :trac:`11935`, the element/parent classes of a join
        category over base only depend on the element/parent class of
        its super categories.

        .. SEEALSO::

            - :meth:`CategoryWithParameters`
            - :meth:`CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Modules(ZZ)._make_named_class_key('element_class')
            Category of euclidean domains
            sage: Modules(QQ)._make_named_class_key('parent_class')
            Category of quotient fields
            sage: Schemes(Spec(ZZ))._make_named_class_key('parent_class')
            Category of Schemes
            sage: ModularAbelianVarieties(QQ)._make_named_class_key('parent_class')
            Category of quotient fields
        """
        return tuple(getattr(cat, name) for cat in self._super_categories)

    def super_categories(self):
        """
        Returns the immediate super categories, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: from sage.categories.category import JoinCategory
            sage: JoinCategory((Semigroups(), FiniteEnumeratedSets())).super_categories()
            [Category of semigroups, Category of finite enumerated sets]
        """
        return self.__super_categories

    def _subcategory_hook_(self, category):
        """
        Returns whether ``category`` is a subcategory of this join category

        INPUT:

        - ``category`` -- a category.

        .. note::

            ``category`` is a sub-category of this join category if
            and only if it is a sub-category of all super categories
            of this join category.

        EXAMPLE::

            sage: QQ['x'].category().is_subcategory(Category.join([Rings(), VectorSpaces(QQ)]))  # indirect doctest
            True
        """
        return all(category.is_subcategory(X) for X in self._super_categories)

    def is_subcategory(self, C):
        """
        Check whether this join category is subcategory of another
        category ``C``.

        EXAMPLES::

            sage: Category.join([Rings(),Modules(QQ)]).is_subcategory(Category.join([Rngs(),Bimodules(QQ,QQ)]))
            True
        """
        if C is self:
            return True
        hook = C._subcategory_hook_(self)
        if hook is Unknown:
            return any(X.is_subcategory(C) for X in self._super_categories)
        return hook

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: Category.join((Groups(), CommutativeAdditiveMonoids())) #indirect doctest
            Join of Category of groups and Category of commutative additive monoids

        TODO: find a better place to implement this representation improvement:

            sage: Category.join((Groups(), Sets().Facades()))
            Category of facade groups
            sage: Category.join((Sets().Facades(), Groups()))
            Category of facade groups
        """
        categories = self._super_categories
        from sets_cat import Sets
        if len(categories) == 2 and Sets().Facades() in categories:
            categories = list(categories)
            categories.remove(Sets().Facades())
            return "Category of facade %s"%(categories[0]._repr_object_names())
        return "Join of %s"%(" and ".join(str(cat) for cat in categories))
