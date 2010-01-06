"""
Categories

AUTHORS:

- David Kohel, William Stein and Nicolas M. Thiery

Every Sage object lies in a category. Categories in Sage are
modeled on the mathematical idea of category, and are distinct from
Python classes, which are a programming construct.

In most cases, typing ``x.category()`` returns the
category to which `x` belongs. If `C` is a category
and `x` is any object, `C(x)` tries to make an
object in `C` from `x`.

See :class:`Category` and :mod:`sage.categories.primer` for more details.

EXAMPLES:

We create a couple of categories::

    sage: Sets()
    Category of sets
    sage: GSets(AbelianGroup([2,4,9]))
    Category of G-sets for Multiplicative Abelian Group isomorphic to C2 x C4 x C9
    sage: Semigroups()
    Category of semigroups
    sage: VectorSpaces(FiniteField(11))
    Category of vector spaces over Finite Field of size 11
    sage: Ideals(IntegerRing())
    Category of ring ideals in Integer Ring

The default category for elements `x` of an object `O` is the category
of all objects of `O`. For example::

    sage: V = VectorSpace(RationalField(), 3)
    sage: x = V.gen(1)
    sage: x.category()
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

from sage.misc.abstract_method import abstract_method, abstract_methods_of_class
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
#from sage.misc.misc import attrcall

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.dynamic_class import dynamic_class

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

        sage: from sage.categories.all import Category
        sage: from sage.misc.lazy_attribute import lazy_attribute
        sage: class As (Category):
        ...       @cached_method
        ...       def super_categories(self):
        ...           return []
        ...
        ...       class ParentMethods:
        ...           def fA(self):
        ...               return "A"
        ...           f = fA
        ...
        sage: class Bs (Category):
        ...       @cached_method
        ...       def super_categories(self):
        ...           return [As()]
        ...
        ...       class ParentMethods:
        ...           def fB(self):
        ...               return "B"
        ...
        sage: class Cs (Category):
        ...       @cached_method
        ...       def super_categories(self):
        ...           return [As()]
        ...
        ...       class ParentMethods:
        ...           def fC(self):
        ...               return "C"
        ...           f = fC
        ...
        sage: class Ds (Category):
        ...       @cached_method
        ...       def super_categories(self):
        ...           return [Bs(),Cs()]
        ...
        ...       class ParentMethods:
        ...           def fD(self):
        ...               return "D"
        ...

    Categories should always have uniq representation. We check
    this before proceeding:

        sage: id(As()) == id(As())
        True
        sage: As().parent_class == As().parent_class
        True

    We construct a parent in the category Ds() (that is an instance of
    Ds().parent_class, and check that it has access to all the
    methods provided by all the categories, with the appropriate
    inheritance order.

        sage: D = Ds().parent_class()
        sage: [ D.fA(), D.fB(), D.fC(), D.fD() ]
        ['A', 'B', 'C', 'D']
        sage: D.f()
        'C'

        sage: C = Cs().parent_class()
        sage: [ C.fA(), C.fC() ]
        ['A', 'C']
        sage: C.f()
        'C'

    Here is the parallel hierarchy of classes which has been built
    automatically, together with the method resolution order (.mro())::

        sage: As().parent_class
        <class '__main__.As.parent_class'>
        sage: As().parent_class.__bases__
        (<type 'object'>,)
        sage: As().parent_class.mro()
        [<class '__main__.As.parent_class'>, <type 'object'>]

        sage: Bs().parent_class
        <class '__main__.Bs.parent_class'>
        sage: Bs().parent_class.__bases__
        (<class '__main__.As.parent_class'>,)
        sage: Bs().parent_class.mro()
        [<class '__main__.Bs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>]

        sage: Cs().parent_class
        <class '__main__.Cs.parent_class'>
        sage: Cs().parent_class.__bases__
        (<class '__main__.As.parent_class'>,)
        sage: Cs().parent_class.__mro__
        (<class '__main__.Cs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>)

        sage: Ds().parent_class
        <class '__main__.Ds.parent_class'>
        sage: Ds().parent_class.__bases__
        (<class '__main__.Bs.parent_class'>, <class '__main__.Cs.parent_class'>)
        sage: Ds().parent_class.mro()
        [<class '__main__.Ds.parent_class'>, <class '__main__.Bs.parent_class'>, <class '__main__.Cs.parent_class'>, <class '__main__.As.parent_class'>, <type 'object'>]

    Note that that two categories in the same class need not have the
    same super_categories. For example, Algebras(QQ) has
    VectorSpaces(QQ) as super category, whereas Algebras(ZZ) only has
    Modules(ZZ) as super category. In particular, the constructed
    parent_class and element_class will differ (inheriting, or not,
    methods specific for vector spaces). On the other hand, caching
    ensures that two identical hierarchy of classes are built only
    once::

        # TODO: redo the same with Algebras
        # and show the mro for Algebras(QQ) w.r.t Algebras(ZZ)
        # 2009/03/11: this feature is temporarily broken, due to the current work around for pickling
        sage: Coalgebras(QQ).parent_class is Coalgebras(FractionField(QQ[x])).parent_class # todo: not implemented
        True

    We now construct a parent in the usual way:

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
        <class '__main__.Bs.parent_class'>,
        <class '__main__.Cs.parent_class'>,
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

        sage: D.element_class
        <class '__main__.myparent_with_category.element_class'>
        sage: D.element_class.mro()
        [<class '__main__.myparent_with_category.element_class'>,
        <class __main__.Element at ...>,
        <class 'sage.categories.category.Ds.element_class'>,
        <class 'sage.categories.category.Bs.element_class'>,
        <class 'sage.categories.category.Cs.element_class'>,
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
    def __init__(self, s=None):
        """
        Initializes this category.

        INPUT:
        - s -- A string giving the name of this category.  If None,
          the name is determined from the name of the class.

        EXAMPLES::

            sage: class SemiprimitiveRings(Category):
            ...       @cached_method
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
        """
        if s is None:  # figure out from the type name!
            t = str(type(self))
            t = t[t.rfind('.')+1:]
            s = t[:t.rfind("'")]
            self.__label = s
            i = -1
            while i < len(s)-1:
                for i in range(len(s)):
                    if s[i].isupper():
                        s = s[:i] + " " + s[i].lower() + s[i+1:]
                        break
            s = s.lstrip()
        elif isinstance(s, str):
            self.__label = s
        else:
            raise TypeError, "Argument string must be a string."
        self.__category = s

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: Rings.an_instance()
            Category of rings

        Parametrized categories should overload this default
        implementation to provide appropriate arguments:

            sage: Algebras.an_instance()
            Category of algebras over Rational Field
            sage: Bimodules.an_instance()
            Category of bimodules over Rational Field on the left and Real Field with 53 bits of precision on the right
            sage: AlgebraIdeals.an_instance()
            Category of algebra ideals in Univariate Polynomial Ring in x over Rational Field
        """
        return cls()

    def __call__(self, x):
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
        return self._call_(x)

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
        return "Category of %s"%self.__category

    def _latex_(self):
        """
        Returns the latex representation of this category.

        EXAMPLES::

            sage: latex(Sets()) #indirect doctest
            \mathbf{Sets}
        """
        return "\\mathbf{%s}"%self.__label

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

#     This method preexisted in the original category code and was barely used (only in functors)
#     Conventions for short names should be discussed at the level of Sage,
#     and only then applied accordingly here.
#     def short_name(self):
#         """
#         A CamelCase representation of this category.

#         EXAMPLES::

#             sage: C = Algebras(ZZ)
#             sage: C.short_name()
#             'Algebras'
#         """
#         return ''.join([x.capitalize() for x in self.__category.split()])

    def __contains__(self, x):
        """
        Returns whether ``x`` is an object in this category.

        More specifically, returns True if and only if ``x`` has a
        category which is a subcategory of this one.

        EXAMPLES::

            sage: ZZ in Sets()
            True
        """
        try:
            c = x.category()
        except AttributeError:
            return False
        return c.is_subcategory(self)

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

    @cached_method
    def all_super_categories(self, proper = False):
        r"""
        Returns a linear extension (topological sort) of all the
        (proper) super categories of this category, and cache the
        result.

        INPUT:

         - ``proper``: a boolean; defaults to False.  Whether to exclude this category.

        FIXME:

        - make sure that this is compatible with the python algorithm
          for method resolution and make it O(n+m)

        EXAMPLES::

            sage: C = GradedHopfAlgebrasWithBasis(QQ).abstract_category(); C
            Category of abstract graded hopf algebras with basis over Rational Field
            sage: C.all_super_categories()
            [Category of abstract graded hopf algebras with basis over Rational Field,
             Category of graded hopf algebras over Rational Field,
             Category of graded bialgebras over Rational Field,
             Category of graded algebras over Rational Field,
             Category of graded coalgebras over Rational Field,
             Category of graded modules over Rational Field,
             Category of hopf algebras over Rational Field,
             Category of bialgebras over Rational Field,
             Category of algebras over Rational Field,
             ...]
        """
#            done = set()
        all_categories = []
        def add_successors(cat): # Factor this out!
#                if cat in done:
#                    return;
#                done.add(cat)
            for super in cat.super_categories():
                all_categories.append(super)
                add_successors(super)
        add_successors(self)
        all_categories.reverse()
        done = set()
        linear_extension = []
        for cat in all_categories:
            if not cat in done:
                done.add(cat)
                linear_extension.append(cat)
        linear_extension.reverse()
        if proper:
            return linear_extension
        else:
            return [self] + linear_extension

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

    @lazy_attribute
    def parent_class(self):
        """
        A common super class for all parents in this category.

        EXAMPLES::

            sage: C = Algebras(QQ).parent_class; C
            <class 'sage.categories.algebras.Algebras.parent_class'>
            sage: type(C)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>
        """
        return dynamic_class("%s.parent_class"%self.__class__.__name__,
                             tuple(cat.parent_class for cat in self.super_categories()),
                             self.ParentMethods,
                             reduction = (getattr, (self, "parent_class")))

    @lazy_attribute
    def element_class(self):
        """
        A common super class for all elements of parents in this category.

        EXAMPLES::

            sage: C = Algebras(QQ).element_class; C
            <class 'sage.categories.algebras.Algebras.element_class'>
            sage: type(C)
            <class 'sage.structure.dynamic_class.DynamicMetaclass'>
        """
        return dynamic_class("%s.element_class"%self.__class__.__name__,
                             (cat.element_class for cat in self.super_categories()),
                             self.ElementMethods,
                             reduction = (getattr, (self, "element_class"))
                             )

    def required_methods(self):
        """
        Returns the methods that are required and optional for parents
        in this category and their elements.

        EXAMPLES::

            sage: Algebras(QQ).required_methods()
            {'parent': {'required': ['__contains__'], 'optional': []}, 'element': {'required': [], 'optional': ['_add_', '_mul_']}}

        """
        return { "parent"  : abstract_methods_of_class(self.parent_class),
                 "element" : abstract_methods_of_class(self.element_class) }


    # Operations on the lattice of categories

    def is_subcategory(self, c):
        """
        Returns True if self is naturally embedded as a subcategory of c.

        EXAMPLES::

            sage: Rings  = Rings()
            sage: AbGrps = CommutativeAdditiveGroups()
            sage: Rings.is_subcategory(AbGrps)
            True
            sage: AbGrps.is_subcategory(Rings)
            False

        The ``is_subcategory`` function takes into account the
        base.

        ::

            sage: M3 = VectorSpaces(FiniteField(3))
            sage: M9 = VectorSpaces(FiniteField(9, 'a'))
            sage: M3.is_subcategory(M9)
            False

        TODO: handle join categories properly::

            sage: Rings().is_subcategory(Category.join((CommutativeAdditiveGroups(), SemiGroups()))) # todo: not implemented
            True


        """
        assert(isinstance(c, Category))
#        print "%s in %s: %s"%(self, c, c in self.all_super_categories())
        return c in self.all_super_categories()

    def _is_subclass(self, c,):
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
            return any(isinstance(cat, c) for cat in self.all_super_categories())

    def meet(self, other):
        """
        Returns the largest common subcategory of self and other:

        EXAMPLES::

            sage: Monoids().meet(Monoids())
            Category of monoids
            sage: Rings().meet(Rings())
            Category of rings
            sage: Rings().meet(Monoids())
            Category of monoids
            sage: Monoids().meet(Rings())
            Category of monoids

            sage: VectorSpaces(QQ).meet(Modules(ZZ))
            Category of commutative additive groups
            sage: Algebras(ZZ).meet(Algebras(QQ))
            Category of rings
            sage: Groups().meet(Rings())
            Category of monoids

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
           appropriate convention for A<B.  using subcategory calls
           for A<B, but the current meet and join call for A>B.

         - :meth:`meet` should be consistent with :meth:`join` and
           take an iterable of categories
        """
        if self is other: # useful? fast pathway
            return self
        elif self.is_subcategory(other):
            return other
        else:
            return Category.join(self.meet(sup) for sup in other.super_categories())

    @staticmethod
    def join(categories, **options):
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
            sage: J.all_super_categories(proper = True)
            [Category of groups,
             Category of monoids,
             Category of semigroups,
             Category of commutative additive monoids,
             Category of commutative additive semigroups,
             Category of sets,
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

        If the optional parameter as_list is True, then just return
        the super categories of the join as a list, without
        constructing the join category itself::

            sage: Category.join((Groups(), CommutativeAdditiveMonoids()), as_list=True)
            [Category of groups, Category of commutative additive monoids]
            sage: Category.join((Sets(), Rings(), Monoids()), as_list=True)
            [Category of rings]

        """
        categories = tuple(categories)
        # Since Objects() is the top category, it is the neutral element of join
        if len(categories) == 0:
            from objects import Objects
            return Objects()

        # Ensure associativity by flattening JoinCategory's
        # Invariant: the super categories of a JoinCategory are not JoinCategories themselves
        categories = sum( (tuple(category.super_categories()) if isinstance(category, JoinCategory) else (category,)
                           for category in categories), ())

        # remove redundant categories which are super categories of others
        result = ()
        for category in categories:
            if any(cat.is_subcategory(category) for cat in result):
                continue
            result = tuple( cat for cat in result if not category.is_subcategory(cat) ) + (category,)
        if "as_list" in options and options["as_list"]:
            return list(result)
        if len(result) == 1:
            return result[0]
        else:
            return JoinCategory(result)

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
        if hasattr(self, "HomCategory"):
            return self.HomCategory(self)
        else:
            return Category.join((category.hom_category() for category in self.super_categories()))

    def abstract_category(self):
        r"""
        An abstract parent is a parent which models an abstract
        algebraic structure which has several concrete representations.

        This returns a mostly technical category which provides
        support tools for handling the different representation, and
        in particular the coercions between them.

        It can be manually specified by defining a class
        AbstractCategory as a member of this category.

        Typically, ``FiniteDimensionalModulesWithBasis(QQ).abstract_category()``
        will be in charge, whenever a coercion `\phi: A\mapsto B` is
        registered, to register `\phi^{-1}` as coercion `B \mapsto A`
        if there is none defined yet.

        This is the analog of the `*WithSeveralBases` categories in MuPAD-Combinat.

        TODO: find a better name!

        The hierarchy of all abstract categories is built in parallel
        to that of their base categories, optimizing away those
        categories which do not have an AbstractCategory.

        Design question: currently ``self.abstract_category()`` is a
        subcategory of self by default. Is this desirable? For example,
        ``Algebras(QQ).abstract_category()`` should definitely be a
        subcategory of ``Algebras(QQ)``. On the other hand,
        ``AlgebrasWithBasis(QQ).abstract_category()`` should be a
        subcategory of ``Algebras(QQ)``, but not of
        ``AlgebrasWithBasis(QQ)``. This is because
        ``AlgebrasWithBasis(QQ)`` is specifying something about the
        concrete representation.

        EXAMPLES::

            sage: Semigroups().abstract_category()
            Category of semigroups
            sage: C = GradedHopfAlgebrasWithBasis(QQ).abstract_category(); C
            Category of abstract graded hopf algebras with basis over Rational Field
            sage: C.all_super_categories()
            [Category of abstract graded hopf algebras with basis over Rational Field,
             Category of graded hopf algebras over Rational Field,
             Category of graded bialgebras over Rational Field,
             Category of graded algebras over Rational Field,
             Category of graded coalgebras over Rational Field,
             Category of graded modules over Rational Field,
             Category of hopf algebras over Rational Field,
             Category of bialgebras over Rational Field,
             Category of algebras over Rational Field,
             ...]

        """
        if hasattr(self, "AbstractCategory"):
            return self.AbstractCategory(self)
        else:
            return Category.join(([self]+[category.abstract_category() for category in self.super_categories()]))

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

        Technical note: by default FooBar(...).example() is
        constructed by looking up
        sage.categories.examples.foo_bar.Example and calling it as
        ``Example(category = FooBar)``. Extra positional or named
        parameters are also passed down. Categories are welcome to
        override this.

        EXAMPLES::

            sage: Semigroups().example()
            An example of a semigroup: the left zero semigroup
        """
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
        ['commutative additive groups',
        'commutative additive monoids',
        'commutative additive semigroups',
        'monoids',
        'objects',
        'rings',
        'rngs',
        'semigroups',
        'sets']
        sage: G.plot()

        sage: sage.categories.category.category_graph().plot()
    """
    from sage import graphs
    import sage.categories.all
    if categories is None:
        import all
        abstract_classes_for_categories = [Category, HomCategory, AbstractCategory, all.CartesianProductCategory, all.TensorCategory, all.CategoryWithCartesianProduct, all.CategoryWithTensorProduct, all.DualityCategory]
        categories = [ cat.an_instance() for cat in sage.categories.all.__dict__.values() if isinstance(cat, type) and issubclass(cat, Category) and cat not in abstract_classes_for_categories ]
    cats = set()
    for category in categories:
        for cat in category.all_super_categories():
            cats.add(cat)
    def name(cat):
        return repr(cat)[12:]
    categories = cats
    g = graphs.digraph.DiGraph()
    for cat in categories:
        g.add_vertex(name(cat))
        for source in categories:
            for target in source.super_categories():
                g.add_edge([name(source), name(target)])
    return g

#############################################################
# Homsets categories
#############################################################

class HomCategory(Category):
    """
    An abstract base class for all categories of homsets

    The hierarchy of homset categories is built in parallel to that of
    their base categories (which is plain wrong!!!)

    """
    def __init__(self, category, name=None):
        """
        Initializes this HomCategory

        INPUT:
         - ``category`` -- the category whose Homsets are the objects of this category.
         - ``name`` -- An optional name for this category.

        EXAMPLES::

            sage: C = sage.categories.category.HomCategory(Rings()); C
            Category of hom sets in Category of rings
            sage: TestSuite(C).run()
        """
        Category.__init__(self, name)
        self.base_category = category

    def _repr_(self): # improve?
        """
        Print representation.

        EXAMPLES::

            sage: Sets().hom_category() #indirect doctest
            Category of hom sets in Category of sets
        """
        return "Category of hom sets in %s"%self.base_category

#    def construction(self):
#        return (attrcall("hom_category"), self.base_category)

    @cached_method
    def super_categories(self):
        """
        Returns the immediate super categories, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: HomCategory(Sets()).super_categories()
            [Category of hom sets in Category of objects]
        """
        return Category.join(self.extra_super_categories() +
                             [category.hom_category()
                              for category in self.base_category.super_categories()],
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


#############################################################
# Categories of abstract parents
#############################################################

class AbstractCategory(Category):
    """
    An abstract base class for all categories of abstract parents

    See Category.abstract_category.

    Caveat: specifications subject to change shortly.
    """
    def __init__(self, category, name=None):
        """
        Initializes this AbstractCategory

        INPUT:
        - ``category`` -- the category of which this category forms the abstract category.
        - ``name`` -- An optional name for this category.

            sage: C = sage.categories.category.AbstractCategory(Sets())
            sage: C
            Category of abstract sets
            sage: TestSuite(C).run()
        """
        Category.__init__(self, name)
        self.base_category = category

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: C = GradedHopfAlgebrasWithBasis(QQ).abstract_category(); C #indirect doctest
            Category of abstract graded hopf algebras with basis over Rational Field
        """
        s = repr(self.base_category)
        return s[:11]+" abstract"+s[11:]

#    def construction(self):
#        return (attrcall("abstract_category"), self.base_category)

    @cached_method
    def super_categories(self):
        """
        Returns the immediate super categories, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: C = GradedHopfAlgebrasWithBasis(QQ).abstract_category()
            sage: C.super_categories()
            [Category of graded hopf algebras over Rational Field]
        """
        return Category.join(self.extra_super_categories() +
                             [category.abstract_category()
                              for category in self.base_category.super_categories()],
                             as_list=True)
    @cached_method
    def extra_super_categories(self):
        """
        The super categories of self that are not derived from the
        inheritance diagram of the base category, as a list.

        EXAMPLES::

            sage: C = GradedHopfAlgebrasWithBasis(QQ).abstract_category()
            sage: C.extra_super_categories()
            [Category of graded hopf algebras with basis over Rational Field]
        """
        return [self.base_category]



#############################################################
# Join of several categories
#############################################################

class JoinCategory(Category):
    """
    A class for joins of several categories. Do not use directly;
    see Category.join instead.

    EXAMPLES::

        sage: from sage.categories.category import JoinCategory
        sage: J = JoinCategory((Groups(), CommutativeAdditiveMonoids())); J
        Join of Category of groups and Category of commutative additive monoids
        sage: J.super_categories()
        [Category of groups, Category of commutative additive monoids]
        sage: J.all_super_categories(proper = True)
        [Category of groups, Category of monoids, Category of semigroups, Category of commutative additive monoids, Category of commutative additive semigroups, Category of sets, Category of objects]

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
        if kwds.has_key('name'):
            Category.__init__(self, kwds['name'])
        else:
            Category.__init__(self)
        self._super_categories = list(super_categories)

    def super_categories(self):
        """
        Returns the immediate super categories, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: from sage.categories.category import JoinCategory
            sage: JoinCategory((Semigroups(), FiniteEnumeratedSets())).super_categories()
            [Category of semigroups, Category of finite enumerated sets]
        """
        return self._super_categories

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: J = Category.join((Groups(), CommutativeAdditiveMonoids())); J #indirect doctest
            Join of Category of groups and Category of commutative additive monoids
        """
        return "Join of %s"%(" and ".join(str(cat) for cat in self.super_categories()))

#############################################################
# Functorial Constructions
#############################################################

class CovariantFunctorialConstruction(SageObject):
    """
    An abstract class for construction functors F (eg F = cartesian product,
    tensor product, cartesian product, ...) such that:

     - Each category Cat (eg Cat=Groups()) can provide a category
       FOf(Cat) for parents constructed via this functor (eg
       FOf(Cat) = CartesianProductsOf(Groups())).

     - For parents A, B, C respectively in the categories CatA, CatB,
       CatC, the category of F(A,B,C) is defined by taking the join of
       FOf(CatD) for every common super category of CatA, CatB, CatC.

    Note: CartesianProductsOf(Groups) needs not to specify that it is
    a subcategory of CartesianProductsOf(Monoids). This part of the
    hierarchy is taken care automatically.

    FIXME: Rework entirely the internal mechanism. In particular,
    CartesianProductsOf(Monoids) should be added automatically to the
    super_categories of CartesianProductsOf(Groups); currently the
    resulting category is built as a join of both.

    TODO: What syntax to use to get FOf(Cat)? For example, for the
    tensor product construction, which one of the followings do we
    want (see chat on IRC, on 07/12/2009):

     - tensor(Cat)
     - tensor((Cat, Cat))
     - tensor.of((Cat, Cat))
     - tensor.category_from_categories((Cat, Cat, Cat))
     - Cat.tensor_category()
     - tensor_category(Cat)

    The syntax Cat.tensor_category() does not supports well situations
    like tensor.of([Algebras(), HopfAlgebras(), ...]). Also it forces
    every category to be (somehow) aware of all the tensorial
    construction that could apply to it, even those which are only
    induced from super categories.
    """

    #Each subclass of this class should provide the following data:

    # functor_name (string, required).  This should match a method on
    # parents that applies the construction

    # functor_category (string, required).  This should match a method
    # on categories that gives the category that the result of this
    # construction applied to objects in that category (ie FOf(Cat) in
    # the docstring)

    # FunctorialCategory (class, required).  This is a subclass of
    # Category from which categories supporting this operation should
    # inherit.

    def category_from_parents(self, parents):
        """
        Returns the category of `F(parents)`.

        INPUT:
         - self: a functor F
         - parents: a list (or iterable) of parents.

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"])
            sage: tensor.category_from_parents((E, E, E)) # todo: not implemented (see upcoming category patch #5985)

        """
        from sage.structure.parent import Parent
        assert(all(isinstance(parent, Parent) for parent in parents))
        return self.category_from_categories(tuple(set(parent.category() for parent in parents)))

    @cached_method
    def category_from_categories(self, categories):
        """
        Returns the category of `F(A,B,C)` for `A,B,C` parents in the given categories

        INPUT:
         - self: a functor `F`
         - categories: a tuple of categories

        EXAMPLES::

            sage: Cat = ModulesWithBasis(QQ)
            sage: tensor.category_from_categories((Cat, Cat, Cat))
            Category of tensor products of modules with basis over Rational Field
        """
        assert(len(categories) > 0)

        super_categories_of_first  = categories[0].all_super_categories()
        super_categories_of_others = [set(category.all_super_categories())
                                      for category in categories[1:len(categories)] ]

        return Category.join((getattr(category, self.functor_category)()
                              for category in super_categories_of_first
                              if isinstance(category, self.FunctorialCategory) and
                              all(category in categories for categories in super_categories_of_others)))

    def __call__(self, args):
        """
        Functorial construction application

        INPUT:
         - self: a covariant functorial construction `F`
         - args: a tuple (or iterable) of parents or elements

        Returns `F(args)`

        EXAMPLES::

            sage: E = CombinatorialFreeModule(QQ, ["a", "b", "c"])
            sage: tensor((E, E, E))           # todo: not implemented (see upcoming category patch #5985)
        """
        args = tuple(args) # a bit brute force; let's see if this becomes a bottleneck later
        assert(all( hasattr(arg, self.functor_name) for arg in args))
        assert(len(args) > 0)
        return getattr(args[0].__class__, self.functor_name)(*args)
