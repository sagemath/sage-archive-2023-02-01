"""
Categories with axiom

EXAMPLES::

    sage: Magmas()
    Category of magmas
    sage: Magmas().Unital()
    Category of unital magmas
    sage: Magmas().Commutative().Unital()
    Category of commutative unital magmas
    sage: Magmas().Associative()
    Category of semigroups
    sage: Magmas().Associative() & Magmas().Unital().Inverse() & Sets().Finite()
    Category of finite groups

    sage: _ is Groups().Finite()
    True

    sage: AdditiveMagmas().AdditiveAssociative().AdditiveCommutative()
    Category of commutative additive semigroups

    sage: from sage.categories.distributive_magmas_and_additive_magmas import DistributiveMagmasAndAdditiveMagmas
    sage: C = CommutativeAdditiveMonoids() & Monoids() & DistributiveMagmasAndAdditiveMagmas(); C
    Category of semirings
    sage: C.AdditiveInverse()
    Category of rings
    sage: Rings().axioms()
    frozenset([...])
    sage: sorted(Rings().axioms())
    ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative', 'Unital']


TESTS::

    sage: from sage.categories.objects import Objects
    sage: Objects()
    Category of objects
    sage: from sage.categories.magmas import Magmas
    sage: Magmas()
    Category of magmas
    sage: Magmas().Finite()
    Category of finite magmas

    sage: from sage.categories.semigroups import Semigroups
    sage: Semigroups()
    Category of semigroups
    sage: Semigroups().Finite()
    Category of finite semigroups

    sage: from sage.categories.modules_with_basis import ModulesWithBasis
    sage: ModulesWithBasis(QQ) is Modules(QQ).WithBasis()
    True
    sage: ModulesWithBasis(ZZ) is Modules(ZZ).WithBasis()
    True

    sage: Semigroups().Unital()
    Category of monoids
    sage: Semigroups().Unital().Commutative() # CHECK
    Category of commutative monoids
    sage: Semigroups().Commutative()
    Category of commutative semigroups
    sage: Semigroups().Commutative().Unital() # oops!
    Category of commutative monoids
    sage: Semigroups().Commutative().Unital().super_categories()
    [Category of monoids, Category of commutative magmas]

    sage: Domains().Commutative()
    Category of integral domains

Wedderburn's theorem::

    sage: DivisionRings().Finite()
    Category of finite fields

    sage: FiniteMonoids().Algebras(QQ)
    Join of Category of finite dimensional algebras with basis over Rational Field and Category of monoid algebras over Rational Field and Category of finite set algebras over Rational Field
    sage: FiniteGroups().Algebras(QQ)
    Join of Category of finite dimensional hopf algebras with basis over Rational Field and Category of group algebras over Rational Field and Category of finite set algebras over Rational Field


.. TODO:

    - primer

    - Implement compatibility axiom / functorial constructions

      E.g. join(A.CartesianProducts(), B.CartesianProducts()) = join(A,B).CartesianProducts()

    - Should an axiom category of a singleton category be
      systematically a singleton category? Same thing for category
      with base ring?

    - Once full subcategories are implemented (see :trac:`10668`),
      make category with axioms be such. Should all full subcategories
      be implemented in term of axioms?
"""
#*****************************************************************************
#  Copyright (C) 2011-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import importlib
import re
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.misc import call_method
from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_types import Category_over_base_ring
from sage.structure.dynamic_class import DynamicMetaclass
# modifier / modified
# specifier / specified
#             variant
# qualifier / qualified category
# axiom / constraint

all_axioms = ("Flying", "Blue",
              "Facade", "Finite", "Infinite",
              "FiniteDimensional", "Connected", "WithBasis",
              "Irreducible",
              "Commutative", "Associative", "Inverse", "Unital", "Division", "NoZeroDivisors",
              "AdditiveCommutative", "AdditiveAssociative", "AdditiveInverse", "AdditiveUnital",
              )

# The order of the axioms implies that
# Magmas().Commutative().Unital() is printed as
# ``Category of commutative unital magmas''

@cached_function
def axioms_rank(axiom):
    """
    Internal function to get the index of an axiom.

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import axioms_rank
        sage: axioms_rank("Finite")
        3
        sage: axioms_rank("FiniteDimensional")
        5

    This is mostly used by :meth:`canonicalize_axioms`

    """
    return all_axioms.index(axiom)

def canonicalize_axioms(axioms):
    r"""
    Canonicalize a set of axioms

    INPUT:

     - ``axioms`` -- a set (or iterable) of axioms

    OUTPUT: the same, as a tuple sorted according to the order of
        sage.categories.category_with_axiom.all_axioms

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import canonicalize_axioms
        sage: canonicalize_axioms(["Commutative", "Connected", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
        sage: canonicalize_axioms(["Commutative", "Connected", "Commutative", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
    """
    return tuple(sorted(set(axioms), key = axioms_rank))

def uncamelcase(s,separator=" "):
    """
    EXAMPLES::

        sage: sage.categories.category_with_axiom.uncamelcase("FiniteDimensionalAlgebras")
        'finite dimensional algebras'
        sage: sage.categories.category_with_axiom.uncamelcase("FiniteDimensionalAlgebras", "_")
        'finite_dimensional_algebras'
    """
    return re.sub("[a-z][A-Z]", lambda match: match.group()[0]+separator+match.group()[1], s).lower()

def base_category_class_and_axiom(cls):
    """
    Try to guess the axioms from the class name and, if the
    class is not a nested class, the plain category.

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import base_category_class_and_axiom
        sage: base_category_class_and_axiom(FiniteSets)
        (<class 'sage.categories.sets_cat.Sets'>, 'Finite')
        sage: Sets.Finite
        <class 'sage.categories.finite_sets.FiniteSets'>
        sage: base_category_class_and_axiom(Sets.Finite)
        (<class 'sage.categories.sets_cat.Sets'>, 'Finite')

        sage: base_category_class_and_axiom(FiniteDimensionalHopfAlgebrasWithBasis)
        (<class 'sage.categories.hopf_algebras_with_basis.HopfAlgebrasWithBasis'>, 'FiniteDimensional')

        sage: Sets.Infinite
        <class 'sage.categories.sets_cat.Sets.Infinite'>
        sage: base_category_class_and_axiom(Sets.Infinite)
        Traceback (most recent call last):
        ...
        TypeError: Could not retrieve base category class for <class 'sage.categories.sets_cat.Sets.Infinite'>

    .. NOTE::

        In the later case, we could possibly retrieve ``Sets`` from
        the class name; however we ca notdo it robustly until #9107 is
        fixed, we are stuck, and anyway we haven't needed it so far.
    """
    if "." in cls.__name__:
        # Case 1: class name of the form Sets.Infinite
        #axiom = cls.__name__.split(".")[-1]
        raise TypeError("Could not retrieve base category class for %s"%cls)

    # Case 2: class name of the form FiniteSets or AlgebrasWithBasis
    # The following checks that the class is defined in the module
    # with the corresponding name (e.g. sage.categories.finite_sets)
    # Otherwise the category should specify explicitly _base_category_class_and_axiom
    # TODO: document this!
    name = cls.__name__
    module = uncamelcase(name, "_")
    assert cls.__module__ == "sage.categories."+module,\
        "%s should be implemented in `sage.categories.%s`"%(cls, module)
    for axiom in all_axioms:
        if axiom == "WithBasis" and name.endswith(axiom):
            base_name = name[:-len(axiom)]
            base_module_name = module[:-len(uncamelcase(axiom))-1]
        elif name.startswith(axiom):
            base_name = name[len(axiom):]
            base_module_name = module[len(uncamelcase(axiom))+1:]
        else:
            continue
        if base_module_name == "sets": # Special case for Sets which is in sets_cat
            base_module_name = "sets_cat"
        try:
            base_module = importlib.import_module("sage.categories."+base_module_name)
            base_category_class = getattr(base_module, base_name)
            assert getattr(base_category_class, axiom, None) is cls, \
                "Missing (lazy import) link for %s to %s for axiom %s?"%(base_category_class, cls, axiom)
            return base_category_class, axiom
        except (ImportError,AttributeError):
            pass
    raise TypeError("Could not retrieve base category class for %s"%cls)

@cached_function
def axiom_of_nested_class(cls, nested_cls):
    """
    Given a class and a nested axiom class, return the axiom.

    EXAMPLES:

    This uses some heuristic like checking if the nested_cls carries
    the name of the axiom, or is built by appending or prepending the
    name of the axiom to that of the class.


        sage: from sage.categories.category_with_axiom import TestObjects, axiom_of_nested_class
        sage: axiom_of_nested_class(TestObjects, TestObjects.FiniteDimensional)
        'FiniteDimensional'
        sage: axiom_of_nested_class(TestObjects.FiniteDimensional, TestObjects.FiniteDimensional.Finite)
        'Finite'
        sage: axiom_of_nested_class(Sets, FiniteSets)
        'Finite'
        sage: axiom_of_nested_class(Algebras, AlgebrasWithBasis)
        'WithBasis'

    In all other cases, the nested class should provide an attribute
    _base_category_class_and_axiom::

        sage: Semigroups._base_category_class_and_axiom
        [<class 'sage.categories.magmas.Magmas'>, 'Associative']
        sage: axiom_of_nested_class(Magmas, Semigroups)
        'Associative'
    """
    try:
        axiom =nested_cls.__dict__["_base_category_class_and_axiom"][1]
    except KeyError:
        assert not isinstance(cls, DynamicMetaclass)
        nested_cls_name = nested_cls.__name__.split(".")[-1]
        if nested_cls_name in all_axioms:
            axiom = nested_cls_name
        else:
            cls_name = cls.__name__.split(".")[-1]
            if nested_cls_name.startswith(cls_name):
                axiom = nested_cls_name[len(cls_name):]
            elif nested_cls_name.endswith(cls_name):
                axiom = nested_cls_name[:-len(cls_name)]
            else:
                raise ValueError, "could not infer axiom for the nested class %s of %s"%(nested_cls, cls)
    assert axiom in all_axioms, \
        "Incorrect guess (%s) for the name of the axiom for the nested class %s of %s"%(axiom, nested_cls, cls)
    assert axiom in cls.__dict__ and cls.__dict__[axiom] == nested_cls, \
        "%s not a nested axiom class of %s for axiom %s"%(nested_cls, cls, axiom)
    return axiom

class CategoryWithAxiom(Category):
    r"""
    An abstract class for categories obtained by adding an axiom to a base category.

    .. TODO::

        Explanations on how they work + link to the tutorial above

    Design goals:

     - Freedom to implement the category of finite sets either in a
       separate file (in a class FiniteSets), or within the Sets
       category (in a nested class Sets.Finite).

    .. NOTE::

        The constructor for instances of this class takes as input the
        base category. Hence, they should in principle be constructed
        as::

            sage: FiniteSets(Sets())    # todo: not tested
            Category of finite sets

            sage: Sets.Finite(Sets())   # todo: not tested
            Category of finite sets

        None of those syntaxes are really practical for the user. So instead,
        this object is to be constructed using any of the following idioms::

            sage: Sets()._with_axiom('Finite')
            Category of finite sets
            sage: FiniteSets()
            Category of finite sets
            sage: Sets().Finite()
            Category of finite sets

        The later two are implemented using respectively
        :meth:`__classcall__` and :meth:`__classget__` which see.
    """

    @lazy_class_attribute
    def _base_category_class_and_axiom(cls):
        r"""
        The class of the base category and the axiom for this class.

        By default, this attribute is guessed from the name of this
        category (see :func:`base_category_class_and_axiom`). It
        should be overriden whenever it cannot be guessed properly.

        When this attribute has been guessed, the attribute
        ``_base_category_class_and_axiom_was_guessed`` is set to
        ``True``.

        .. SEEALSO:: :meth:`_axiom`

        EXAMPLES::

            sage: FiniteSets()._base_category_class_and_axiom
            (<class 'sage.categories.sets_cat.Sets'>, 'Finite')
            sage: FiniteSets()._base_category_class_and_axiom_was_guessed
            True

            sage: Fields()._base_category_class_and_axiom
            (<class 'sage.categories.division_rings.DivisionRings'>, 'Commutative')
            sage: Fields()._base_category_class_and_axiom_was_guessed
            Traceback (most recent call last):
            ...
            AttributeError: 'Fields_with_category' object has no attribute '_base_category_class_and_axiom_was_guessed'

        .. NOTE::

            The base category class is often another category with
            axiom, therefore having a special ``__classget__`` method.
            Storing the base category class and the axiom in a single
            tuple attribute -- instead of two separate attributes --
            has the advantage of not trigerring, for example,
            ``Semigroups.__classget__`` upon
            ``Monoids._base_category_class``.
        """
        base_category_class, axiom = base_category_class_and_axiom(cls)
        cls._base_category_class_and_axiom_was_guessed = True
        return (base_category_class, axiom)

    @lazy_class_attribute
    def _axiom(cls):
        r"""
        The axiom for this category with axiom.

        .. SEEALSO:: :meth:`_base_category_class_and_axiom`

        EXAMPLES::

            sage: FiniteSets._axiom
            'Finite'
            sage: Sets.Finite._axiom
            'Finite'
            sage: Algebras.Commutative._axiom
            'Commutative'

        The result can be less obvious::

            sage: Semigroups._axiom
            'Associative'
            sage: Rings._axiom
            'Unital'
            sage: Fields._axiom
            'Commutative'
        """
        return cls._base_category_class_and_axiom[1]

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Make ``FooBars(**)`` an alias for ``Foos(**)._with_axiom("Bar")``.

        EXAMPLES::


            sage: FiniteGroups()
            Category of finite groups
            sage: ModulesWithBasis(ZZ)
            Category of modules with basis over Integer Ring
            sage: AlgebrasWithBasis(QQ)
            Category of algebras with basis over Rational Field

        This is relevant when e.g. ``Foos(**)`` does some non trivial
        transformations::

            sage: Modules(QQ) is VectorSpaces(QQ)
            True
            sage: type(Modules(QQ))
            <class 'sage.categories.vector_spaces.VectorSpaces_with_category'>

            sage: ModulesWithBasis(QQ) is VectorSpaces(QQ).WithBasis()
            True
            sage: type(ModulesWithBasis(QQ))
            <class 'sage.categories.vector_spaces.VectorSpaces.WithBasis_with_category'>
        """
        (base_category_class, axiom) = cls._base_category_class_and_axiom
        if len(args) == 1 and not options and isinstance(args[0], base_category_class):
            return super(CategoryWithAxiom, cls).__classcall__(cls, args[0])
        else:
            # The following fails with Modules(QQ), as the later returns
            # VectorSpaces(QQ) which is not an instance of the
            # base_category_class of ModulesWithBasis
            # return cls(base_category_class(*args, **options))
            return base_category_class(*args, **options)._with_axiom(axiom)

    @staticmethod
    def __classget__(cls, base_category, base_category_class):
        """
        A bit of black magic to support the following syntax.

        EXAMPLES::

            sage: Sets().Finite()
            Category of finite sets

        When a category does not provide code or properties for its
        finite objects, we get a join category::

            sage: Rings().Finite()
            Category of finite rings

        Well, we have to open the hood to actually see it::

            sage: Rings().Finite()._repr_(as_join=True)
            'Join of Category of rings and Category of finite monoids'

        .. NOTE:: the above example is subject to change

        Thanks to this magic, the documentation obtained by
        introspection on ``Sets().Infinite`` is that of
        :func:`sage.categories.category_with_axiom.Infinite`.

            sage: Sets().Infinite
            Cached version of <function Infinite at ...>

        TESTS::

            sage: Sets().Infinite.f == Sets.SubcategoryMethods.Infinite.f
            True

        We check that this also works when the class is implemented in
        a separate file, and lazy imported::

            sage: Sets().Finite
            Cached version of <function Finite at ...>

        There is no black magic when accessing ``Finite`` or
        ``Infinite`` from the class of the category instead of the
        category::

            sage: Sets.Finite
            <class 'sage.categories.finite_sets.FiniteSets'>
            sage: Sets.Infinite
            <class 'sage.categories.sets_cat.Sets.Infinite'>
        """
        # TODO: this is super paranoid; see if this can be simplified a bit
        if base_category is not None:
            assert base_category.__class__ is base_category_class
            assert isinstance(base_category_class, DynamicMetaclass)
        if isinstance(base_category_class, DynamicMetaclass):
            base_category_class = base_category_class.__base__
        if "_base_category_class_and_axiom" not in cls.__dict__:
            cls._base_category_class_and_axiom = (base_category_class, axiom_of_nested_class(base_category_class, cls))
            cls._base_category_class_and_axiom_was_guessed = False
        else:
            assert cls._base_category_class_and_axiom[0] is base_category_class, \
                "base category class for %s mismatch; expected %s, got %s"%(cls, cls._base_category_class_and_axiom[0], base_category_class)
        if base_category is None:
             return cls
        # For Rings().Finite(), this returns the method
        # Category.Finite, with its first argument bound to Rings()
        return getattr(super(base_category.__class__.__base__, base_category), cls._axiom)

    def __init__(self, base_category):
        """
        TESTS::

            sage: C = Sets.Finite(); C
            Category of finite sets
            sage: type(C)
            <class 'sage.categories.finite_sets.FiniteSets_with_category'>
            sage: type(C).__base__.__base__
            <class 'sage.categories.category_with_axiom.CategoryWithAxiom_singleton'>

            sage: TestSuite(C).run()
        """
        # A hack to upgrade axiom categories of singleton categories
        # to be singleton categories themselves
        if isinstance(base_category, Category_singleton) and not isinstance(self, CategoryWithAxiom_singleton):
            cls = self.__class__
            assert cls.__base__ == CategoryWithAxiom
            cls.__bases__ = (CategoryWithAxiom_singleton,)+cls.__bases__[1:]

        self._base_category = base_category
        Category.__init__(self)

    def _test_category_with_axiom(self, **options):
        r"""
        Run generic tests on this category with axioms

        .. SEEALSO:: :class:`TestSuite`.

        This check that an axiom category of a
        :class:`Category_singleton` is a singleton category, and
        similarwise for :class`Category_over_base_ring`.

        EXAMPLES::

            sage: Sets().Finite()._test_category_with_axiom()
            sage: Modules(ZZ).FiniteDimensional()._test_category_with_axiom()
        """
        tester = self._tester(**options)
        base = self.base_category()
        if isinstance(base, Category_singleton):
            tester.assertIsInstance(self, CategoryWithAxiom_singleton)
        if isinstance(base, Category_over_base_ring):
            tester.assertIsInstance(self, CategoryWithAxiom_over_base_ring)

    def extra_super_categories(self):
        """
        Returns the extra super categories of a category with axiom

        Default implementation which returns ``[]``

        EXAMPLES::

            sage: FiniteSets().extra_super_categories()
            []
        """
        return []

    @cached_method
    def super_categories(self):
        """
        Returns a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        This implements the property that if ``As`` is a subcategory
        of ``Bs``, then the intersection of As with ``FiniteSets()``
        is a subcategory of ``As`` and of the intersection of ``Bs``
        with ``FiniteSets()``.

        EXAMPLES::

            sage: FiniteSets().super_categories()
            [Category of sets]

            sage: FiniteSemigroups().super_categories()
            [Category of semigroups, Category of finite enumerated sets]

        EXAMPLES:

        A finite magma is both a magma and a finite set::

            sage: Magmas().Finite().super_categories()
            [Category of magmas, Category of finite sets]

        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects
            sage: C = TestObjects().FiniteDimensional().Unital().Commutative().Finite()
            sage: sorted(C.super_categories(), key=str)
            [Category of finite commutative test objects,
             Category of finite dimensional commutative unital test objects,
             Category of finite finite dimensional test objects]
        """
        base_category = self._base_category
        axiom = self._axiom
        extra = self.extra_super_categories()
        return Category.join((self._base_category,) +
                             tuple(base_category.super_categories()) +
                             tuple(extra),
                             axioms = (axiom,),
                             uniq=False,
                             ignore_axioms = ((base_category, axiom),),
                             as_list = True)

    @staticmethod
    def _repr_object_names_static(category, axioms):
        r"""
        INPUT:

        - ``base_category`` -- a category
        - ``axioms`` -- a list or iterable of strings

        EXAMPLES::

            sage: from sage.categories.category_with_axiom import CategoryWithAxiom
            sage: CategoryWithAxiom._repr_object_names_static(Semigroups(), ["Flying", "Blue"])
            'flying blue semigroups'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ), ["Flying", "WithBasis", "Blue"])
            'flying blue algebras with basis over Rational Field'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ), ["WithBasis"])
            'algebras with basis over Rational Field'
            sage: CategoryWithAxiom._repr_object_names_static(Sets().Finite().Subquotients(), ["Finite"])
            'subquotients of finite sets'
            sage: CategoryWithAxiom._repr_object_names_static(Monoids(), ["Unital"])
            'monoids'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ['x']['y']), ["Flying", "WithBasis", "Blue"])
            'flying blue algebras with basis over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field'

        If the axioms is a set or frozen set, then they are first
        sorted using :func:`canonicalize_axioms`::

            sage: CategoryWithAxiom._repr_object_names_static(Semigroups(), set(["Finite", "Commutative", "Facade"]))
            'facade finite commutative semigroups'

        .. SEEALSO:: :meth:`_repr_object_names`

        .. NOTE:: The logic here is shared between :meth:`_repr_object_names`
            and :meth:`.category.JoinCategory._repr_object_names`
        """
        axioms = canonicalize_axioms(axioms)
        base_category = category._without_axioms(named=True)
        if isinstance(base_category, CategoryWithAxiom): # Smelly runtime type checking
            result = super(CategoryWithAxiom, base_category)._repr_object_names()
        else:
            result = base_category._repr_object_names()
        for axiom in reversed(axioms):
            # TODO: find a more generic way to handle the special cases below
            if axiom in base_category.axioms():
                # If the base category already has this axiom, we
                # need not repeat it here. See the example with
                # Sets().Finite().Subquotients() or Monoids()
                continue
            if axiom == "WithBasis":
                result = result.replace(" over ", " with basis over ", 1)
            elif axiom == "Connected" and "graded " in result:
                result = result.replace("graded ", "graded connected ", 1)
            else:
                result = uncamelcase(axiom) + " " + result
        return result

    def _repr_object_names(self):
        r"""
        The names of the objects of this category, as used by `_repr_`

        .. SEEALSO:: :meth:`Category._repr_object_names`

        EXAMPLES::

            sage: FiniteSets()._repr_object_names()
            'finite sets'
            sage: AlgebrasWithBasis(QQ).FiniteDimensional()._repr_object_names()
            'finite dimensional algebras with basis over Rational Field'
            sage: Monoids()._repr_object_names()
            'monoids'
            sage: Semigroups().Unital().Finite()._repr_object_names()
            'finite monoids'
            sage: Algebras(QQ).Commutative()._repr_object_names()
            'commutative algebras over Rational Field'

        .. NOTE::

            This is implemented by taking _repr_object_names from
            self._without_axioms(named=True), and adding the names
            of the relevant axioms in appropriate order.
        """
        return CategoryWithAxiom._repr_object_names_static(self, self.axioms())

    def base_category(self):
        r"""
        Returns the plain category of ``self``

        EXAMPLES::

            sage: C = Sets.Finite(); C
            Category of finite sets
            sage: C.base_category()
            Category of sets
            sage: C._without_axioms()
            Category of sets

        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects, CategoryWithAxiom
            sage: C = TestObjects().Commutative().Facade()
            sage: assert isinstance(C, CategoryWithAxiom)
            sage: C._without_axioms()
            Category of test objects
        """
        return self._base_category

    def __reduce__(self):
        r"""
        Implement the pickle protocol.

        This overides the implementation in
        :meth:`UniqueRepresentation.__reduce__` in order to not
        exposes the implementation detail that, for example, the
        category of magmas which distribute over an associative
        additive magma is implemented as
        ``DistributiveMagmasAndAdditiveMagmas.AdditiveAssociative.AdditiveCommutative``
        and not
        ``DistributiveMagmasAndAdditiveMagmas.AdditiveCommutative.AdditiveAssociative``::

        EXAMPLES::

            sage: C = Semigroups()
            sage: reduction = C.__reduce__(); reduction
            (<function call_method at ...>, (Category of magmas, '_with_axiom', 'Associative'))
            sage: loads(dumps(C)) is C
            True
            sage: FiniteSets().__reduce__()
            (<function call_method at ...>, (Category of sets, '_with_axiom', 'Finite'))

            sage: C = DistributiveMagmasAndAdditiveMagmas().AdditiveAssociative().AdditiveCommutative()
            sage: C.__class__
            <class 'sage.categories.distributive_magmas_and_additive_magmas.AdditiveAssociative.AdditiveCommutative_with_category'>
            sage: C.__reduce__()
            (<function call_method at ...>, (Category of additive associative distributive magmas and additive magmas, '_with_axiom', 'AdditiveCommutative'))
        """
        return (call_method, (self._base_category, "_with_axiom", self._axiom))

    @cached_method
    def _without_axiom(self, axiom):
        r"""
        Return this category with axiom ``axiom`` removed.

        OUTPUT: a category ``C`` which does not have axiom ``axiom``
        and such that either ``C`` is ``self``, or adding back all the
        axioms of ``self`` gives back ``self``.

        .. SEEALSO:: :meth:`Category._without_axiom`

        .. WARNING:: This is not guaranteed to be robust.

        EXAMPLES::

            sage: Groups()._without_axiom("Unital")
            Category of semigroups
            sage: Groups()._without_axiom("Associative")
            Category of inverse unital magmas
            sage: Groups().Commutative()._without_axiom("Unital")
            Category of commutative semigroups
        """
        axioms = self.axioms().difference([axiom])
        return self._without_axioms()._with_axioms(axioms)

    def _without_axioms(self, named=False):
        """
        Returns the category without the axioms that have been added to create it

        EXAMPLES::

            sage: Sets().Finite()._without_axioms()
            Category of sets
            sage: Monoids().Finite()._without_axioms()
            Category of magmas

        This is because::

            sage: Semigroups().Unital() is Monoids()
            True

        If ``named`` is True, then `_without_axioms` stops at the
        first category that has a explicit name of its own::

            sage: Sets().Finite()._without_axioms(named=True)
            Category of sets
            sage: Monoids().Finite()._without_axioms(named=True)
            Category of monoids

        Technically we tests this by checking if the class has a
        predefined attribute ``_base_category_class_and_axiom`` (if
        instead this attribute is computed, then this computation
        should set the attribute
        ``_base_category_class_and_axiom_was_guessed``.

        Some more examples::

            sage: Algebras(QQ).Commutative()._without_axioms()
            Category of magmatic algebras over Rational Field
            sage: Algebras(QQ).Commutative()._without_axioms(named=True)
            Category of algebras over Rational Field
        """
        if named:
            base_category = self
            axioms = []
            while isinstance(base_category, CategoryWithAxiom) and hasattr(base_category.__class__, "_base_category_class_and_axiom_was_guessed"): # and "." in Sets().Infinite().__class__.__name__:
                axioms.append(base_category._axiom)
                base_category = base_category._base_category
            return base_category
        else:
            return self._base_category._without_axioms()

    @cached_method
    def axioms(self):
        r"""
        Returns the axioms of ``self``

        EXAMPLES::

            sage: C = Sets.Finite(); C
            Category of finite sets
            sage: C.axioms()
            frozenset(['Finite'])

            sage: C = Modules(GF(5)).FiniteDimensional(); C
            Category of finite finite dimensional vector spaces over Finite Field of size 5
            sage: sorted(C.axioms())
            ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Finite', 'FiniteDimensional']

            sage: sorted(FiniteMonoids().Algebras(QQ).axioms())
            ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative', 'FiniteDimensional', 'Unital', 'WithBasis']
            sage: sorted(FiniteMonoids().Algebras(GF(3)).axioms())
            ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative', 'Finite', 'FiniteDimensional', 'Unital', 'WithBasis']

            sage: DistributiveMagmasAndAdditiveMagmas().Unital().axioms()
            frozenset(['Unital'])

            sage: D = DistributiveMagmasAndAdditiveMagmas()
            sage: X = D.AdditiveAssociative().AdditiveCommutative().Associative()
            sage: X.Unital().super_categories()[1]
            Category of monoids
            sage: X.Unital().super_categories()[1] is Monoids()
            True
        """
        # We would want to write the following line:
        #     return super(CategoryWithAxiom, self).axioms() | {self._axiom}
        # However one currently can't use super to call a cached
        # method in a super class. So we dup the code from there ...
        return frozenset(axiom
                         for category in self.super_categories()
                         for axiom in category.axioms()) | {self._axiom}

class CategoryWithAxiom_over_base_ring(CategoryWithAxiom, Category_over_base_ring):

    def __init__(self, base_category):
        """
        TESTS::

            sage: C = Modules(ZZ).FiniteDimensional(); C
            Category of finite dimensional modules over Integer Ring
            sage: type(C)
            <class 'sage.categories.modules.Modules.FiniteDimensional_with_category'>
            sage: type(C).__base__.__base__
            <class 'sage.categories.category_with_axiom.CategoryWithAxiom_over_base_ring'>

            sage: TestSuite(C).run()
        """
        # FIXME: this is basically a duplicates the code from
        # CategoryWithAxiom.__init__; but we can't call the later
        # without calling twice Category.__init__; or maybe playing
        # around with super?
        self._base_category = base_category
        Category_over_base_ring.__init__(self, base_category.base_ring())

class CategoryWithAxiom_singleton(Category_singleton, CategoryWithAxiom):#, Category_singleton, FastHashable_class):
    pass

"""
The following work around is needed until any CategoryWithAxiom of a
Category_over_base_ring becomes automatically a
CategoryWithAxiom_over_base_ring::

    sage: from sage.categories.category_with_axiom import TestObjectsOverBaseRing, Category_over_base_ring
    sage: from sage.categories.category import JoinCategory
    sage: isinstance(TestObjectsOverBaseRing(QQ), Category_over_base_ring)
    True
    sage: C = TestObjectsOverBaseRing(QQ).Commutative()
    sage: isinstance(C, Category_over_base_ring)          # todo: not implemented
    True
    sage: C.FiniteDimensional()
    Category of finite dimensional commutative test objects over base ring over Rational Field
    sage: C.Commutative()
    Category of commutative test objects over base ring over Rational Field
    sage: C.Unital()
    Category of commutative unital test objects over base ring over Rational Field

    sage: C = TestObjectsOverBaseRing(IntegerModRing(2)).Connected()
    sage: isinstance(C, JoinCategory)
    True
    sage: isinstance(C, Category_over_base_ring)          # todo: not implemented
    True
    sage: C.FiniteDimensional()
    Category of finite dimensional connected test objects over base ring over Ring of integers modulo 2
    sage: C.Connected()
    Category of connected test objects over base ring over Ring of integers modulo 2
"""

##############################################################################
# Utilities and tests tools

class SmallTestObjects(Category):
    class Finite(CategoryWithAxiom):
        pass
    # Those should be ignored
    Connected = 1
    class Commutative:
        pass

def axiom(axiom):
    """
    Return a function/method ``self -> self._with_axiom(axiom)``.

    This can used as a shorthand to define axioms, in particular in
    the tests below. Usually one will want to attach documentation to
    an axiom, so the need for such a shorthand in real life might not
    be that clear, unless we start creating lots of axioms.

    In the long run maybe this could evolve into an @axiom decorator.

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import axiom
        sage: axiom("Finite")(Semigroups())
        Category of finite semigroups

    Upon assigning the result to a class this becomes a method::

        sage: class As:
        ....:     def _with_axiom(self, axiom): return self, axiom
        ....:     Finite = axiom("Finite")
        sage: As().Finite()
        (<__main__.As instance at ...>, 'Finite')
    """
    def with_axiom(self):
        return self._with_axiom(axiom)
    with_axiom.__name__ = axiom
    return with_axiom

class Blahs(Category):

    def super_categories(self):
        """
        TESTS::

             sage: from sage.categories.category_with_axiom import Blahs
             sage: Blahs().super_categories()
             [Category of sets]
        """
        from sage.categories.sets_cat import Sets
        return [Sets()]

    class SubcategoryMethods:
        FiniteDimensional = axiom("FiniteDimensional")
        Commutative       = axiom("Commutative")
        Unital            = axiom("Unital")
        Connected         = axiom("Connected")
        Flying            = axiom("Flying")
        Blue              = axiom("Blue")

    class FiniteDimensional(CategoryWithAxiom):
        pass
    class Commutative(CategoryWithAxiom):
        pass
    class Connected(CategoryWithAxiom):
        pass
    class Unital(CategoryWithAxiom):
        class Blue(CategoryWithAxiom):
            pass
    class Flying(CategoryWithAxiom):
        def extra_super_categories(self):
            """
            This illustrates a way to have an axiom imply another one.

            TESTS::

                sage: from sage.categories.category_with_axiom import Blahs, TestObjects, Bars
                sage: Blahs().Flying().extra_super_categories()
                [Category of unital blahs]
                sage: Blahs().Flying()
                Category of flying unital blahs
            """
            return [Blahs().Unital()]
    def Blue_extra_super_categories(self):
        """
        Tries to illustrate another way to have an axiom imply another one.

        This currently fails because there is no base axiom category
        Blahs.Blue, and thus somewhere during the join calculation the
        axiom is lost.

        TESTS::

            sage: from sage.categories.category_with_axiom import Blahs, TestObjects, Bars
            sage: Blahs().Blue_extra_super_categories()
            [Category of unital blahs]
            sage: Blahs().Blue()                          # todo: not implemented
            Category of blue unital blahs
        """
        return [Blahs().Unital()]

class Bars(Category):
    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.category_with_axiom import Bars
            sage: Bars().super_categories()
            [Category of blahs]
        """
        return [Blahs()]

    def Unital_extra_super_categories(self):
        """
        Return extraneous super categories for the unital objects of ``self``.

        This method specifies that a unital bar is a test object.
        Thus, the categories of unital bars and of unital test objects
        coincide.

        EXAMPLES::

            sage: from sage.categories.category_with_axiom import Bars, TestObjects
            sage: Bars().Unital_extra_super_categories()
            [Category of test objects]
            sage: Bars().Unital()
            Category of unital test objects
            sage: TestObjects().Unital().all_super_categories()
            [Category of unital test objects,
             Category of unital blahs,
             Category of test objects,
             Category of bars,
             Category of blahs,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """
        return [TestObjects()]

class TestObjects(Category):

    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects
            sage: TestObjects().super_categories()
            [Category of bars]
        """
        return [Bars()]

    class FiniteDimensional(CategoryWithAxiom):
         class Finite(CategoryWithAxiom):
              pass
         class Unital(CategoryWithAxiom):
              class Commutative(CategoryWithAxiom):
                   pass

    class Commutative(CategoryWithAxiom):
         class Facade(CategoryWithAxiom):
             pass
         class FiniteDimensional(CategoryWithAxiom):
             pass
         class Finite(CategoryWithAxiom):
             pass

    class Unital(CategoryWithAxiom):
        pass

class TestObjectsOverBaseRing(Category_over_base_ring):
    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjectsOverBaseRing
            sage: TestObjectsOverBaseRing(QQ).super_categories()
            [Category of test objects]
        """
        return [TestObjects()]

    class FiniteDimensional(CategoryWithAxiom):
         class Finite(CategoryWithAxiom):
              pass
         class Unital(CategoryWithAxiom):
              class Commutative(CategoryWithAxiom):
                   pass

    class Commutative(CategoryWithAxiom):
         class Facade(CategoryWithAxiom):
             pass
         class FiniteDimensional(CategoryWithAxiom):
             pass
         class Finite(CategoryWithAxiom):
             pass

    class Unital(CategoryWithAxiom):
        pass

class BrokenTestObjects(Category):
    class Commutative(CategoryWithAxiom):
         class Finite(CategoryWithAxiom):
              pass
    class Finite(CategoryWithAxiom):
         class Commutative(CategoryWithAxiom):
              pass

