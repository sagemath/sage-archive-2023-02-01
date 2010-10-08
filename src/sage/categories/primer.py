r"""
Elements, parents, and categories in Sage: a (draft of) primer

Abstract
--------

    The purpose of categories in Sage is to translate the mathematical
    concept of categories (category of groups, of vector spaces, ...)
    into a concrete software engineering design pattern for:

    - organizing and promoting generic code
    - fostering consistency across the Sage library (naming conventions, doc, tests)
    - embedding more mathematical knowledge into the system

    This design pattern is largely inspired from Axiom and its
    followers (Aldor, Fricas, MuPAD, ...). It differs from those by:

    - blending in the Magma inspired concept of Parent/Element

    - being built on top of (and not into) the standard Python object
      oriented and class hierarchy mechanism. This did not require
      changing the language, and could in principle be implemented in
      any language allowing for dynamically creating new classes.


Introduction
------------

To motivate the design, we need to come back to the basics: the
purpose of Sage is to do mathematical computations, and this requires
modeling mathematical objects and operations on the computer. Let us
review some examples of such objects and operations:

- An integer:               `+`, `*`, factorization

- A matrix:                 `+`, `*`, determinant

- A word:                   pattern matching, ...

- The permutations of 5, the rational points of an elliptic curve: counting, listing, random generation

- A language (set of words): rationality testing, counting elements, generating series

- A finite semigroup: left/right ideals, center, representation theory

- A vector space, an algebra: cartesian product, tensor product, quotient

- The category of algebras: what's its initial object? its super categories? its dual category?


Object oriented programming paradigm
------------------------------------

The fundamental principle of object oriented programming (OOP) is that
any object that a program is to manipulate should be modeled by an
*instance* of a *class*. The class implements:

 - a *data structure*: which describes how the object is stored
 - *methods*: which describe the operations on the object

The instance itself contains the data for the given object, according
to the specified data structure.

Hence, all the objects mentionned above should be instances of some
classes. For example, an integer in Sage is an instance of the class
Integer (and it knows about it!)::

    sage: i = 12
    sage: type(i)
    <type 'sage.rings.integer.Integer'>

Applying an operation is generally done by *calling a method*::

    sage: i.factor()
    2^2 * 3

    sage: p = 6*x^2 + 12*x + 6
    sage: type(p)
    <type 'sage.symbolic.expression.Expression'>
    sage: p.factor()
    6*(x + 1)^2

    sage: R.<x> = PolynomialRing(QQ, sparse=True)
    sage: pQ = R ( p )
    sage: type(pQ)
    <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_sparse_field'>
    sage: pQ.factor()
    (6) * (x + 1)^2

    sage: pZ = ZZ[x] ( p )
    sage: type(pZ)
    <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    sage: pZ.factor()
    2 * 3 * (x + 1)^2

Factoring integers, expressions, or polynomials are distinct tasks,
with completely different algorithms. Yet, from a user (or caller)
point of view, all those objects can be manipulated alike. This
illustrates the OOP concepts of *polymorphism*, *data abstraction*,
and *encapsulation*.

Let us be curious, and see where some methods are defined. This can be
done by introspection::

    sage: i._mul_??                   # not tested

For plain Python methods, one can also just ask in which module they
are implemented::

    sage: i._pow_.__module__
    'sage.categories.semigroups'

    sage: pQ._mul_.__module__
    'sage.rings.polynomial.polynomial_element_generic'
    sage: pQ._pow_.__module__
    'sage.categories.semigroups'

We see that integers and polynomials have each their own
multiplication method: the multiplication algorithms are indeed
unrelated and deeply tied to their respective datastructures.

On the other hand, they share the same powering method.
Mathematically speaking the sets `\ZZ` of integers, and the set
`\QQ[x]` of polynomials are both semigroups: i.e. sets endowed with an
associative binary internal lay `*`. Once we know how to compute
products, the powering operation is the same for all semigroups; hence
it can be implemented once for all for all semigroups.

This illustrates the use of *hierarchy of classes* to share common
code between classes having common behaviour. Namely, the class for
integers and the class for polynomials both derive from an *abstract
class* for elements of a semigroups, which factors out the *generic*
methods like ``_pow_``.

OOP design is all about isolating the objects that one wants to model
and their operations, and building up an appropriate hierarchy of
classes for organizing the code. Luckily for mathematicians, a lot of
this work is given for free by abstract algebra.  Namely, the
hierarchy of classes can be modeled upon the hierarchy of abstract
algebraic structures (i.e. *categories*): semigroups, monoids, groups,
rings, algebras, etc. For example, a group is a special case of a
monoid (i.e. the category of groups is a subcategory of the category
of monoids). Therefore, the abstract class for elements of a group
should derive from the abstract classes for elements of monoid.

Before going further, let us recall that we do not just want to compute
with elements of mathematical sets, but with the sets themselves::

    sage: ZZ.one()
    1

    sage: R = QQ['x,y']
    sage: R.krull_dimension()
    2
    sage: A = R.quotient( R.ideal(x^2 - 2) )
    sage: A.krull_dimension() # todo: not implemented

Hence there should be another hierarchy of classes for those sets.

Elements, Parents, Categories
-----------------------------

Philosophy: *Building mathematical information into the system yields
more expressive, more conceptual and, at the end, easier to maintain
and faster code*.

(within a programming realm; this would not necessarily apply to
specialized libraries like gmp!)

Here is some mathematical information that Sage is aware of::

    sage: 1.parent()
    Integer Ring

    sage: ZZ
    Integer Ring

    sage: ZZ.category()
    Category of euclidean domains

    sage: ZZ.categories()
    [Category of euclidean domains,
     Category of principal ideal domains,
     Category of gcd domains,
     Category of integral domains,
     Category of commutative rings,
     Category of domains,
     Category of rings,
     Category of rngs,
     Category of commutative additive groups,
     Category of semirings,
     Category of commutative additive monoids,
     Category of commutative additive semigroups,
     Category of additive magmas,
     Category of monoids,
     Category of semigroups,
     Category of magmas,
     Category of sets,
     Category of sets with partial maps,
     Category of objects]

    sage: EuclideanDomains().category_graph().plot(talk = True)

This illustrates the following relations between mathematical objects:

- elements    <-> parent
- parents     <-> category
- subcategory <-> super category

**Parent**
    A parent is a Python instance representing a set of mathematical
    elements together with its additional (algebraic) strcuture.

    Examples include the ring of integers, the group `S_3`, the set of
    prime numbers, the set of linear maps between a given two vector
    spaces, and a given finite semigroup.

    These sets are often equipped with additional structure: the set
    of all integers forms a ring.  The main way of encoding this
    information is specifying which categories a parent belongs to.

    It is completely possible to have different Python instances
    representing the same set of elements.  For example, one might
    want to consider the ring of integers, or the poset of integers
    under their standard order, or the poset of integers under
    divisibility or the semiring of integers under the operations of
    addition and maximum.  Each of these would be a different
    instance, belonging to different categories.

    For a given model, there should be a unique instance in Sage
    representing that parent::

        sage: IntegerRing() is IntegerRing()
        True

**Element**
    An element is a Python instance representing a mathematical
    element of a set.

    Examples of element include 5 in the integer ring, `x^3 - x` in
    the polynomial ring in `x` over the rationals, `4 + O(3^3)` in the
    3-adics, the transposition `(1 2)` in `S_3`, and the identity
    morphism in the set of linear maps from `\QQ^3` to `\QQ^3`.

    Every element in Sage has a parent.  The standard mechanism in
    Sage for creating elements is to create their parent, and then
    provide enough data to define the element::

        sage: R = PolynomialRing(ZZ, name='x')
        sage: R([1,2,3])
        3*x^2 + 2*x + 1

    One can also create elements using various methods on the parent
    and arithmetic of elements::

        sage: x = R.gen()
        sage: 1 + 2*x + 3*x^2
        3*x^2 + 2*x + 1

    Unlike parents, elements in Sage are not necessarily unique::

        sage: ZZ(5040) is ZZ(5040)
        False

    Many parents represent algebraic structures, and their elements
    support arithmetic operations.  One often further wants to do
    arithmetic by combining elements from different parents: adding
    together integers and rationals for example. Sage supports this
    feature using coercion (see :mod:`sage.structure.coerce` for more
    details).

    It is possible for a parent to also have simultaneously the
    structure of an element. Consider for example the monoid of all
    finite groups, endowed with the cartesian product operation.
    Then, every finite group (which is a parent) is also an element of
    this monoid. This is not yet implemented, and the design details
    are not yet fixed but experiments are underway in this direction.

    TODO: give a concrete example, typically using :class:`ElementWrapper`.

**Category**
    A category is a Python instance representing a mathematical category.

    Examples of categories include the category of finite semigroups,
    the category of `\ZZ`-algebras, the category of bigraded
    `\QQ`-algebras, the category of all (Python) objects, and the
    category of cartesian products of `\ZZ`-algebras::

        sage: FiniteSemigroups()
        Category of finite semigroups
        sage: Algebras(ZZ)
        Category of algebras over Integer Ring
        sage: Objects()
        Category of objects
        sage: Algebras(ZZ).CartesianProducts()
        Category of Cartesian products of algebras over Integer Ring

    (mind the 's' in the names of the categories above;
    ``GroupAlgebra`` and ``GroupAlgebras`` are distinct things).

    Every parent belongs to a collection of categories.  Moreover,
    these categories are related by the relation of being
    super categories.  For example, the category of rings is a
    super category of the category of fields, because every field is
    also a ring.

    A category serves two roles:

    - to provide a model for the mathematical concept of a category
      and the associated structures (homsets, morphisms, functorial
      constructions)

    - to organize and promote generic code, naming conventions,
      documentation, and tests across similar mathematical structures.

**CategoryObjects**
    Objects of a mathematical category are not necessarily parents.
    Parent has a superclass that provides a means of modeling such.

    For example, schemes do not have a faithful forgetful functor to
    the category of Sets, so it does not make sense to talk about
    Schemes as Parents.

(Todo: include a picture!)

Let us look at a semigroup::

    sage: S = FiniteSemigroups().example()
    sage: S?                    # not tested

Now would be a good time to go through the
:mod:`sage.categories.tutorial`!

Where do all the operations on ``S`` and its elements come from?

::

    sage: x = S('a')

``_repr_`` is a technical method which comes with the data structure (ElementWrapper)::

    sage: x._repr_.__module__
    'sage.structure.element_wrapper'

``__pow__`` is a generic method for all finite semigroups::

    sage: x.__pow__.__module__
    'sage.categories.semigroups'

``_mul_`` is a default implementation from the :class:`Magmas`
category (a *magma* is a set with an inner law `*`, not necessarily
associative)::

    sage: x.__mul__.__module__
    'sage.categories.magmas'

It delegates the work to the parent (if you do not know what to do,
ask your parent)::

    sage: x.__mul__??           # not tested

``product`` is a mathematical method implemented by the parent::

    sage: S.product.__module__
    'sage.categories.examples.finite_semigroups'

``cayley_graph`` is a generic method on the parent, provided by the
:class:`FiniteSemigroups` category::

    sage: S.cayley_graph.__module__
    'sage.categories.semigroups'

``multiplication_table`` is a generic method on the parent, provided
by the :class:`Magmas` category (it does not require associativity)::

    sage: S.multiplication_table.__module__
    'sage.categories.magmas'

Consider now the implementation of the semigroup::

    sage: S??                   # not tested

This implementation specifies a data structure for the parents and the
elements, and makes a promise: the implemented parent is a
semigroup. Then it fullfills the promise by implementing the basic
operations ``product`` and ``semigroup_generators``. In exchange of
this promise, `S` and its elements receive generic implementations of
all the other operations. `S` may override any of the operations by
more efficient ones. It could typically implement ``is_idempotent`` to
always return ``True``.

A (not yet complete) list of mandatory and optional methods to be
implemented can be found by introspection with::

    sage: from sage.misc.abstract_method import abstract_methods_of_class
    sage: abstract_methods_of_class(FiniteSemigroups().element_class)
    {'required': [], 'optional': ['_mul_']}
    sage: abstract_methods_of_class(FiniteSemigroups().parent_class)
    {'required': ['__contains__'], 'optional': []}

Documentation about those methods can be obtained with::

    sage: C = FiniteSemigroups().element_class
    sage: C._mul_?                               # not tested

See also the :class:`abstract_methods` decorator.

Here is the code for the finite semigroups category::

    sage: FiniteSemigroups??    # not tested


Wrapup: the mathematical relations between elements, parents, and
categories translates into *inheritance* between classes::

    sage: FiniteSemigroups().all_super_categories()
    [Category of finite semigroups,
     Category of semigroups,
     Category of magmas,
     Category of finite enumerated sets,
     Category of enumerated sets,
     Category of sets,
     Category of sets with partial maps,
     Category of objects]
    sage: S.__class__.mro()
    [<class 'sage.categories.examples.finite_semigroups.LeftRegularBand_with_category'>,
     <class 'sage.categories.examples.finite_semigroups.LeftRegularBand'>,
     <class 'sage.structure.unique_representation.UniqueRepresentation'>,
     <type 'sage.structure.parent.Parent'>,
     <type 'sage.structure.category_object.CategoryObject'>,
     <type 'sage.structure.sage_object.SageObject'>,
     <class 'sage.categories.finite_semigroups.FiniteSemigroups.parent_class'>,
     <class 'sage.categories.semigroups.Semigroups.parent_class'>,
     <class 'sage.categories.magmas.Magmas.parent_class'>,
     <class 'sage.categories.finite_enumerated_sets.FiniteEnumeratedSets.parent_class'>,
     <class 'sage.categories.enumerated_sets.EnumeratedSets.parent_class'>,
     <class 'sage.categories.sets_cat.Sets.parent_class'>,
     <class 'sage.categories.category.SetsWithPartialMaps.parent_class'>,
     <class 'sage.categories.objects.Objects.parent_class'>,
     <type 'object'>]
    sage: x.__class__.mro()
    [<class 'sage.categories.examples.finite_semigroups.LeftRegularBand_with_category.element_class'>,
     <class 'sage.categories.examples.finite_semigroups.LeftRegularBand.Element'>,
     <class 'sage.structure.element_wrapper.ElementWrapper'>,
     <type 'sage.structure.element.Element'>,
     <type 'sage.structure.sage_object.SageObject'>,
     <class 'sage.categories.category.FiniteSemigroups.element_class'>,
     <class 'sage.categories.semigroups.Semigroups.element_class'>,
     <class 'sage.categories.magmas.Magmas.element_class'>,
     <class 'sage.categories.category.FiniteEnumeratedSets.element_class'>,
     <class 'sage.categories.enumerated_sets.EnumeratedSets.element_class'>,
     <class 'sage.categories.sets_cat.Sets.element_class'>,
     <class 'sage.categories.category.SetsWithPartialMaps.element_class'>,
     <class 'sage.categories.objects.Objects.element_class'>,
     <type 'object'>]

Generic tests
-------------

Another feature that a parent receive from its categories are generic
tests; their purpose is to check (at least to some extent) that the
parent has the required mathematical properties (is my semigroup
indeed associative?) and is implemented according to the
specifications (does the method ``an_element`` indeed return an
element of the parent?)::

    sage: S = FiniteSemigroups().example(alphabet=('a', 'b'))
    sage: TestSuite(S).run(verbose = True)
    running ._test_an_element() . . . pass
    running ._test_associativity() . . . pass
    running ._test_category() . . . pass
    running ._test_elements() . . .
      Running the test suite of self.an_element()
      running ._test_category() . . . pass
      running ._test_eq() . . . pass
      running ._test_not_implemented_methods() . . . pass
      running ._test_pickling() . . . pass
      pass
    running ._test_elements_eq() . . . pass
    running ._test_enumerated_set_contains() . . . pass
    running ._test_enumerated_set_iter_cardinality() . . . pass
    running ._test_enumerated_set_iter_list() . . . pass
    running ._test_eq() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_pickling() . . . pass
    running ._test_some_elements() . . . pass
    sage: S._test_associativity?? # not tested
    sage: S._test_associativity(elements=S)

See :class:`TestSuite` for more information.

Let us see what happens when the test fails. Here we redefine the
product of `S` to something definitely not associative::

    sage: %pdb                  # not tested
    sage: S.product = lambda x, y: S("("+x.value +y.value+")")

Combined with the use of the debugger and introspection, those tests
gives instantly a counter example to the associativity of our broken
semigroup::

    sage: S._test_associativity()
    Traceback (most recent call last):
      ...
      File ".../sage/categories/semigroups.py", line ..., in _test_associativity
        tester.assert_((x * y) * z == x * (y * z))
      ...
    AssertionError


Wrapup:

 - Categories bring not only code, but also testing tools
 - This enforces robustness and consistency (despite using an interpreted language)

Advanced algebraic structures
-----------------------------

::

    sage: HopfAlgebrasWithBasis(QQ)?? # not tested
    sage: HopfAlgebrasWithBasis(QQ).category_graph().plot()

Functorial constructions
------------------------

::

    sage: A = AlgebrasWithBasis(QQ).example(); A.rename("A")
    sage: B = HopfAlgebrasWithBasis(QQ).example(); B.rename("B")
    sage: C = cartesian_product([A, B, B]); C
    A (+) B (+) B
    sage: C.one()
    B[(0, word: )] + B[(1, ())] + B[(2, ())]
    sage: cartesian_product([A.one(), B.one(), B.one()])
    B[(0, word: )] + B[(1, ())] + B[(2, ())]
    sage: C.one??                                                # not tested
    sage: C.product??                                            # not tested
    sage: C.categories()
    [Category of Cartesian products of algebras with basis over Rational Field,
     Category of algebras with basis over Rational Field,
     Category of Cartesian products of modules with basis over Rational Field,
     Category of modules with basis over Rational Field,
     Category of Cartesian products of algebras over Rational Field,
     Category of algebras over Rational Field,
     Category of rings,
     Category of rngs,
     Category of semirings,
     Category of vector spaces over Rational Field,
     Category of modules over Rational Field,
     Category of bimodules over Rational Field on the left and Rational Field on the right,
     Category of left modules over Rational Field,
     Category of right modules over Rational Field,
     Category of commutative additive groups,
     Category of commutative additive monoids,
     Category of commutative additive semigroups,
     Category of additive magmas,
     Category of Cartesian products of monoids,
     Category of monoids,
     Category of Cartesian products of semigroups,
     Category of semigroups,
     Category of Cartesian products of magmas,
     Category of magmas,
     Category of Cartesian products of sets,
     Category of sets,
     Category of sets with partial maps,
     Category of objects]

Wrapup:

 - All the mathematical information about algebras with basis
   (properties of elements, parents, morphisms, cartesian products,
   tensor products) is gathered in AlgebrasWithBasis

Dynamic hierarchies for categories and classes
----------------------------------------------

Let `A` and `B` be two parents, and let us construct their cartesian
product `A\times B`. In which category should this new parent be?  For
example, if both `A` and `B` are monoids, then `A\times B` is
naturally endowed with a monoid structure. Similarly, it could
potentially be an algebra, a coalgebra, a differential module, and
possibly finite dimensional, graded, or ...

This can only be decided at runtime, by introspection into the
properties of `A` and `B`. Furthermore the number of possible
combinations (e.g. finite dimensional differential algebra) grows
exponnentially with the number of properties.

Functorial constructions:
 - Cartesian products (see :const:`~sage.categories.cartesian_product.cartesian_product`)
 - Tensor products (see :const:`~sage.categories.tensor.tensor`)
 - Subquotients / quotients / subobjects / isomorphic objects
   (see ``Sets().Subquotients()``, ``Sets().Quotients(), ``Sets().Subobjects()``,
   ``Sets().IsomorphicObjects()``
 - Dual
 - Algebras (see e.g. ``Sets().Algebras(QQ)``)

 - Morphisms

Full subcategories (work in progress):
 - finite dimensional / graded? / graded connected
 - finite / infinite
 - commutative

Wrapup:
 - There is a combinatorial explosion of potential classes
 - This explosion can be controlled by implementing "few" building
   blocks, and using dynamic classes to *compose* them together,
   lazily, according to the needs.

Writing a new category
----------------------

Each category should come with a good example, in sage.categories.examples.

The order between super categories should not be mathematically
relevant (otherwise this usually means the category hierarchy is
wrong). On the other hand, it should be consistent, to help Python
build the method resolution order for the generated classes (it always
respects inheritance, and tries to respect the order of the bases).

The current convention is to order them lexicographically w.r.t. the
following criterions:

 - Graded... or Finite dimensional... first
 - ...WithBasis first
 - Algebras before Coalgebras
 - Modules first

This gives the following order::

    sage: GradedHopfAlgebrasWithBasis(QQ).all_super_categories()
    [Category of graded hopf algebras with basis over Rational Field,
     Category of graded bialgebras with basis over Rational Field,
     Category of graded algebras with basis over Rational Field,
     Category of graded coalgebras with basis over Rational Field,
     Category of graded modules with basis over Rational Field,
     Category of graded hopf algebras over Rational Field,
     Category of graded bialgebras over Rational Field,
     Category of graded algebras over Rational Field,
     Category of graded coalgebras over Rational Field,
     Category of graded modules over Rational Field,
     Category of hopf algebras with basis over Rational Field,
     Category of bialgebras with basis over Rational Field,
     Category of algebras with basis over Rational Field,
     Category of coalgebras with basis over Rational Field,
     Category of modules with basis over Rational Field,
     Category of hopf algebras over Rational Field,
     Category of bialgebras over Rational Field,
     Category of algebras over Rational Field,
     Category of rings,
     Category of rngs,
     Category of semirings,
     Category of monoids,
     Category of semigroups,
     Category of magmas,
     Category of coalgebras over Rational Field,
     Category of vector spaces over Rational Field,
     Category of modules over Rational Field,
     Category of bimodules over Rational Field on the left and Rational Field on the right,
     Category of left modules over Rational Field,
     Category of right modules over Rational Field,
     Category of commutative additive groups,
     Category of commutative additive monoids,
     Category of commutative additive semigroups,
     Category of additive magmas,
     Category of sets,
     Category of sets with partial maps,
     Category of objects]

Todo: any better convention? Maybe we should further specify that subcategories of Modules() go first?

Caveats
-------

See :mod:`sage.misc.nested_class_test`

"""
