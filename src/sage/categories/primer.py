r"""
Elements, parents, and categories in Sage: a (draft of) primer

Abstract
========

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

    The general philosophy is that *Building mathematical information
    into the system yields more expressive, more conceptual and, at
    the end, easier to maintain and faster code* (within a programming
    realm; this would not necessarily apply to specialized libraries
    like gmp!).

Introduction: Sage as a library of objects and algorithms
=========================================================

The Sage library, with more than one million lines of code,
documentation, and tests, implements:

- Thousands of different kinds of objects (classes):

  Integers, polynomials, matrices, groups, number fields, elliptic
  curves, permutations, morphisms, languages, ... and a few racoons ...

- Tens of thousands methods and functions:

  Arithmetic, integer and polynomial factorization, pattern matching
  on words, ...

Some challenges
---------------

- How to organize this library?

  One needs some bookshelves to group together related objects and algorithms.

- How to ensure consistency?

  Similar objects should behave similarly::

    sage: Permutations(5).cardinality()
    120

    sage: GL(2,2).cardinality()
    6

    sage: A=random_matrix(ZZ,3,6,x=7)
    sage: L=LatticePolytope(A)
    sage: L.npoints()                # oops!   # random
    37

- How to ensure robustness?

- How to reduce duplication?

  Example: binary powering::

      sage: m = 3
      sage: m^8 == m*m*m*m*m*m*m*m == ((m^2)^2)^2
      True

      sage: m=random_matrix(QQ, 4, algorithm='echelonizable', rank=3, upper_bound=60); m # randomm = GL(2,QQ).random_element(); m  # random

  ::

      sage: m^8 == m*m*m*m*m*m*m*m == ((m^2)^2)^2
      True

  We want to implement binary powering only once, as *generic* code
  that will apply in all cases.


A bit of help from abstract algebra
===================================

The hierarchy of categories
---------------------------

What makes binary powering work in the above examples? In both cases,
we have *a set* endowed with a *multiplicative binary operation* which
is *associative*. Such a set is called a *semigroup*, and binary
powering works generally for any semigroup.

Sage knows about semigroups::

    sage: Semigroups()
    Category of semigroups

and sure enough, binary powering is defined there::

    sage: m._pow_.__module__
    'sage.categories.semigroups'

That's our bookshelf! And it's used in many places::

    sage: GL(2,ZZ) in Semigroups()
    True
    sage: NN in Semigroups()
    True

What else do we have on that bookshelf?::

    sage: Semigroups
    <class 'sage.categories.semigroups.Semigroups'>

For a less trivial bookshelf we can consider euclidean rings: once we
know how to do euclidean division in some set `R`, we can compute
gcd's in `R` generically using the Euclidean algorithm.

We are in fact very lucky: abstract algebra provides us right away
with a large and robust set of bookshelves which is the result of
centuries of work of mathematicians to identify the important
concepts. This includes for example::

    sage: Sets()
    Category of sets

    sage: Groups()
    Category of groups

    sage: Rings()
    Category of rings

    sage: Fields()
    Category of fields

    sage: HopfAlgebras(QQ)
    Category of hopf algebras over Rational Field

Each of the above is called a *category*. It typically specifies what
are the operations on the elements, as well as the axioms satisfied by
those operations. For example the category of groups specifies that a
group is a set endowed with a binary operation (the multiplication)
which is associative and admits a unit and inverses.

Each set in Sage knows which bookshelf of generic algorithms it can
use, that is to which category it belongs::

    sage: G = GL(2,ZZ)
    sage: G.category()
    Category of groups

In fact a group is a semigroup, and Sage knows about this::

    sage: Groups().is_subcategory(Semigroups())
    True
    sage: G in Semigroups()
    True

Altogether, our group gets algorithms from a bunch of bookshelves::

    sage: G.categories()
    [Category of groups, Category of monoids, Category of semigroups,
     ...,
     Category of magmas,
     Category of sets, ...]

Those can be viewed graphically::

    sage: g = Groups().category_graph()
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

In case ``dot2tex`` is not available, you can use instead::

    sage: g.show(vertex_shape=None, figsize=20)

And here is an overview of all categories in Sage::

    sage: g = sage.categories.category.category_graph()
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

Wrap-up: generic algorithms in Sage are organized in a hierarchy of
bookshelves modelled upon the usual hierarchy of categories provided
by abstract algebra.

Elements, Parents, Categories
-----------------------------

**Parent**
    A *parent* is a Python instance modelling a set of mathematical
    elements together with its additional (algebraic) structure.

    Examples include the ring of integers, the group `S_3`, the set of
    prime numbers, the set of linear maps between two given vector
    spaces, and a given finite semigroup.

    These sets are often equipped with additional structure: the set
    of all integers forms a ring. The main way of encoding this
    information is specifying which categories a parent belongs to.

    It is completely possible to have different Python instances
    modelling the same set of elements.  For example, one might want
    to consider the ring of integers, or the poset of integers under
    their standard order, or the poset of integers under divisibility,
    or the semiring of integers under the operations of addition and
    maximum.  Each of these would be a different instance, belonging
    to different categories.

    For a given model, there should be a unique instance in Sage
    representing that parent::

        sage: IntegerRing() is IntegerRing()
        True

**Element**
    An *element* is a Python instance modelling a mathematical element
    of a set.

    Examples of element include `5` in the integer ring, `x^3 - x` in
    the polynomial ring in `x` over the rationals, `4 + O(3^3)` in the
    3-adics, the transposition `(1 2)` in `S_3`, and the identity
    morphism in the set of linear maps from `\QQ^3` to `\QQ^3`.

    Every element in Sage has a parent.  The standard idiom in Sage
    for creating elements is to create their parent, and then provide
    enough data to define the element::

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

    Many parents model algebraic structures, and their elements
    support arithmetic operations. One often further wants to do
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
    A *category* is a Python instance modelling a mathematical category.

    Examples of categories include the category of finite semigroups,
    the category of all (Python) objects, the category of
    `\ZZ`-algebras, and the category of cartesian products of
    `\ZZ`-algebras::

        sage: FiniteSemigroups()
        Category of finite semigroups
        sage: Objects()
        Category of objects
        sage: Algebras(ZZ)
        Category of algebras over Integer Ring
        sage: Algebras(ZZ).CartesianProducts()
        Category of Cartesian products of algebras over Integer Ring

    Mind the 's' in the names of the categories above;
    ``GroupAlgebra`` and ``GroupAlgebras`` are distinct things.

    Every parent belongs to a collection of categories. Moreover,
    categories are interrelated by the *super categories*
    relation. For example, the category of rings is a super category
    of the category of fields, because every field is also a ring.

    A category serves two roles:

    - to provide a model for the mathematical concept of a category
      and the associated structures: homsets, morphisms, functorial
      constructions, axioms.

    - to organize and promote generic code, naming conventions,
      documentation, and tests across similar mathematical structures.

**CategoryObject**
    Objects of a mathematical category are not necessarily parents.
    Parent has a superclass that provides a means of modeling such.

    For example, the category of schemes does not have a faithful
    forgetful functor to the category of sets, so it does not make
    sense to talk about schemes as parents.

**Morphisms, Homsets**
    As category theorist will expect, *Morphisms* and *Homsets* will
    play an ever more important role, as support for them will
    improve.

Much of the mathematical information in Sage is encoded as relations
between elements and their parents, parents and their categories, and
categories and their super categories::

    sage: 1.parent()
    Integer Ring

    sage: ZZ
    Integer Ring

    sage: ZZ.category()
    Category of euclidean domains

    sage: ZZ.categories()
    [Category of euclidean domains, Category of principal ideal domains,
     Category of unique factorization domains, Category of gcd domains,
     Category of integral domains, Category of domains,
     Category of commutative rings, Category of rings, ...
     Category of distributive magmas and additive magmas,
     Category of monoids, Category of semigroups,
     Category of commutative magmas, Category of unital magmas, Category of magmas,
     Category of commutative additive groups, ..., Category of additive magmas,
     Category of sets,
     Category of sets with partial maps,
     Category of objects]

    sage: g = EuclideanDomains().category_graph()
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

A bit of help from computer science
===================================

Hierarchy of classes
--------------------

How are the bookshelves implemented in practice?

Sage uses the classical design paradigm of Object Oriented Programming
(OOP). It's fundamental principle is that any object that a program is
to manipulate should be modeled by an *instance* of a *class*. The
class implements:

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

    sage: x = var('x')
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
unrelated and deeply tied to their respective datastructures. On the
other hand, as we have seen above, they share the same powering method
because the sets `\ZZ` of integers, and the set `\QQ[x]` of
polynomials are both semigroups. Namely, the class for integers and
the class for polynomials both derive from an *abstract class* for
elements of a semigroups, which factors out the *generic* methods like
``_pow_``. This illustrates the use of *hierarchy of classes* to share
common code between classes having common behaviour.

OOP design is all about isolating the objects that one wants to model
and their operations, and designing an appropriate hierarchy of
classes for organizing the code. As we have seen above, the design of
the class hierarchy is easy since it can be modelled upon the
hierarchy of categories (bookshelves). Here is for example a piece of
the hierarchy of classes for an element of a group of matrices::

    sage: G = GL(2,ZZ)
    sage: m = G.an_element()
    sage: for cls in m.__class__.mro(): print cls
    <class 'sage.groups.matrix_gps.group_element.LinearMatrixGroup_gap_with_category.element_class'>
    <class 'sage.groups.matrix_gps.group_element.MatrixGroupElement_gap'>
    ...
    <class 'sage.categories.groups.Groups.element_class'>
    <class 'sage.categories.monoids.Monoids.element_class'>
    <class 'sage.categories.semigroups.Semigroups.element_class'>
    ...

On the top, we see concrete classes that describe the data structure
for matrices and provide the operations that are tied to this data
structure. Then follows abstract classes that are attached to the
hierarchy of categories and provide generic algorithms.

The full hierarchy is best viewed graphically:

    sage: g = class_graph(m.__class__)
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

Parallel hierarchy of classes for parents
-----------------------------------------

Let us recall that we do not just want to compute with elements of
mathematical sets, but with the sets themselves::

    sage: ZZ.one()
    1

    sage: R = QQ['x,y']
    sage: R.krull_dimension()
    2
    sage: A = R.quotient( R.ideal(x^2 - 2) )
    sage: A.krull_dimension() # todo: not implemented

Here are some typical operations that one may want to carry on various
kinds of sets:

    - The permutations of 5, the rational points of an elliptic curve: counting, listing, random generation

    - A language (set of words): rationality testing, counting elements, generating series

    - A finite semigroup: left/right ideals, center, representation theory

    - A vector space, an algebra: cartesian product, tensor product, quotient

Hence, following the OOP fundamental principle, parents should also be
modelled by instances of some (hierarchy of) classes. For example, our
group `G` is an instance of the following class:

    sage: G = GL(2,ZZ)
    sage: type(G)
    <class 'sage.groups.matrix_gps.linear.LinearMatrixGroup_gap_with_category'>

Here is a piece of the hierarchy of classes above it::

    sage: for cls in G.__class__.mro(): print cls
    <class 'sage.groups.matrix_gps.linear.LinearMatrixGroup_gap_with_category'>
    ...
    <class 'sage.categories.groups.Groups.parent_class'>
    <class 'sage.categories.monoids.Monoids.parent_class'>
    <class 'sage.categories.semigroups.Semigroups.parent_class'>
    ...

Note that the hierarchy of abstract classes is again attached to
categories and parallel to that we had seen for the elements. This is
best viewed graphically:

    sage: g = class_graph(m.__class__)
    sage: g.relabel(lambda x: x.replace("_","\_"))
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

.. NOTE::

    This is a progress upon systems like Axiom or MuPAD where a parent
    is modelled by the class of its elements; this oversimplification
    leads to confusion between methods on parents and elements, and
    makes parents special; in particular it prevents potentially
    interesting constructions like "groups of groups".

Sage categories
===============

Why this business of categories? And to start with, why don't we just
have a good old hierarchy of classes Group, Semigroup, Magma, ... ?

Dynamic hierarchy of classes
----------------------------

As we have just seen, when we manipulate groups, we actually
manipulate several kinds of objects:

- groups
- group elements
- morphisms between groups
- and even the category of groups itself!

Thus, on the group bookshelf, we want to put generic code for each of
the above. We therefore need three parallel hierarchies of abstract
classes:

- Group, Monoid, Semigroup, Magma, ...
- GroupElement, MonoidElement, SemigroupElement, MagmaElement, ...
- GroupMorphism, SemigroupElement, SemigroupMorphism, MagmaMorphism, ...

(and in fact many more as we will see).

We could implement the above hierarchies as usual::

    class Group(Monoid):
        # generic methods that apply to all groups

    class GroupElement(MonoidElement):
        # generic methods that apply to all group elements

    class GroupMorphism(MonoidMorphism):
        # generic methods that apply to all group morphisms

And indeed that's how it was done in Sage before 2009, and there are
still many traces of this. The drawback of this approach is
duplication: the fact that a group is a monoid is repeated three times
above!

Instead, Sage now uses the following syntax, where the :class:`Groups`
bookshelf is structured into units with *nested classes*::

    class Groups(Category):

        def super_categories(self):
            return [Monoids(), ...]

        class ParentMethods:
            # generic methods that apply to all groups

        class ElementMethods:
            # generic methods that apply to all group elements

        class MorphismMethods:
            # generic methods that apply to all group morphisms (not yet implemented)

        class SubcategoryMethods:
            # generic methods that apply to all subcategories of Groups()

With this syntax, the information that a group is a monoid is
specified only once, in the :meth:`Category.super_categories`
method. And indeed, when the category of unital inverse magmas was
introduced, there was a *single point of truth* to update::

    sage: Groups().super_categories()
    [Category of monoids, Category of inverse unital magmas]

The price to pay (there is no free lunch) is that some magic is
required to construct the actual hierarchy of classes for parents,
elements, and morphisms.  Namely, ``Groups.ElementMethods`` should be
seen as just a bag of methods, and the actual class
``Groups().element_class`` is constructed from it by adding the
appropriate super classes according to
``Groups().super_categories()``::

    sage: Groups().element_class
    <class 'sage.categories.groups.Groups.element_class'>

    sage: Groups().element_class.__bases__
    (<class 'sage.categories.monoids.Monoids.element_class'>,
     <class 'sage.categories.magmas.Unital.Inverse.element_class'>)

We now see that the hierarchy of classes for parents and elements is
parallel to the hierarchy of categories::

    sage: Groups().all_super_categories()
    [Category of groups,
     Category of monoids,
     Category of semigroups,
     ...
     Category of magmas,
     ...]

    sage: for cls in Groups().element_class.mro(): print cls
    <class 'sage.categories.groups.Groups.element_class'>
    <class 'sage.categories.monoids.Monoids.element_class'>
    <class 'sage.categories.semigroups.Semigroups.element_class'>
    ...
    <class 'sage.categories.magmas.Magmas.element_class'>
    ...
    sage: for cls in Groups().parent_class.mro(): print cls
    <class 'sage.categories.groups.Groups.parent_class'>
    <class 'sage.categories.monoids.Monoids.parent_class'>
    <class 'sage.categories.semigroups.Semigroups.parent_class'>
    ...
    <class 'sage.categories.magmas.Magmas.parent_class'>
    ...

Another advantage of building the hierarchy of classes dynamically is
that, for parametrized categories, the hierarchy may depend on the
parameters. For example an algebra over `\QQ` is a `\QQ`-vector space,
but not an algebra over `\ZZ` which is just a `\ZZ`-module!

.. NOTE::

    At this point this whole infrastructure may feel like
    overdesigning, right? We felt like this too! But we will see later
    that, once one gets used to it, this approach scales very
    naturally.

Categories are instances and have operations
--------------------------------------------

Note that categories themselves are naturally modelled by instances
because they can have operations of their own. An important one is::

    sage: Groups().example()
    General Linear Group of degree 4 over Rational Field

which gives an example of object of the category. Besides illustrating
the category, the example provides a minimal template for implementing
a new object in the category::

    sage: S = Semigroups().example(); S
    An example of a semigroup: the left zero semigroup

Its source code can be obtained by introspection::

    sage: S??                                     # not tested

This example is also typically used for testing generic methods. See
:meth:`Category.example` for more.

Other operations on categories include querying the super categories
or the axioms satisfied by the operations of a category::

    sage: Groups().super_categories()
    [Category of monoids, Category of inverse unital magmas]
    sage: Groups().axioms()
    frozenset(['Inverse', 'Associative', 'Unital'])

or constructing the intersection of two categories, or the smallest
category containing them::

    sage: Groups() & FiniteSets()
    Category of finite groups
    sage: Algebras(QQ) | Groups()
    Category of monoids

Specifications and generic documentation
----------------------------------------

Categories do not only contain code but also the specifications of the
operations. In particular a list of mandatory and optional methods to
be implemented can be found by introspection with::

    sage: Groups().required_methods()
    {'parent': {'required': ['__contains__'], 'optional': []},
     'element': {'required': [], 'optional': ['_mul_']}}

Documentation about those methods can be obtained with::

    sage: G = Groups()
    sage: G.element_class._mul_?        # not tested
    sage: G.parent_class.one?           # not tested

See also the :meth:`abstract_method` decorator.

.. WARNING::

    Well, more precisely, that's how things should be, but there is
    still some work to do in this direction. For example, the inverse
    operation is not specified above. Also, we are still missing a
    good programmatic syntax to specify the input and output types of
    the methods. Finally, in many cases the implementer must provide
    at least one of two methods, each having a default implementation
    using the other one (e.g. listing or iterating for a finite
    enumerated set); there is currently no good programmatic way to
    specify this.

Generic tests
-------------

Another feature that parents and elements receives from categories is
generic tests; their purpose is to check (at least to some extent)
that the parent satisfies the required mathematical properties (is my
semigroup indeed associative?) and is implemented according to the
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
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
    running ._test_enumerated_set_contains() . . . pass
    running ._test_enumerated_set_iter_cardinality() . . . pass
    running ._test_enumerated_set_iter_list() . . . pass
    running ._test_eq() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_pickling() . . . pass
    running ._test_some_elements() . . . pass

Tests can be run individually::

    sage: S._test_associativity()

Here is the code of this test::

    sage: S._test_associativity?? # not tested

Here is how to run the test on all elements::

    sage: L = S.list()
    sage: S._test_associativity(elements=L)

See :class:`TestSuite` for more information.

Let us see what happens when a test fails. Here we redefine the
product of `S` to something definitely not associative::

    sage: S.product = lambda x, y: S("("+x.value +y.value+")")

And rerun the test::

    sage: S._test_associativity(elements=L)
    Traceback (most recent call last):
    ...
      File ".../sage/categories/semigroups.py", line ..., in _test_associativity
        tester.assert_((x * y) * z == x * (y * z))
    ...
    AssertionError: False is not true

We can recover instantly the actual values of x,y,z, that is a counter
example to the associativity of our broken semigroup, using post
mortem introspection with the Python debugger ``pdb`` (this does not
work yet in the notebook)::

    sage: import pdb
    sage: pdb.pm()                       # not tested
    > /opt/sage-5.11.rc1/local/lib/python/unittest/case.py(424)assertTrue()
    -> raise self.failureException(msg)
    (Pdb) u
    > /opt/sage-5.11.rc1/local/lib/python2.7/site-packages/sage/categories/semigroups.py(145)_test_associativity()
    -> tester.assert_((x * y) * z == x * (y * z))
    (Pdb) p x, y, z
    ('a', 'a', 'a')
    (Pdb) p (x * y) * z
    '((aa)a)'
    (Pdb) p x * (y * z)
    '(a(aa))'

Wrap-up
-------

- Categories provide a natural hierarchy of bookshelves to organize
  not only code, but also specifications and testing tools.

- Everything about, say, algebras with a distinguished basis is
  gathered in :class:`AlgebrasWithBasis` or its super categories.
  This includes properties and algorithms for elements, parents,
  morphisms, but also, as we will see about constructions like
  cartesian products or quotients.

- The mathematical relations between elements, parents, and categories
  translates dynamically into a traditional hierarchy of classes.

- This design enforces robustness and consistency, which is
  particularly welcome given that Python is an interpreted language
  without static type checking.

Case study
==========

In this section, we study an existing parent in detail; a good
followup is to go through the :mod:`sage.categories.tutorial` to learn
how to implement a new one!

We consider the example of finite semigroup provided by the category::

    sage: S = FiniteSemigroups().example(); S
    An example of a finite semigroup: the left regular band generated by ('a', 'b', 'c', 'd')
    sage: S?                    # not tested

Where do all the operations on ``S`` and its elements come from?

::

    sage: x = S('a')

``_repr_`` is a technical method which comes with the data structure
(:class:`ElementWrapper`); since it's implemented in Cython, we need
to use Sage's introspection tools to recover where it's implemented::

    sage: x._repr_.__module__
    sage: sage.misc.sageinspect.sage_getfile(x._repr_)
    '.../sage/structure/element_wrapper.pyx'

``__pow__`` is a generic method for all finite semigroups::

    sage: x.__pow__.__module__
    'sage.categories.semigroups'

``_mul_`` is a default implementation from the :class:`Magmas`
category (a *magma* is a set with an inner law `*`, not necessarily
associative)::

    sage: x.__mul__.__module__
    'sage.categories.magmas'

It delegates the work to the parent (following the advice: if you do
not know what to do, ask your parent)::

    sage: x.__mul__??                             # not tested

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

    sage: S??                                     # not tested

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

    sage: FiniteSemigroups().required_methods()
    {'parent': {'required': ['__contains__'], 'optional': []},
     'element': {'required': [], 'optional': ['_mul_']}}

Documentation about those methods can be obtained with::

    sage: C = FiniteSemigroups().element_class
    sage: C._mul_?                                # not tested

See also the :meth:`abstract_method` decorator.

Here is the code for the finite semigroups category::

    sage: FiniteSemigroups??                      # not tested

.. TODO::

    - Specifying the category of a parent
    - Adding code to an existing parent

Scaling further: functorial constructions, axioms, ...
======================================================

In this section, we explore more advanced features of categories.
Along the way, we illustrate that a large hierarchy of categories is
desirable to model complicated mathematics, and that scaling to
support such a large hierarchy is the driving motivation for the
design of the category infrastructure.

Functorial constructions
------------------------

Sage has support for a certain number of so-called *covariant
functorial constructions* which can be used to construct new sets from
existing ones while carrying over as much as possible of their
algebraic structure. This includes:

 - Cartesian products
   (see :const:`~sage.categories.cartesian_product.cartesian_product`)
 - Tensor products
   (see :const:`~sage.categories.tensor.tensor`)
 - Subquotients / quotients / subobjects / isomorphic objects
   (see ``Sets().Subquotients()``, ``Sets().Quotients()``,
   ``Sets().Subobjects()``, ``Sets().IsomorphicObjects()``)
 - Dual
 - Algebras, as in group algebras, monoid algebras, ...
   (see e.g. ``Sets().Algebras(QQ)``)

Let for example `A` and `B` be two parents, and let us construct their
cartesian product `A\times B\times B`::

    sage: A = AlgebrasWithBasis(QQ).example();     A.rename("A")
    sage: B = HopfAlgebrasWithBasis(QQ).example(); B.rename("B")
    sage: C = cartesian_product([A, B, B]); C
    A (+) B (+) B

In which category should this new parent be? Since `A` and `B` are
vector spaces, the result is, as a vector space, the direct sum
`A\oplus B\oplus B`, hence the notation. Also, since both `A` and `B`
are monoids, `A\times B\times B` is naturally endowed with a monoid
structure for pointwise multiplication::

    sage: C in Monoids()
    True

the unit being the cartesian product of the units of the operands::

    sage: C.one()
    B[(0, word: )] + B[(1, ())] + B[(2, ())]
    sage: cartesian_product([A.one(), B.one(), B.one()])
    B[(0, word: )] + B[(1, ())] + B[(2, ())]

The pointwise product can be implemented generically for all magmas
(i.e. sets endowed with a multiplicative operation) that are
constructed as cartesian products. It's thus implemented in the
:class:`Magmas` category::

    sage: C.product.__module__
    'sage.categories.magmas'

More specifically, keeping on using nested classes to structure the
code, the product method is put in the nested class
:class:`Magmas.CartesianProducts.ParentMethods`::

    class Magmas(Category):
        class ParentMethods:
            # methods for magmas
        class ElementMethods:
            # methods for elements of magmas
        class CartesianProduct(CartesianProductCategory):
            class ParentMethods:
                # methods for magmas that are constructed as cartesian products
                def product(self, x, y):
                    # ...
            class ElementMethods:
                # ...

.. NOTE::

    The support for nested classes in Python is relatively
    recent. Their intensive use for the category infrastructure did
    reveal some glitches in their implementation, in particular around
    class naming and introspection. Sage currently works around the
    more anoying ones but some remain visible. See
    e.g. :mod:`sage.misc.nested_class_test`.


Let us now look at the categories of ``C``::

    sage: C.categories()
    [Category of Cartesian products of algebras with basis over Rational Field, ...
     Category of Cartesian products of monoids, Category of monoids,
     Category of Cartesian products of semigroups, Category of semigroups,
     Category of Cartesian products of magmas, ..., Category of magmas, ...
     Category of Cartesian products of sets, Category of sets, ...]

This reveals the parallel hierarchy of categories for cartesian
products of monoids, semigroups, ... We are thus glad that Sage uses
its knowledge that a monoid is a semigroup to automatically deduce
that a cartesian product of monoids is a cartesian product of
semigroups, and build the hierarchy of classes for parents and
elements accordingly.

In general, the cartesian product of `A` and `B` can potentially be an
algebra, a coalgebra, a differential module, and be finite
dimensional, or graded, or ...  This can only be decided at runtime,
by introspection into the properties of `A` and `B`; furthermore, the
number of possible combinations (e.g. finite dimensional differential
algebra) grows exponnentially with the number of properties.

Axioms
------

We have seen that several categories can be defined by specifying the
axioms that are satisfied by the operations of its super
categories. For example, starting from the category of magmas we can
build all the following categories just by specifying the axioms
satisfied by the multiplication::

    sage: Magmas()
    Category of magmas
    sage: Magmas().Unital()
    Category of unital magmas

::

    sage: Magmas().Commutative().Unital()
    Category of commutative unital magmas
    sage: Magmas().Unital().Commutative()
    Category of commutative unital magmas

::

    sage: Magmas().Associative()
    Category of semigroups

::

    sage: Magmas().Associative().Unital()
    Category of monoids

::

    sage: Magmas().Associative().Unital().Commutative()
    Category of commutative monoids

Note by the way that intersecting categories plays nicely::

    sage: Magmas().Associative() & Magmas().Unital().Inverse() & Sets().Finite()
    Category of finite groups

Even for a very simple category like ``Magmas``, there are about `2^5`
potential combinations of the axioms! Think about what this becomes
for a category with two operations ``+`` and ``*``::

    sage: DistributiveMagmasAndAdditiveMagmas()
    Category of distributive magmas and additive magmas

    sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse()
    Category of rngs

    sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().Unital()
    Category of semirings

    sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse().Unital()
    Category of rings

    sage: Rings().Division()
    Category of division rings

    sage: Rings().Division().Commutative()
    Category of fields

    sage: Rings().Division() & Sets().Finite()    # Wedderburn theorem!
    Category of finite fields

or for more advanced categories::

    sage: g = GradedHopfAlgebrasWithBasis(QQ).category_graph()
    sage: g.set_latex_options(format="dot2tex")
    sage: view(g, tightpage=True)                 # not tested

For more about axioms, see :mod:`sage.categories.category_with_axiom`.

Wrap-up
-------

As we have seen, there is a combinatorial explosion of possible
classes. Constructing by hand the full class hierarchy would not scale
unless one would restrict to a very rigid subset. Even if it was
possible to construct automatically the full hierarchy, this would not
scale with respect to system ressources.

When designing software systems with large hierarchies of abstract
classes for business objects, the difficulty is usually to identify a
proper set of key concepts. Here we are lucky, as the key concepts
have been long identified and are relatively few:

- Operations (`+`, `*`, ...)
- Axioms on those operations (associativity, ...)
- Constructions (cartesian products, ...)

Better, those concepts are sufficiently well known so that a user can
reasonably be expected to be familiar with the concepts that are
involved for his own needs.

Instead, the difficulty is concentrated in the huge number of possible
combinations, an unpredictable large subset of which being potentially
of interest; at the same time, only a small -- but moving -- subset
has code naturally attached to it.

This has led to the current design where one focuses on writing the
relatively few classes for which there is actual code or mathematical
information, and let Sage *compose dynamically and lazily* those
building blocks to construct the minimal hierarchy of classes needed
for the computation at hand. This allows for the infrastructure to
scale smoothly as bookshelves are added, extended, or reorganized.


Writing a new category
======================

Each category `C` **must** be provided with a method
``C.super_categories()`` and *can* be provided with a method
``C._subcategory_hook_(D)``. Also, it may be needed to insert `C` into
the output of the ``super_categories()`` method of some other
category. This determines the position of `C` in the category graph.

A category *may* provide methods that can be used by all its objects,
respectively by all elements of its objects.

Each category *should* come with a good example, in
:mod:`sage.categories.examples`.

Inserting the new category into the category graph
--------------------------------------------------

``C.super_categories()`` must return a list of categories, namely the
*immediate* super categories of `C`.  Of course, if you know that your
new category `C` is an immediate super category of some existing
category `D`, then you should also update `D`'s ``super_categories``
to include `C`.

In several cases, the category `C` is directly provided with a generic
implementation of ``super_categories``; a typical example is when `C`
implements an axiom or a functorial construction; in such a case, `C`
may implement ``C.extra_super_categories()`` to complement the super
categories discovered by the generic implementation. This method needs
not return immediate super categories; instead it's usually best to
specify the largest super category providing the desired mathematical
information. For example, the category
:class:`Magmas.Commutative.Algebras` just states that the algebra of a
commutative magma is a commutative magma. This is sufficient to let
Sage deduce that it's in fact a commutative algebra.

Methods for objects and elements
--------------------------------

Different objects of the same category share some algebraic features, and
very often these features can be encoded in a method, in a generic way.
For example, for every commutative additive monoid, it makes sense to ask
for the sum of a list of elements. Sage's category framework allows to
provide a generic implementation for all objects of a category.

If you want to provide your new category with generic methods for
objects (or elements of objects), then you simply add a nested class
called ``ParentMethods`` (or ``ElementMethods``). The methods of that
class will automatically become methods of the objects (or the
elements). For instance::

    sage: P.<x,y> = ZZ[]
    sage: P.prod([x,y,2])
    2*x*y
    sage: P.prod.__module__
    'sage.categories.monoids'
    sage: P.prod.__func__ is Monoids().ParentMethods.prod.__func__
    True

We recommend to study the code of one example::

    sage: C = CommutativeAdditiveMonoids()
    sage: C??                               # not tested

On the order of super categories
--------------------------------

The generic method ``C.all_super_categories()`` determine recursively
the list of *all* super categories of `C`.

The order of the categories in this list does influences the
inheritance of methods for parents and elements. Namely, if `P` is an
object in the category `C` and if `C_1` and `C_2` are both super
categories of `C` defining some method ``foo`` in ``ParentMethods``,
then `P` will use `C_1`'s version of ``foo`` if and only if `C_1`
appears in ``C.all_super_categories()`` before `C_2`.

However this must be considered as an *implementation detail*: if
`C_1` and `C_2` are incomparable categories, then the order in which
they appear must be mathematically irrelevant: in particular, the
methods ``foo`` in `C_1` and `C_2` must have the same semantic. Code
should not rely on any specific order, as it it subject to later
change. In case one of the implementations is prefered in a common
subcategory of `C_1` and `C_2`, for example for efficiency reasons,
then the ambiguity should be resolved explicitly by definining a
method ``foo`` in this category. See the method ``some_elements`` in
the code of the category :class:`FiniteCoxeterGroups` for an example.

Since :trac:`11943`, ``C.all_super_categories()`` is computed by the
so-called ``C3`` algorithm used by Python to compute Method Resolution
Order of new-style classes. Thus the order in
``C.all_super_categories()``, ``C.parent_class.mro()`` and
``C.element_class.mro()`` are guaranteed to be consistent.

Since :trac:`13589`, the ``C3`` algorithm is put under control of some
total order on categories. This order is not necessarily meaningful,
but it guarantees that ``C3`` always finds a consistent Method
Resolution Order. For background, see
:mod:`sage.misc.c3_controlled`. A visible effect is that the order in
which categories are specified in ``C.super_categories()``, or in a
join category, no longer influences the result of
``C.all_super_categories()``.

Subcategory hook (advanced optimization feature)
------------------------------------------------

The default implementation of the method ``C.is_subcategory(D)`` is to
look up whether `D` appears in ``C.all_super_categories()``. However,
building the list of all the super categories of `C` is an expensive
operation that is sometimes best avoided. For example, if both `C` and
`D` are categories defined over a base, but the bases differ, then one
knows right away that they can not be subcategories of each other.

When such a short-path is known, one can implement a method
``_subcategory_hook_``. Then, ``C.is_subcategory(D)`` first calls
``D._subcategory_hook_(C)``. If this returns ``Unknown``, then
``C.is_subcategory(D)`` tries to find ``D`` in
``C.all_super_categories()``. Otherwise, ``C.is_subcategory(D)``
returns the result of ``D._subcategory_hook_(C)``.

By default, ``D._subcategory_hook_(C)`` tests whether
``issubclass(C.parent_class,D.parent_class)``, which is very often
giving the right answer::

    sage: Rings()._subcategory_hook_(Algebras(QQ))
    True
    sage: HopfAlgebras(QQ)._subcategory_hook_(Algebras(QQ))
    False
    sage: Algebras(QQ)._subcategory_hook_(HopfAlgebras(QQ))
    True
"""
