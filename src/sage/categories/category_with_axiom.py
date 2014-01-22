"""
Axioms

This documentation covers how to implement axioms and proceeds with an
overview of the implementation of the axiom infrastructure. It assumes
that the reader is familiar with the :mod:`category primer
<sage.categories.primer>`, and in particular its :ref:`section about
axioms <category-primer-axioms>`.

Implementing axioms
===================

Simple case involving a single predefined axiom
-----------------------------------------------

Suppose that one wants to provide code (and documentation, tests, ...)
for the objects of some existing category ``Cs()`` that satisfy some
predefined axiom ``A``.

The first step is to open the hood and check whether there already
exists a class implementing the category ``Cs().A()``. For example,
taking ``Cs=Semigroups`` and the ``Finite`` axiom, there already
exists a class for the category of finite semigroups::

    sage: Semigroups().Finite()
    Category of finite semigroups
    sage: type(Semigroups().Finite())
    <class 'sage.categories.finite_semigroups.FiniteSemigroups_with_category'>

In this case, we say that the category of semigroups *implements* the
axiom ``Finite``, and code about finite semigroups should go in the
class :class:`FiniteSemigroups` (or, as usual, in its nested classes
ParentMethods, ElementsMethods and so on).

On the other hand, there is no class for the category of infinite
semigroups::

    sage: Semigroups().Infinite()
    Category of infinite semigroups
    sage: type(Semigroups().Infinite())
    <class 'sage.categories.category.JoinCategory_with_category'>

This category is indeed just constructed as the intersection of the
categories of semigroups and of infinite sets respectively::

    sage: Semigroups().Infinite().super_categories()
    [Category of semigroups, Category of infinite sets]

In this later case, one needs to create a new class to implement the
axiom for this category. This boils down to adding a nested class
inheriting from :class:`CategoryWithAxiom`. Here we implement a
category ``Cs``, with a subcategory for the objects satisfying the
``Finite`` axiom defined in the super category ``Sets`` (we will see
later on how to *define* new axioms)::

    sage: from sage.categories.category_with_axiom import CategoryWithAxiom
    sage: class Cs(Category):
    ....:     def super_categories(self):
    ....:         return [Sets()]
    ....:     class Finite(CategoryWithAxiom):
    ....:         class ParentMethods:
    ....:             def foo(self):
    ....:                 print "I am a method on finite C's"

::

    sage: Cs().Finite()
    Category of finite cs
    sage: Cs().Finite().super_categories()
    [Category of finite sets, Category of cs]
    sage: Cs().Finite().all_super_categories()
    [Category of finite cs, Category of finite sets,
     Category of cs, Category of sets, ...]
    sage: Cs().Finite().axioms()
    frozenset(['Finite'])

Now a parent declared in the category ``Cs().Finite()`` inherits from
all the methods of finite sets and of finite `C`'s, as desired::

    sage: P = Parent(category=Cs().Finite())
    sage: P.is_finite()             # Provided by Sets.Finite.ParentMethods
    True
    sage: P.foo()                   # Provided by Cs.Finite.ParentMethods
    I am a method on finite C's

.. NOTE::

    - This follows the same idiom as for
      :mod:`covariant functorial constructions
      <sage.categories.covariant_functorial_construction>`.

      From an object oriented point of view, any subcategory ``Cs`` of
      Sets inherits from a ``Finite`` method, and it can complement
      this method with extra data (here a mixin class) in the form of
      a class attribute ``Cs.Finite``. The fact that ``Cs.Finite``
      *complements* the ``Finite`` method rather than *overriding* is
      not unnatural: that's what sequences of ``super`` calls are
      usually about.

    - Under the hood, the category ``Cs().Finite()`` is aware that it
      has been constructed from the category ``Cs()`` and adding the
      axiom ``Finite``::

        sage: Cs().Finite().base_category()
        Category of cs
        sage: Cs().Finite()._axiom
        'Finite'

Over time, the nested class ``Cs.Finite`` may become large and too
cumbersome to keep as a nested class of ``Cs``. Or the category with
axiom may have a name of its own in the litterature, like *semigroups*
rather than *associative magmas*, or *fields* rather than *commutative
division rings*. In this case, the category with axiom can be put
elsewhere, typically in a separate file, with just links to and from
``Cs``::

    sage: class Cs(Category):
    ....:     def super_categories(self):
    ....:         return [Sets()]
    sage: class FiniteCs(CategoryWithAxiom):
    ....:     _base_category_class_and_axiom = (Cs, 'Finite')
    ....:     class ParentMethods:
    ....:         def foo(self):
    ....:             print "I am a method on finite C's"
    sage: Cs.Finite = FiniteCs

For a real example, look up the code of the class
:class:`FiniteGroups` and the link to it in :class:`Groups`. Note that
the link is implemented using :class:`LazyImport`; this is highly
recommended: it makes sure that :class:`FiniteGroups` is imported
after :class:`Groups` it depends upon, and makes it explicit that the
class :class:`Groups` can be imported and is fully functional without
importing :class:`FiniteGroups`.

.. NOTE::

    Some categories with axioms are created upon Sage's startup. In
    such a case, one needs to pass the ``at_startup=True`` option to
    :class:`LazyImport`, in order to quiet the warning about that lazy
    import being resolved upon startup. See for example ``Sets.Finite``.

    This is undoubtly a code smell. Nethertheless it should be kept as
    is, first to resolve the import order properly, and more
    importantly as a reminder that the category would be best not
    constructed upon Sage's startup. This to entice developpers to
    reduce the number of parents (and therefore categories) that are
    constructed upon startup. Each "at_startup=True" that will be
    removed will be a measure of progress in this direction.

Note that the category with axiom can be created either directly or
through its base category::

    sage: FiniteGroups()
    Category of finite groups
    sage: Groups().Finite()
    Category of finite groups
    sage: Groups().Finite() is FiniteGroups()
    True

For the former idiom to work, and with the current implementation of
axioms, the class :class:`FiniteGroups` needs to be aware of the base
category class (here, :class:`Groups`) and of the axiom (here,
``Finite``)::

    sage: FiniteGroups._base_category_class_and_axiom
    (<class 'sage.categories.groups.Groups'>, 'Finite')
    sage: Semigroups._base_category_class_and_axiom
    (<class 'sage.categories.magmas.Magmas'>, 'Associative')
    sage: Fields._base_category_class_and_axiom
    (<class 'sage.categories.division_rings.DivisionRings'>, 'Commutative')
    sage: FiniteDimensionalAlgebrasWithBasis._base_category_class_and_axiom
    (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>, 'FiniteDimensional')

As a syntactic sugar, Sage tries some obvious heuristics to guess
those from the name of the category with axiom (see
:func:`base_category_class_and_axiom` for the details). When this
fails, typically because the category has a name of its own like
:class:`Fields`, the attribute ``_base_category_class_and_axiom``
should be set explicitly. For more examples, see the code of the
classes :class:`Semigroups` or :class:`Fields`.

In our example ``FiniteCs``, Sage failed to guess automatically the
base category class and axiom because the class ``Cs`` is not in the
standard location ``sage.categories.cs``.

Defining a new axiom
--------------------

We describe now how to define a new axiom. The first step is to figure
out the largest category where the axiom makes sense. For example
``Sets`` for ``Finite``, ``Magmas`` for ``Associative``, or
``Modules`` for ``FiniteDimensional``. Here we define the axiom
``Green`` for the category ``Cs`` and its subcategories::

    sage: from sage.categories.category_with_axiom import CategoryWithAxiom
    sage: class Cs(Category):
    ....:     def super_categories(self):
    ....:         return [Sets()]
    ....:     class SubcategoryMethods:
    ....:         def Green(self):
    ....:             '<documentation of the axiom Green>'
    ....:             return self._with_axiom("Green")
    ....:     class Green(CategoryWithAxiom):
    ....:         class ParentMethods:
    ....:             def foo(self):
    ....:                 print "I am a method on green C's"

With the current implementation, the name of the axiom must also be
added to a global tuple::

    sage: sage.categories.category_with_axiom.all_axioms += ("Green",)

We can now use the axiom as usual::

    sage: Cs().Green()
    Category of green cs

    sage: P = Parent(category=Cs().Green())
    sage: P.foo()
    I am a method on green C's

Compared with our first example, the only newcomer is the method
``.Green()`` that can be used by any subcategory ``Ds()`` of ``Cs()``
to add the axiom ``Green``. Due to some magic, ``Ds().Green`` always
gives back this method, regardless of whether ``Ds`` has a nested
class ``Ds.Green`` or not (an implementation detail)::

    sage: Cs().Green
    <bound method Cs_with_category.Green of Category of cs>

Thanks to this feature, the user is systematically referred to the
documentation of this method when doing introspection on
``Ds().Green``::

    sage: C = Cs()
    sage: C.Green?             # not tested
    sage: Cs().Green.__doc__
    '<documentation of the axiom Green>'

It is therefore the natural spot for the documentation of the axiom.

.. NOTE::

    The presence of the nested class ``Green`` in ``Cs`` is mandatory
    even if it is empty.

.. TODO::

    Specify whether or not one should systematically use
    @cached_method in the definition of the axiom. And make sure all
    the definition of axioms in Sage are consistent in this respect!

.. TODO::

    We could possibly define an @axiom decorator? This could hide two
    little implementation details: whether or not to make the method a
    cached method, and the call to _with_axiom(...) under the hood. It
    could do possibly do some more magic. The gain is not obvious though.

.. NOTE::

    ``all_axioms`` is only used marginally, for sanity checks and when
    trying to guess the base category class. The order of the axioms
    in this tuple also controls the order in which they appear when
    printing categories with axioms (see
    :meth:`CategoryWithAxiom._repr_object_names_static`).

    During a Sage session, new axioms should only be added at the end
    of ``all_axioms`` as above, so as to not break the cache of
    :func:`axioms_rank`. Otherwise, they can be inserted statically
    anywhere in the tuple. For axioms defined within the Sage library,
    the name is best inserted by editing directly the definition of
    ``all_axioms`` in :mod:`category_with_axiom`.

.. NOTE::

    Let us state again that, unlike what the existence of
    ``all_axioms`` might suggests, the definition of an axiom is local
    to a category and its subcategories. In particular, two
    independent categories ``Cs()`` and ``Ds()`` can very well define
    axioms with the same name and different semantic. As long as the
    two hierarchies of subcategories don't intersect, this is not a
    problem. And if they do intersect naturally (that is if one is
    likely to create a parent belonging to both categories), this
    probably means that ``Cs`` and ``Ds`` are about related enough
    areas of mathematics that one should clear the ambiguity by having
    either the same semantic or different names.

    This caveat is no different from that of name clash issues in
    hierarchy of classes involving multiple inheritance.

.. TODO::

    Explore ways to get rid of this global ``all_axioms`` tuple,
    and/or have automatic registration there, and/or having a
    register_axiom(...) method.


Handling multiple axioms, tree structure of the code
----------------------------------------------------

Prelude
^^^^^^^

Let us consider the category of magmas, together with two of its
axioms, namely ``Associative`` and ``Unital``. An associative magma is
a semigroup and a unital semigroup is a monoid; we have also seen that
axioms commute::

    sage: Magmas().Unital()
    Category of unital magmas
    sage: Magmas().Associative()
    Category of semigroups
    sage: Magmas().Associative().Unital()
    Category of monoids
    sage: Magmas().Unital().Associative()
    Category of monoids

At the level of the classes implementing those categories, the
following comes as a general naturalization of the previous section::

    sage: Magmas.Unital
    <class 'sage.categories.magmas.Magmas.Unital'>
    sage: Magmas.Associative
    <class 'sage.categories.semigroups.Semigroups'>
    sage: Magmas.Associative.Unital
    <class 'sage.categories.monoids.Monoids'>

However, the following may look suspicous at first::

    sage: Magmas.Unital.Associative
    Traceback (most recent call last):
    ...
    AttributeError: type object 'Magmas.Unital' has no attribute 'Associative'

The purpose of this section is to explain the design of the code
layout and the rationale for this mismatch.

Abstract model
^^^^^^^^^^^^^^

As we have seen in the :ref:`Primer <category-primer-axioms-explosion>`,
the objects of a category ``Cs()`` can usually satisfy, or not, many
different axioms. Out of all combinations of axioms, only a small
number are relevant in practice; in the sense that we actually want to
provide features for the objects satisfying those axioms.

Therefore, in the context of the category class `Cs`, we want to
provide the system with a collection `(D_S)_{S\in \mathcal S}` were
`S` is a subset of the axioms and `D_S` is a class for the subcategory
of the object of ``Cs()`` satisfying the axioms in `S`. For example,
if ``Cs()`` is the category of magmas, the pairs would include::

    {Associative}                 : Semigroups
    {Associative, Unital}         : Monoids
    {Associative, Unital, Inverse}: Groups
    {Associative, Commutative}    : Commutative Semigroups
    {Unital,      Inverse}        : Loops

Then, given a subset `T` of axioms, we want the system to be able to
select automatically the relevant classes
`(D_S)_{S\in \mathcal S, S\subset T}`,
and build from them a category for the objects of ``Cs`` satisfying
the axioms in `T`, together with its hierarchy of super categories. If
`T` is in `\mathcal S`, then the class of the resulting category is
directly `D_T`::

    sage: C = Magmas().Unital().Inverse().Associative(); C
    Category of groups
    sage: type(C)
    <class 'sage.categories.groups.Groups_with_category'>

Otherwise, we get a join category::

    sage: C = Magmas().Infinite().Unital().Associative(); C
    Category of infinite monoids
    sage: type(C)
    <class 'sage.categories.category.JoinCategory_with_category'>
    sage: C.super_categories()
    [Category of monoids, Category of infinite sets]

Concrete model as a tree of nested classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We further want the construction to be efficient and amenable to
lazyness. This led us to the following design decision: the collection
`(D_S)_{S\in \mathcal S}` of classes should be structured as a rooted
tree. The root is ``Cs``, corresponding to `S=\emptyset`. Any other
class `D_S` should be the child of a single class `D_{S'}` where `S'`
is obtained from `S` by removing exactly one axiom `A`. Of course,
`D_{S'}` and `A` are respectively the base category class and axiom of
the category with axiom `D_S` that we have met in the first section.

At this point, we urge the reader to explore the code of
:class:`Magmas` and :class:`DistributiveMagmasAndAdditiveMagmas` and
see how the tree structure on the categories with axioms is reflected
by the nesting of category classes.

Discussion of the design
^^^^^^^^^^^^^^^^^^^^^^^^

Asymmetry
~~~~~~~~~

As we have seen at the beginning of this section this design
introduces an asymmetry. It's not so bad in practice as in most
practical cases, we want to work incrementally. It's for example more
natural to describe :class:`FiniteFields` as :class:`Fields` with the
axiom ``Finite`` rather than :class:`Magmas` and
:class:`AdditiveMagmas` with all (or at least sufficiently many) of
the following axioms (not counting ``Distributive`` which is not yet
implemented)::

    sage: Fields().axioms()
    frozenset(['Division', 'AdditiveUnital', 'NoZeroDivisors', 'Commutative',
               'AdditiveInverse', 'AdditiveAssociative', 'Unital',
	       'AdditiveCommutative', 'Associative'])

The main limitation is that the infrastructure currently imposes to be
incremental by steps of one axiom.

In practice, among the roughly 60 categories with axioms that are
currently implemented in Sage, most admitted a (rather) natural choice
of a base category and single axiom to add. For example, one usually
thinks more naturally of a monoid as a semigroup which is unital
rather than as a unital magma which is associative. Only in a few
cases is a choice made that feels mathematically arbitrary. This is
essentially in the chain of nested classes
:class:`DistributiveMagmasAndAdditiveMagmas.AdditiveAssociative.AdditiveCommutative.AdditiveUnital`.

Placeholder classes
~~~~~~~~~~~~~~~~~~~

Given that we can only add a single axiom at a time when implementing
a :class:`CategoryWithAxiom`, we need to create a few category classes
that are just placeholders. See for example the chain of nested
classes
:class:`DistributiveMagmasAndAdditiveMagmas.AdditiveAssociative.AdditiveCommutative.AdditiveUnital`.

Yet, this is well within the scope of the axiom infrastructure which
is to go from a potentially exponential number of placeholder
categories to just a couple.

Mismatch between the tree of nested classes and the hierarchy of categories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This may sound suspicious at first! However, as mentionned in the
primer, this is actually a big selling point of the axioms
infrastructure: by calculating automatically the hierarchy relation
between categories with axioms one avoids the nightmare of maintaining
it by hand.  Instead, only a rather minimal number of links needs to
be maintainted in the code (one per category with axiom).

Besides, with the flexibility introduced by runtime deduction rules
(see below) the hierarchy of categories may depend on the parameters
of the categories and not just their class. So it's best to make it
clear from the onset that the two relations do not match.

Flexibility
~~~~~~~~~~~

This design also brings in quite some flexibility, with the
possibility to support features such as defining new axioms depending
on other axioms. See below.

Limitations
~~~~~~~~~~~

As a small variant of what we have seen above, defining the
``Distributive`` axiom within this infrastructure would require to
implement a placeholder category ``MagmasAndAdditiveMagmas``. What is
more annoying is that Sage can't currently be taught that
``MagmasAndAdditiveMagmas`` is the join of :class:`Magmas` and
:class:`AdditiveMagmas`, which is why ``MagmasAndAdditiveMagmas``
together with its ``Distributive`` axiom are not yet implemented.

.. TODO::

   Make :meth:`Category.join` more powerful, to handle such
   joins. Another typical use case is::

       sage: C = Algebras(QQ) & Coalgebras(QQ); C   # todo: not implemented
       Category of bialgebras

       sage: C.Hopf()                               # todo: not implemented
       Category of hopf algebras

    One could also have the category :class:`Bimodules` be the join of
    :class:`LeftModules` and :class:`RightModules` with an axiom
    stating the compatibility between the left and right action.

    This will require some thought on appropriate idioms and
    algorithms. Potential directions include refactoring the whole
    implementation with a database mapping sets of axioms to
    categories together with other deduction rules, or just adding
    some hook or registration process to specify what the join of two
    categories is.


Axioms defined upon other axioms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes an axiom can only be defined when some other axiom
holds. For example, the axiom ``NoZeroDivisors`` only makes sense if
there is a zero, that is if the axiom ``AdditiveUnital`` holds. Hence,
for the category :class:`DistributiveMagmasAndAdditiveMagmas`, we
consider in the abstract model only those subsets of axioms where the
presence of ``NoZeroDivisors`` implies that of ``AdditiveUnital``.  We
also want the axiom to be only available if meaningful::

    sage: Rings().NoZeroDivisors()
    Category of domains
    sage: Rings().Commutative().NoZeroDivisors()
    Category of integral domains
    sage: DistributiveMagmasAndAdditiveMagmas().NoZeroDivisors()
    Traceback (most recent call last):
    ...
    AttributeError: 'DistributiveMagmasAndAdditiveMagmas_with_category' object has no attribute 'NoZeroDivisors'

Concretely, this is to be implemented by defining the new axiom in the
(``SubcategoryMethods`` nested class of the) appropriate category with
axiom. For example the axiom ``NoZeroDivisors`` would be naturally
defined in :class:`DistributiveMagmasAndAdditiveMagmas.AdditiveUnital`.

.. NOTE::

    The axiom ``NoZeroDivisors`` is currently defined in
    :class:`Rings`, by simple lack of need for the feature; it should
    be lifted up as soon as relevant, that is when some code will be
    available for parents with no zero divisors that are not
    necessarily rings.

.. _axioms-deduction-rules:

Deduction rules
^^^^^^^^^^^^^^^

A similar situation is when an axiom ``A`` of a category ``Cs``
implies some other axiom B, with the same consequence as above on the
subsets of axioms appearing in the abstract model. For example, a
division ring necessarily has no zero divisors::

    sage: 'NoZeroDivisors' in Rings().Division().axioms()
    True
    sage: 'NoZeroDivisors' in Rings().axioms()
    False

This deduction rule is implemented by the method
:meth:`Rings.Division.extra_super_categories`::

    sage: Rings().Division().extra_super_categories()
    (Category of domains,)

In general, this is to be implemented by a method
``Cs.A.extra_super_categories`` returning a tuple ``(Cs().B(),)``, or
preferably ``(Ds().B(),)`` where ``Ds`` is the category defining the
axiom ``B``.

This follows the same idiom as for deduction rules about functorial
constructions (see :meth:`.covariant_functorial_constructions.CovariantConstructionCategory.extra_super_categories`).
For example, the fact that a cartesian product of associative magmas
(i.e. of semigroups) is associative is implemented in
:meth:`Semigroups.Algebras.extra_super_categories`::

    sage: Magmas().Associative()
    Category of semigroups
    sage: Magmas().Associative().CartesianProducts().extra_super_categories()
    [Category of semigroups]

Similarly, the fact that the algebra of a commutative magma is
commutative is implemented in
:meth:`Magmas.Commutative.Algebras.extra_super_categories`::

    sage: Magmas().Commutative().Algebras(QQ).extra_super_categories()
    [Category of commutative magmas]

.. WARNING::

    In some situations this idiom is inapplicable as it would require
    to implement two classes for the same category. This is the
    purpose of the next section.

Special case
~~~~~~~~~~~~

In the previous examples, the deduction rule only had an influence on
the super categories of the category being constructed. For example,
when constructing ``Rings().Division()``, the rule
:meth:`Rings.Division.extra_super_categories` simply adds
``Rings().NoZeroDivisors()`` as a super category thereof.

.. TODO::

    rephrase that as: sometimes we do not want to create a class for
    the category with axiom just to specify some properties of this
    class. This is typically because a class for this category already
    exists elsewhere.

In some situations this idiom is inapplicable. Take for example
Wedderburn's theorem: any finite division ring is commutative, i.e. is
a finite field. The new feature of this situation is that this is not
just about the super categories of the category of finite division
rings. Instead it's stating that the category of finite division rings
*coincides* with that of finite fields::

        sage: DivisionRings().Finite()
        Category of finite fields

Therefore, we can't have two separate classes ``DivisionRings.Finite``
(with an ``extra_super_categories`` method as above) and
``Fields.Finite`` for this object.

A natural idiom would be to have ``DivisionRings.Finite`` be a
link to ``Fields.Finite`` (locally introducing a cycle in the tree
of nested classes). It would be a bit tricky to implement though,
since one would need to detect, upon constructing
``DivisionRings().Finite()``, that ``DivisionRings.Finite`` is
actually ``Fields.Finite``, in order to construct appropriately
``Fields().Finite()``; and reciprocally, upon computing the super
categories of ``Fields().Finite()``, to not try to add
``DivisionRings().Finite()`` as super category.

Instead the current idiom is to have a method
``DivisionRings.Finite_extra_super_categories`` which mimicks the
behavior of the would be
``DivisionRings.Finite.extra_super_categories``::

    sage: DivisionRings().Finite_extra_super_categories()
    (Category of commutative magmas,)

This idiom is admittedly rudimentary, but consistent with how
mathematical facts specifying non trivial inclusion relations between
categories are implemented elsewhere in the various
``extra_super_categories`` methods of axiom categories and covariant
functorial constructions. Besides, it gives a natural spot to document
and test the modeling of the mathematical fact.


In general, if several categories ``C1s(), C2s(), ... are mapped to
the same category when applying some axiom ``A`` (that is ``C1s().A()
== C2s().A() == ...``), then that category with axiom should be
implemented in ``Cs.A`` where ``Cs`` is the smallest of them, and the
other categories should implement a method
``A_extra_super_categories`` returning ``(Cs(),)``.

.. NOTE::

    An open question is whether there will always be a natural
    smallest category ``Cs``.

    Let's put a bit of formalism to explore the question. Let `A` be
    an axiom, and consider the application `\phi_A` which maps a
    category to its category of objects satisfying `A`. Equivalently,
    `\phi_A` is computing the intersection with the defining category
    with axiom of `A`. It follows immediately from the latter that
    `\phi_A` is a regressive endomorphism of the lattice of
    categories.

    The set `L` of categories that can be *constructed* in Sage is
    also a lattice (with joins = intersections that coincides with
    that of categories, but not meets). The application `\phi_A`
    restricts to a regressive endomorphism ``Cs() -> Cs().A()``.

    Consider ``C1s(), C2s(), ... `` as above. They form the
    intersection `S` of some fiber of `\phi_A` with the upper set
    `I_A` of categories that do not satisfy ``A``. The fiber itself is
    a sublattice. However `I_A` is not guaranteed to be stable under
    intersection (though exceptions should be rare). Therefore, there
    is a priori no guarantee to `S` would be stable under
    intersection. Also it's presumably finite, in fact small, but this
    is not guaranteed either.

    Altogether I (Nicolas) believe that this is not an issue in
    practice, but I would need to see a whole array of use cases to
    check that this is indeed not an issue in practice and, if at all
    possible, to specify a bit more the framework to guarantee the
    existence of ``Cs``.

Supporting such deduction rules will be an important feature in the
future, with quite a few occurences already implemented in upcoming
tickets. For the time being though there is a single occurence of this
idiom outside of the tests. So this would be an easy thing to
refactor after #10963 if a better idiom is found.


Specifications
==============

After fixing some vocabulary, we sumarize here some specifications
about categories and axioms.

The lattice of constructible categories
---------------------------------------

A mathematical category `C` is *implemented* if there is a class in
Sage modelling it; it is *constructible* if it is either implemented,
or is the intersection of *implemented* categories; in the later it is
modelled by a :class:`JoinCategory`. The comparison of two
constructible categories with the :meth:`Category.is_super_category`
method is supposed to model the comparison of the corresponding
mathematical categories for inclusion of the objects. For example::

    sage: Fields().is_subcategory(Rings())

However this modelling may be incomplete. It can happen that a
mathematical fact implying that a category `A` is a subcategory of a
category `B` is not implemented. Still, the comparison should endow
the set of constructible categories with a poset structure and in fact
a lattice structure.

In this lattice, the join of two categories (:meth:`Category.join`) is
supposed to model their intersection. Given that we compare categories
for inclusion, it would be more natural to call this operation the
*meet*; blames go to Nicolas who was originally comparing categories
by *amount of structure* rather than by *inclusion*. In practice, the
join of two categories may be a strict super category of their
intersection; first because this intersection might not be
constructible; second because Sage might miss some mathematical
information to recover the smallest constructible super category of
the intersection.

Axioms
------

We say that an axiom ``A`` is *defined by* a category ``Cs()`` if
``Cs`` defines an appropriate method ``Cs.SubcategoryMethods.A``, with
the semantic of the axiom specified in the documentation; for all
subcategory ``Ds()``, ``Ds().A()`` models the subcategory of the
objects of ``Ds()`` satisfying ``A``. In this case, we say that the
axiom ``A`` is *defined for* the category ``Ds()``. Furthermore,
``Ds`` *implements the axiom* ``A`` if ``Ds`` has a category with
axiom as nested class ``Ds.A``. The category ``Ds()`` *satisfies* the
axiom if ``Ds()`` is a subcategory of ``Cs().A()`` (meaning that all
the objects of ``Ds()`` are known to satisfy the axiom ``A``).

Specifications
--------------

- Any constructible category ``C`` should admit a finite number of
  larger constructible categories.

- The base category of a :class:`CategoryWithAxiom` should be an
  implemented category (i.e. not a :class:`JoinCategory`). This is
  checked by :meth:`CategoryWithAxiom._test_category_with_axiom`.

  ??? The base category of a :class:`CategoryWithAxiom` should be a
  direct super category of it.

- Tree structure: Let ``Cs()`` be a category, and `S` be some set of
  axioms defined in some super categories of ``Cs()`` but not
  satisfied by ``Cs()``. Suppose we want to provide a category with
  axiom for the elements of ``Cs()`` satisfying the axioms in
  `S`. Then, there should be a single enumeration ``A1, A2, ..., Ak``
  without repetition of the axioms in `S` such that
  ``Cs.A1.A2....Ak=Ds``. Furthermore, every intermediate step, like
  ``Cs.A1.A2....Ai`` with `i\leq k` should be a category with axiom
  having ``Ai`` as axiom and ``Cs.A1.A2.. Ai-1`` as base category
  class.

- If ``Cs()`` is a category that satisfies some axiom ``A`` (e.g. from
  one of its super categories), then it should not implement that
  axiom. In particular, a category class ``Cs`` can never have a
  nested class ``Cs.A.A``. Similarly, applying the specification
  recursively a category satisfying ``A`` cannot have a nested class
  ``Cs.A1.A2.A3.A`` where ``A1``, ``A2``, ``A3`` are axioms.

- A category can only implement an axiom if this axiom is defined by
  some super category. The code has not been systematically checked to
  support having two super categories defining the same axiom (which
  should of course have the same semantic). You are welcome to try, at
  your own risk :-)

- When a category defines an axiom or functorial construction ``A``,
  this fixes the semantic of ``A`` for all the subcategories. In
  particular, if two categories define ``A``, then those categories
  should be independent, and either the semantic of ``A`` should be
  the same, or there should be no natural intersection between the two
  hierarchies of subcategories.

- Any super category of a :class:`CategoryWithParameters` should
  either be a :class:`CategoryWithParameters` or a
  :class:`CategorySingleton`.

- A :class:`CategoryWithAxiom` having a :class:`Category_singleton` as
  base category should be a :class:`CategoryWithAxiom_singleton`. This
  is handled automatically by :meth:`CategoryWithAxiom.__init__` and
  checked in :meth:`CategoryWithAxiom._test_category_with_axiom`.

- A :class:`CategoryWithAxiom` having a
  :class:`Category_over_base_ring` as base category should be a
  :class:`Category_over_base_ring`. This currently has to be handled
  by hand, using :class:`CategoryWithAxiom_over_base_ring`. This is
  checked in :meth:`CategoryWithAxiom._test_category_with_axiom`.

.. TODO::

    The following specifications would be desirable but are not yet
    implemented:

    - A functorial construction category (Graded, CartesianProducts,
      ...) having a :class:`Category_singleton` as base category
      should be a :class:`CategoryWithAxiom_singleton`.

      Nothing difficult to implement, but this will need to rework the
      current "no subclass of a concrete class" assertion test of
      :meth:`Category_singleton.__classcall__`.

    - Similarly, a covariant functorial construction category having a
      :class:`Category_over_base_ring` as base category should be a
      :class:`Category_over_base_ring`.

    The following specification might be desirable, or not:

    - A join category involving a :class:`Category_over_base_ring`
      should a be a :class:`Category_over_base_ring`. In the mean
      time, a ``base_ring`` method is automatically provided for most
      of those by :meth:`Modules.SubcategoryMethods.base_ring`.


Design goals
============

As pointed out in the primer, the main design goal of the axioms
infrastructure is to subdue the potential combinatorial explosion of
the category hierarchy by letting the developper focus on implementing
a few bookshelves for which there is actual code or mathematical
information, and let Sage *compose dynamically and lazily* those
building blocks to construct the minimal hierarchy of classes needed
for the computation at hand. This to allows for the infrastructure to
scale smoothly as bookshelves are added, extended, or reorganized.

Other design goals include:

 - Flexibility in the code layout: the category of, say, finite sets
   can be implemented either within the Sets category (in a nested
   class ``Sets.Finite``), or in a separate file (typically in a class
   ``FiniteSets`` in a lazily imported module
   sage.categories.finite_sets).

 - Single point of truth: a theorem, like Wedderburn's, should be
   implemented in a single spot.

 - Single entry point: for example, from the entry :class:`Rings`, one
   can explore a whole range of related categories just by applying
   axioms and constructions::

       sage: Rings().Commutative().Finite().NoZeroDivisors()
       Category of finite integral domains
       sage: Rings().Finite().Division()
       Category of finite fields

   This will allow for progressively getting rid of all the
   :class:`GradedHopfAlgebrasWithBasis` which are polluting the global
   name space

 - Concise idioms for the users (adding axioms, ...)

 - Concise idioms and well highlighted hierarchy of bookshelves for
   the developper (especially with code folding)

 - Introspection friendly (listing the axioms, recovering the mixins)

.. NOTE::

    The constructor for instances of this class takes as input the
    base category. Hence, they should in principle be constructed
    as::

        sage: FiniteSets(Sets())
        Category of finite sets

        sage: Sets.Finite(Sets())
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

Upcoming features
=================

.. TODO:

    - Implement compatibility axiom / functorial constructions. For
      example, one would want to have::

          A.CartesianProducts() & B.CartesianProducts()) = (A&B).CartesianProducts()

    - Once full subcategories are implemented (see :trac:`10668`),
      make category with axioms be such. Should all full subcategories
      be implemented in term of axioms?

Description of the algorithmic
==============================

Computing joins
---------------

The workhorse of the axiom infrastructure is the algorithm for
computing the join `J` of a set `C_1,\dots,C_k` of categories (see
:meth:`Category.join`). Formally, `J` is defined as the largest
constructible category such that `J\subset C_i` for all `i`, and
`J\subset C.A()` for every constructible category `C\supset J` and any
axiom `A` satisfied by `J`.

The join `J` is naturally computed as a closure in the lattice of
constructible categories: it starts with the `C_i`'s, gathers the set
`S` of all the axioms satisfied by theim, and repeteadly adds each
axiom `A` to those categories that do not yet satisfy `A` using
:meth:`Category._with_axiom`. Due to deduction rules or (extra) super
categories, new categories or new axioms may appear in the
process. The process stops when each remaining categories has been
combined with each axiom. In practice, only the smallest categories
are kept along the way; this is correct because adding an axiom is
covariant: ``C.A()`` is a subcategory of ``D.A()`` whenever ``D`` is a
subcategory of ``C``.

As usual in such closure computations, the result does not depend on
the order of execution. Futhermore, given that adding an axiom is an
idempotent and regressive operation, the process is guaranteed to stop
in a number of steps which is bounded by the number of super
categories of `J`. In particular, it is a finite process.

.. TODO::

    Detail this a bit. What could typically go wrong is a situation
    where, for some category ``C1`, ``C1.A()`` specifies a category
    ``C2`` as super category such that ``C2.A()`` specifies ``C3`` as
    super category such that ...; this would clearly cause an infinite
    execution. Note that this situation violates the specifications
    since C1.A() is supposed to be a subcategory of C2.A(), ... so we
    would have an infinite increasing chain of constructible
    categories.

    It's reasonnable to assume that there is a finite number of axioms
    defined in the code. There remains to use this assumptiono to
    argue that any infinite execution of the algorithm would give rise
    to such an infinite sequence.

Adding an axiom
---------------

Let ``Cs`` be a category and ``A`` an axiom defined for this
category. To compute ``Cs().A()``, there are two cases.

Adding an axiom ``A`` to a category ``Cs()`` not implementing it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, ``Cs().A()`` returns the join of:

- ``Cs()``
- ``Bs().A()`` for every super category ``Bs()`` of ``Cs()``
- the categories appearing in ``Cs().A_extra_super_categories()``

This is a highly recursive process. In fact, as such, it would run
right away into an infinite loop! Indeed, the join of ``Cs()`` with
``Bs().A()`` would trigger the construction of ``Cs().A()`` and
reciprocally. To avoid this, the :meth:`Category.join` method itself
does not use :meth:`Category._with_axiom` to add axioms, but its
sister :meth:`Category._with_axiom_as_tuple`; the latter builds a
tuple of categories that should be joined together but leaves the
computation of the join to its caller the master join calculation.

Adding an axiom ``A`` to a category ``Cs()`` implementing it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case ``Cs().A()`` simply constructs an instance `D` of
``Cs.A`` which models the desired category. The non trivial part is
the construction of the the super categories of `D`. Very much like
above, this includes:

- ``Cs()``
- ``Bs().A()`` for every super category ``Bs()`` of ``Cs()``
- the categories appearing in ``D.extra_super_categories()``

This by itself may not be sufficient, due in particular to deduction
rules. On may for example discover a new axiom ``A1`` satisfied by
`D`, imposing to add ``A1`` to all of the above categories. Therefore
the super categories are computed as the join of the above categories.
Up to one twist: as is, the computation of this join would trigger
recursively a recalculation of ``Cs().A()``! To avoid this,
:meth:`Category.join` is given an optional argument to specify that
the axiom ``A`` should *not* be applied to ``Cs()``.

Examples where things go wrong if we are not careful
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As we have seen, this is a highly recursive process. In particular,
some appropriate steps need to be taken to not run in an infinite
recursion, in particular in case of deduction rule.

.. TODO:: Examples

Conclusion
==========

This is the end of the axioms documentation. Congratulations on
having read that far!


Tests
=====


TESTS:

.. NOTE::

    Quite a few categories with axioms are constructed early on during
    Sage's startup. Therefore, when playing around with the
    implementation of the axiom infrastructure, it is easy to break
    Sage. The following sequence of tests is designed to test the
    infrastructure from the ground up even in a partially broken
    Sage. Don't remove the imports!

::

    sage: from sage.categories.magmas import Magmas
    sage: Magmas()
    Category of magmas
    sage: Magmas().Finite()
    Category of finite magmas

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
    sage: Semigroups().Unital().Commutative()
    Category of commutative monoids
    sage: Semigroups().Commutative()
    Category of commutative semigroups
    sage: Semigroups().Commutative().Unital()
    Category of commutative monoids
    sage: Semigroups().Commutative().Unital().super_categories()
    [Category of monoids, Category of commutative magmas]

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

    sage: Domains().Commutative()
    Category of integral domains

    sage: DivisionRings().Finite() # Wedderburn's theorem
    Category of finite fields

    sage: FiniteMonoids().Algebras(QQ)
    Join of Category of finite dimensional algebras with basis over Rational Field and Category of monoid algebras over Rational Field and Category of finite set algebras over Rational Field
    sage: FiniteGroups().Algebras(QQ)
    Join of Category of finite dimensional hopf algebras with basis over Rational Field and Category of group algebras over Rational Field and Category of finite set algebras over Rational Field
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
from sage.misc.lazy_import import LazyImport
from sage.misc.misc import call_method
from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_types import Category_over_base_ring
from sage.structure.dynamic_class import DynamicMetaclass

# The order of the axioms in this lists implies that
# Magmas().Commutative().Unital() is printed as
# ``Category of commutative unital magmas''

all_axioms = ("Flying", "Blue",
              "Facade", "Finite", "Infinite",
              "FiniteDimensional", "Connected", "WithBasis",
              "Irreducible",
              "Commutative", "Associative", "Inverse", "Unital", "Division", "NoZeroDivisors",
              "AdditiveCommutative", "AdditiveAssociative", "AdditiveInverse", "AdditiveUnital",
              )

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
    Try to guess the base category and the axiom from the name of ``cls``.

    The heuristic is to try to decompose the name as the concatenation
    of the name of a category and the name of an axiom, and looking up
    that category in the standard location (i.e. in
    :mod:`sage.categories.hopf_algebras` for :class:`HopfAlgebras`,
    and in :mod:`sage.categories.sets_cat` as a special case for :class:`Sets`).

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import base_category_class_and_axiom, CategoryWithAxiom
        sage: base_category_class_and_axiom(FiniteSets)
        (<class 'sage.categories.sets_cat.Sets'>, 'Finite')
        sage: Sets.Finite
        <class 'sage.categories.finite_sets.FiniteSets'>
        sage: base_category_class_and_axiom(Sets.Finite)
        (<class 'sage.categories.sets_cat.Sets'>, 'Finite')

        sage: base_category_class_and_axiom(FiniteDimensionalHopfAlgebrasWithBasis)
        (<class 'sage.categories.hopf_algebras_with_basis.HopfAlgebrasWithBasis'>, 'FiniteDimensional')

        sage: base_category_class_and_axiom(HopfAlgebrasWithBasis)
        (<class 'sage.categories.hopf_algebras.HopfAlgebras'>, 'WithBasis')

    Along the way, this does some sanity checks::

        sage: class FacadeSemigroups(CategoryWithAxiom):
        ....:     pass
        sage: base_category_class_and_axiom(FacadeSemigroups)
        Traceback (most recent call last):
        ...
        AssertionError: Missing (lazy import) link for <class 'sage.categories.semigroups.Semigroups'> to <class '__main__.FacadeSemigroups'> for axiom Facade?

        sage: Semigroups.Facade = FacadeSemigroups
        sage: base_category_class_and_axiom(FacadeSemigroups)
        (<class 'sage.categories.semigroups.Semigroups'>, 'Facade')

    .. NOTE::

        In the following example, we could possibly retrieve ``Sets``
        from the class name. However this cannot be implemented
        robustly until #9107 is fixed. Anyway this feature has not
        been needed so far::

            sage: Sets.Infinite
            <class 'sage.categories.sets_cat.Sets.Infinite'>
            sage: base_category_class_and_axiom(Sets.Infinite)
            Traceback (most recent call last):
            ...
            TypeError: Could not retrieve the base category class and axiom for <class 'sage.categories.sets_cat.Sets.Infinite'>
            ...
    """
    if "." in cls.__name__:
        # Case 1: class name of the form Sets.Infinite
        # Start of implementation when #9107 will be fixed:
        # axiom = cls.__name__.split(".")[-1]
        # ...
        pass
    else:
        # Case 2: class name of the form FiniteSets or AlgebrasWithBasis,
        # with the base class (say Algebras) being implemented in the
        # standard location (sage.categories.algebras)
        name = cls.__name__
        for axiom in all_axioms:
            if axiom == "WithBasis" and name.endswith(axiom):
                base_name = name[:-len(axiom)]
            elif name.startswith(axiom):
                base_name = name[len(axiom):]
            else:
                continue
            if base_name == "Sets": # Special case for Sets which is in sets_cat
                base_module_name = "sets_cat"
            else:
                base_module_name = uncamelcase(base_name, "_")
            try:
                base_module = importlib.import_module("sage.categories."+base_module_name)
                base_category_class = getattr(base_module, base_name)
                assert getattr(base_category_class, axiom, None) is cls, \
                    "Missing (lazy import) link for %s to %s for axiom %s?"%(base_category_class, cls, axiom)
                return base_category_class, axiom
            except (ImportError,AttributeError):
                pass
    raise TypeError("Could not retrieve the base category class and axiom for %s\nPlease specify it explictly using the attribute _base_category_class_and_axiom.\nSee CategoryWithAxiom for details."%cls)

@cached_function
def axiom_of_nested_class(cls, nested_cls):
    """
    Given a class and a nested axiom class, return the axiom.

    EXAMPLES:

    This uses some heuristics like checking if the nested_cls carries
    the name of the axiom, or is built by appending or prepending the
    name of the axiom to that of the class::

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
    ``_base_category_class_and_axiom``::

        sage: Semigroups._base_category_class_and_axiom
        (<class 'sage.categories.magmas.Magmas'>, 'Associative')
        sage: axiom_of_nested_class(Magmas, Semigroups)
        'Associative'
    """
    try:
        axiom = nested_cls.__dict__["_base_category_class_and_axiom"][1]
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

    See the :mod:`category primer <sage.categories.primer>`, and in
    particular its :ref:`section about axioms <category-primer-axioms>`
    for an introduction to axioms, and :class:`CategoryWithAxiom` for
    how to implement axioms and the documentation of the axiom
    infrastructure.
    """

    @lazy_class_attribute
    def _base_category_class_and_axiom(cls):
        r"""
        The class of the base category and the axiom for this class.

        By default, this attribute is guessed from the name of this
        class (see :func:`base_category_class_and_axiom`). For a
        nested class, when the category is first created from its base
        category, as in e.g. ``Sets().Infinite()``, this attribute is
        instead set explicitly by :meth:``__classget__``.

        When this is not sufficient, that is when ``cls`` is not
        implemented as a nested class and the base category and the
        axiom cannot be guessed from the name of ``cls``, this
        attribute should be set explicitly by ``cls``.

        The origin of the attribute is stored in the attribute
        ``_base_category_class_and_axiom_origin``.

        .. SEEALSO:: :meth:`_axiom`

        EXAMPLES:

        ``CommutativeRings`` is not a nested class, but the name of
        the base category and the axiom can be guessed::

            sage: CommutativeRings()._base_category_class_and_axiom
            (<class 'sage.categories.rings.Rings'>, 'Commutative')
            sage: CommutativeRings()._base_category_class_and_axiom_origin
            'guessed by base_category_class_and_axiom'

        ``Sets.Infinite`` is a nested class, so the attribute is set
        by :meth:`CategoryWithAxiom.__classget__` the first time
        ``Sets().Infinite()`` is called::

            sage: Sets().Infinite()
            Category of infinite sets
            sage: Sets.Infinite._base_category_class_and_axiom
            (<class 'sage.categories.sets_cat.Sets'>, 'Infinite')
            sage: Sets.Infinite._base_category_class_and_axiom_origin
            'set by __classget__'

        ``Fields`` is not a nested class, and the name of the base
        category and axioms cannot be guessed; so this attributes
        needs to be set explicitly in the ``Fields`` class::

            sage: Fields()._base_category_class_and_axiom
            (<class 'sage.categories.division_rings.DivisionRings'>, 'Commutative')
            sage: Fields()._base_category_class_and_axiom_origin
            'hardcoded'

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
        cls._base_category_class_and_axiom_origin = "guessed by base_category_class_and_axiom"
        return (base_category_class, axiom)

    _base_category_class_and_axiom_origin = "hardcoded"

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
            cls._base_category_class_and_axiom_origin = "set by __classget__"
        else:
            assert cls._base_category_class_and_axiom[0] is base_category_class, \
                "base category class for %s mismatch; expected %s, got %s"%(cls, cls._base_category_class_and_axiom[0], base_category_class)

        # Workaround #15648: if Rings.Finite is a LazyImport object,
        # this forces the substitution of the object back into Rings
        # to avoid resolving the lazy import over and over
        if isinstance(base_category_class.__dict__[cls._axiom], LazyImport):
            setattr(base_category_class, cls._axiom, cls)

        if base_category is None:
             return cls
        # For Rings().Finite, this returns the method
        # Sets.SubcategoryMethods.Finite, with its first argument bound to Rings()
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
        Run generic tests on this category with axioms.

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
        Return the extra super categories of a category with axiom.

        Default implementation which returns ``[]``

        EXAMPLES::

            sage: FiniteSets().extra_super_categories()
            []
        """
        return []

    @cached_method
    def super_categories(self):
        """
        Return a list of the (immediate) super categories of
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
        return Category.join((base_category,) +
                             tuple(cat
                                   for category in base_category._super_categories
                                   for cat in category._with_axiom_as_tuple(axiom)) +
                             tuple(self.extra_super_categories()),
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

        .. NOTE::

            The logic here is shared between :meth:`_repr_object_names`
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
        Return the base category of ``self``.

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
        Return the category without the axioms that have been added to create it.

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

        Technically we tests this by checking if the class specifies
        explicitly the attribute ``_base_category_class_and_axiom``
        by looking up ``_base_category_class_and_axiom_origin``.

        Some more examples::

            sage: Algebras(QQ).Commutative()._without_axioms()
            Category of magmatic algebras over Rational Field
            sage: Algebras(QQ).Commutative()._without_axioms(named=True)
            Category of algebras over Rational Field
        """
        if named:
            base_category = self
            axioms = []
            while isinstance(base_category, CategoryWithAxiom) and base_category._base_category_class_and_axiom_origin != "hardcoded":
                axioms.append(base_category._axiom)
                base_category = base_category._base_category
            return base_category
        else:
            return self._base_category._without_axioms()

    @cached_method
    def axioms(self):
        r"""
        Return the axioms known to be satisfied by all the objects of ``self``.

        .. SEEALSO:: :meth:`Category.axioms`

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
                         for category in self._super_categories
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

class Blahs(Category_singleton):

    def super_categories(self):
        """
        TESTS::

             sage: from sage.categories.category_with_axiom import Blahs
             sage: Blahs().super_categories()
             [Category of sets]
            sage: TestSuite(Blahs()).run()
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

            Here, we want ``Flying`` to imply ``Unital``, and to put
            the class for the category of unital flying blahs in
            ``Blahs.Flying`` rather than ``Blahs.Unital.Flying``.

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
        Illustrates a current limitation in the way to have an axiom imply another one.

        Here, we would want ``Blue`` to imply ``Unital``, and to put
        the class for the category of unital blue blahs in
        ``Blahs.Unital.Blue`` rather than ``Blahs.Blue``.

        This currently fails because ``Blahs`` is the category where
        the axiom ``Blue`` is defined, and the specifications
        currently impose that a category defining an axiom should also
        implement it (here in an category with axiom
        ``Blahs.Blue``). In practice, due to this violation of the
        specifications, the axiom is lost during the join calculation.

        .. TODO::

            Decide whether we care about this feature. In such a
            situation, we are not really defining a new axiom, but
            just defining an axiom as an alias for a couple others,
            which might not be that useful.

        .. TODO::

            Improve the infrastructure to detect and report this
            violation of the specifications, if this is
            easy. Otherwise, it's not so bad: when defining an axiom A
            in a category ``Cs`` the first thing one is supposed to
            doctest is that ``Cs().A()`` works. So the problem should
            not go unnoticed.

        TESTS::

            sage: from sage.categories.category_with_axiom import Blahs, TestObjects, Bars
            sage: Blahs().Blue_extra_super_categories()
            [Category of unital blahs]
            sage: Blahs().Blue()                          # todo: not implemented
            Category of blue unital blahs
        """
        return [Blahs().Unital()]

class Bars(Category_singleton):
    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.category_with_axiom import Bars
            sage: Bars().super_categories()
            [Category of blahs]
            sage: TestSuite(Bars()).run()
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

class TestObjects(Category_singleton):

    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects
            sage: TestObjects().super_categories()
            [Category of bars]
            sage: TestSuite(TestObjects()).run()
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
            sage: TestObjectsOverBaseRing.Unital.an_instance()
            Category of unital test objects over base ring over Rational Field
            sage: TestObjectsOverBaseRing.FiniteDimensional.Unital.an_instance()
            Category of finite dimensional unital test objects over base ring over Rational Field
            sage: TestSuite(TestObjectsOverBaseRing(QQ).FiniteDimensional().Unital().Commutative()).run()
        """
        return [TestObjects()]

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
         class Finite(CategoryWithAxiom_over_base_ring):
              pass
         class Unital(CategoryWithAxiom_over_base_ring):
              class Commutative(CategoryWithAxiom_over_base_ring):
                   pass

    class Commutative(CategoryWithAxiom_over_base_ring):
         class Facade(CategoryWithAxiom_over_base_ring):
             pass
         class FiniteDimensional(CategoryWithAxiom_over_base_ring):
             pass
         class Finite(CategoryWithAxiom_over_base_ring):
             pass

    class Unital(CategoryWithAxiom_over_base_ring):
        pass

class BrokenTestObjects(Category):
    class Commutative(CategoryWithAxiom):
         class Finite(CategoryWithAxiom):
              pass
    class Finite(CategoryWithAxiom):
         class Commutative(CategoryWithAxiom):
              pass

