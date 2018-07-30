---------------------------------------------------------
Affine Root Systems and Integrable Highest Weight modules
---------------------------------------------------------

Among infinite-dimensional Lie algebras, *Kac-Moody Lie algebras*
are generalizations of finite-dimensional simple Lie algebras.
They include finite-dimensional simple Lie algebras as a special
case but are usually infinite-dimensional.  many concepts and
results from the representation theory of finite-dimensional Lie groups
and Lie algebras extend to Kac-Moody Lie algebras.  This includes the root
system, Weyl group, weight lattice, the parametrization of an important
representations (the integrable highest weight ones) by dominant weights
and the Weyl character formula for these representations.

Among Kac-Moody Lie algebras, *affine Lie algebras* are an important
infinite-dimensional class, and their infinite-dimensional
integrable highest-weight representations are an important class
among their representations.  In this section many of the concepts are
applicable to general Kac-Moody Lie algebras. But the code that we will
discuss is is mainly for the affine case. This is sufficiently different
from the general case (and important) to merit special attention.

In this section we will review some of the Kac-Moody theory,
taking [Kac]_ as our primary reference. We will also
explain what facilities there are in Sage for computing
with these. We will often restrict ourselves to the case
of affine Lie algebras.

One realization of affine Lie algebras, described in Chapter 7
of [Kac]_ begins with a
finite-dimensional isimple Lie algebra `\mathfrak{g}^\circ`,
with Cartan type `X_\ell` (``['X',l]`` in Sage). Tensoring with the
Laurent polynomial ring gives the loop Lie algebra
`\mathfrak{g}^\circ\otimes\CC[t,t^{-1}]`. This is the Lie algebra of
vector fields in `\mathfrak{g}^\circ` on the circle. Then one may make a
central extension:

.. MATH::

   0 \rightarrow \CC\cdot K\rightarrow {\mathfrak{g}}'
   \rightarrow\mathfrak{g}^\circ\otimes\CC[t,t^{-1}]\rightarrow 0.

After that, it is convenient to adjoin another basis element,
which acts on `\mathfrak{g}'` as a derivation `d`. If `\mathfrak{h}^\circ`
is a Cartan subalgebra of `\mathfrak{g}^\circ` then we obtain a
Cartan subalgebra `\mathfrak{h}'` of `\mathfrak{g}'` by adjoining
the central element `K`, and then a Cartan subalgebra `\mathfrak{h}`
by further adjoining the derivation `d`.

The resulting Lie algebra `{\mathfrak{g}}` is the *untwisted affine
Lie algebra*.  The Cartan type is designated to be `X_\ell^{(1)}`
in Kac' notation, which is rendered as ``['X',l,1]`` or ``"Xl~"``
in Sage. The Dynkin diagram of this
Cartan type is the extended Dykin-diagram of `\mathfrak{g}^\circ`::

    sage: CartanType("E6~").dynkin_diagram()

                  O 0
                  |
                  |
                  O 2
                  |
                  |
          O---O---O---O---O
          1   3   4   5   6

From the Dynkin diagram, we can read off generators and relations
for the affine Weyl group, which is a Coxeter group with generators
`s_i` of order 2, that commute if `i` and `j` are not adjacent in
the Dynkin diagram, and othrwise are subject to a braid relation.
We can infer some Levi subalgebras of `\mathfrak{g}`, obtained by
omitting one node from the Dynkin diagram; particularly omitting
the "affine node" `0` gives `E_6`, that is `\mathfrak{g}^\circ`.

There are two versions of the weight lattice, depending on
whether we are working with `\mathfrak{g}` or `\mathfrak{g}'`.
The larger Lie algebra of `\mathfrak{g}` is called the
*extended* weight lattice. We may create these as follows::

    sage: RootSystem("A2~").weight_lattice()
    Weight lattice of the Root system of type ['A', 2, 1]
    sage: RootSystem("A2~").weight_lattice(extended=True)
    Extended weight lattice of the Root system of type ['A', 2, 1]

Usually there is an advantage to working with `\mathfrak{g}` instead of
`\mathfrak{g}'`. (Thus we prefer the extended weight lattice,
though this is not the default.) The reason for this is as 
follows. If `V` is a representation of `\mathfrak{g}` then
usually the weight spaces `V_\lambda`, in a decomposition
with respect to characters (weights) of `\mathfrak{h}` are
finite-dimensional; but the corresponding weight spaces for
`\mathfrak{h}'` would not be.

There are exceptions to this rule of preferring the extended
weight lattice in certain finite-dimensional representions of
`\mathfrak{g}'` that cannot be extended to `\mathfrak{g}` (although they
do have infinite-dimensional analogs). These finite-dimensional
representations have crystal bases, including the Kirillov-Reshetikhin
crystals. Thus for Kirillov-Reshetikhin crystals we prefer tue See
:ref:`AffineFinite`.

Twisted Types
-------------

There are also *twisted* types with Cartan type `X_\ell^{(m)}` or
``['X',l,m]`` where `m` is the order of an
automorphism of the Dynkin diagram of `\mathfrak{g}^\circ`. These are
described in [Kac]_ Chapter 8.  Alternative descriptions of the twisted
types may be found in [Macdonald2003]_. Examining the tables Aff1, Aff2
and Aff3 in Chapter 4 of Kac, you will see that each twisted type is dual
to an untwisted type. For example the twisted type `['E',6,2]` in Aff2 is
dual to the untwisted type `['F',4,1]`.

Referring to the above Dynkin diagram for `['E',6,1]`, if
we collapse the nodes 1 and 6 together, and the nodes 3 and 5,
we obtain the Dynkin diagram for `['E',6,2]`::

     sage: CartanType(['E',6,2]).dynkin_diagram()
     O---O---O=<=O---O
     0   1   2   3   4
     F4~*

We must explain why Sage calls this Cartan type `F4~*`.
The Cartan type `['F',4,1]` is obtained by adding one
Dynkin node to the Cartan type "F4"::

    sage: CartanType(['F',4,1]).dynkin_diagram()
    O---O---O=>=O---O
    0   1   2   3   4
    F4~

The Cartan types `['E',6,2]` and `['F',4,1]` (abbreviated ``F4~``) are dual
in the sense that long roots of one correspond to short roots of the other.
(Thus 0,1 and 2 are short roots of `['E',6,2]`, they are long roots of
`['F',4,1]`.) More generally, every twisted affine type is dual to a
unique untwisted type, and the Macdonald convention is to refer to
the Cartan type as the dual of the corresponding untwisted type::

    sage: CartanType(['F',4,1]).dual()==CartanType(['E',6,2])
    True

.. _roots_and_weights:

Roots and Weights
-----------------

The Lie algebra `\mathfrak{g}` has a triangular decomposition

.. MATH::

    \mathfrak{g} = \mathfrak{h} \oplus \mathfrak{n}_+ \oplus \mathfrak{n}_-

where `\mathfrak{n}_-` and `\mathfrak{n}_+` are nilpotent Lie algebras.

If `V` is a `\mathfrak{g}`-module then we often have
a *weight space decomposition*

.. MATH::

    V = \bigoplus_{\lambda\in\mathfrak{h}^*} V_\lambda

where `V_\lambda` is finite-dimensional, and where `\mathfrak{h}`
acts by `X\,v=\lambda(X)v` for `X\in\mathfrak{h}`, `v\in V_\lambda`.
The linear functional `\lambda` is called a *weight*.
The space `V_\lambda` is called the *weight space* and its
dimension is the *multiplicity* of the weight `\lambda`.

As a special case, `\mathfrak{g}` is a module over itself
under the adjoint representation, and it has a weight
decomposition.

The nonzero weights in the adjoint representation of `\mathcal{g}`
on itself are called *roots*. In contrast with the finite-dimensional
case, if `\mathcal{g}` is an infinte Kac-Moody Lie algebra there are two
types of roots, called *real* and imaginary. The real roots have
multiplicity 1, while the imaginary roots can have multiplicity
`>1`. In the case of the affine Kac-Moody Lie algebra the
imaginary roots have bounded multiplicity, while in non-affine
cases the multiplicities of the imaginary roots is somewhat
mysterious.

The roots may be divided into those in the adjoint
representation of `\mathfrak{h}` on `\mathbf{n}_+`,
called *positive*, and those for `\mathbf{n}_-`,
called *negative*. 

Returning to the general module `V` with a weight space
decomposition, a vector in the module `V` that is annihilated by
`\mathfrak{n}_+` is called a *highest weight vector*. If the space of
highest weight vectors is one-dimensional, and if `V` is generated by a
highest weight vector `v` then `\CC\,v=V_\lambda` for a weight
`\lambda`, called the *highest weight*, and `v` is called a *highest weight vector*.

If `\lambda` is any linear functional on `\mathfrak{h}` then there
is a *universal highest weight module* `M(\lambda)` such that any
highest weight module with highest weight `\lambda` is a quotient
of `M(\lambda)`. In particular `M(\lambda)` (which is also called
a *Verma module*) has a unique irreducible quotient denoted `L(\lambda)`.
Looking ahead to crystal bases, the infinity crystal `\mathcal{B}(\infty)`
is a crystal base of the Verma module `M(0)`.

Affine Root System and Weyl Group
---------------------------------

We now specialize to affine Kac-Moody Lie algebras and their
root systems. The basic reference for the affine root system and Weyl
group is [Kac]_ Chapter 6.

There is a minimal imaginary root `\delta`. The imaginary roots
are the vectors `n\delta` where `n` is a nonzero integer. This
root is positive if and only if `n>0`. The root system `\Delta`
contains a copy of the finite root system `Delta^\circ` of
`\mathfrak{g}^\circ`. In the untwisted case, the real roots
are `\alpha+n\delta` where `n` is an integer; the root is
positive if `n>0` or if `n=0` and `\alpha` is positive. For
a description of the real roots in the twisted case, see
[Kac]_ Proposition 6.3.

The multiplicity `m(\alpha)` is the dimension of `\mathfrak{g}_\alpha`.
It is 1 if `\alpha` is a real root. For the affine Lie algebras
that concern us now, the multiplicity of an imaginary root is
the rank `\ell` of `\mathfrak{g}^\circ`.

In Sage, many important things such as the roots, and Weyl group and are methods
of the ambient space::

    sage: V=RootSystem(['A',2,1]).ambient_space()
    sage: V.positive_roots()
    Disjoint union of Family (Positive real roots of type ['A', 2, 1], Positive imaginary roots of type ['A', 2, 1])
    sage: V.simple_roots()
    Finite family {0: -e[0] + e[2] + e['delta'], 1: e[0] - e[1], 2: e[1] - e[2]}
    sage: V.weyl_group()
    Weyl Group of type ['A', 2, 1] (as a matrix group acting on the ambient space)
    sage: V.basic_imaginary_roots()[0]
    e['delta']

However it may be better for weights to have their parents to be
the weight lattice instead of its ambient vector space. Therefore
you may create all of the above vectors as parents of the weight
lattice. In this case we recommend creating the lattice with the
option ``extended=True``::

    sage: WL = RootSystem(['A',2,1]).weight_lattice(extended=True); WL
    Extended weight lattice of the Root system of type ['A', 2, 1]
    sage: WL.positive_roots()
    Disjoint union of Family (Positive real roots of type ['A', 2, 1], Positive imaginary roots of type ['A', 2, 1])
    sage: WL.simple_roots()
    Finite family {0: 2*Lambda[0] - Lambda[1] - Lambda[2] + delta, 1: -Lambda[0] + 2*Lambda[1] - Lambda[2], 2: -Lambda[0] - Lambda[1] + 2*Lambda[2]}
    sage: WL.weyl_group()
    Weyl Group of type ['A', 2, 1] (as a matrix group acting on the extended weight lattice)
    sage: WL.basic_imaginary_roots()[0]
    delta

Certain constants `a_i` label the vertices `i=0,\cdots,\ell` in
the tables Aff1, Aff2 and Aff3 in [Kac]_ Chapter 4. They
play an important role in the theory. In Sage they are available
as follows::

    sage: CartanType(['B',5,1]).a()
    Finite family {0: 1, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2}

Be aware that for the exceptional groups, the ordering of the indices
are different from those in [Kac]_. This is because Sage uses the Bourbaki
ordering of the roots, and Kac does not. Thus in Bourbaki (and in Sage)
the `G_2` short root is `\alpha_1`::

    sage: CartanType(['G',2,1]).dynkin_diagram()
      3
    O=<=O---O
    1   2   0
    G2~
  
By contrast in Kac, `\alpha_2` is the short root.

The Weyl Group and extended Affine Weyl Group
---------------------------------------------

The ambient space of the root system comes with an
(indefinite) inner product. The real roots have
nonzero length but the imaginary roots are isotropic.
If `\alpha` is a real root we may define a reflection `r_\alpha`
in the hyperplane orthogonal to `\alpha`. In particular
the `\ell+1` reflections `s_i` with respect to the *simple positive roots*
`\alpha_i` (`i=0,1,2,\cdots,\ell`) generate a Coxeter group.
This is the *Weyl group* `W`.

The subgroup `W^\circ` generated by `s_1,\cdots,s_\ell`
is a finite Coxeter group that may be identified with
the Weyl group of the finite-dimensional simple
Lie algebra `\mathfrak{g}^\circ`.

Geometrically, `W` may be interpreted as the semidirect product
of the finite Weyl group `W^\circ` by a discrete group of
translations `Q^\vee`; this group is isomorphic to the coroot
lattice. A larger *extended affine Weyl group* is the semidirect
product of `W^\circ` by the coweight lattice `P^\vee`. If
`P^\vee` is strictly larger than `Q^\vee` this
is not a Coxeter group but arises naturally in many problems.
It may be constructed in Sage as follows::

    sage: E = ExtendedAffineWeylGroup(["A",2,1]); E
    Extended affine Weyl group of type ['A', 2, 1]

See the documentation in
:file:`~sage.combinat.root_system.extended_affine_weyl_group` if you need this.

Weight Lattice
--------------

The rank of the weight lattice of `{\mathfrak{g}}` is larger
by 2 than the weight lattice of `\mathfrak{g}^\circ`. It contains
fundamental weights `\Lambda_1,\cdots,\Lambda_l`
corresponding to the fundamental weights of `\mathfrak{g}`
and one more, the *affine* fundamental weight `\Lambda_0`.

A finite linear combination with nonnegative integer
coefficients of `\Lambda_0,\cdots,\Lambda_l` is a
*dominant weight*.

If `\Lambda` is a dominant weight then `\mathfrak{g}` has
an infinite-dimensional irreducible representation with highest
weight `\Lambda`. We can study these using the
``IntegrableRepresentation`` class of Sage.

Now there is a distinction between the weight lattices
of `\mathfrak{g}` and `\mathfrak{g}'`.

Integrable Highest Weight Representations
-----------------------------------------

In this section `\mathfrak{g}` can be an arbitrary
Kac-Moody Lie Algebra.

Suppose that `V` is a representation with a weight
decomposition as in :ref:`roots_and_weights`.
Let `\alpha` be a real root, and let `\mathfrak{g}_\alpha`
be the corresponding weight space, called a *root space*.
Then `-\alpha` is also a root. The two
one-dimensional spaces `\mathfrak{g}_\alpha` and
`\mathfrak{g}_{-\alpha}` generate a Lie algebra
isomorphic to `\mathfrak{sl}_2`. The module `V`
is called *integrable* if for each such `\alpha`
the representation of `\mathfrak{sl}_2` obtained this
way integrates to a representation of the Lie group
`\text{SL}_2`.

If the Kac-Moody Lie algebra `\mathfrak{g}` is finite-dimensional
then the integrable highest weight representations are
just the irreducible finite-dimensional ones. For a general
Kac-Moody Lie algebra the integrable highest weight representations
are the analogs of the finite-dimensional ones,
discussed in :file:`weyl_character_ring`, and their
theory has many aspects in common with the finite-dimensional
representations of finite-dimensional simple Lie algebras,
such as the parametrization by dominant weights, and
generalizations of the Weyl denominator and character
formulas, due to Macdonald and Kac respectively.

Within `\mathfrak{h}^\ast` there is a lattice `\Lambda`,
called the *weight* lattice such that if `V` is an
integrable highest weight representation, the weights
in the weight space decomposition (:ref:`roots_and_weights`)
are in `\Lambda`. Moreover, there exists a cone `\Lambda^+`
of *dominant weights* such that `\lambda\in\mathfrak{h}^\ast` is
the highest weight of a (unique) integrable highest
weight module if and only if `\lambda\in\Lambda^+`. See
[Kac]_ Chapters 9 and 10 for the theory of integrable
highest weight representations.

There exists a basis `\Lambda_i` of the lattice `\Lambda` such that the
dominant weights are the nonnegative linear combinations of the
`\Lambda_i`. In the affine cases we label the weights
`i=0,1,\cdots,r-1`. If `\mathfrak{g}` is the untwisted affine Lie
algebra of Cartan type `X_\ell^{(1)}`` then `r=\ell+1`. The labels correspond
to the nodes in the Dynkin diagram.

:class:`~sage.combinat.root_system.integrable_representations.IntegrableRepresentation`