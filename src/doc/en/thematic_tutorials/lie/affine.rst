-------------------------------
Affine Root Systems and methods
-------------------------------

Among infinite-dimensional Lie algebras, *Kac-Moody Lie algebras*
are generalizations of finite-dimensional semisimple Lie algebras.
Among Kac-Moody Lie algebras, *affine Lie algebras* are an important
infinite-dimensional class, and their infinite-dimensional
integrable highest-weight representations are an important class.

One realization of affine Lie algebras, described in Chapter 7
of [Kac]_ begins with a
finite-dimensional semisimple Lie algebra `\mathfrak{g}^\circ`,
with Cartan type `['X',\ell]`. Tensoring with the Laurent polynomial ring gives
the loop Lie algebra `\mathfrak{g}^\circ\otimes\CC[t,t^{-1}]`. This
is the Lie algebra of vector fields in `\mathfrak{g}` on
the circle. Then one may make a central extension:

.. MATH::

   0 \rightarrow \CC\cdot K\rightarrow \widehat{\mathfrak{g}}'
   \rightarrow\mathfrak{g}^\circ\otimes\CC[t,t^{-1}]\rightarrow 0.

After that, it is convenient to adjoin another basis element,
which acts on `\mathfrak{g}'` as a derivation `d`. If `\mathfrak{h}^\circ`
is a Cartan subalgebra of `\mathfrak{g}^\circ` then we obtain a
Cartan subalgebra `\mathfrak{h}'` of `\mathfrak{g}'` by adjoining
the central element `K`, and then a Cartan subalgebra `\mathfrak{h}`
by further adjoining the derivation `d`.

Concepts from the representation theory of finite-dimensional
Lie groups extend to affine (and more general Kac-Moody) Lie algebras.
This includes the root system, Weyl group, weight lattice and
the parametrization of important representations by dominant weights.
All of these objects have realizations in Sage.

The resulting Lie algebra `\widehat{\mathfrak{g}}` is the *untwisted affine
Lie algebra*.  The Cartan type is `['X',\ell,1]`, which
we can abbreviate as `"Xl~".The Dynkin diagram of this
Cartan type is the extended Dykin-diagram of `\mathfrak{g}`::

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

Usually there is an advantage to working with `\mathfrak{g}` instead of
`\mathfrak{g}'`. If `V` is a representation of `\mathfrak{g}` then
usually the weight spaces `V_\lambda`, in a decomposition
with respect to characters (weights) of `\mathfrak{h}` are
finite-dimensional; but the corresponding weight spaces for
`\mathfrak{h}'` would not be. There are exceptions in certain
finite-dimensional representions of `\mathfrak{g}'`
cannot be extended to `\mathfrak{g}` (though they do have
infinite-dimensional analogs). These representations have
crystal bases, including the Kirillov-Reshetikhin crystals.
See :ref:`AffineFinite`.

Twisted Types
-------------

There are also *twisted* types `['X',\ell,r]` where `r` is the order of an
automorphism of the Dynkin diagram. These are described in [Kac]_ Chapter 8.
Alternative descriptions of the twisted types may be found in
[Macdonald2003]_. Examining the tables Aff1, Aff2 and Aff3
in Chapter 4 of Kac, you will see that each twisted type is
dual to an untwisted type. For example the twisted type `['E',6,2]` in Aff2 is
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

The Cartan types `['E',6,2]` and `['F',4,1]` (abbreviated `F4~`) are dual
in the sense that long roots of one correspond to short roots of the other.
(Thus 0,1 and 2 are short roots of `['E',6,2]`, they are long roots of
`['F',4,1]`.) More generally, every twisted affine type is dual to a
unique untwisted type, and the Macdonald convention is to refer to
the Cartan type as the dual of the corresponding untwisted type::

    sage: CartanType(['F',4,1]).dual()==CartanType(['E',6,2])
    True

Affine Root System and Weyl Group
---------------------------------

The basic reference for the affine root system and Weyl group
is [Kac]_ Chapter 6. In particular the reader who wishes to

As with finite-dimensional reductive Lie algebras, described in
:ref:`LieBasics`, the root system of `\mathfrak{g}` is a discrete
(but now infinite) subset of `\mathfrak{h}^\ast`, and the Weyl
group acts on `\mathfrak{h}` permuting the roots. Now however
the root system is infinite. Moreover, it is partitioned into
*real roots* and *imaginary roots*. Let `\Delta` be the root
system, `\Delta_{\text{re}}` the set of real roots and
`\Delta_{\text{im}}` the set of imaginary roots. Then we
may decompose `\mathfrak{g}` into root eigenspaces `\mathfrak{g}_\alpha`
with `\alpha\in\Delta`.

.. MATH::

    \mathfrak{g} = \mathfrak{h} \oplus \bigoplus_{\alpha\in\Delta} \mathfrak{g}_\alpha.

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
    sage: V.simple_roots()
    Finite family {0: -e[0] + e[2] + e['delta'], 1: e[0] - e[1], 2: e[1] - e[2]}
    sage: V.basic_imaginary_roots()[0]
    e['delta']

The Weyl Group and extended Affine Weyl Group
---------------------------------------------

The ambient space of the root system comes with an
inner product.


Weight Lattice
--------------

The rank of the weight lattice of `\widehat{\mathfrak{g}}` is larger
by 2 than the weight lattice of `\mathfrak{g}`. It contains
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

Integrable Highest Weight Representations
-----------------------------------------