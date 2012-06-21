-------------------------------------
Maximal Subgroups and Branching Rules
-------------------------------------


Branching rules
---------------

If `G` is a Lie group and `H` is a subgroup, one often needs to know
how representations of `G` restrict to `H`. Irreducibles usually do
not restrict to irreducibles. In some cases the restriction is regular
and predictable, in other cases it is chaotic. In some cases it
follows a rule that can be described combinatorially, but the
combinatorial description is subtle. The description of how
irreducibles decompose into irreducibles is called a *branching rule*.

References for this topic:

- [FauserEtAl2006]_

- [King1975]_

- [HoweEtAl2005]_

- [McKayPatera1981]_

Sage can compute how a character of `G` restricts to `H`. It does
so not by memorizing a combinatorial rule, but by computing the
character and restricting the character to a maximal torus of `H`.
What Sage has memorized (in a series of built-in encoded rules)
are the various embeddings of maximal tori of maximal subgroups of `G`.
This approach to computing branching rules has a limitation: the
character must fit into memory and be computable by Sage's
internal code (based on the Freudenthal multiplicity formula)
in real time.

Sage has enough built-in branching rules to handle all cases where `G`
is a classical group, that is, type A, B, C or D. It also has many
built-in cases where `G` is an exceptional group.

Clearly it is sufficient to consider the case where `H` is a maximal
subgroup of `G`, since if this is known then one may branch down
successively through a series of subgroups, each maximal in its
predecessors. A problem is therefore to understand the maximal
subgroups in a Lie group, and to give branching rules for each.

For convenience Sage includes some branching rules to non-maximal
subgroups, but strictly speaking these are not necessary, since
one could do any branching rule to a non-maximal subgroup using
only branching rules to maximal subgroups. The goal is
to give a sufficient set of built-in branching rules for all maximal
subgroups, and this is accomplished for classical groups (types A, B,
C or D) at least up to rank 8, and for many maximal subgroups of
exceptional groups.


Levi subgroups
--------------

A Levi subgroup may or may not be maximal. They are easily
classified. If one starts with a Dynkin diagram for `G` and removes a
single node, one obtains a smaller Dynkin diagram, which is the Dynkin
diagram of a smaller subgroup `H`.

For example, here is the A3 Dynkin diagram:

::

    sage: A3 = WeylCharacterRing("A3")
    sage: A3.dynkin_diagram()
    O---O---O
    1   2   3
    A3

We see that we may remove the node 3 and obtain A2, or the node 2 and
obtain A1xA1. These correspond to the Levi subgroups `GL(3)` and
`GL(2) \times GL(2)` of `GL(4)`. Let us construct the irreducible
representations of `GL(4)` and branch them down to these down to
`GL(3)` and `GL(2) \times GL(2)`:

.. link

::

    sage: reps = [A3(v) for v in A3.fundamental_weights()]; reps
    [A3(1,0,0,0), A3(1,1,0,0), A3(1,1,1,0)]
    sage: A2 = WeylCharacterRing("A2")
    sage: A1xA1 = WeylCharacterRing("A1xA1")
    sage: [pi.branch(A2, rule="levi") for pi in reps]
    [A2(0,0,0) + A2(1,0,0), A2(1,0,0) + A2(1,1,0), A2(1,1,0) + A2(1,1,1)]
    sage: [pi.branch(A1xA1, rule="levi") for pi in reps]
    [A1xA1(1,0,0,0) + A1xA1(0,0,1,0),
     A1xA1(1,1,0,0) + A1xA1(1,0,1,0) + A1xA1(0,0,1,1),
      A1xA1(1,1,1,0) + A1xA1(1,0,1,1)]

Let us redo this calculation in coroot notation. As we have explained,
coroot notation does not distinguish between representations of
`GL(4)` that have the same restriction to `SL(4)`, so in effect we are
now working with the groups `SL(4)` and its Levi subgroups `SL(3)` and
`SL(2) \times SL(2)`::

    sage: A3 = WeylCharacterRing("A3", style="coroots")
    sage: reps = [A3(v) for v in A3.fundamental_weights()]; reps
    [A3(1,0,0), A3(0,1,0), A3(0,0,1)]
    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
    sage: [pi.branch(A2, rule="levi") for pi in reps]
    [A2(0,0) + A2(1,0), A2(0,1) + A2(1,0), A2(0,0) + A2(0,1)]
    sage: [pi.branch(A1xA1, rule="levi") for pi in reps]
    [A1xA1(1,0) + A1xA1(0,1), 2*A1xA1(0,0) + A1xA1(1,1), A1xA1(1,0) + A1xA1(0,1)]

Now we may observe a distinction difference in branching from

.. MATH::

    GL(4) \to GL(2) \times GL(2)

versus

.. MATH::

    SL(4) \to SL(2) \times SL(2).

Consider the representation `A3(0,1,0)`, which is the six dimensional exterior
square. In the coroot notation, the restriction contained two copies of the
trivial representation, ``2*A1xA1(0,0)``. The other way, we had instead three
distinct representations in the restriction, namely ``A1xA1(1,1,0,0)`` and
``A1xA1(0,0,1,1)``, that is, `\det \otimes 1` and `1 \otimes \det`.

The Levi subgroup ``A1xA1`` is actually not maximal. Indeed, we may
factor the embedding:

.. MATH::

    SL(2) \times SL(2) \to Sp(4) \to SL(4).

Therfore there are branching rules ``A3 -> C2`` and ``C2 -> A2``, and
we could accomplish the branching in two steps, thus::

    sage: A3 = WeylCharacterRing("A3", style="coroots")
    sage: C2 = WeylCharacterRing("C2", style="coroots")
    sage: B2 = WeylCharacterRing("B2", style="coroots")
    sage: D2 = WeylCharacterRing("D2", style="coroots")
    sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
    sage: reps = [A3(fw) for fw in A3.fundamental_weights()]
    sage: [pi.branch(C2, rule="symmetric").branch(B2, rule="isomorphic"). \
             branch(D2, rule="extended").branch(A1xA1, rule="isomorphic") for pi in reps]
    [A1xA1(1,0) + A1xA1(0,1), 2*A1xA1(0,0) + A1xA1(1,1), A1xA1(1,0) + A1xA1(0,1)]

As you can see, we've redone the branching rather circuitously this
way, making use of the branching rules ``A3->C2`` and ``B2->D2``, and
two accidental isomorphisms ``C2=B2`` and ``D2=A1xA1``. It is much
easier to go in one step using ``rule="levi"``, but reassuring that we
get the same answer!


Subgroups classified by the extended Dynkin diagram
---------------------------------------------------

It is also true that if we remove one node from the extended Dynkin
diagram that we obtain the Dynkin diagram of a subgroup. For example::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: G2.extended_dynkin_diagram()
      3
    O=<=O---O
    1   2   0
    G2~

Observe that by removing the 1 node that we obtain an A2 Dynkin
diagram. Therefore the exceptional group G2 contains a copy of
`SL(3)`. We branch the two representations of G2 corresponding to the
fundamental weights to this copy of A2::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: A2 = WeylCharacterRing("A2", style="coroots")
    sage: [G2(f).degree() for f in G2.fundamental_weights()]
    [7, 14]
    sage: [G2(f).branch(A2, rule="extended") for f in G2.fundamental_weights()]
    [A2(0,0) + A2(0,1) + A2(1,0), A2(0,1) + A2(1,0) + A2(1,1)]

The two representations of G2, of degrees 7 and 14 respectively, are
the action on the octonions of trace zero and the adjoint
representation.

For embeddings of this type, the rank of the subgroup `H` is the same
as the rank of `G`. This is in contrast with embeddings of Levi type,
where `H` has rank one less than `G`.


Orthogonal and symplectic subgroups of orthogonal and symplectic groups
-----------------------------------------------------------------------

If `G = SO(r+s)` then `G` has a subgroup `SO(r) \times SO(s)`. This
lifts to an embedding of the universal covering groups

.. MATH::

    \hbox{spin}(r) \times \hbox{spin}(s) \to \hbox{spin}(r+s).

Sometimes this embedding is of extended type, and sometimes it is
not. It is of extended type unless `r` and `s` are both odd. If it is
of extended type then you may use ``rule="extended"``. In any case you
may use ``rule="orthogonal_sum"``. The name refer to the origin of the
embedding `SO(r) \times SO(s) \to SO(r+s)` from the decomposition of
the underlying quadratic space as a direct sum of two orthogonal
subspaces.

There are four cases depending on the parity of `r` and `s`. For
example, if `r = 2k` and `s = 2l` we have an embedding::

    ['D',k] x ['D',l] --> ['D',k+l]

This is of extended type. Thus consider the embedding
``D4xD3 -> D7``. Here is the extended Dynkin diagram::

      0 O           O 7
        |           |
        |           |
    O---O---O---O---O---O
    1   2   3   4   5   6

Removing the 4 vertex results in a disconnected Dynkin diagram::

      0 O           O 7
        |           |
        |           |
    O---O---O       O---O
    1   2   3       5   6

This is ``D4xD3``.  Therefore use the "extended" branching rule:

::

    sage: D7 = WeylCharacterRing("D7", style="coroots")
    sage: D4xD3 = WeylCharacterRing("D4xD3", style="coroots")
    sage: spin = D7(D7.fundamental_weights()[7]); spin
    D7(0,0,0,0,0,0,1)
    sage: spin.branch(D4xD3, rule="extended")
    D4xD3(0,0,1,0,0,1,0) + D4xD3(0,0,0,1,0,0,1)

But we could equally well use the "orthogonal_sum" rule.

.. link

::

    sage: spin.branch(D4xD3, rule="orthogonal_sum")
    D4xD3(0,0,1,0,0,1,0) + D4xD3(0,0,0,1,0,0,1)

Similarly we have embeddings::

    ['D',k] x ['B',l] --> ['B',k+l]

These are also of extended type. For example consider the embedding of
``D3xB2->B5``. Here is the ``B5`` extended Dynkin diagram::

        O 0
        |
        |
    O---O---O---O=>=O
    1   2   3   4   5

Removing the 3 node gives::

        O 0
        |
    O---O       O=>=O
    1   2       4   5

and this is the Dynkin diagram or ``D3xB2``. For such branchings we
again use either ``rule="extended"`` or ``rule="orthogonal_sum"``.

Finally, there is an embedding ::

    ['B',k] x ['B',l] --> ['D',k+l+1]

This is *not* of extended type, so you may not use ``rule="extended"``.
You *must* use ``rule="orthogonal_sum"``::

     sage: D5=WeylCharacterRing("D5",style="coroots")
     sage: B2xB2=WeylCharacterRing("B2xB2",style="coroots")
     sage: [D5(v).branch(B2xB2,rule="orthogonal_sum") for v in D5.fundamental_weights()]
     [B2xB2(1,0,0,0) + B2xB2(0,0,1,0),
      B2xB2(0,2,0,0) + B2xB2(1,0,1,0) + B2xB2(0,0,0,2),
      B2xB2(0,2,0,0) + B2xB2(0,2,1,0) + B2xB2(1,0,0,2) + B2xB2(0,0,0,2),
      B2xB2(0,1,0,1), B2xB2(0,1,0,1)]

Symmetric subgroups
-------------------

If `G` admits an outer automorphism (usually of order two) then we may
try to find the branching rule to the fixed subgroup `H`. It can be
arranged that this automorphism maps the maximal torus `T` to itself
and that a maximal torus `U` of `H` is contained in `T`.

Suppose that the Dynkin diagram of `G` admits an automorphism. Then
`G` itself admits an outer automorphism. The Dynkin diagram of the
group `H` of invariants may be obtained by "folding" the Dynkin
diagram of `G` along the automorphism. The exception is the branching
rule `GL(2r) \to SO(2r)`.

Here are the branching rules that can be obtained using
``rule="symmetric"``.

+------------+-------------+---------------------------+
| `G`        | `H`         | Cartan Types              |
+============+=============+===========================+
| `GL(2r)`   | `Sp(2r)`    | ``['A',2r-1] => ['C',r]`` |
+------------+-------------+---------------------------+
| `GL(2r+1)` | `SO(2r+1)`  | ``['A',2r] => ['B',r]``   |
+------------+-------------+---------------------------+
| `GL(2r)`   | `SO(2r)`    | ``['A',2r-1] => ['D',r]`` |
+------------+-------------+---------------------------+
| `SO(2r)`   | `SO(2r-1)`  | ``['D',r] => ['B',r-1]``  |
+------------+-------------+---------------------------+
| `E_6`      | `F_4`       | ``['E',6] => ['F',4]``    |
+------------+-------------+---------------------------+


Tensor products
---------------

If `G_1` and `G_2` are Lie groups, and we have representations
`\pi_1: G_1 \to GL(n)` and `\pi_2: G_2 \to GL(m)` then the tensor
product is a representation of `G_1 \times G_2`. It has its image
in `GL(nm)` but sometimes this is conjugate to a subgroup of `SO(nm)`
or `Sp(nm)`. In particular we have the following cases.

+-------------------+---------------------------+------------------------------------------+
| Group             | Subgroup                  | Cartan Types                             |
+===================+===========================+==========================================+
| `GL(rs)`          | `GL(r)\times GL(s)`       | ``['A', rs-1] => ['A',r-1] x ['A',s-1]`` |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs+2r+2s+1)` | `SO(2r+1)\times SO(2s+1)` | ``['B',2rs+r+s] => ['B',r] x ['B',s]``   |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs+2s)`      | `SO(2r+1)\times SO(2s)`   | ``['D',2rs+s] => ['B',r] x ['D',s]``     |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs)`         | `SO(2r)\times SO(2s)`     | ``['D',2rs] => ['D',r] x ['D',s]``       |
+-------------------+---------------------------+------------------------------------------+
| `SO(4rs)`         | `Sp(2r)\times Sp(2s)`     | ``['D',2rs] => ['C',r] x ['C',s]``       |
+-------------------+---------------------------+------------------------------------------+
| `Sp(4rs+2s)`      | `SO(2r+1)\times Sp(2s)`   | ``['C',2rs+s] => ['B',r] x ['C',s]``     |
+-------------------+---------------------------+------------------------------------------+
| `Sp(4rs)`         | `Sp(2r)\times SO(2s)`     | ``['C',2rs] => ['C',r] x ['D',s]``       |
+-------------------+---------------------------+------------------------------------------+

These branching rules are obtained using ``rule="tensor"``.


Symmetric powers
----------------

The `k`-th symmetric and exterior power homomorphisms map
`GL(n) \to GL \left({n+k-1 \choose k} \right)` and
`GL \left({n \choose k} \right)`. The corresponding branching rules
are not implemented but a special case is. The `k`-th symmetric power
homomorphism `SL(2) \to GL(k+1)` has its image inside of `SO(2r+1)` if
`k = 2r` and inside of `Sp(2r)` if `k = 2r-1`. Hence there are
branching rules::

    ['B',r] => A1
    ['C',r] => A1

and these may be obtained using ``rule="symmetric_power"``.


Plethysms
---------

The above branching rules are sufficient for most cases, but a few
fall between the cracks. Mostly these involve maximal subgroups of
fairly small rank.

The rule ``rule="plethysm"`` is a powerful rule that includes any
branching rule from types A, B, C or D as a special case. Thus it
could be used in place of the above rules and would give the same
results. However, it is most useful when branching from `G` to a
maximal subgroup `H` such that `rank(H) < rank(G)-1`.

We consider a homomorphism `H \to G` where `G` is one of `SL(r+1)`,
`SO(2r+1)`, `Sp(2r)` or `SO(2r)`. The function
``branching_rule_from_plethysm`` produces the corresponding branching
rule. The main ingredient is the character `\chi` of the
representation of `H` that is the homomorphism to `GL(r+1)`,
`GL(2r+1)` or `GL(2r)`.

Let us consider the symmetric fifth power representation of
`SL(2)`. This is implemented above by ``rule="symmetric_power"``, but
suppose we want to use ``rule="plethysm"``. First we construct the
homomorphism by invoking its character, to be called ``chi``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: chi = A1([5])
    sage: chi.degree()
    6
    sage: chi.frobenius_schur_indicator()
    -1

This confirms that the character has degree 6 and is symplectic, so it
corresponds to a homomorphism `SL(2) \to Sp(6)`, and there is a
corresponding branching rule ``C3 => A1``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: C3 = WeylCharacterRing("C3", style="coroots")
    sage: chi = A1([5])
    sage: sym5rule = branching_rule_from_plethysm(chi, "C3")
    sage: [C3(hwv).branch(A1, rule=sym5rule) for hwv in C3.fundamental_weights()]
    [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

This is identical to the results we would obtain using
``rule="symmetric_power"``::

    sage: A1 = WeylCharacterRing("A1", style="coroots")
    sage: C3 = WeylCharacterRing("C3", style="coroots")
    sage: [C3(v).branch(A1, rule="symmetric_power") for v in C3.fundamental_weights()]
    [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

But the next example of plethysm gives a branching rule not available
by other methods::

    sage: G2 = WeylCharacterRing("G2", style="coroots")
    sage: D7 = WeylCharacterRing("D7", style="coroots")
    sage: ad = G2(0,1); ad.degree()
    14
    sage: ad.frobenius_schur_indicator()
    1
    sage: for r in D7.fundamental_weights():  # long time (26s on sage.math, 2012)
    ...      print D7(r).branch(G2, rule=branching_rule_from_plethysm(ad, "D7"))
    ...
    G2(0,1)
    G2(0,1) + G2(3,0)
    G2(0,0) + G2(2,0) + G2(3,0) + G2(0,2) + G2(4,0)
    G2(0,1) + G2(2,0) + G2(1,1) + G2(0,2) + G2(2,1) + G2(4,0) + G2(3,1)
    G2(1,0) + G2(0,1) + G2(1,1) + 2*G2(3,0) + 2*G2(2,1) + G2(1,2) + G2(3,1) + G2(5,0) + G2(0,3)
    G2(1,1)
    G2(1,1)


In this example, `ad` is the 14-dimensional adjoint representation of the
exceptional group `G_2`. Since the Frobenius-Schur indicator is 1, the
representation is orthogonal, and factors through `SO(14)`, that is, `D7`.

Miscellaneous other subgroups
-----------------------------

Use ``rule="miscellaneous"`` for the branching rule ``B3 => G2``. This
may also be obtained using a plethysm but for convenience this one is
hand-coded.


Nuts and bolts of branching rules
---------------------------------

Sage has many built-in branching rules, enough to handle most
cases. However, if you find a case where there is no existing rule,
you may code it by hand. Moreover, it may be useful to understand how
the built-in rules work.

Suppose you want to branch from a group `G` to a subgroup `H`.
Arrange the embedding so that a Cartan subalgebra `U` of `H` is
contained in a Cartan subalgebra `T` of `G`. There is thus a mapping
from the weight spaces `\hbox{Lie}(T)^* \to \hbox{Lie}(U)^*`.  Two
embeddings will produce identical branching rules if they differ by an
element of the Weyl group of `H`. The *rule* is this map
`\hbox{Lie}(T)^*` ``= G.space()`` to
`\hbox{Lie}(U)^*` ``= H.space()``, which you may implement as a
function.

As an example, let us consider how to implement the branching rule
``A3 -> C2``.  Here ``H = C2 = Sp(4)`` embedded as a subgroup in
``A3 = GL(4).`` The Cartan subalgebra `\hbox{Lie}(U)` consists of
diagonal matrices with eigenvalues ``u1, u2, -u2, -u1``. Then
``C2.space()`` is the two dimensional vector spaces consisting of the
linear functionals ``u1`` and ``u2`` on ``U``. On the other hand
`Lie(T) = \mathbf{R}^4`. A convenient way to see the restriction is to
think of it as the adjoint of the map ``[u1,u2] -> [u1,u2,-u2,-u1]``,
that is, ``[x0,x1,x2,x3] -> [x0-x3,x1-x2].`` Hence we may encode the
rule::

    def brule(x):
        return [x[0]-x[3], x[1]-x[2]]

or simply::

    brule = lambda x: [x[0]-x[3], x[1]-x[2]]

Let us check that this agrees with the built-in rule::

    sage: A3 = WeylCharacterRing(['A', 3])
    sage: C2 = WeylCharacterRing(['C', 2])
    sage: brule = lambda x: [x[0]-x[3], x[1]-x[2]]
    sage: A3(1,1,0,0).branch(C2, rule=brule)
    C2(0,0) + C2(1,1)
    sage: A3(1,1,0,0).branch(C2, rule="symmetric")
    C2(0,0) + C2(1,1)


Automorphisms and triality
--------------------------

The case where `G=H` can be treated as a special case of a branching
rule. In most cases if `G` has a nontrivial outer automorphism, it
has order two, corresponding to the symmetry of the Dynkin diagram.
Such an involution exists in the cases `A_r`, `D_r`, `E_6`.

So the automorphism acts on the representations of `G`, and its
effect may be computed using the branching rule code::

    sage: A4=WeylCharacterRing("A4",style="coroots")
    sage: A4(1,0,1,0).degree()
    45
    sage: A4(0,1,0,1).degree()
    45
    sage: A4(1,0,1,0).branch(A4,rule="automorphic")
    A4(0,1,0,1)

In the special case where `G=D4`, the Dynkin diagram has
extra symmetries::


    sage: D4=WeylCharacterRing("D4",style="coroots")
    sage: D4.dynkin_diagram()
        O 4
        |
        |
    O---O---O
    1   2   3
    D4

The automorphism group of the Dynkin diagram has order 6,
corresponds to any permutation of the three outer nodes.
Therefore the group has extra outer automorphisms.  One
where an additional automorphism of order three can be
obtained using ``rule="triality"``. This is not an
automorphisms of `SO(8)`, but of its double cover
`spin(8)`. Note that `spin(8)` has three representations
of degree 8, namely the standard representation of
`SO(8)` and the two eight-dimensional spin
representations. These are permuted by triality:

.. link:

::

    sage: D4(0,0,0,1).branch(D4,rule="triality")
    D4(1,0,0,0)
    sage: D4(0,0,0,1).branch(D4,rule="triality").branch(D4,rule="triality")
    D4(0,0,1,0)
    sage: D4(0,0,0,1).branch(D4,rule="triality").branch(D4,rule="triality").branch(D4,rule="triality")
    D4(0,0,0,1)

By contrast, ``rule="automorphic"`` simply interchanges the two
spin representations, as it always does in Type D:

.. link:

::

    sage: D4(0,0,0,1).branch(D4,rule="automorphic")
    D4(0,0,1,0)
    sage: D4(0,0,1,0).branch(D4,rule="automorphic")
    D4(0,0,0,1)
