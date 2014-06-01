.. linkall

=====================================
Young's lattice and the RSK algorithm
=====================================

This section provides some examples on the RSK algorithm explained in
Chapter 8 of Stanley's book [Stanley2013]_.

We begin by creating the first few levels of Young's lattice `Y`. For
we need to define the elements and the order relation for the poset,
which is containment of partitions::

    sage: level = 6
    sage: elements = [b for n in range(level) for b in Partitions(n)]
    sage: ord = lambda x,y: y.contains(x)
    sage: Y = Poset((elements,ord), facade=True)
    sage: H = Y.hasse_diagram()
    sage: view(H)      # optional

.. image:: ../media/young_lattice.png
   :scale: 60
   :align: center

We can now define the up and down operators on `\QQ Y`. First we do
so on partitions::

    sage: QQY = CombinatorialFreeModule(QQ,elements)

    sage: def up_operator_on_partitions(la):
    ....:     covers = Y.upper_covers(la)
    ....:     return sum(QQY(c) for c in covers)

    sage: def down_operator_on_partitions(la):
    ....:     covers = Y.lower_covers(la)
    ....:     return sum(QQY(c) for c in covers)

Here is the result when we apply the operators to the partition `(2,1)`::

    sage: la = Partition([2,1])
    sage: up_operator_on_partitions(la)
    B[[2, 1, 1]] + B[[2, 2]] + B[[3, 1]]
    sage: down_operator_on_partitions(la)
    B[[1, 1]] + B[[2]]

Now we define the up and down operator on `\QQ Y`::

    sage: def up_operator(b):
    ....:     return sum(b.coefficient(p)*up_operator_on_partitions(p) for p in b.support())

    sage: def down_operator(b):
    ....:     return sum(b.coefficient(p)*down_operator_on_partitions(p) for p in b.support())

We can check the identity `D_{i+1} U_i - U_{i-1} D_i = I_i` explicitly on
all partitions of 3 (so that `i=3`)::

    sage: for p in Partitions(3):
    ....:     b = QQY(p)
    ....:     assert down_operator(up_operator(b)) - up_operator(down_operator(b)) == b
    ....:

We can also check that the coefficient of `\lambda \vdash n` in `U^n(\emptyset)` is equal
to the number of standard Young tableaux of shape `\lambda`::

    sage: u = QQY(Partition([]))
    sage: for i in range(4):
    ....:     u = up_operator(u)
    ....:
    sage: u
    B[[1, 1, 1, 1]] + 3*B[[2, 1, 1]] + 2*B[[2, 2]] + 3*B[[3, 1]] + B[[4]]

For example, the number of standard Young tableaux of shape `(2,1,1)` is 3::

    sage: StandardTableaux([2,1,1]).cardinality()
    3

We can test this in general::

    sage: for la in u.support():
    ....:     assert u[la] == StandardTableaux(la).cardinality()
    ....:

Let us now turn to the RSK algorithm. If we want to verify Example 8.12, we can
do so as follows::

    sage: p = Permutation([4,2,7,3,6,1,5])
    sage: p.robinson_schensted()
    [[[1, 3, 5], [2, 6], [4, 7]], [[1, 3, 5], [2, 4], [6, 7]]]

The tableaux can also be displayed as tableaux::

    sage: P,Q = p.robinson_schensted()
    sage: P.pp()
    1  3  5
    2  6
    4  7
    sage: Q.pp()
    1  3  5
    2  4
    6  7

The inverse RSK algorithm is also implemented::

    sage: RSK_inverse(P,Q, output='word')
    word: 4273615


