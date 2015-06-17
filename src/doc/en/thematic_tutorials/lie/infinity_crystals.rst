-----------------
Infinity Crystals
-----------------

Infinity crystals are the crystal analogue of Verma modules with highest weight
`0` associated to a symmetrizable Kac-Moody algebra.  As such, they are
infinite-dimensional and any irreducible highest weight crystal may be obtained
from an infinity crystal via some cutting procedure.  On the other hand, the
crystal `B(\infty)` is the direct limit of all irreducible highest weight
crystals `B(\lambda)`, so there are natural embeddings of each `B(\lambda)` in
`B(\infty)`.  Below, we outline the various implementations of the crystal
`B(\infty)` in Sage and give examples of how each `B(\lambda)` interacts with
`B(\infty)`.

All infinity crystals that are currently implemented in Sage can be accessed
by typing ``crystals.infinity.<tab>``.


Marginally large tableaux
-------------------------

Marginally large tableaux were introduced by J. Hong and H. Lee as a realization
of the crystal `B(\infty)` in types `A_n`, `B_n`, `C_n`, `D_{n+1}`, and `G_2`.
The marginally large condition guarantees that all tableau have exactly `n`
rows and that the number of `i`-boxes in the `i`-th row from the top (in
the English convention) is exactly one more than the total number of boxes in
the `(i+1)`-st row.  Other specific conditions on the tableaux vary by type.
See [HongLee2008]_ for more information.

Here is an example in type `A_2`::

    sage: B = crystals.infinity.Tableaux(['A',2])
    sage: b = B.highest_weight_vector()
    sage: b.pp()
      1  1
      2
    sage: b.f_string([1,2,2,1,2,1,2,2,2,2,2]).pp()
      1  1  1  1  1  1  1  1  1  2  2  3
      2  3  3  3  3  3  3  3

Since the crystal is infinite, we must specify the subcrystal we would like to
view::

    sage: S = B.subcrystal(max_depth=4)
    sage: G = B.digraph(subset=S)
    sage: view(G, tightpage=True) # optional - dot2tex graphviz, not tested (opens external window)

.. image:: ../media/BinfA2.png
   :scale: 50
   :align: center

Using elementary crystals, we can cut irreducible highest weight crystals from
`B(\infty)`.  In this example, we cut out `B(\rho)` from `B(\infty)` in type
`A_2`::

    sage: T = crystals.elementary.T(['A',2], B.Lambda()[1] + B.Lambda()[2])
    sage: t = T[0]
    sage: C = crystals.elementary.Component(['A',2])
    sage: c = C[0]
    sage: TP = crystals.TensorProduct(C,T,B)
    sage: t0 = TP(c,t,b)
    sage: STP = TP.subcrystal(generators=[t0])
    sage: GTP = TP.digraph(subset=STP)
    sage: view(GTP, tightpage=True) # optional - dot2tex graphviz, not tested (opens external window)

.. image:: ../media/BinfTCrhoA2.png
   :scale: 50
   :align: center

Note that the above code can be simplified using the R-crystal::

    sage: R = crystals.elementary.R(['A',2], B.Lambda()[1] + B.Lambda()[2])
    sage: r = R[0]
    sage: TP2 = crystals.TensorProduct(R,B)
    sage: t2 = TP2(r,b)
    sage: STP2 = TP2.subcrystal(generators=[t2])
    sage: GTP2 = TP2.digraph(subset=STP2)
    sage: view(GTP2, tightpage=True) # optional - dot2tex graphviz, not tested (opens external window)

On the other hand, we can embed the irreducible highest weight crystal
`B(\rho)` into `B(\infty)`::

    sage: Brho = crystals.Tableaux(['A',2],shape=[2,1])
    sage: brho = Brho.highest_weight_vector()
    sage: B = crystals.infinity.Tableaux(['A',2])
    sage: binf = B.highest_weight_vector()
    sage: Psi = Brho.crystal_morphism({brho : b})
    sage: BG = B.digraph(subset=[Psi(x) for x in Brho])
    sage: view(BG, tightpage=True) # optional - dot2tex graphviz, not tested (opens external window)


Generalized Young walls
-----------------------

Generalized Young walls were introduced by J.-A. Kim and D.-U. Shin as a model
for `B(\infty)` and each `B(\lambda)` solely in affine type `A_n^{(1)}`. See
[KimShin2010]_ for more information on the construction of generalized Young
walls.

Since this model is only valid for one Cartan type, the input to initialize the
crystal is simply the rank of the underlying type::

    sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
    sage: y = Y.highest_weight_vector()
    sage: y.f_string([0,1,2,2,2,1,0,0,1,2]).pp()
             2|
              |
              |
           1|2|
           0|1|
     2|0|1|2|0|

In the ``weight`` method for this model, we can choose whether to view weights
in the root lattice or in the extended weight lattice::




Modified Nakajima monomials
---------------------------



Rigged configurations
---------------------



