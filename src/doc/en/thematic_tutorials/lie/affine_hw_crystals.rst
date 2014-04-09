==============================
Affine Highest Weight Crystals
==============================

Affine highest weight crystals are infinite-dimensional. Hence, to work
with them in Sage, we need some further tools.

Littelmann path model
---------------------

The Littelmann path model for highest weight crystals is implemented
in Sage. It models finite highest crystals as well as affine highest weight
crystals which are infinite dimensional. The elements of the crystal are
piecewise linear maps in the weight space. For more information on the
Littelmann path model, see [L1995]_.

Since the affine highest weight crystals are infinite, it is not possible
to list all elements or draw the entire crystal graph. However, if the user
is only interested in the crystal up to a certain distance or depth from the highest
weight element, then one can work with the corresponding subcrystal.
To view the corresponding upper part of the crystal, one can build the
associated digraph::

    sage: R = RootSystem(['C',3,1])
    sage: La = R.weight_space().basis()
    sage: LS = CrystalOfLSPaths(2*La[1]); LS
    The crystal of LS paths of type ['C', 3, 1] and weight 2*Lambda[1]
    sage: C = LS.subcrystal(max_depth=3)
    sage: sorted(C, key=str)
    [(-Lambda[0] + 2*Lambda[1] - Lambda[2] + Lambda[3], Lambda[1]),
     (-Lambda[0] + Lambda[1] + Lambda[2], Lambda[0] - Lambda[1] + Lambda[2]),
     (-Lambda[0] + Lambda[1] + Lambda[2], Lambda[1]),
     (-Lambda[1] + 2*Lambda[2], Lambda[1]),
     (2*Lambda[0] - 2*Lambda[1] + 2*Lambda[2],),
     (2*Lambda[1],),
     (Lambda[0] + Lambda[2] - Lambda[3], Lambda[1]),
     (Lambda[0] - Lambda[1] + Lambda[2], Lambda[1]),
     (Lambda[0] - Lambda[2] + Lambda[3], Lambda[0] - Lambda[1] + Lambda[2]),
     (Lambda[0] - Lambda[2] + Lambda[3], Lambda[1])]

    sage: G = LS.digraph(subset = C)
    sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

.. image:: ../media/LScrystal.png
   :scale: 50
   :align: center

REFERENCES:

.. [L1995] P. Littelmann, Paths and root operators in representation theory. Ann. of Math. (2) 142 (1995), no. 3, 499-525.

The Littelmann path model also lends itself as a model for to level zero crystals which are bi-infinite.
To cut out a slice of these crystals, one can use the direction option in subcrystal::

    sage: R = RootSystem(['A',2,1])
    sage: La = R.weight_space(extended = True).basis()
    sage: LS = CrystalOfLSPaths(La[1]-La[0]); LS
    The crystal of LS paths of type ['A', 2, 1] and weight -Lambda[0] + Lambda[1]
    sage: C = LS.subcrystal(max_depth=2, direction = 'both')
    sage: G = LS.digraph(subset = C)
    sage: G.set_latex_options(edge_options = lambda (u,v,label): ({}))
    sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

.. image:: ../media/level_zero_crystal.png
   :scale: 50
   :align: center
