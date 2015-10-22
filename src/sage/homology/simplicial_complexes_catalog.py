# -*- coding: utf-8 -*-
"""
Examples of simplicial complexes

There are two main types: manifolds and examples related to graph
theory.

For manifolds, there are functions defining the `n`-sphere for any
`n`, the torus, `n`-dimensional real projective space for any `n`, the
complex projective plane, surfaces of arbitrary genus, and some other
manifolds, all as simplicial complexes.

Aside from surfaces, this file also provides functions for
constructing some other simplicial complexes: the simplicial complex
of not-`i`-connected graphs on `n` vertices, the matching complex on n
vertices, the chessboard complex for an `n` by `i` chessboard, and
others.  These provide examples of large simplicial complexes; for
example, ``simplicial_complexes.NotIConnectedGraphs(7,2)`` has over a
million simplices.

All of these examples are accessible by typing
``simplicial_complexes.NAME``, where ``NAME`` is the name of the example.

- :func:`BarnetteSphere`
- :func:`BrucknerGrunbaumSphere`
- :func:`ChessboardComplex`
- :func:`ComplexProjectivePlane`
- :func:`K3Surface`
- :func:`KleinBottle`
- :func:`MatchingComplex`
- :func:`MooreSpace`
- :func:`NotIConnectedGraphs`
- :func:`PoincareHomologyThreeSphere`
- :func:`PseudoQuaternionicProjectivePlane`
- :func:`RandomComplex`
- :func:`RealProjectivePlane`
- :func:`RealProjectiveSpace`
- :func:`Simplex`
- :func:`Sphere`
- :func:`SumComplex`
- :func:`SurfaceOfGenus`
- :func:`Torus`

You can also get a list by typing ``simplicial_complexes.`` and hitting the
TAB key.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
    sage: S.homology()
    {0: 0, 1: 0, 2: Z}
    sage: simplicial_complexes.SurfaceOfGenus(3)
    Simplicial complex with 15 vertices and 38 facets
    sage: M4 = simplicial_complexes.MooreSpace(4)
    sage: M4.homology()
    {0: 0, 1: C4, 2: 0}
    sage: simplicial_complexes.MatchingComplex(6).homology()
    {0: 0, 1: Z^16, 2: 0}
"""

from sage.homology.examples import (Sphere, Simplex, Torus, ProjectivePlane,
        RealProjectivePlane, KleinBottle, SurfaceOfGenus, MooreSpace,
        ComplexProjectivePlane, PseudoQuaternionicProjectivePlane,
        PoincareHomologyThreeSphere, RealProjectiveSpace, K3Surface,
        BarnetteSphere, BrucknerGrunbaumSphere, NotIConnectedGraphs,
        MatchingComplex, ChessboardComplex, RandomComplex, SumComplex)
