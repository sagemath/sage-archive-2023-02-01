"""
Homomorphisms Between Matrix Groups

AUTHORS:

- David Joyner and William Stein (2006-03): initial version

- David Joyner (2006-05): examples

- Simon King (2011-01): cleaning and improving code

- Volker Braun (2013-1) port to new Parent, libGAP.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.categories.morphism import Morphism
from sage.misc.latex import latex
from sage.groups.libgap_morphism import GroupMorphism_libgap


def to_libgap(x):
    """
    Helper to convert ``x`` to a LibGAP matrix or matrix group
    element.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.morphism import to_libgap
        sage: to_libgap(GL(2,3).gen(0))
        [ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        sage: to_libgap(matrix(QQ, [[1,2],[3,4]]))
        [ [ 1, 2 ], [ 3, 4 ] ]
    """
    try:
        return x.gap()
    except AttributeError:
        from sage.libs.gap.libgap import libgap
        return libgap(x)


class MatrixGroupMap(Morphism):
    r"""
    """
    def _repr_type(self):
        """
        Part of the implementation of :meth:`_repr_`

        EXAMPLES::

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))._repr_type()
            'MatrixGroup'
        """
        return "MatrixGroup"


class MatrixGroupMorphism(MatrixGroupMap):
    pass


class MatrixGroupMorphism_im_gens(MatrixGroupMorphism, GroupMorphism_libgap):
    """
    Group morphism specified by images of generators.

    EXAMPLES::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi
        Group morphism:
        From: Matrix group over Finite Field of size 5 with 1 generators (
        [1 1]
        [0 1]
        )
        To:   Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [1 1]
        )
        sage: phi(MS([1,1,0,1]))
        [1 0]
        [1 1]
        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])

    TESTS:

    Check that :trac:`19406` is fixed::

        sage: G = GL(2, GF(3))
        sage: H = GL(3, GF(2))
        sage: mat1 = H([[-1,0,0],[0,0,-1],[0,-1,0]])
        sage: mat2 = H([[1,1,1],[0,0,-1],[-1,0,0]])
        sage: phi = G.hom([mat1, mat2])
        Traceback (most recent call last):
        ...
        ValueError: images do not define a group homomorphism

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi.gap()
        CompositionMapping( [ (6,7,8,10,9)(11,13,14,12,15)(16,19,20,18,17)(21,25,22,24,23) ]
        -> [ [ [ Z(5)^0, 0*Z(5) ], [ Z(5)^0, Z(5)^0 ] ] ], <action isomorphism> )
        sage: type(_)
        <type 'sage.libs.gap.element.GapElement'>

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi.kernel()
        Matrix group over Finite Field of size 7 with 1 generators (
        [6 0]
        [0 1]
        )

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi
        Group endomorphism of Matrix group over Finite Field of size 7 with 1 generators (
        [3 0]
        [0 1]
        )
        sage: g = G.gens()[0]
        sage: phi(g)
        [2 0]
        [0 1]
        sage: H = MatrixGroup([MS(a.list())])
        sage: H
        Matrix group over Finite Field of size 7 with 1 generators (
        [2 0]
        [0 1]
        )

    The following tests against :trac:`10659`::

        sage: phi(H)   # indirect doctest
        Matrix group over Finite Field of size 7 with 1 generators (
        [4 0]
        [0 1]
        )
        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: g = MS([1,1,0,1])
        sage: G = MatrixGroup([g])
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 1]
        [0 1]
        sage: phi(G(g^2))
        [1 2]
        [0 1]

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([1,2,  -1,1]),MS([1,1,  0,1])]
        sage: G = MatrixGroup(gens)
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 2]
        [4 1]
        sage: phi(G.1)
        [1 1]
        [0 1]

    The following tests that the call method was successfully
    improved in :trac:`10659`::

        sage: O = WeylGroup(['D',6])
        sage: r = prod(O.gens())
        sage: r_ = r^-1
        sage: f = O.hom([r*x*r_ for x in O.gens()])  # long time (19s on sage.math, 2011)
        sage: [f(x) for x in O.gens()]  # long time
        [
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
        [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
        [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
        <BLANKLINE>
        [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
        [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
        [1 0 0 0 0 0], [-1  0  0  0  0  0]
        ]
        sage: f(O)  # long time
        Matrix group over Rational Field with 6 generators
        sage: f(O).gens()   # long time
        (
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
        [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
        [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
        <BLANKLINE>
        [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
        [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
        [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
        [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
        [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
        [1 0 0 0 0 0], [-1  0  0  0  0  0]
        )

    We check that :trac:`19780` is fixed::

        sage: G = groups.matrix.SO(3, 3)
        sage: H = groups.matrix.GL(3, 3)
        sage: phi = G.hom([H(x) for x in G.gens()])
        sage: phi(G.one()).parent()
        General Linear Group of degree 3 over Finite Field of size 3
    """
    pass

