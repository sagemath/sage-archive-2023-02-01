"""
Homomorphisms Between Matrix Groups

AUTHORS:

- David Joyner and William Stein (2006-03): initial version

- David Joyner (2006-05): examples

- Simon King (2011-01): cleaning and improving code

- Volker Braun (2013-1): port to new Parent, libGAP.

- Simon Brandhorst (2018-5): deprecation
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
from sage.misc.superseded import deprecation

def to_libgap(x):
    """
    Helper to convert ``x`` to a LibGAP matrix or matrix group
    element.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.morphism import to_libgap
        doctest:...: DeprecationWarning: MatrixGroupMorphism_im_gens is deprecated. Use GroupMorphism_libgap instead.
        See https://trac.sagemath.org/25444 for details.
        sage: to_libgap(GL(2,3).gen(0))
        doctest:...: DeprecationWarning: this function is deprecated.Try x.gap() or libgap(x) instead
        See https://trac.sagemath.org/25444 for details.
        [ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        sage: to_libgap(matrix(QQ, [[1,2],[3,4]]))
        [ [ 1, 2 ], [ 3, 4 ] ]
    """
    deprecation(25444, "this function is deprecated."
                "Try x.gap() or libgap(x) instead")
    try:
        return x.gap()
    except AttributeError:
        from sage.libs.gap.libgap import libgap
        return libgap(x)


class MatrixGroupMap(Morphism):
    r"""
    This class is deprecated; use :class:`sage.groups.libgap_morphism.GroupMorphism_libgap` instead.
    """
    def _repr_type(self):
        """
        Part of the implementation of :meth:`_repr_`

        EXAMPLES::

            sage: from sage.groups.matrix_gps.morphism import MatrixGroupMap
            sage: MatrixGroupMap(ZZ.Hom(ZZ))._repr_type()
            ...
            doctest:...: DeprecationWarning: MatrixGroupMap is deprecated.Use GroupMorphism_libgap instead.
            See https://trac.sagemath.org/25444 for details.
            'MatrixGroup'

            'MatrixGroup'
        """
        deprecation(25444, "MatrixGroupMap is deprecated."
                "Use GroupMorphism_libgap instead.")
        return "MatrixGroup"

class MatrixGroupMorphism(MatrixGroupMap):
    r"""
    This class is deprecated; use :class:`sage.groups.libgap_morphism.GroupMorphism_libgap` instead.
    """
    deprecation(25444, "MatrixGroupMorphism is deprecated."
                "Use GroupMorphism_libgap instead.")


class MatrixGroupMorphism_im_gens(GroupMorphism_libgap):
    """
    This class is deprecated; use :class:`sage.groups.libgap_morphism.GroupMorphism_libgap` instead.

    Group morphism specified by images of generators.

    EXAMPLES::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: G.Hom(H)
        Set of Morphisms from Matrix group over Finite Field of size 5 with 1 generators (
        [1 1]
        [0 1]
        ) to Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [1 1]
        ) in Category of finite groups
    """
    deprecation(25444, "MatrixGroupMorphism_im_gens is deprecated. "
                "Use GroupMorphism_libgap instead.")

