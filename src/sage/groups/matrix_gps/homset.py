"""
Matrix Group Homsets

AUTHORS:

- William Stein (2006-05-07): initial version

- Volker Braun (2013-1) port to new Parent, libGAP
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.groups.libgap_morphism import GroupHomset_libgap
from sage.misc.superseded import deprecation

def is_MatrixGroupHomset(x):
    r"""
    Test whether ``x`` is a matrix group homset.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.homset import is_MatrixGroupHomset
        doctest:...: DeprecationWarning: MatrixGroupHomset is deprecated. Use GroupHomset_libgap instead.
        See https://trac.sagemath.org/25444 for details.
        sage: is_MatrixGroupHomset(4)
        doctest:...: DeprecationWarning: MatrixGroupHomset is deprecated. Use GroupHomset_libgap instead.
        See https://trac.sagemath.org/25444 for details.
        False

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens)
        sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
        sage: M = MatrixGroupHomset(G, G)
        sage: is_MatrixGroupHomset(M)
        True
    """
    deprecation(25444, "MatrixGroupHomset is deprecated. "
                "Use GroupHomset_libgap instead.")
    return isinstance(x, MatrixGroupHomset)


class MatrixGroupHomset(GroupHomset_libgap):
    r"""
    Homsets of matrix groups.

    EXAMPLES::

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens)
        sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
        sage: MatrixGroupHomset(G, G)
        Set of Morphisms from Matrix group over Finite Field of size 5 with 2 generators (
        [1 2]  [1 1]
        [4 1], [0 1]
        ) to Matrix group over Finite Field of size 5 with 2 generators (
        [1 2]  [1 1]
        [4 1], [0 1]
        ) in Category of finite groups
        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens)
        sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
        sage: M = MatrixGroupHomset(G, G)
        sage: M(gens)
        Group endomorphism of Matrix group over Finite Field of size 5 with 2 generators (
        [1 2]  [1 1]
        [4 1], [0 1]
        )
    """
    deprecation(25444, "MatrixGroupHomset is deprecated. "
                "Use GroupHomset_libgap instead.")
