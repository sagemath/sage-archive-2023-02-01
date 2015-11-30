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

from sage.groups.group_homset import GroupHomset_generic


def is_MatrixGroupHomset(x):
    r"""
    Test whether ``x`` is a homset.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.homset import is_MatrixGroupHomset
        sage: is_MatrixGroupHomset(4)
        False

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens)
        sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
        sage: M = MatrixGroupHomset(G, G)
        sage: is_MatrixGroupHomset(M)
        True
    """
    return isinstance(x, MatrixGroupHomset)


class MatrixGroupHomset(GroupHomset_generic):

    def __init__(self, G, H, category=None):
        r"""
        Return the homset of two matrix groups.

        INPUT:

        - ``G`` -- a matrix group

        - ``H`` -- a matrix group

        OUTPUT:

        The homset of two matrix groups.

        EXAMPLES::

            sage: F = GF(5)
            sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
            sage: MatrixGroupHomset(G, G)
            Set of Homomorphisms from
            Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            ) to Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            )
        """
        from sage.categories.homset import HomsetWithBase
        HomsetWithBase.__init__(self, G, H, category, G.base_ring())

    def __call__(self, im_gens, check=True):
        """
        Return the homomorphism defined by images of generators.

        INPUT:

        - ``im_gens`` -- iterable, the list of images of the
          generators of the domain

        - ``check`` -- bool (optional, default: ``True``), whether to
          check if images define a valid homomorphism

        OUTPUT:

        Group homomorphism.

        EXAMPLES::

            sage: F = GF(5)
            sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
            sage: M = MatrixGroupHomset(G, G)
            sage: M(gens)
            Homomorphism : Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            ) --> Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            )
        """
        from sage.groups.matrix_gps.morphism import MatrixGroupMorphism_im_gens
        try:
            return MatrixGroupMorphism_im_gens(self, im_gens, check=check)
        except (NotImplementedError, ValueError) as err:
            raise TypeError('images do not define a group homomorphism')


