"""
Symplectic Linear Groups

EXAMPLES::

    sage: G = Sp(4,GF(7));  G
    Symplectic Group of degree 4 over Finite Field of size 7
    sage: g = prod(G.gens());  g
    [3 0 3 0]
    [1 0 0 0]
    [0 1 0 1]
    [0 2 0 0]
    sage: m = g.matrix()
    sage: m * G.invariant_form() * m.transpose() == G.invariant_form()
    True
    sage: G.order()
    276595200

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, NamedMatrixGroup_generic, NamedMatrixGroup_gap )



###############################################################################
# Symplectic Group
###############################################################################

def Sp(n, R, var='a'):
    r"""
    Return the symplectic group.

    The special linear group `GL( d, R )` consists of all `d \times d`
    matrices that are invertible over the ring `R` with determinant
    one.

    .. note::

        This group is also available via ``groups.matrix.Sp()``.

    INPUT:

    - ``n`` -- a positive integer.

    - ``R`` -- ring or an integer. If an integer is specified, the
      corresponding finite field is used.

    - ``var`` -- variable used to represent generator of the finite
      field, if needed.

    EXAMPLES::

        sage: Sp(4, 5)
        Symplectic Group of degree 4 over Finite Field of size 5

        sage: Sp(4, IntegerModRing(15))
        Symplectic Group of degree 4 over Ring of integers modulo 15

        sage: Sp(3, GF(7))
        Traceback (most recent call last):
        ...
        ValueError: the degree must be even

    TESTS::

        sage: groups.matrix.Sp(2, 3)
        Symplectic Group of degree 2 over Finite Field of size 3

        sage: G = Sp(4,5)
        sage: TestSuite(G).run()
    """
    degree, ring = normalize_args_vectorspace(n, R, var=var)
    if degree % 2 != 0:
        raise ValueError('the degree must be even')
    name = 'Symplectic Group of degree {0} over {1}'.format(degree, ring)
    ltx  = r'\text{{Sp}}_{{{0}}}({1})'.format(degree, latex(ring))
    from sage.libs.gap.libgap import libgap
    try:
        cmd  = 'Sp({0}, {1})'.format(degree, ring._gap_init_())
        return SymplecticMatrixGroup_gap(degree, ring, True, name, ltx, cmd)
    except ValueError:
        return SymplecticMatrixGroup_generic(degree, ring, True, name, ltx)



class SymplecticMatrixGroup_generic(NamedMatrixGroup_generic):

    @cached_method
    def invariant_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: Sp(4, QQ).invariant_form()
            [0 0 0 1]
            [0 0 1 0]
            [0 1 0 0]
            [1 0 0 0]
        """
        from sage.matrix.constructor import zero_matrix
        m = zero_matrix(self.base_ring(), self.degree())
        for i in range(self.degree()):
            m[i, self.degree()-i-1] = 1
        m.set_immutable()
        return m

    def _check_matrix(self, x, *args):
        """
        Check whether the matrix ``x`` is symplectic.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = Sp(4,GF(5))
            sage: G._check_matrix(G.an_element().matrix())
        """
        F = self.invariant_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be symplectic')


class SymplecticMatrixGroup_gap(SymplecticMatrixGroup_generic, NamedMatrixGroup_gap):
    r"""
    Symplectic group in GAP

    EXAMPLES::

        sage: Sp(2,4)
        Symplectic Group of degree 2 over Finite Field in a of size 2^2

        sage: latex(Sp(4,5))
        \text{Sp}_{4}(\Bold{F}_{5})
    """

    @cached_method
    def invariant_form(self):
        """
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: Sp(4, GF(3)).invariant_form()
            [0 0 0 1]
            [0 0 1 0]
            [0 2 0 0]
            [2 0 0 0]
        """
        m = self.gap().InvariantBilinearForm()['matrix'].matrix()
        m.set_immutable()
        return m





