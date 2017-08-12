"""
Heisenberg Group

AUTHORS:

- Hilder Vitor Lima Pereira (2017-08): initial version
"""

#*****************************************************************************
#    Copyright (C) 2017 Hilder Vitor Lima Pereira <hilder.vitor at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.latex import latex
from sage.matrix.matrix_space import MatrixSpace
from sage.categories.groups import Groups
from sage.rings.all import ZZ

class HeisenbergGroup(UniqueRepresentation, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The Heisenberg group of degree `n`.

    Let `R` be the ring `\ZZ` or `\ZZ / p\ZZ` and let `n` be a positive
    integer. The Heisenberg group of degree `n` over `R` is
    a multiplicative group whose elements are matrices with the following form:

    .. MATH::

        \begin{pmatrix}
        1 & x & z \\
        0 & I_n & y \\
        0 & 0   & 1
        \end{pmatrix}

    where `x` is a row vector in `R^n`, `y` is a column vector in `R^n`,
    `z` is a scalar in `R`, and `I_n` is the identity matrix of size `n`.

    REFERENCES:

    - :wikipedia:`Heisenberg_group`
    """
    def __init__(self, p = 0, n = 1):
        """
        Initialize ``self``.

        INPUT:

        - ``p`` -- If it is a positive integer, then the underlying group is
        set to be the ring ``\ZZ / p\ZZ``. Otherwise, the set of integers
        is used as the underlying ring. The default value of ``p`` is 0.

        -``n`` -- Elements of the Heisenberg group are (n+2)x(n+2) matrices.
        Notice that if ``n`` equals 1 (which is the default value), the usal
        Heseinberg group of 3x3 matrices is constructed.


        EXAMPLES::
            sage: H = groups.matrix.Heisenberg(); H
            Heisenberg group of degree 1 (3x3 matrices) over Integer Ring
            sage: H.gens()
            (
            [1 1 0]  [1 0 0]  [1 0 1]
            [0 1 0]  [0 1 1]  [0 1 0]
            [0 0 1], [0 0 1], [0 0 1]
            )
            sage: X, Y, Z = H.gens()
            sage: Z*X*Y**-1
            [ 1  1  0]
            [ 0  1 -1]
            [ 0  0  1]
            sage: X * Y * X**-1 * Y**-1 == Z
            True
            sage: H = groups.matrix.Heisenberg(p=5); H
            Heisenberg group of degree 1 (3x3 matrices) over Ring of integers modulo 5
            sage: H = groups.matrix.Heisenberg(p=13, n=3); H
            Heisenberg group of degree 3 (5x5 matrices) over Ring of integers modulo 13
        """

        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ

        def elementary_matrix(i, j, dim, base_ring):
            elm = matrix.identity(base_ring, dim)
            elm[i,j] = elm[0,0] # put a one at position (i,j)
            return elm

        if n in ZZ and n > 0:
            self._n = n
        else:
            raise TypeError(
                    "Degree of Heisenberg group must be a positive integer.")

        if p in ZZ and p > 1:
            self._ring = ZZ.quo(p)
            self._p = p
        else:
            self._ring = ZZ
            self._p = 1

        R = self._ring
        self._dim = n+2
        MS = MatrixSpace(R, self._dim)
        gens_x = [elementary_matrix(0, j, self._dim, R) for j\
                                in range(1, self._dim-1) ]
        gens_y = [elementary_matrix(i, self._dim-1, self._dim, R) for i\
                                in range(1, self._dim-1)]
        gen_z = elementary_matrix(0, self._dim-1, self._dim, R)
        gens = gens_x + gens_y + [ gen_z ]

        from sage.libs.gap.libgap import libgap
        gap_gens = [libgap(single_gen) for single_gen in gens]
        gap_group = libgap.Group(gap_gens)

        if self._ring is ZZ:
            FinitelyGeneratedMatrixGroup_gap.__init__(self, ZZ(self._dim), R,
                    gap_group, category=Groups().FinitelyGenerated())
        else:
            FinitelyGeneratedMatrixGroup_gap.__init__(self, ZZ(self._dim), R,
                    gap_group, category=Groups().Finite())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: groups.matrix.Heisenberg()
            Heisenberg group of degree 1 (3x3 matrices) over Integer Ring
        """
        return "Heisenberg group of degree {} ({}x{} matrices) over {}"\
                    .format(self._n, self._dim, self._dim, self._ring)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg();
            sage: latex(H)
            H_{1}({\Bold{Z}})
        """
        return "H_{{{}}}({{{}}})".format(self._n, latex(self._ring))

    def order(self):
        """
        Return the order of ``self``.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg();
            sage: H.order()
            +Infinity
            sage: H = groups.matrix.Heisenberg(n=4);
            sage: H.order()
            +Infinity
            sage: H = groups.matrix.Heisenberg(p=3);
            sage: H.order()
            27
            sage: H = groups.matrix.Heisenberg(n=2, p=3); H.order()
            243
        """
        if self._ring is ZZ:
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            return ZZ(self._p ** (2*self._n + 1))

    cardinality = order

