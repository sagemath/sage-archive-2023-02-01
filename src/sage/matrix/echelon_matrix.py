r"""
Echelon matrices
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def echelon_matrix_iterator(K,n,k,sparse=True):
    r"""
    An iterator over `(k,n)` echelon matrices over `K`.

    INPUT:

    - ``K`` -- a finite ring

    - ``n`` -- number of rows

    - ``k`` -- number of columns

    - ``sparse`` -- boolean

    EXAMPLES::

        sage: from sage.matrix.echelon_matrix import echelon_matrix_iterator
        sage: it = echelon_matrix_iterator(GF(2),3,2)
        sage: for m in it: print m,"\n*******"
        [1 0 0]
        [0 1 0]
        *******
        [1 0 0]
        [0 1 1]
        *******
        [1 0 1]
        [0 1 0]
        *******
        [1 0 1]
        [0 1 1]
        *******
        [1 0 0]
        [0 0 1]
        *******
        [1 1 0]
        [0 0 1]
        *******
        [0 1 0]
        [0 0 1]
        *******
    """
    assert n >= k >= 0

    from sage.matrix.matrix_space import MatrixSpace
    from itertools import combinations,product

    M = MatrixSpace(K, k, n, sparse=sparse)
    Klist = list(K)

    # First, we select which columns will be pivots:
    for pivots in combinations(range(n),k):
        m = M()
        free_positions = []
        for i in range(k):
            m[i, pivots[i]] = 1
            for j in range(pivots[i]+1,n):
                if j not in pivots:
                    free_positions.append((i,j))
        # Next, we fill in those entries that are not
        # determined by the echelon form alone:
        num_free_pos = len(free_positions)
        for v in product(Klist, repeat=num_free_pos):
            for i in range(num_free_pos):
                m[free_positions[i]] = v[i]
            # Finally, we have to multiply by the basis matrix
            # to take corresponding linear combinations of the basis
            yield m
