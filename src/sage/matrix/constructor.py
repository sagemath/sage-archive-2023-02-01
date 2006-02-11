"""
Constructor for making matrices
"""


#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.all as rings
import sage.matrix.matrix_space as matrix_space

def Matrix(R, nrows, ncols, entries=0, sparse=False):
    """
    Create a matrix.

    INPUT:
        R -- ring
        nrows -- int; number of rows
        ncols -- int; number of columns
        entries -- list; entries of the matrix
        sparse -- bool (default: False); whether or not to store matrices as sparse
    OUTPUT:
        a matrix

    EXAMPLES:
        sage: Matrix(RationalField(), 2, 2, [1,2,3,4])
        [1 2]
        [3 4]

        sage: Matrix(FiniteField(5), 2, 3, range(6))
        [0 1 2]
        [3 4 0]

        sage: Matrix(IntegerRing(), 10, 10, range(100)).parent()
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring

        sage: Matrix(IntegerRing(), 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring
    """
    if not rings.is_Ring(R):
        raise TypeError, "R (=%s) must be a ring."%R
    (nrows, ncols) = (int(nrows), int(ncols))
    return matrix_space.MatrixSpace(R, nrows, ncols, sparse=sparse)(entries)


