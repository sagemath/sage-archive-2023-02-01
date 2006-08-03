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

        sage: R = PolynomialRing(QQ,9,'x')
        sage: A = matrix(R,3,R.gens())
        sage: A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
    """
    if not rings.is_Ring(R):
        return matrix_over_R(R,nrows)
    if isinstance(ncols, (list, tuple)):
        entries = ncols
        ncols = len(entries) / nrows
    (nrows, ncols) = (rings.Integer(nrows), rings.Integer(ncols))
    return matrix_space.MatrixSpace(R, nrows, ncols, sparse=sparse)(list(entries))


matrix = Matrix


def matrix_over_R(x, R):
    """
    Return the \sage matrix over $R$ obtained from x, if possible.
    """
    try:
        return x._matrix_(R)
    except AttributeError:
        raise TypeError, "No known way to create a matrix from %s"%x

