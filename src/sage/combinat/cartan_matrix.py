#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from dynkin_diagram import dynkin_diagram, dynkin_diagram_as_function
import cartan_type
from sage.matrix.all import matrix, MatrixSpace
from sage.rings.all import ZZ

def cartan_matrix_as_function(t):
    """
    Returns a function that represents the Cartan matrix
    of type t.

    EXAMPLES:

    """
    ct = cartan_type.CartanType(t)
    f = dynkin_diagram_as_function(ct)
    s = ct.index_set()

    def cmf(i, j):
        if i == j:
            return 2
        else:
            return -f(j,i)

    return cmf

def cartan_matrix(t):
    """
    Returns the Cartan matrix corresponding to type t.

    EXAMPLES:
        sage: cartan_matrix(['A', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        sage: cartan_matrix(['B', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1  0]
        [ 0  0  0 -1  2 -2]
        [ 0  0  0  0 -1  2]
        sage: cartan_matrix(['C', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -2  2]
        sage: cartan_matrix(['D', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1 -1]
        [ 0  0  0 -1  2  0]
        [ 0  0  0 -1  0  2]
        sage: cartan_matrix(['E', 8])
        [ 2  0 -1  0  0  0  0  0]
        [ 0  2  0 -1  0  0  0  0]
        [-1  0  2 -1  0  0  0  0]
        [ 0 -1 -1  2 -1  0  0  0]
        [ 0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: cartan_matrix(['F', 4])
        [ 2 -1  0  0]
        [-1  2 -2  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        sage: cartan_matrix(['G', 2])
        [ 2 -1]
        [-3  2]

    """
    ct = cartan_type.CartanType(t)
    index_set = ct.index_set()
    cmf = cartan_matrix_as_function(ct)
    MS = MatrixSpace(ZZ, len(index_set), sparse=True)
    m = MS(0)
    for i in range(len(index_set)):
        for j in range(len(index_set)):
            m[i,j] = cmf(index_set[i],index_set[j])
    return m
