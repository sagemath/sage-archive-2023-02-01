"""
Coxeter matrices
"""
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
from cartan_type import CartanType
from dynkin_diagram import dynkin_diagram_as_function
from sage.matrix.all import MatrixSpace
from sage.rings.all import ZZ

def coxeter_matrix_as_function(t):
    """
    Returns the coxeter matrix associated to the Cartan type t.

    EXAMPLES:
        sage: from sage.combinat.root_system.coxeter_matrix import coxeter_matrix_as_function
        sage: f = coxeter_matrix_as_function(['A',4])
        sage: matrix([[f(i,j) for j in range(1,5)] for i in range(1,5)])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
    """
    ct = CartanType(t)
    f = dynkin_diagram_as_function(ct)

    if ct.letter == "G":
        coxeter = lambda i,j: 1 if i == j else max(f(j,i), f(i,j))+4
    else:
        coxeter = lambda i,j: 1 if i == j else max(f(j,i), f(i,j))+2

    return coxeter


def coxeter_matrix(t):
    """
    Returns the Coxeter matrix of type t.

    EXAMPLES:
        sage: coxeter_matrix(['A', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
        sage: coxeter_matrix(['B', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: coxeter_matrix(['C', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: coxeter_matrix(['D', 4])
        [1 3 2 2]
        [3 1 3 3]
        [2 3 1 2]
        [2 3 2 1]

        sage: coxeter_matrix(['E', 6])
        [1 2 3 2 2 2]
        [2 1 2 3 2 2]
        [3 2 1 3 2 2]
        [2 3 3 1 3 2]
        [2 2 2 3 1 3]
        [2 2 2 2 3 1]

        sage: coxeter_matrix(['F', 4])
        [1 3 2 2]
        [3 1 4 2]
        [2 4 1 3]
        [2 2 3 1]

        sage: coxeter_matrix(['G', 2])
        [1 7]
        [7 1]

    """
    ct = CartanType(t)
    cf = coxeter_matrix_as_function(ct)
    index_set = ct.index_set()
    MS = MatrixSpace(ZZ, len(index_set))
    m = MS(0)
    for i in range(len(index_set)):
        for j in range(len(index_set)):
            m[i,j] = cf(index_set[i],index_set[j])
    return m
