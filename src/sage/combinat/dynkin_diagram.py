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

from sage.graphs.all import DiGraph, Graph

def precheck(t, letter=None, length=None, affine=None, n_ge=None, n=None):
    if letter is not None:
        if t[0] != letter:
            raise ValueError, "t[0] must be = '%s'"%letter

    if length is not None:
        if len(t) != length:
            raise ValueError, "len(t) must be = %s"%length


    if affine is not None:
        if t[2] != affine:
            raise ValueError, "t[2] must be = %s"%affine


    if n_ge is not None:
        if t[1] < n_ge:
            raise ValueError, "t[1] must be >= %s"%n_ge

    if n is not None:
        if t[1] != n:
            raise ValueError, "t[1] must be = %s"%n

def type_a(t):
    """
    Returns the graph corresponding to the Dynkin diagram
    of type A.


    """
    precheck(t, letter="A", length=2)
    n = t[1]
    g = DiGraph(multiedges=True)
    for i in range(1, n):
        g.add_edge(i, i+1)
        g.add_edge(i+1, i)

    return g

def type_a_affine(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of affine type A.

    EXAMPLES:
        sage: T = ['A', 1, 1]
        sage: cartan_matrix(T)
        [ 2 -2]
        [-2  2]

        sage: T = ['A', 2, 1]
        sage: cartan_matrix(T)
        [ 2 -1 -1]
        [-1  2 -1]
        [-1 -1  2]
    """
    precheck(t, letter="A", length=3, affine=1)
    n = t[1]
    g = type_a(['A', n])

    g.add_edge(0, 1)
    g.add_edge(1, 0)
    g.add_edge(0, n)
    g.add_edge(n, 0)

    return g

def type_b(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type B.

    TESTS:
        sage: T = ['B', 3]
        sage: cartan_matrix(T)
        [ 2 -1  0]
        [-1  2 -2]
        [ 0 -1  2]


    """
    precheck(t, letter='B', length=2, n_ge=2)
    n = t[1]
    g = type_a(['A', n])
    g.add_edge(n, n-1)
    return g


def type_c(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type C.

    TESTS:
        sage: T = ['C', 3]
        sage: cartan_matrix(T)
        [ 2 -1  0]
        [-1  2 -1]
        [ 0 -2  2]
    """
    precheck(t, letter='C', length=2, n_ge=2)
    n = t[1]
    g = type_a(['A', n])
    g.add_edge(n-1,n)
    return g

def type_c_affine(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of affine type C.

    TESTS:
        sage: T = ['C', 3, 1]
        sage: cartan_matrix(T)
        [ 2 -1  0  0]
        [-2  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -2  2]
    """
    precheck(t, letter='C', length=3, affine=1)
    n = t[1]
    g = type_c(['C', n])
    g.add_edge(0,1)
    g.add_edge(0,1)
    g.add_edge(1,0)
    return g

def type_d(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type D.

    TESTS:
        sage: cartan_matrix(['D', 3])
        [ 2 -1 -1]
        [-1  2  0]
        [-1  0  2]
        sage: cartan_matrix(['D', 4])
        [ 2 -1  0  0]
        [-1  2 -1 -1]
        [ 0 -1  2  0]
        [ 0 -1  0  2]

    """
    precheck(t, letter="D", length=2, n_ge=3)
    n = t[1]
    g = type_a(['A', n-1])
    g.add_edge(n-2,n)
    g.add_edge(n, n-2)
    return g

def type_d_affine(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type D.

    TESTS:
        sage: cartan_matrix(['D',3,1])
        [ 2  0 -1  0]
        [ 0  2 -1 -1]
        [-1 -1  2  0]
        [ 0 -1  0  2]

        sage: cartan_matrix(['D',4,1])
        [ 2  0 -1  0  0]
        [ 0  2 -1  0  0]
        [-1 -1  2 -1 -1]
        [ 0  0 -1  2  0]
        [ 0  0 -1  0  2]
    """
    precheck(t, letter="D", length=3, affine=1)
    n = t[1]
    g = type_d(['D', n])
    g.add_edge(0,2)
    g.add_edge(2,0)
    return g

def type_e(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type E.

    TESTS:
        sage: cartan_matrix(['E',6])
        [ 2  0 -1  0  0  0]
        [ 0  2  0 -1  0  0]
        [-1  0  2 -1  0  0]
        [ 0 -1 -1  2 -1  0]
        [ 0  0  0 -1  2 -1]
        [ 0  0  0  0 -1  2]
    """
    precheck(t, letter="E", length=2, n_ge=3)
    n = t[1]
    g = DiGraph(multiedges=True)
    g.add_edge(1,3)
    g.add_edge(3,1)
    g.add_edge(2,4)
    g.add_edge(4,2)
    for i in range(3,n):
        g.add_edge(i, i+1)
        g.add_edge(i+1, i)
    return g

def type_f(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type F.

    TESTS:
        sage: cartan_matrix(['F',4])
        [ 2 -1  0  0]
        [-1  2 -2  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]

    """
    precheck(t, letter='F', length=2, n=4)
    g = DiGraph(multiedges=True)
    for i in range(1, 4):
        g.add_edge(i, i+1)
        g.add_edge(i+1, i)
    g.add_edge(3,2)
    return g

def type_g(t):
    """
    Returns the DiGraph corresponding to the Dynkin diagram
    of type G.

    TESTS:
        sage: cartan_matrix(['G',2])
        [ 2 -1]
        [-3  2]
    """
    precheck(t, letter='G', length=2, n=2)
    g = DiGraph(multiedges=True)
    for i in range(3):
        g.add_edge(1,2)
    g.add_edge(2,1)
    return g

def dynkin_diagram(t):
    """
    Returns a DiGraph corresponding to the Dynkin diagram of
    type t.

    EXAMPLES:

    """
    import cartan_type
    f = "type_"
    letter = t[0].lower()
    affine = ""
    ct = cartan_type.CartanType(t)
    if ct.is_affine():
        affine = "_affine"

    function = eval(f+letter+affine)
    try:
        return function(t)
    except RuntimeError:
        raise TypeError, "Dynkin diagram data not yet hardcoded for type %s"%t



def dynkin_diagram_as_function(t):
    """
    Returns a function corresponding to the Dynkin diagram of
    type t.

    EXAMPLES:

    """
    d = dynkin_diagram(t)
    def f(i, j):
        if i == j:
            return -2
        elif d.has_edge(i,j):
            return len(filter(lambda x: x[1] == j, d.outgoing_edges(i)))
        else:
            return 0

    return f
