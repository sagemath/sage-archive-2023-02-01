"""
Dynkin diagrams
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
from sage.graphs.all import DiGraph

def precheck(t, letter=None, length=None, affine=None, n_ge=None, n=None):
    """
    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import precheck
        sage: ct = CartanType(['A',4])
        sage: precheck(ct, letter='C')
        Traceback (most recent call last):
        ...
        ValueError: t[0] must be = 'C'
        sage: precheck(ct, affine=1)
        Traceback (most recent call last):
        ...
        ValueError: t[2] must be = 1
        sage: precheck(ct, length=3)
        Traceback (most recent call last):
        ...
        ValueError: len(t) must be = 3
        sage: precheck(ct, n=3)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be = 3
        sage: precheck(ct, n_ge=5)
        Traceback (most recent call last):
        ...
        ValueError: t[1] must be >= 5

    """
    if letter is not None:
        if t[0] != letter:
            raise ValueError, "t[0] must be = '%s'"%letter

    if length is not None:
        if len(t) != length:
            raise ValueError, "len(t) must be = %s"%length


    if affine is not None:
        try:
            if t[2] != affine:
                raise ValueError, "t[2] must be = %s"%affine
        except IndexError:
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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_a
        sage: ct = CartanType(['A',3])
        sage: a = type_a(ct);a
        Multi-digraph on 3 vertices
        sage: e = a.edges(); e.sort(); e
        [(1, 2, None), (2, 1, None), (2, 3, None), (3, 2, None)]

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
        sage: from sage.combinat.root_system.dynkin_diagram import type_a_affine
        sage: ct = CartanType(['A',3,1])
        sage: a = type_a_affine(ct); a
        Multi-digraph on 4 vertices
        sage: e = a.edges(); e.sort(); e
        [(0, 1, None),
         (0, 3, None),
         (1, 0, None),
         (1, 2, None),
         (2, 1, None),
         (2, 3, None),
         (3, 0, None),
         (3, 2, None)]
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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_b
        sage: ct = CartanType(['B',3])
        sage: b = type_b(ct);b
        Multi-digraph on 3 vertices
        sage: e = b.edges(); e.sort(); e
        [(1, 2, None), (2, 1, None), (2, 3, None), (3, 2, None), (3, 2, None)]


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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_c
        sage: ct = CartanType(['C',3])
        sage: c = type_c(ct);c
        Multi-digraph on 3 vertices
        sage: e = c.edges(); e.sort(); e
        [(1, 2, None), (2, 1, None), (2, 3, None), (2, 3, None), (3, 2, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_c_affine
        sage: ct = CartanType(['C',3,1])
        sage: c = type_c_affine(ct);c
        Multi-digraph on 4 vertices
        sage: e = c.edges(); e.sort(); e
        [(0, 1, None),
         (0, 1, None),
         (1, 0, None),
         (1, 2, None),
         (2, 1, None),
         (2, 3, None),
         (2, 3, None),
         (3, 2, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_d
        sage: ct = CartanType(['D',3])
        sage: d = type_d(ct);d
        Multi-digraph on 3 vertices
        sage: e = d.edges(); e.sort(); e
        [(1, 2, None), (1, 3, None), (2, 1, None), (3, 1, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_d_affine
        sage: ct = CartanType(['D',3,1])
        sage: d = type_d_affine(ct);d
        Multi-digraph on 4 vertices
        sage: e = d.edges(); e.sort(); e
        [(0, 2, None),
         (1, 2, None),
         (1, 3, None),
         (2, 0, None),
         (2, 1, None),
         (3, 1, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_e
        sage: ct = CartanType(['E',6])
        sage: e = type_e(ct);e
        Multi-digraph on 6 vertices
        sage: edges = e.edges(); edges.sort(); edges
        [(1, 3, None),
         (2, 4, None),
         (3, 1, None),
         (3, 4, None),
         (4, 2, None),
         (4, 3, None),
         (4, 5, None),
         (5, 4, None),
         (5, 6, None),
         (6, 5, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_f
        sage: ct = CartanType(['F',4])
        sage: f = type_f(ct);f
        Multi-digraph on 4 vertices
        sage: e = f.edges(); e.sort(); e
        [(1, 2, None),
         (2, 1, None),
         (2, 3, None),
         (3, 2, None),
         (3, 2, None),
         (3, 4, None),
         (4, 3, None)]

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

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_g
        sage: ct = CartanType(['G',2])
        sage: g = type_g(ct);g
        Multi-digraph on 2 vertices
        sage: e = g.edges(); e.sort(); e
        [(1, 2, None), (1, 2, None), (1, 2, None), (2, 1, None)]

    """
    precheck(t, letter='G', length=2, n=2)
    g = DiGraph(multiedges=True)
    for _ in range(3):
        g.add_edge(1,2)
    g.add_edge(2,1)
    return g

def dynkin_diagram(t):
    """
    Returns a DiGraph corresponding to the Dynkin diagram of
    type t.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import dynkin_diagram
        sage: dynkin_diagram(['A', 4])
        Multi-digraph on 4 vertices

    """
    import cartan_type
    f = "type_"
    letter = t[0].lower()
    affine = ""
    ct = cartan_type.CartanType(t)
    if ct.is_affine():
        affine = "_affine"

    function = globals()[f+letter+affine]
    try:
        return function(t)
    except KeyError:
        raise TypeError, "Dynkin diagram data not yet hardcoded for type %s"%t



def dynkin_diagram_as_function(t):
    """
    Returns a function corresponding to the Dynkin diagram of
    type t.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import dynkin_diagram_as_function
        sage: f = dynkin_diagram_as_function(['A',4])
        sage: f(1,1)
        -2
    """
    d = dynkin_diagram(t)
    f = lambda i,j: -2 if i == j else ( len(filter(lambda x: x[1] == j, d.outgoing_edges(i))) if d.has_edge(i,j) else 0)

    return f
