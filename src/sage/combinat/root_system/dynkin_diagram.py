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
from cartan_type import CartanType, CartanType_abstract
from cartan_matrix import cartan_matrix
import cartan_type
from root_system import RootSystem

def DynkinDiagram(*args):
    """
    INPUT:
       ct -- A Cartan Type
    Returns a Dynkin diagram for type ct.

    The edge multiplicities are encoded as edge labels. This uses the
    convention in Kac / Fulton Harris Representation theory wikipedia
    http://en.wikipedia.org/wiki/Dynkin_diagram, that is for i != j:

    j -k-> i  <==>  a_ij = -k  <==>  -scalar(coroot[i], root[j]) = k
    <==>  multiple arrows point from the longer root to the shorter one

    TODO: say something about the node labelling conventions.

    EXAMPLES:
        sage: DynkinDiagram(['A', 4])
        Dynkin diagram of type ['A', 4]
        sage: DynkinDiagram(['A',1],['A',1])
        Dynkin diagram of type A1xA1
        sage: R = RootSystem("A2xB2xF4")
        sage: DynkinDiagram(R)
        Dynkin diagram of type A2xB2xF4
    """
    if len(args) == 1:
        t = args[0]
    else:
        t = CartanType(args)

    if isinstance(t, RootSystem):
        ct = t.cartan_type()
    else:
        ct = CartanType(t)
    if ct.is_reducible():
        function = globals()["type_reducible"]
    else:
        letter = ct[0].lower()
        affine = ""
        ct = CartanType(ct)
        if ct.is_affine():
            affine = "_affine"
        function = globals()["type_"+letter+affine]
    try:
        return function(ct)
    except KeyError:
        raise TypeError, "Dynkin diagram data not yet hardcoded for type %s"%t

def dynkin_diagram(t):
    """
    Deprecated; please use DynkinDiagram
    """
    return DynkinDiagram(t)


class DynkinDiagram_class(DiGraph, CartanType_abstract):
    def __init__(self, t):
        DiGraph.__init__(self)
        self._cartan_type = t

    def __repr__(self):
        """
        EXAMPLES:
          sage: DynkinDiagram(['A',3])
          Dynkin diagram of type ['A', 3]
        """
        if self._cartan_type is None:
            return "Dynkin diagram of rank %s"%self.rank()
        else:
            return "Dynkin diagram of type %s"%self._cartan_type

    def add_edge(self, i, j, label=1):
        DiGraph.add_edge(self, i, j, label)
        if not self.has_edge(j,i):
            self.add_edge(j,i,1)

    def index_set(self):
        """
        EXAMPLES:
          sage: DynkinDiagram(['C',3]).index_set()
          [1, 2, 3]
          sage: DynkinDiagram("A2","B2","F4").index_set()
          [1, 2, 3, 4, 5, 6, 7, 8]
        """
        return self.vertices()

    def cartan_type(self):
        """
        EXAMPLES:
          sage: DynkinDiagram("A2","B2","F4").cartan_type()
          A2xB2xF4
        """
        return self._cartan_type

    def rank(self):
        r"""
        returns the index set for this Dynkin diagram

        EXAMPLES:
        """
        return self.num_verts()

    def dynkin_diagram(self):
        return self

    def cartan_matrix(self):
        r"""
        returns the Cartan matrix for this Dynkin diagram

        EXAMPLES:
            sage: DynkinDiagram(['C',3]).cartan_matrix()
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return cartan_matrix(self)

    def dual(self):
        r"""
        Returns the dual Dynkin diagram, obtained by reversing all edges.

        EXAMPLES:
            sage: D = DynkinDiagram(['C',3])
            sage: D.edges()
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]
            sage: D.dual()
            Dynkin diagram of type ['B', 3]
            sage: D.dual().edges()
            [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]
            sage: D.dual() == DynkinDiagram(['B',3])
            True

        TESTS:
            sage: D = DynkinDiagram(['A',0]); D
            Dynkin diagram of type ['A', 0]
            sage: D.edges()
            []
            sage: D.dual()
            Dynkin diagram of type ['A', 0]
            sage: D.dual().edges()
            []
            sage: D = DynkinDiagram(['A',1])
            sage: D.edges()
            []
            sage: D.dual()
            Dynkin diagram of type ['A', 1]
            sage: D.dual().edges()
            []

        """
        result = DynkinDiagram_class(None)
        result.add_vertices(self.vertices())
        for source, target, label in self.edges():
            result.add_edge(target, source, label)
        result._cartan_type = self._cartan_type.dual() if not self._cartan_type is None else None
        return result

    def __getitem__(self, i):
        """
        With a tuple (i,j) as argument, returns the scalar product $\langle
        \alpha^\vee_i, \alpha_j\rangle$.

        Otherwise, behaves as the usual DiGraph.__getitem__

        EXAMPLES:
        We use the $C_4$ dynkin diagram as a cartan matrix:
            sage: g = DynkinDiagram(['C',4])
            sage: matrix([[g[i,j] for j in range(1,5)] for i in range(1,5)])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]

        The neighbors of a node can still be obtained in the usual way:
            sage: [g[i] for i in range(1,5)]
            [[2], [1, 3], [2, 4], [3]]
        """
        if not isinstance(i, tuple):
            return DiGraph.__getitem__(self,i)
        [i,j] = i
        if i == j:
            return 2
        elif self.has_edge(j, i):
            return -self.edge_label(j, i)
        else:
            return 0

    def column(self, j):
        """
        Returns the $j$-th column $(a_{i,j})_i$ of the Cartan matrix
        corresponding to this Dynkin diagram, as a container (or iterator)
        of tuples (i, a_{i,j})

        EXAMPLES:
            sage: g = DynkinDiagram(["B",4])
            sage: [ (i,a) for (i,a) in g.column(3) ]
            [(3, 2), (2, -1), (4, -2)]

        Caveat: broken in sage < 3.0.3
        """
        return [(j,2)] + [(i,-m) for (j1, i, m) in self.outgoing_edges(j)]

    def row(self, i):
        """
        Returns the $i$-th row $(a_{i,j})_j$ of the Cartan matrix
        corresponding to this Dynkin diagram, as a container (or iterator)
        of tuples (j, a_{i,j})

        EXAMPLES:
            sage: g = DynkinDiagram(["C",4])
            sage: [ (i,a) for (i,a) in g.row(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return [(i,2)] + [(j,-m) for (j, i1, m) in self.incoming_edges(i)]

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

##############################################################################
# Everything below is type by type hardcoded data. It probably should be moved
# into the type_... files
##############################################################################


def type_a(t):
    """
    Returns the graph corresponding to the Dynkin diagram
    of type A.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_a
        sage: ct = CartanType(['A',3])
        sage: a = type_a(ct); a
        Dynkin diagram of type ['A', 3]
        sage: e = a.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 1)]

    TEST:
        sage: a = DynkinDiagram(['A',1])
        sage: a
        Dynkin diagram of type ['A', 1]
        sage: a.vertices(), a.edges()
        ([1], [])
    """
    precheck(t, letter="A", length=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    return g

def type_a_affine(t):
    """
    Returns the extended Dynkin diagram for affine type A.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_a_affine
        sage: ct = CartanType(['A',3,1])
        sage: a = type_a_affine(ct); a
        Dynkin diagram of type ['A', 3, 1]
        sage: e = a.edges(); e.sort(); e
        [(0, 1, 1),
         (0, 3, 1),
         (1, 0, 1),
         (1, 2, 1),
         (2, 1, 1),
         (2, 3, 1),
         (3, 0, 1),
         (3, 2, 1)]
    """
    precheck(t, letter="A", length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.add_edge(0, 1)
    g.add_edge(0, n)

    return g

def type_b(t):
    """
    Returns a Dynkin diagram for type B.

    EXAMPLES:
         sage: from sage.combinat.root_system.dynkin_diagram import type_b
         sage: ct = CartanType(['B',3])
         sage: b = type_b(ct);b
         Dynkin diagram of type ['B', 3]
         sage: e = b.edges(); e.sort(); e
         [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

    """
    precheck(t, letter='B', length=2, n_ge=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n-1, n, 2)
    return g

def type_b_affine(t):
    """
    Returns the extended Dynkin diagram for affine type B.

    EXAMPLES:
       sage: DynkinDiagram(['B',3,1]).edges()
       [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

    """
    precheck(t, letter='B', length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n-1, n, 2)
    g.add_edge(0,2)
    return g

def type_c(t):
    """
    Returns a Dynkin diagram for type C.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_c
        sage: ct = CartanType(['C',3])
        sage: c = type_c(ct); c
        Dynkin diagram of type ['C', 3]
        sage: e = c.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

    """
    precheck(t, letter='C', length=2, n_ge=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n,n-1,2)
    return g

def type_c_affine(t):
    """
    Returns the extended Dynkin diagram for affine type C.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_c_affine
        sage: ct = CartanType(['C',3,1])
        sage: c = type_c_affine(ct);c
        Dynkin diagram of type ['C', 3, 1]
        sage: e = c.edges(); e.sort(); e
        [(0, 1, 2), (1, 0, 1), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

    """
    precheck(t, letter='C', length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n,n-1,2)
    g.add_edge(0,1,1)
    g.add_edge(0,1,2)
    return g

def type_d(t):
    """
    Returns a Dynkin diagram for type D.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_d
        sage: ct = CartanType(['D',4])
        sage: d = type_d(ct);d
        Dynkin diagram of type ['D', 4]
        sage: e = d.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 1), (2, 4, 1), (3, 2, 1), (4, 2, 1)]

    """
    precheck(t, letter="D", length=2, n_ge=3)
    n = t[1]

    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n-1):
        g.add_edge(i, i+1)
    g.add_edge(n-2,n)
    return g

def type_d_affine(t):
    """
    Returns the extended Dynkin diagram for affine type D.

    EXAMPLES:
       sage: DynkinDiagram(CartanType(['D', 4, 1]))
       Dynkin diagram of type ['D', 4, 1]
       sage: DynkinDiagram(CartanType(['D', 4, 1])).edges()
       [(0, 2, 1),
        (1, 2, 1),
        (2, 0, 1),
        (2, 1, 1),
        (2, 3, 1),
        (2, 4, 1),
        (3, 2, 1),
        (4, 2, 1)]

    """
    precheck(t, letter="D", length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n-1):
        g.add_edge(i, i+1)
    g.add_edge(n-2,n)
    g.add_edge(0,2)
    return g

def type_e(t):
    """
    Returns a Dynkin diagram for type E.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_e
        sage: ct = CartanType(['E',6])
        sage: e = type_e(ct);e
        Dynkin diagram of type ['E', 6]
        sage: edges = e.edges(); edges.sort(); edges
        [(1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1), (4, 2, 1), (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1), (6, 5, 1)]

    """
    precheck(t, letter="E", length=2, n_ge=3)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_edge(1,3)
    g.add_edge(2,4)
    for i in range(3,n):
        g.add_edge(i, i+1)
    return g

def type_e_affine(t):
    """
    Returns the extended Dynkin diagram for affine type E.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_e_affine
        sage: e = DynkinDiagram(['E', 6, 1])
        sage: edges = e.edges(); edges.sort(); edges
        [(0, 2, 1),
         (1, 3, 1),
         (2, 0, 1),
         (2, 4, 1),
         (3, 1, 1),
         (3, 4, 1),
         (4, 2, 1),
         (4, 3, 1),
         (4, 5, 1),
         (5, 4, 1),
         (5, 6, 1),
         (6, 5, 1)]

    """
    precheck(t, letter="E", length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_edge(1,3)
    g.add_edge(2,4)
    for i in range(3,n):
        g.add_edge(i, i+1)
    if n == 6:
        g.add_edge(0, 2)
    elif n == 7:
        g.add_edge(0, 1)
    elif n == 8:
        g.add_edge(0, 8)
    else:
        raise ValueError, "Invalid Cartan Type for Type E affine"
    return g

def type_f(t):
    """
    Returns a Dynkin diagram for type F.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_f
        sage: ct = CartanType(['F',4])
        sage: f = type_f(ct);f
        Dynkin diagram of type ['F', 4]
        sage: e = f.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 3, 1)]

    """
    precheck(t, letter='F', length=2, n=4)
    g = DynkinDiagram_class(t)
    for i in range(1, 4):
        g.add_edge(i, i+1)
    g.set_edge_label(2,3,2)
    return g

def type_f_affine(t):
    """
    Returns the extended Dynkin diagram for affine type F.

    EXAMPLES:
        sage: f = DynkinDiagram(['F', 4, 1])
        sage: edges = f.edges(); edges.sort(); edges
        [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 3, 1)]

    """
    precheck(t, letter="F", length=3, affine=1)
    g = DynkinDiagram_class(t)
    for i in range(1, 4):
        g.add_edge(i, i+1)
    g.set_edge_label(2,3,2)
    g.add_edge(0, 1)
    return g

def type_g(t):
    """
    Returns a Dynkin diagram for type G.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_g
        sage: ct = CartanType(['G',2])
        sage: g = type_g(ct);g
        Dynkin diagram of type ['G', 2]
        sage: e = g.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 3)]

    """
    precheck(t, letter='G', length=2, n=2)
    g = DynkinDiagram_class(t)
    g.add_edge(1,2)
    g.set_edge_label(2,1,3)
    return g

def type_g_affine(t):
    """
    Returns the extended Dynkin diagram for type G.

    EXAMPLES:
        sage: from sage.combinat.root_system.dynkin_diagram import type_g_affine
        sage: ct = CartanType(['G',2,1])
        sage: g = type_g_affine(ct); g
        Dynkin diagram of type ['G', 2, 1]
        sage: e = g.edges(); e.sort(); e
        [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 3)]

    """
    precheck(t, letter="G", length=3, affine=1)
    g = DynkinDiagram_class(t)
    g.add_edge(1, 2)
    g.set_edge_label(2,1,3)
    g.add_edge(0, 2)
    return g

def type_reducible(t):
    """
    Returns a Dynkin diagram for type reducible.
    EXAMPLES:
      sage: t = CartanType("A2xB2xF4")
      sage: dd = DynkinDiagram(t); dd
      Dynkin diagram of type A2xB2xF4
      sage: dd.edges()
      [(1, 2, 1), (2, 1, 1), (3, 4, 2), (4, 3, 1), (5, 6, 1), (6, 5, 1), (6, 7, 2), (7, 6, 1), (7, 8, 1), (8, 7, 1)]

    """
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(len(t._types)):
        for [e1, e2, l] in DynkinDiagram(t._types[i]).edges():
            shift = t._rshifts[i]
            g.add_edge(e1+shift, e2+shift, label=l)
    return g


