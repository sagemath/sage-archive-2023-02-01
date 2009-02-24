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
from root_system import RootSystem

def DynkinDiagram(*args):
    """
    INPUT:


    -  ``ct`` - A Cartan Type Returns a Dynkin diagram for
       type ct.


    The edge multiplicities are encoded as edge labels. This uses the
    convention in Kac / Fulton Harris Representation theory wikipedia
    http://en.wikipedia.org/wiki/Dynkin_diagram, that is for i != j::

       j --k--> i <==> a_ij = -k
                  <==> -scalar(coroot[i], root[j]) = k
                  <==> multiple arrows point from the longer root
                       to the shorter one

    TODO: say something about the node labelling conventions.

    EXAMPLES::

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

    if ct.is_affine():
        function = ct.tools.affine_dynkin_diagram
    else:
        function = ct.tools.dynkin_diagram

    try:
        return function(ct)
    except KeyError:
        raise TypeError, "Dynkin diagram data not yet hardcoded for type %s"%t

def dynkin_diagram(t):
    """
    Returns the Dynkin diagram of type t.

    Note that this function is deprecated, and that you should use
    DynkinDiagram instead as this will be disappearing in the near
    future.

    EXAMPLES::

        sage: dynkin_diagram(["A", 3])
        doctest:1: DeprecationWarning: dynkin_diagram is deprecated, use DynkinDiagram instead!
        Dynkin diagram of type ['A', 3]
    """
    from sage.misc.misc import deprecation
    deprecation("dynkin_diagram is deprecated, use DynkinDiagram instead!")
    return DynkinDiagram(t)


class DynkinDiagram_class(DiGraph, CartanType_abstract):
    def __init__(self, t):
        """
        EXAMPLES::

            sage: d = DynkinDiagram(["A", 3])
            sage: d == loads(dumps(d))
            True
        """
        DiGraph.__init__(self)
        self._cartan_type = t

    def __repr__(self):
        """
        EXAMPLES::

            sage: DynkinDiagram(['A',3])
            Dynkin diagram of type ['A', 3]
        """
        if self._cartan_type is None:
            return "Dynkin diagram of rank %s"%self.rank()
        else:
            return "Dynkin diagram of type %s"%self._cartan_type

    def add_edge(self, i, j, label=1):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
            sage: d = DynkinDiagram_class(CartanType(['A',3]))
            sage: list(sorted(d.edges()))
            []
            sage: d.add_edge(2, 3)
            sage: list(sorted(d.edges()))
            [(2, 3, 1), (3, 2, 1)]
        """
        DiGraph.add_edge(self, i, j, label)
        if not self.has_edge(j,i):
            self.add_edge(j,i,1)

    def index_set(self):
        """
        EXAMPLES::

            sage: DynkinDiagram(['C',3]).index_set()
            [1, 2, 3]
            sage: DynkinDiagram("A2","B2","F4").index_set()
            [1, 2, 3, 4, 5, 6, 7, 8]
        """
        return self.vertices()

    def cartan_type(self):
        """
        EXAMPLES::

            sage: DynkinDiagram("A2","B2","F4").cartan_type()
            A2xB2xF4
        """
        return self._cartan_type

    def rank(self):
        r"""
        Returns the index set for this Dynkin diagram

        EXAMPLES::

            sage: DynkinDiagram(['C',3]).rank()
            3
            sage: DynkinDiagram("A2","B2","F4").rank()
            8
        """
        return self.num_verts()

    def dynkin_diagram(self):
        """
        EXAMPLES::

            sage: DynkinDiagram(['C',3]).dynkin_diagram()
            Dynkin diagram of type ['C', 3]
        """
        return self

    def cartan_matrix(self):
        r"""
        returns the Cartan matrix for this Dynkin diagram

        EXAMPLES::

            sage: DynkinDiagram(['C',3]).cartan_matrix()
            [ 2 -1  0]
            [-1  2 -2]
            [ 0 -1  2]
        """
        return cartan_matrix(self)

    def dual(self):
        r"""
        Returns the dual Dynkin diagram, obtained by reversing all edges.

        EXAMPLES::

            sage: D = DynkinDiagram(['C',3])
            sage: D.edges()
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]
            sage: D.dual()
            Dynkin diagram of type ['B', 3]
            sage: D.dual().edges()
            [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]
            sage: D.dual() == DynkinDiagram(['B',3])
            True

        TESTS::

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
        r"""
        With a tuple (i,j) as argument, returns the scalar product
        `\langle
                \alpha^\vee_i, \alpha_j\rangle`.

        Otherwise, behaves as the usual DiGraph.__getitem__

        EXAMPLES: We use the `C_4` dynkin diagram as a cartan
        matrix::

            sage: g = DynkinDiagram(['C',4])
            sage: matrix([[g[i,j] for j in range(1,5)] for i in range(1,5)])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]

        The neighbors of a node can still be obtained in the usual way::

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
        Returns the `j`-th column `(a_{i,j})_i` of the
        Cartan matrix corresponding to this Dynkin diagram, as a container
        (or iterator) of tuples (i, a_i,j)

        EXAMPLES::

            sage: g = DynkinDiagram(["B",4])
            sage: [ (i,a) for (i,a) in g.column(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return [(j,2)] + [(i,-m) for (j1, i, m) in self.outgoing_edges(j)]

    def row(self, i):
        """
        Returns the `i`-th row `(a_{i,j})_j` of the
        Cartan matrix corresponding to this Dynkin diagram, as a container
        (or iterator) of tuples (j, a_i,j)

        EXAMPLES::

            sage: g = DynkinDiagram(["C",4])
            sage: [ (i,a) for (i,a) in g.row(3) ]
            [(3, 2), (2, -1), (4, -2)]
        """
        return [(i,2)] + [(j,-m) for (j, i1, m) in self.incoming_edges(i)]

def precheck(t, letter=None, length=None, affine=None, n_ge=None, n=None):
    """
    EXAMPLES::

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
