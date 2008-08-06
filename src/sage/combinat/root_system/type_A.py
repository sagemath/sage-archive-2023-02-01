from ambient_space import AmbientSpace
from sage.rings.all import ZZ
from sage.combinat.family import Family

class ambient_space(AmbientSpace):
    r"""
    EXAMPLES:
        sage: R = RootSystem(["A",3])
        sage: e = R.ambient_space(); e
        Ambient space of the Root system of type ['A', 3]
        sage: e == loads(dumps(e))
        True
    """
    @classmethod
    def smallest_base_ring(cls):
        """
        EXAMPLES:
            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.smallest_base_ring()
            Integer Ring
        """
        return ZZ

    def dimension(self):
        """
        EXAMPLES:
            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.dimension()
            4
        """
        return self.root_system.cartan_type().rank()+1

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.root(0,1)
            (1, -1, 0, 0)
        """
        return self._term(i) - self._term(j)

    def simple_root(self, i):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}
        """
        return self.root(i-1, i)

    def negative_roots(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 1, 0),
             (0, -1, 0, 1),
             (0, 0, -1, 1)]
        """
        res = []
        for j in range(self.n-1):
            for i in range(j+1,self.n):
                res.append(  self.root(i,j) )
        return res

    def positive_roots(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.positive_roots()
            [(1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]

        """
        res = []
        for j in range(self.n):
            for i in range(j):
                res.append(  self.root(i,j) )
        return res

    def highest_root(self):
        """
        EXAMPLE:
           sage: e = RootSystem(['A',3]).ambient_lattice()
           sage: e.highest_root()
           (1, 0, 0, -1)
        """
        return self.root(0,self.n-1)

    def fundamental_weight(self, i):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0)}

        """
        return self.sum(self.term(j) for j in range(i))

def dynkin_diagram(t):
    """
    Returns the graph corresponding to the Dynkin diagram
    of type A.

    EXAMPLES:
        sage: from sage.combinat.root_system.type_A import dynkin_diagram
        sage: ct = CartanType(['A',3])
        sage: a = dynkin_diagram(ct); a
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
    from dynkin_diagram import precheck, DynkinDiagram_class
    precheck(t, letter="A", length=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    return g

def affine_dynkin_diagram(t):
    """
    Returns the extended Dynkin diagram for affine type A.

    EXAMPLES:
        sage: from sage.combinat.root_system.type_A import affine_dynkin_diagram
        sage: ct = CartanType(['A',3,1])
        sage: a = affine_dynkin_diagram(ct); a
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
    from dynkin_diagram import precheck, DynkinDiagram_class
    precheck(t, letter="A", length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.add_edge(0, 1)
    g.add_edge(0, n)

    return g
