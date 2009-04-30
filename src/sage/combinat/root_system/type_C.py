from ambient_space import AmbientSpace
from sage.rings.all import ZZ

class ambient_space(AmbientSpace):

# The coroots can't be defined with integer coefficients!
#    @classmethod
#    def smallest_base_ring(cls):
#        return ZZ;

    def dimension(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
        """
        return (-1)**p1 * self._term(i) + (-1)**p2 * self._term(j)

    def simple_root(self, i):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_space().simple_roots()
            Finite family {1: (1, -1, 0), 2: (0, 1, -1), 3: (0, 0, 2)}
        """
        assert(i in self.index_set())
        return self.root(i-1, i,0,1) if i < self.n else self.root(self.n-1, self.n-1, 0, 0)

    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_space().positive_roots()
            [(1, 1, 0),
             (1, 0, 1),
             (0, 1, 1),
             (1, -1, 0),
             (1, 0, -1),
             (0, 1, -1),
             (2, 0, 0),
             (0, 2, 0),
             (0, 0, 2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,0,p) for i in range(j)])
        res.extend([self.root(i,i,0,0) for i in range(self.n)])
        return res

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_space().negative_roots()
            [(-1, 1, 0),
             (-1, 0, 1),
             (0, -1, 1),
             (-1, -1, 0),
             (-1, 0, -1),
             (0, -1, -1),
             (-2, 0, 0),
             (0, -2, 0),
             (0, 0, -2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend( [self.root(i,j,1,p) for i in range(j) ] )
        res.extend( [ self.root(i,i,1,1) for i in range(self.n) ] )
        return res


    def fundamental_weight(self, i):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1)}
        """
        return self.sum(self._term(j) for j in range(i))


def dynkin_diagram(t):
    """
    Returns a Dynkin diagram for type C.

    EXAMPLES:
        sage: from sage.combinat.root_system.type_C import dynkin_diagram
        sage: ct = CartanType(['C',3])
        sage: c = dynkin_diagram(ct); c
        O---O=<=O
        1   2   3
        C3
        sage: e = c.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

    """
    from dynkin_diagram import precheck , DynkinDiagram_class
    precheck(t, letter='C', length=2, n_ge=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n,n-1,2)
    return g

def affine_dynkin_diagram(t):
    """
    Returns the extended Dynkin diagram for affine type C.

    EXAMPLES:
        sage: from sage.combinat.root_system.type_C import affine_dynkin_diagram
        sage: ct = CartanType(['C',3,1])
        sage: c = affine_dynkin_diagram(ct);c
         O=>=O---O=<=O
         0   1   2   3
         C3~
        sage: e = c.edges(); e.sort(); e
        [(0, 1, 2), (1, 0, 1), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

    """
    from dynkin_diagram import precheck , DynkinDiagram_class
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

