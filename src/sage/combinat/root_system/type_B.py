from ambient_space import AmbientSpace

class ambient_space(AmbientSpace):
    def dimension(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['B',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['B',3]).ambient_space()
            sage: e.root(0,1)
            (1, -1, 0)

        """
        return self._term(i) - self._term(j)

    def simple_root(self, i):
        """
        EXAMPLES:
            sage: e = RootSystem(['B',4]).ambient_space()
            sage: e.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1), 4: (0, 0, 0, 1)}
            sage: e.positive_roots()
            [(1, -1, 0, 0),
            (1, 1, 0, 0),
            (1, 0, -1, 0),
            (1, 0, 1, 0),
            (1, 0, 0, -1),
            (1, 0, 0, 1),
            (0, 1, -1, 0),
            (0, 1, 1, 0),
            (0, 1, 0, -1),
            (0, 1, 0, 1),
            (0, 0, 1, -1),
            (0, 0, 1, 1),
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, 1)]
            sage: e.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0), 4: (1/2, 1/2, 1/2, 1/2)}
        """
        assert(i in self.index_set())
        return self.root(i-1,i) if i < self.n else self._term(self.n-1)

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_space().negative_roots()
            [(-1, 1, 0),
             (-1, -1, 0),
             (-1, 0, 1),
             (-1, 0, -1),
             (0, -1, 1),
             (0, -1, -1),
             (-1, 0, 0),
             (0, -1, 0),
             (0, 0, -1)]

        """
        return [ -a for a in self.positive_roots()]


    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_space().positive_roots()
            [(1, -1, 0),
             (1, 1, 0),
             (1, 0, -1),
             (1, 0, 1),
             (0, 1, -1),
             (0, 1, 1),
             (1, 0, 0),
             (0, 1, 0),
             (0, 0, 1)]

        """
        res = []
        for i in range(self.n-1):
            for j in range(i+1,self.n):
                res.append(self._term(i) - self._term(j))
                res.append(self._term(i) + self._term(j))
        for i in range(self.n):
            res.append(self._term(i))
        return res

    def fundamental_weight(self, i):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1/2, 1/2, 1/2)}
        """
        assert(i in self.index_set())
        n = self.dimension()
        if i == n:
            return self.sum( self.term(j) for j in range(n) ) / 2
        else:
            return self.sum(self.term(j) for j in range(i))



def dynkin_diagram(t):
    """
    Returns a Dynkin diagram for type B.

    EXAMPLES:
         sage: from sage.combinat.root_system.type_B import dynkin_diagram
         sage: ct = CartanType(['B',3])
         sage: b = dynkin_diagram(ct);b
         O---O=>=O
         1   2   3
         B3
         sage: e = b.edges(); e.sort(); e
         [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]
    """
    from dynkin_diagram import precheck, DynkinDiagram_class
    precheck(t, letter='B', length=2, n_ge=2)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n-1, n, 2)
    return g

def affine_dynkin_diagram(t):
    """
    Returns the extended Dynkin diagram for affine type B.

    EXAMPLES:
       sage: DynkinDiagram(['B',3, 1]).edges()
       [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

    """
    from dynkin_diagram import precheck, DynkinDiagram_class
    precheck(t, letter='B', length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n):
        g.add_edge(i, i+1)
    g.set_edge_label(n-1, n, 2)
    g.add_edge(0,2)
    return g

