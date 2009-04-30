from ambient_space import AmbientSpace

class ambient_space(AmbientSpace):
    def dimension(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['D',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['D',3]).ambient_space()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
            sage: e.root(0, 0, 1, 1)
            (-1, 0, 0)

        """
        if i != j:
            return (-1)**p1 * self._term(i) + (-1)**p2 * self._term(j)
        else:
            return (-1)**p1 * self._term(i)

    def simple_root(self, i):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_space().simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1), 4: (0, 0, 1, 1)}

        """
        assert(i in self.index_set())
        return self.root(i-1, i, 0, 1) if i < self.n else self.root(self.n-2, self.n-1, 0, 0)

    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_space().positive_roots()
            [(1, 1, 0, 0),
             (1, 0, 1, 0),
             (0, 1, 1, 0),
             (1, 0, 0, 1),
             (0, 1, 0, 1),
             (0, 0, 1, 1),
             (1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]

        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,0,p) for i in range(j)])
        return res

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_space().negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (0, -1, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 0, 1),
             (0, 0, -1, 1),
             (-1, -1, 0, 0),
             (-1, 0, -1, 0),
             (0, -1, -1, 0),
             (-1, 0, 0, -1),
             (0, -1, 0, -1),
             (0, 0, -1, -1)]

        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,1,p) for i in range(j)])
        return res


    def fundamental_weight(self, i):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1/2, 1/2, 1/2, -1/2), 4: (1/2, 1/2, 1/2, 1/2)}
        """
        assert(i in self.index_set())
        n = self.dimension()
        if i == n:
            return  self.sum(self.term(j) for j in range(n)) / 2
        elif i == n-1:
            return (self.sum(self.term(j) for j in range(n-1)) - self._term(n-1)) / 2
        else:
            return  self.sum(self.term(j) for j in range(i))


def dynkin_diagram(t):
    """
    Returns a Dynkin diagram for type D.

    EXAMPLES:
        sage: from sage.combinat.root_system.type_D import dynkin_diagram
        sage: ct = CartanType(['D',4])
        sage: d = dynkin_diagram(ct);d
            O 3
            |
            |
        O---O---O
        1   2   4
        D4
        sage: e = d.edges(); e.sort(); e
        [(1, 2, 1), (2, 1, 1), (2, 3, 1), (2, 4, 1), (3, 2, 1), (4, 2, 1)]

    """
    from dynkin_diagram import precheck , DynkinDiagram_class
    precheck(t, letter="D", length=2, n_ge=3)
    n = t[1]

    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n-1):
        g.add_edge(i, i+1)
    g.add_edge(n-2,n)
    return g

def affine_dynkin_diagram(t):
    """
    Returns the extended Dynkin diagram for affine type D.

    EXAMPLES:

       sage: DynkinDiagram(CartanType(['D', 4, 1]))
            O 4
            |
            |
        O---O---O
        1   |2  3
            |
            O 0
        D4~
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
    from dynkin_diagram import precheck , DynkinDiagram_class
    precheck(t, letter="D", length=3, affine=1)
    n = t[1]
    g = DynkinDiagram_class(t)
    g.add_vertices(t.index_set())
    for i in range(1, n-1):
        g.add_edge(i, i+1)
    g.add_edge(n-2,n)
    g.add_edge(0,2)
    return g

