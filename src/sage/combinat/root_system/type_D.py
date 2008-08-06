from ambient_space import AmbientSpace

class ambient_space(AmbientSpace):
    def dimension(self):
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
