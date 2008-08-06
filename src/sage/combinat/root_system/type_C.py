from ambient_space import AmbientSpace
from sage.rings.all import ZZ

class ambient_space(AmbientSpace):

# The coroots can't be defined with integer coefficients!
#    @classmethod
#    def smallest_base_ring(cls):
#        return ZZ;

    def dimension(self):
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
