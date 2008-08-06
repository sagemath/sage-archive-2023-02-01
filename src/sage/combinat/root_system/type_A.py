from ambient_space import AmbientSpace
from sage.rings.all import ZZ
from sage.combinat.family import Family

class ambient_space(AmbientSpace):
    r"""
        sage: R = RootSystem(["A",3])
        sage: R.ambient_space ()
        Ambient space for the Root system of type ['A', 3]
    """

    @classmethod
    def smallest_base_ring(cls):
        return ZZ;

    def dimension(self):
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
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}
        """
        return self.root(i-1, i)

    def negative_roots(self):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
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
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
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
           sage: e = CartanType(['A',3]).root_system().ambient_lattice()
           sage: e.highest_root()
           (1, 0, 0, -1)
        """
        return self.root(0,self.n-1)

    def fundamental_weight(self, i):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0)}

        """
        return self.sum(self.term(j) for j in range(i))
