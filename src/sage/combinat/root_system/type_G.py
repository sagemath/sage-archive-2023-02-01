from ambient_space import AmbientSpace
from sage.combinat.family import Family

# TODO: check whether this can be defined over ZZ

class ambient_space(AmbientSpace):
    """
    TESTS:
        sage: [WeylDim(['G',2],[a,b]) for a,b in [[0,0], [1,0], [0,1], [1,1]]]
        [1, 7, 14, 64]

    EXAMPLES:
        sage: e = RootSystem(['G',2]).ambient_space()
        sage: e == loads(dumps(e))
        True
    """

    def dimension(self):
        return 3

    def simple_root(self, i):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_space().simple_roots()
	    Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
         """
        return self._term(1)-self._term(2) if i == 1 else self._term(0)-2*self._term(1)+self._term(2)

    def positive_roots(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_space().positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return [ self(v) for v in
                 [[0,1,-1],[1,-2,1],[1,-1,0],[1,0,-1],[1,1,-2],[2,-1,-1]]]

    def negative_roots(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_space().negative_roots()
            [(0, -1, 1), (-1, 2, -1), (-1, 1, 0), (-1, 0, 1), (-1, -1, 2), (-2, 1, 1)]
        """
        return [ self(v) for v in
                 [[0,-1,1],[-1,2,-1],[-1,1,0],[-1,0,1],[-1,-1,2],[-2,1,1]]]

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_space().fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return Family({ 1: self([1,0,-1]),
                        2: self([2,-1,-1])})
