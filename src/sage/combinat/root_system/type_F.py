from ambient_space import AmbientSpace
from sage.rings.all import ZZ
from sage.combinat.family import Family

# TODO: double check that this can't be defined over ZZ

class ambient_space(AmbientSpace):
    """
    The lattice behind F4.  The computations are based on Bourbaki, Groupes et Algebres de Lie,
    Ch. 4,5,6 (planche VIII).
    """

    def dimension(self):
        return self.root_system.cartan_type().rank()

    def __init__(self, root_system, base_ring):
        """
        Create the ambient lattice for the root system for F4.
        Specify the Base, i.e., the simple roots w.r. to the canonical
        basis for R^4.

        EXAMPLES:
            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e == loads(dumps(e))
            True
        """
        AmbientSpace.__init__(self, root_system, base_ring)
        v = ZZ(1)/ZZ(2)
        self.Base = [self.root(1,2,p2=1),
                     self.root(2,3,p2=1),
                     self.root(3),
                     v*(self.root(0)-self.root(1)-self.root(2)-self.root(3))]

    def root(self, i, j=None, k=None, l=None, p1=0, p2=0, p3=0, p4=0):
        """
        Compute a root from base elements of the underlying lattice.
        The arguments specify the basis elements and the signs.
        Sadly, the base elements are indexed zero-based.
        We assume that if one of the indices is not given, the rest are not as well.

        EXAMPLES:
            sage: F4 = RootSystem(['F',4])
            sage: LF4 = F4.ambient_space()
            sage: [ LF4.root(i,j,p2=1) for i in xrange(LF4.n) for j in xrange(i+1,LF4.n) ]
            [(1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1), (0, 1, -1, 0), (0, 1, 0, -1), (0, 0, 1, -1)]
        """
        if i == j or j == None:
            return (-1)**p1*self._term(i)
        if k == None:
            return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)
        if l == None:
            return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)+(-1)**p3*self._term(k)
        return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)+(-1)**p3*self._term(k)+(-1)**p4*self._term(l)

    def simple_root(self, i):
        """
        There are computed as what Bourbaki calls the Base:
            a1 = e2-e3, a2 = e3-e4, a3 = e4, a4 = 1/2*(e1-e2-e3-e4)

        EXAMPLES:
            sage: LF4 = RootSystem(['F',4]).ambient_space()
            sage: LF4.simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 1), 4: (1/2, -1/2, -1/2, -1/2)}
        """
        return self.Base[i-1]

    def negative_roots(self):
        """
        Returns the negative roots in self.

        EXAMPLES:
            sage: LF4 =  RootSystem(['F',4]).ambient_space()
            sage: LF4.negative_roots()
            [(-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (-1, -1, 0, 0), (-1, 0, -1, 0), (-1, 0, 0, -1), (0, -1, -1, 0), (0, -1, 0, -1), (0, 0, -1, -1), (-1, 1, 0, 0), (-1, 0, 1, 0), (-1, 0, 0, 1), (0, -1, 1, 0), (0, -1, 0, 1), (0, 0, -1, 1), (-1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2), (-1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, 1/2), (-1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, 1/2)]
        """
        return [ -a for a in self.positive_roots()]

    def positive_roots(self):
        """
        These are the roots positive w.r. to lexicographic ordering of the
        basis elements (e1<...<e4).

        EXAMPLES:
            sage: LF4 = RootSystem(['F',4]).ambient_space()
            sage: LF4.positive_roots()
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 1, 1), (1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1), (0, 1, -1, 0), (0, 1, 0, -1), (0, 0, 1, -1), (1/2, 1/2, 1/2, 1/2), (1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, 1/2), (1/2, 1/2, -1/2, -1/2), (1/2, -1/2, 1/2, 1/2), (1/2, -1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2), (1/2, -1/2, -1/2, -1/2)]
            sage: LF4.rho()
            (11/2, 5/2, 3/2, 1/2)
        """
        v = ZZ(1)/ZZ(2)
        if not hasattr(self, 'PosRoots'):
            self.PosRoots = ([ self._term(i) for i in xrange(self.n) ] +
                            [ self.root(i,j,p2=0) for i in xrange(self.n) for j in xrange(i+1,self.n) ] +
                            [ self.root(i,j,p2=1) for i in xrange(self.n) for j in xrange(i+1,self.n) ] +
                            [ v*self.root(0,1,2,3,0,p2,p3,p4) for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] ])
        return self.PosRoots

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: LF4 =  RootSystem(['F',4]).ambient_space()
            sage: LF4.fundamental_weights()
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0), 3: (3/2, 1/2, 1/2, 1/2), 4: (1, 0, 0, 0)}
        """
        v = ZZ(1)/ZZ(2)
        return Family({ 1: self._term(0)+self._term(1),
                        2: 2*self._term(0)+self._term(1)+self._term(2),
                        3: v*(3*self._term(0)+self._term(1)+self._term(2)+self._term(3)),
                        4: self._term(0)})
