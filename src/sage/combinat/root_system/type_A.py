"""
Root system data for type A
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.all import ZZ
import ambient_space

class AmbientSpace(ambient_space.AmbientSpace):
    r"""
    EXAMPLES::

        sage: R = RootSystem(["A",3])
        sage: e = R.ambient_space(); e
        Ambient space of the Root system of type ['A', 3]
        sage: e == loads(dumps(e))
        True
    """
    @classmethod
    def smallest_base_ring(cls, cartan_type=None):
        """
        Returns the smallest base ring the ambient space can be defined upon

        .. SEEALSO:: :meth:`~sage.combinat.root_system.ambient_space.AmbientSpace.smallest_base_ring`

        EXAMPLES::

            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.smallest_base_ring()
            Integer Ring
        """
        return ZZ

    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(["A",3]).ambient_space()
            sage: e.dimension()
            4
        """
        return self.root_system.cartan_type().rank()+1

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.root(0,1)
            (1, -1, 0, 0)
        """
        return self.monomial(i) - self.monomial(j)

    def simple_root(self, i):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}
        """
        return self.root(i-1, i)

    def negative_roots(self):
        """
        EXAMPLES::

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
        EXAMPLES::

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
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.fundamental_weights()
            Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0)}

        """
        return self.sum(self.monomial(j) for j in range(i))

    def det(self, k=1):
        """
        returns the vector (1, ... ,1) which in the ['A',r]
        weight lattice, interpreted as a weight of GL(r+1,CC)
        is the determinant. If the optional parameter k is
        given, returns (k, ... ,k), the k-th power of the
        determinant.

        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_space()
            sage: e.det(1/2)
            (1/2, 1/2, 1/2, 1/2)
        """
        return self.sum(self.monomial(j)*k for j in range(self.n))


from cartan_type import CartanType_standard_finite, CartanType_simply_laced, CartanType_simple
class CartanType(CartanType_standard_finite, CartanType_simply_laced, CartanType_simple):
    """
    Cartan Type A

    SEE ALSO:func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])
            sage: ct
            ['A', 4]
            sage: ct._repr_(compact = True)
            'A4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystalographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['A', 4, 1]
            sage: ct.dual()
            ['A', 4]

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 0
        CartanType_standard_finite.__init__(self, "A", n)

    AmbientSpace = AmbientSpace

    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram of type A.

        EXAMPLES::

            sage: a = CartanType(['A',3]).dynkin_diagram()
            sage: a
            O---O---O
            1   2   3
            A3
            sage: sorted(a.edges())
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 1)]

        TEST::

            sage: a = DynkinDiagram(['A',1])
            sage: a
            O
            1
            A1
            sage: a.vertices(), a.edges()
            ([1], [])
        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns an ascii art representation of the Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['A',0]).ascii_art()
            sage: print CartanType(['A',1]).ascii_art()
            O
            1
            sage: print CartanType(['A',3]).ascii_art()
            O---O---O
            1   2   3
            sage: print CartanType(['A',5]).ascii_art(label = lambda x: x+2)
            O---O---O---O---O
            3   4   5   6   7
        """
        n = self.n
        if n == 0:
            return ""
        ret  = "---".join("O"           for i in range(1,n+1)) + "\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret

# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_A', 'ambient_space',  AmbientSpace)
