"""
Root system data for type C
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import ambient_space

class  AmbientSpace(ambient_space.AmbientSpace):
    """
    EXAMPLES::

        sage: e = RootSystem(['C',2]).ambient_space(); e
        Ambient space of the Root system of type ['C', 2]

    One cannot construct the ambient lattice because the fundamental
    coweights have rational coefficients::

        sage: e.smallest_base_ring()
        Rational Field

        sage: RootSystem(['B',2]).ambient_space().fundamental_weights()
        Finite family {1: (1, 0), 2: (1/2, 1/2)}

    TESTS::

        sage: TestSuite(e).run()
    """


    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['C',3]).ambient_space()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
        """
        return (-1)**p1 * self.monomial(i) + (-1)**p2 * self.monomial(j)

    def simple_root(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().simple_roots()
            Finite family {1: (1, -1, 0), 2: (0, 1, -1), 3: (0, 0, 2)}
        """
        assert(i in self.index_set())
        return self.root(i-1, i,0,1) if i < self.n else self.root(self.n-1, self.n-1, 0, 0)

    def positive_roots(self):
        """
        EXAMPLES::

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
        EXAMPLES::

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
        EXAMPLES::

            sage: RootSystem(['C',3]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1)}
        """
        return self.sum(self.monomial(j) for j in range(i))

from cartan_type import CartanType_standard_finite, CartanType_simple, CartanType_crystalographic, CartanType_simply_laced
class CartanType(CartanType_standard_finite, CartanType_simple, CartanType_crystalographic):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['C',4])
            sage: ct
            ['C', 4]
            sage: ct._repr_(compact = True)
            'C4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_crystalographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.affine()
            ['C', 4, 1]
            sage: ct.dual()
            ['B', 4]

            sage: ct = CartanType(['C',1])
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['C', 1, 1]

        TESTS:
            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 1
        CartanType_standard_finite.__init__(self, "C", n)
        if n == 1:
            self._add_abstract_superclass(CartanType_simply_laced)

    AmbientSpace = AmbientSpace

    def dual(self):
        """
        Types B and C are in duality:

        EXAMPLES::

            sage: CartanType(["C", 3]).dual()
            ['B', 3]
        """
        import cartan_type
        return cartan_type.CartanType(["B", self.n])

    def dynkin_diagram(self):
        """
        Returns a Dynkin diagram for type C.

        EXAMPLES::

            sage: c = CartanType(['C',3]).dynkin_diagram()
            sage: c
            O---O=<=O
            1   2   3
            C3
            sage: e = c.edges(); e.sort(); e
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

             sage: b = CartanType(['C',1]).dynkin_diagram()
             sage: b
             O
             1
             C1
             sage: sorted(b.edges())
             []
        """
        return self.dual().dynkin_diagram().dual()

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['C',1]).ascii_art()
            O
            1
            sage: print CartanType(['C',2]).ascii_art()
            O=<=O
            1   2
            sage: print CartanType(['C',3]).ascii_art()
            O---O=<=O
            1   2   3
            sage: print CartanType(['C',5]).ascii_art(label = lambda x: x+2)
            O---O---O---O=<=O
            3   4   5   6   7
        """
        return self.dual().ascii_art(label = label).replace("=>=", "=<=")

# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_C', 'ambient_space',  AmbientSpace)
