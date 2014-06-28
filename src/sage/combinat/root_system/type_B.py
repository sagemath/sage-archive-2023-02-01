"""
Root system data for type B
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import ambient_space

class  AmbientSpace(ambient_space.AmbientSpace):
    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['B',3]).ambient_space()
            sage: e.dimension()
            3
        """
        return self.root_system.cartan_type().rank()

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES::

            sage: e = RootSystem(['B',3]).ambient_space()
            sage: e.root(0,1)
            (1, -1, 0)

        """
        return self.monomial(i) - self.monomial(j)

    def simple_root(self, i):
        """
        EXAMPLES::

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
        return self.root(i-1,i) if i < self.n else self.monomial(self.n-1)

    def negative_roots(self):
        """
        EXAMPLES::

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
        EXAMPLES::

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
                res.append(self.monomial(i) - self.monomial(j))
                res.append(self.monomial(i) + self.monomial(j))
        for i in range(self.n):
            res.append(self.monomial(i))
        return res

    def fundamental_weight(self, i):
        """
        EXAMPLES::

            sage: RootSystem(['B',3]).ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1/2, 1/2, 1/2)}
        """
        assert(i in self.index_set())
        n = self.dimension()
        if i == n:
            return self.sum( self.monomial(j) for j in range(n) ) / 2
        else:
            return self.sum(self.monomial(j) for j in range(i))

from cartan_type import CartanType_standard_finite, CartanType_simple, CartanType_crystallographic, CartanType_simply_laced
class CartanType(CartanType_standard_finite, CartanType_simple, CartanType_crystallographic):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['B',4])
            sage: ct
            ['B', 4]
            sage: ct._repr_(compact = True)
            'B4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.affine()
            ['B', 4, 1]
            sage: ct.dual()
            ['C', 4]

            sage: ct = CartanType(['B',1])
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['B', 1, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 1
        CartanType_standard_finite.__init__(self, "B", n)
        if n == 1:
            self._add_abstract_superclass(CartanType_simply_laced)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['B',4]))
            B_{4}
        """
        return "B_{%s}"%self.n

    AmbientSpace = AmbientSpace

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['B',4]).coxeter_number()
            8
        """
        return 2*self.n

    def dual_coxeter_number(self):
        """
        Return the dual Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['B',4]).dual_coxeter_number()
            7
        """
        return 2*self.n - 1

    def dual(self):
        """
        Types B and C are in duality:

        EXAMPLES::

            sage: CartanType(["C", 3]).dual()
            ['B', 3]
        """
        import cartan_type
        return cartan_type.CartanType(["C", self.n])

    def dynkin_diagram(self):
        """
        Returns a Dynkin diagram for type B.

        EXAMPLES::

             sage: b = CartanType(['B',3]).dynkin_diagram()
             sage: b
             O---O=>=O
             1   2   3
             B3
             sage: sorted(b.edges())
             [(1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

             sage: b = CartanType(['B',1]).dynkin_diagram()
             sage: b
             O
             1
             B1
             sage: sorted(b.edges())
             []
        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        if n >= 2:
            g.set_edge_label(n-1, n, 2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType(['B',1]).ascii_art()
            O
            1
            sage: print CartanType(['B',2]).ascii_art()
            O=>=O
            1   2
            sage: print CartanType(['B',5]).ascii_art(label = lambda x: x+2)
            O---O---O---O=>=O
            3   4   5   6   7
        """
        n = self.n
        if n == 1:
            ret = "O\n"
        else:
            ret  = (n-2)*"O---" + "O=>=O\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret

    def _latex_dynkin_diagram(self, label=lambda x: x, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType(['B',4])._latex_dynkin_diagram()
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(5.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0) circle (.25cm) node[below=4pt]{$4$};

        When ``dual=True``, the Dynkin diagram for the dual Cartan
        type `C_n` is returned::

            sage: print CartanType(['B',4])._latex_dynkin_diagram(dual=True)
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(4.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0) circle (.25cm) node[below=4pt]{$4$};

        .. SEEALSO::

            - :meth:`sage.combinat.root_system.type_C.CartanType._latex_dynkin_diagram`
            - :meth:`sage.combinat.root_system.type_BC_affine.CartanType._latex_dynkin_diagram`
        """
        if self.n == 1:
            return "\\draw[fill=white] (0,0) circle (.25cm) node[below=4pt]{$1$};"
        n = self.n
        ret = "\\draw (0 cm,0) -- (%s cm,0);\n"%((n-2)*node_dist)
        ret += "\\draw (%s cm, 0.1 cm) -- +(%s cm,0);\n"%((n-2)*node_dist, node_dist)
        ret += "\\draw (%s cm, -0.1 cm) -- +(%s cm,0);\n"%((n-2)*node_dist, node_dist)
        if dual:
            ret += self._latex_draw_arrow_tip((n-1.5)*node_dist-0.2, 0, 180)
        else:
            ret += self._latex_draw_arrow_tip((n-1.5)*node_dist+0.2, 0, 0)
        for i in range(self.n-1):
            ret += "\\draw[fill=white] (%s cm, 0) circle (.25cm) node[below=4pt]{$%s$};\n"%(i*node_dist, label(i+1))
        ret += "\\draw[fill=white] (%s cm, 0) circle (.25cm) node[below=4pt]{$%s$};"%((n-1)*node_dist, label(n))
        return ret

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['B', 3])._default_folded_cartan_type()
            ['B', 3] as a folding of ['D', 4]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        n = self.n
        return CartanTypeFolded(self, ['D', n+1],
            [[i] for i in range(1, n)] + [[n, n+1]])

# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_B', 'ambient_space',  AmbientSpace)

