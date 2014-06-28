"""
Root system data for type BC affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_affine
from sage.rings.integer_ring import ZZ
class CartanType(CartanType_standard_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['BC',4,2])
            sage: ct
            ['BC', 4, 2]
            sage: ct._repr_(compact = True)
            'BC4~'
            sage: ct.dynkin_diagram()
            O=<=O---O---O=<=O
            0   1   2   3   4
            BC4~

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.classical()
            ['C', 4]

            sage: dual = ct.dual()
            sage: dual.dynkin_diagram()
            O=>=O---O---O=>=O
            0   1   2   3   4
            BC4~*

            sage: dual.special_node()
            0
            sage: dual.classical().dynkin_diagram()
            O---O---O=>=O
            1   2   3   4
            B4

            sage: CartanType(['BC',1,2]).dynkin_diagram()
              4
            O=<=O
            0   1
            BC1~

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n in ZZ and n >= 1
        CartanType_standard_affine.__init__(self, "BC", n, 2)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type BC.

        EXAMPLES::

            sage: c = CartanType(['BC',3,2]).dynkin_diagram()
            sage: c
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

            sage: c = CartanType(["A", 6, 2]).dynkin_diagram() # should be the same as above; did fail at some point!
            sage: c
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

            sage: c = CartanType(['BC',2,2]).dynkin_diagram()
            sage: c
            O=<=O=<=O
            0   1   2
            BC2~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 2)]

            sage: c = CartanType(['BC',1,2]).dynkin_diagram()
            sage: c
              4
            O=<=O
            0   1
            BC1~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 4)]

        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        if n == 1:
            g.add_edge(1,0,4)
            return g
        g.add_edge(1,0,2)
        for i in range(1, n-1):
            g.add_edge(i, i+1)
        g.add_edge(n,n-1,2)
        return g

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['BC',4,2]))
            BC_{4}^{(2)}

            sage: CartanType.global_options['notation'] = 'Kac'
            sage: latex(CartanType(['BC',4,2]))
            A_{8}^{(2)}
            sage: latex(CartanType(['A',8,2]))
            A_{8}^{(2)}
            sage: CartanType.global_options.reset()
        """
        if self.global_options('notation') == "Kac":
            return "A_{%s}^{(2)}"%(2*self.classical().rank())
        else:
            return "BC_{%s}^{(2)}"%self.n

    def _latex_dynkin_diagram(self, label=lambda x: x, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType(['BC',4,2])._latex_dynkin_diagram()
            \draw (0, 0.1 cm) -- +(2 cm,0);
            \draw (0, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(0.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            {
            \pgftransformxshift{2 cm}
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(4.8, 0)}, rotate=180] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0) circle (.25cm) node[below=4pt]{$4$};
            }
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$0$};
            sage: print CartanType(['BC',4,2]).dual()._latex_dynkin_diagram()
            \draw (0, 0.1 cm) -- +(2 cm,0);
            \draw (0, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(1.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            {
            \pgftransformxshift{2 cm}
            \draw (0 cm,0) -- (4 cm,0);
            \draw (4 cm, 0.1 cm) -- +(2 cm,0);
            \draw (4 cm, -0.1 cm) -- +(2 cm,0);
            \draw[shift={(5.2, 0)}, rotate=0] (135 : 0.45cm) -- (0,0) -- (-135 : 0.45cm);
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0) circle (.25cm) node[below=4pt]{$4$};
            }
            \draw[fill=white] (0 cm, 0) circle (.25cm) node[below=4pt]{$0$};
        """
        if self.global_options('mark_special_node') in ['latex', 'both']:
            special_fill = 'black'
        else:
            special_fill = 'white'
        if self.n == 1:
            ret = "\\draw (0, 0.05 cm) -- +(%s cm,0);\n"%node_dist
            ret += "\\draw (0, -0.05 cm) -- +(%s cm,0);\n"%node_dist
            ret += "\\draw (0, 0.15 cm) -- +(%s cm,0);\n"%node_dist
            ret += "\\draw (0, -0.15 cm) -- +(%s cm,0);\n"%node_dist
            if dual:
                ret += self._latex_draw_arrow_tip(0.5*node_dist+0.2, 0, 0)
            else:
                ret += self._latex_draw_arrow_tip(0.5*node_dist-0.2, 0, 180)
            ret += "\\draw[fill=%s] (0,0) circle (.25cm) node[below=4pt]{$%s$};\n"%(special_fill, label(0))
            ret += "\\draw[fill=white] (%s cm,0) circle (.25cm) node[below=4pt]{$%s$};"%(node_dist, label(1))
            return ret
        n = self.n
        ret = "\\draw (0, 0.1 cm) -- +(%s cm,0);\n"%node_dist
        ret += "\\draw (0, -0.1 cm) -- +(%s cm,0);\n"%node_dist
        if dual:
            ret += self._latex_draw_arrow_tip(0.5*node_dist+0.2, 0, 0)
        else:
            ret += self._latex_draw_arrow_tip(0.5*node_dist-0.2, 0, 180)
        ret += "{\n\\pgftransformxshift{%s cm}\n"%node_dist
        classical = self.classical()
        ret += self.classical()._latex_dynkin_diagram(label, node_dist, dual=dual)
        ret += "\n}\n\\draw[fill=%s] (0 cm, 0) circle (.25cm) node[below=4pt]{$%s$};"%(special_fill, label(0))
        return ret

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['BC',2,2]).ascii_art()
            O=<=O=<=O
            0   1   2
            sage: print CartanType(['BC',3,2]).ascii_art()
            O=<=O---O=<=O
            0   1   2   3
            sage: print CartanType(['BC',5,2]).ascii_art(label = lambda x: x+2)
            O=<=O---O---O---O=<=O
            2   3   4   5   6   7

            sage: print CartanType(['BC',1,2]).ascii_art(label = lambda x: x+2)
              4
            O=<=O
            2   3
        """
        n = self.n
        if self.global_options('mark_special_node') in ['printing', 'both']:
            special_str = self.global_options('special_node_str')
        else:
            special_str = 'O'
        if n == 1:
            return "  4\n%s=<=O\n%s   %s"%(special_str, label(0), label(1))
        ret = "%s=<=O"%special_str + (n-2)*"---O"+"=<=O\n%s   "%label(0)
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret

    def classical(self):
        """
        Returns the classical Cartan type associated with self

            sage: CartanType(["BC", 3, 2]).classical()
            ['C', 3]
        """
        import cartan_type
        return cartan_type.CartanType(["C", self.n])

    def basic_untwisted(self):
        r"""
        Return the basic untwisted Cartan type associated with this affine
        Cartan type.

        Given an affine type `X_n^{(r)}`, the basic untwisted type is `X_n`.
        In other words, it is the classical Cartan type that is twisted to
        obtain ``self``.

        EXAMPLES::

            sage: CartanType(['A', 2, 2]).basic_untwisted()
            ['A', 2]
            sage: CartanType(['A', 4, 2]).basic_untwisted()
            ['A', 4]
            sage: CartanType(['BC', 4, 2]).basic_untwisted()
            ['A', 8]
        """
        import cartan_type
        return cartan_type.CartanType(["A", 2*self.n])

    def _default_folded_cartan_type(self):
        """
        Return the default folded Cartan type.

        EXAMPLES::

            sage: CartanType(['BC', 3, 2])._default_folded_cartan_type()
            ['BC', 3, 2] as a folding of ['A', 5, 1]
        """
        from sage.combinat.root_system.type_folded import CartanTypeFolded
        n = self.n
        return CartanTypeFolded(self, ['A', 2*n - 1, 1],
            [[0]] + [[i, 2*n-i] for i in range(1, n)] + [[n]])

