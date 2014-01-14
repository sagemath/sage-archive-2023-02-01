"""
Root system data for (untwisted) type D affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_untwisted_affine, CartanType_simply_laced
class CartanType(CartanType_standard_untwisted_affine, CartanType_simply_laced):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['D',4,1])
            sage: ct
            ['D', 4, 1]
            sage: ct._repr_(compact = True)
            'D4~'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            True
            sage: ct.is_untwisted_affine()
            True
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.classical()
            ['D', 4]
            sage: ct.dual()
            ['D', 4, 1]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 3
        CartanType_standard_untwisted_affine.__init__(self, "D", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type D.

        EXAMPLES::

           sage: d = CartanType(['D', 6, 1]).dynkin_diagram()
           sage: d
              0 O       O 6
                |       |
                |       |
            O---O---O---O---O
            1   2   3   4   5
            D6~
           sage: sorted(d.edges())
           [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 1),
            (3, 2, 1), (3, 4, 1), (4, 3, 1), (4, 5, 1), (4, 6, 1), (5, 4, 1), (6, 4, 1)]

           sage: d = CartanType(['D', 4, 1]).dynkin_diagram()
           sage: d
               O 4
               |
               |
           O---O---O
           1   |2  3
               |
               O 0
           D4~
           sage: sorted(d.edges())
           [(0, 2, 1),
            (1, 2, 1),
            (2, 0, 1),
            (2, 1, 1),
            (2, 3, 1),
            (2, 4, 1),
            (3, 2, 1),
            (4, 2, 1)]

           sage: d = CartanType(['D', 3, 1]).dynkin_diagram()
           sage: d
           0
           O-------+
           |       |
           |       |
           O---O---O
           3   1   2
           D3~
           sage: sorted(d.edges())
           [(0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (2, 0, 1), (2, 1, 1), (3, 0, 1), (3, 1, 1)]

        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        if n == 3:
            import cartan_type
            res = cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).dynkin_diagram()
            res._cartan_type = self
            return res
        g = DynkinDiagram_class(self)
        for i in range(1, n-1):
            g.add_edge(i, i+1)
        g.add_edge(n-2,n)
        g.add_edge(0,2)
        return g

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2, dual=False):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType(['D',4,1])._latex_dynkin_diagram()
            \draw (0,0.7 cm) -- (2 cm,0);
            \draw (0,-0.7 cm) -- (2 cm,0);
            \draw (2 cm,0) -- (2 cm,0);
            \draw (2 cm,0) -- (4 cm,0.7 cm);
            \draw (2 cm,0) -- (4 cm,-0.7 cm);
            \draw[fill=white] (0, 0.7 cm) circle (.25cm) node[left=3pt]{$0$};
            \draw[fill=white] (0, -0.7 cm) circle (.25cm) node[left=3pt]{$1$};
            \draw[fill=white] (2 cm, 0) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0.7 cm) circle (.25cm) node[right=3pt]{$4$};
            \draw[fill=white] (4 cm, -0.7 cm) circle (.25cm) node[right=3pt]{$3$};
        """
        n = self.n
        if n == 3:
            import cartan_type
            relabel = {0:label(0), 1:label(3), 2:label(1), 3:label(2)}
            return cartan_type.CartanType(["A",3,1]).relabel(relabel)._latex_dynkin_diagram(node_dist=node_dist)
        if self.global_options('mark_special_node') in ['latex', 'both']:
            special_fill = 'black'
        else:
            special_fill = 'white'
        rt_most = (n-2)*node_dist
        center_point = rt_most-node_dist
        ret = "\\draw (0,0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (0,-0.7 cm) -- (%s cm,0);\n"%node_dist
        ret += "\\draw (%s cm,0) -- (%s cm,0);\n"%(node_dist, center_point)
        ret += "\\draw (%s cm,0) -- (%s cm,0.7 cm);\n"%(center_point, rt_most)
        ret += "\\draw (%s cm,0) -- (%s cm,-0.7 cm);\n"%(center_point, rt_most)
        ret += "\\draw[fill=%s] (0, 0.7 cm) circle (.25cm) node[left=3pt]{$%s$};\n"%(special_fill, label(0))
        ret += "\\draw[fill=white] (0, -0.7 cm) circle (.25cm) node[left=3pt]{$%s$};\n"%label(1)
        for i in range(1, self.n-2):
            ret += "\\draw[fill=white] (%s cm, 0) circle (.25cm) node[below=4pt]{$%s$};\n"%(i*node_dist, label(i+1))
        ret += "\\draw[fill=white] (%s cm, 0.7 cm) circle (.25cm) node[right=3pt]{$%s$};\n"%(rt_most, label(n))
        ret += "\\draw[fill=white] (%s cm, -0.7 cm) circle (.25cm) node[right=3pt]{$%s$};"%(rt_most, label(n-1))
        return ret

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        TESTS::

            sage: print CartanType(['D',6,1]).ascii_art(label = lambda x: x+2)
              2 O       O 8
                |       |
                |       |
            O---O---O---O---O
            3   4   5   6   7

            sage: print CartanType(['D',4,1]).ascii_art(label = lambda x: x+2)
                O 6
                |
                |
            O---O---O
            3   |4  5
                |
                O 2

            sage: print CartanType(['D',3,1]).ascii_art(label = lambda x: x+2)
            2
            O-------+
            |       |
            |       |
            O---O---O
            5   3   4
        """
        n = self.n
        if n == 3:
            import cartan_type
            return cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).ascii_art(label)
        if self.global_options('mark_special_node') in ['printing', 'both']:
            special_str = self.global_options('special_node_str')
        else:
            special_str = 'O'
        if n == 4:
            return "    " + special_str + " %s\n    |\n    |\nO---O---O\n%s   |%s  %s\n    |\n    O %s"%tuple(label(i) for i in (4,1,2,3,0))
        ret = "  %s O"%label(0) +(4*(n-4)-1)*" "+"O %s"%label(n)+ "\n"
        ret += "    |"          +(4*(n-4)-1)*" "                  +"|\n"
        ret += "    |"          +(4*(n-4)-1)*" "                  +"|\n"
        ret += (n-2)*"O---"+"O\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n))
        return ret

