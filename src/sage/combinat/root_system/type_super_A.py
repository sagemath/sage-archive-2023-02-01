"""
Root system data for super type A
"""
#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from sage.rings.all import ZZ
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from . import ambient_space
from .cartan_type import SuperCartanType_standard

class AmbientSpace(ambient_space.AmbientSpace):
    r"""
    EXAMPLES::

        sage: R = RootSystem(["A",3])
        sage: e = R.ambient_space(); e
        Ambient space of the Root system of type ['A', 3]
        sage: TestSuite(e).run()
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
        EXAMPLES::

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

        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_space()
            sage: e.det(1/2)
            (1/2, 1/2, 1/2, 1/2)
        """
        return self.sum(self.monomial(j)*k for j in range(self.n))

    __doc__ += """
    By default, this ambient space uses the barycentric projection for plotting::

        sage: L = RootSystem(["A",2]).ambient_space()
        sage: e = L.basis()
        sage: L._plot_projection(e[0])
        (1/2, 989/1142)
        sage: L._plot_projection(e[1])
        (-1, 0)
        sage: L._plot_projection(e[2])
        (1/2, -989/1142)
        sage: L = RootSystem(["A",3]).ambient_space()
        sage: l = L.an_element(); l
        (2, 2, 3, 0)
        sage: L._plot_projection(l)
        (0, -1121/1189, 7/3)

    .. SEEALSO::

        - :meth:`sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations.ParentMethods._plot_projection`
    """
    _plot_projection = RootLatticeRealizations.ParentMethods.__dict__['_plot_projection_barycentric']

class CartanType(SuperCartanType_standard):
    """
    Cartan Type `A(m|n)`.

    .. SEEALSO:: :func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """
    def __init__(self, m, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A', [4,2]])
            sage: ct
            ['A', [4, 2]]
            sage: ct._repr_(compact = True)
            'A4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['A', 4, 1]
            sage: ct.dual()
            ['A', 4]

        TESTS::

            sage: TestSuite(ct).run()
        """
        self.m = m
        self.n = n
        self.letter = "A"

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['A',[4,3]]))
            A(4|3)
        """
        return "A(%s|%s)"%(self.m, self.n)

    def index_set(self):
        """
        Return the index set of ``self``.
        """
        return tuple(range(-self.m, self.n+1))

    AmbientSpace = AmbientSpace

    def is_irreducible(self):
        """
        Return whether ``self`` is irreducible, which is ``True``.

        EXAMPLES::

            sage: CartanType(['A', [3,4]]).is_irreducible()
            True
        """
        return True

    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram of super type A.

        EXAMPLES::

            sage: a = CartanType(['A',3]).dynkin_diagram()
            sage: a
            O---O---O
            1   2   3
            A3
            sage: sorted(a.edges())
            [(1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 1)]

        TESTS::

            sage: a = DynkinDiagram(['A',1])
            sage: a
            O
            1
            A1
            sage: a.vertices(), a.edges()
            ([1], [])
        """
        from .dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        return g

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A',4])._latex_dynkin_diagram())
            \draw (0 cm,0) -- (6 cm,0);
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            \draw[fill=white] (2 cm, 0 cm) circle (.25cm) node[below=4pt]{$2$};
            \draw[fill=white] (4 cm, 0 cm) circle (.25cm) node[below=4pt]{$3$};
            \draw[fill=white] (6 cm, 0 cm) circle (.25cm) node[below=4pt]{$4$};
            <BLANKLINE>

            sage: print(CartanType(['A',0])._latex_dynkin_diagram())
            <BLANKLINE>
            sage: print(CartanType(['A',1])._latex_dynkin_diagram())
            \draw[fill=white] (0 cm, 0 cm) circle (.25cm) node[below=4pt]{$1$};
            <BLANKLINE>
        """
        if node is None:
            node = self._latex_draw_node
        if self.n > 1:
            ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format((self.n-1)*node_dist)
        else:
            ret = ""
        return ret + "".join(node((i-1)*node_dist, 0, label(i))
                             for i in self.index_set())

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: t = CartanType(['A', [3,2]])
            sage: print(t.ascii_art())
            O---O---O---X---O---O
            -3  -2  -1  0   1   2
            sage: t = CartanType(['A', [3,7]])
            sage: print(t.ascii_art())
            O---O---O---X---O---O---O---O---O---O---O
            -3  -2  -1  0   1   2   3   4   5   6   7

            sage: t = CartanType(['A', [0,7]])
            sage: print(t.ascii_art())
            X---O---O---O---O---O---O---O
            0   1   2   3   4   5   6   7
            sage: t = CartanType(['A', [0,0]])
            sage: print(t.ascii_art())
            X
            0
            sage: t = CartanType(['A', [5,0]])
            sage: print(t.ascii_art())
            O---O---O---O---O---X
            -5  -4  -3  -2  -1  0
        """
        if node is None:
            node = lambda i: 'O'
        ret  = "---".join(node(label(i)) for i in range(1,self.m+1))
        if self.m == 0:
            if self.n == 0:
                ret = "X"
            else:
                ret += "X---"
        else:
            if self.n == 0:
                ret += "---X"
            else:
                ret += "---X---"
        ret += "---".join(node(label(i)) for i in range(1,self.n+1)) + "\n"
        ret += "".join("{!s:4}".format(-label(i)) for i in reversed(range(1,self.m+1)))
        ret += "{!s:4}".format(label(0))
        ret += "".join("{!s:4}".format(label(i)) for i in range(1,self.n+1))
        return ret

