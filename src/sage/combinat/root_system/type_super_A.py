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
from __future__ import print_function, absolute_import

from sage.rings.all import ZZ
from sage.misc.cachefunc import cached_method
from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
from . import ambient_space
from .cartan_type import SuperCartanType_standard
from six import iteritems

class AmbientSpace(ambient_space.AmbientSpace):
    r"""
    EXAMPLES::

        sage: R = RootSystem(['A', [4,2]])
        sage: e = R.ambient_space(); e
        Ambient space of the Root system of type ['A', [4, 2]]
        sage: TestSuite(e).run()
    """
    def __init__(self, root_system, base_ring):
        ct = root_system.cartan_type()
        I = tuple(range(-ct.m-1,0) + range(1,ct.n+2))
        ambient_space.AmbientSpace.__init__(self, root_system, base_ring,
                                            index_set=I)

    @classmethod
    def smallest_base_ring(cls, cartan_type=None):
        """
        Return the smallest base ring the ambient space can be defined upon

        .. SEEALSO::

            :meth:`~sage.combinat.root_system.ambient_space.AmbientSpace.smallest_base_ring`

        EXAMPLES::

            sage: e = RootSystem(['A',[3,1]]).ambient_space()
            sage: e.smallest_base_ring()
            Integer Ring
        """
        return ZZ

    def dimension(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A', [4,2]]).ambient_space()
            sage: e.dimension()
            8
        """
        ct = self.root_system.cartan_type()
        return ct.m + ct.n + 2

    def simple_root(self, i):
        """
        EXAMPLES::  

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: list(e.simple_roots())
            [(1, -1, 0, 0, 0), (0, 1, -1, 0, 0),
             (0, 0, 1, -1, 0), (0, 0, 0, 1, -1)]
        """
        if i < 0:
            return self.monomial(i-1) - self.monomial(i)
        if i == 0:
            return self.monomial(-1) - self.monomial(1)
        return self.monomial(i) - self.monomial(i+1)

    def positive_roots(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_roots()
            [(0, 1, -1, 0, 0),
             (1, 0, -1, 0, 0),
             (1, -1, 0, 0, 0),
             (0, 0, 0, 1, -1),
             (0, 0, 1, -1, 0),
             (0, 0, 1, 0, -1),
             (0, 1, 0, -1, 0),
             (0, 1, 0, 0, -1),
             (1, 0, 0, -1, 0),
             (1, 0, 0, 0, -1)]
        """
        return self.positive_even_roots() + self.positive_odd_roots()

    def positive_even_roots(self):
        """
        Return the positive even roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_even_roots()
            [(0, 1, -1, 0, 0), (1, 0, -1, 0, 0),
             (1, -1, 0, 0, 0), (0, 0, 0, 1, -1)]
        """
        ct = self.root_system.cartan_type()
        ret = []
        ret += [self.monomial(-j) - self.monomial(-i)
                for i in range(1,ct.m+2)
                for j in range(i+1,ct.m+2)]
        ret += [self.monomial(i) - self.monomial(j)
                for i in range(1,ct.n+2)
                for j in range(i+1,ct.n+2)]
        return ret

    def positive_odd_roots(self):
        """
        Return the positive odd roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.positive_odd_roots()
            [(0, 0, 1, -1, 0),
             (0, 0, 1, 0, -1),
             (0, 1, 0, -1, 0),
             (0, 1, 0, 0, -1),
             (1, 0, 0, -1, 0),
             (1, 0, 0, 0, -1)]
        """
        ct = self.root_system.cartan_type()
        return [self.monomial(-i) - self.monomial(j)
                for i in range(1,ct.m+2)
                for j in range(1,ct.n+2)]

    def highest_root(self):
        """
        EXAMPLES::

           sage: e = RootSystem(['A', [4,2]]).ambient_lattice()
           sage: e.highest_root()
           (1, 0, 0, 0, 0, 0, 0, -1)
        """
        ct = self.root_system.cartan_type()
        return self.monomial(-ct.m-1) - self.monomial(ct.n+1)

    def negative_roots(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_roots()
            [(0, -1, 1, 0, 0),
             (-1, 0, 1, 0, 0),
             (-1, 1, 0, 0, 0),
             (0, 0, 0, -1, 1),
             (0, 0, -1, 1, 0),
             (0, 0, -1, 0, 1),
             (0, -1, 0, 1, 0),
             (0, -1, 0, 0, 1),
             (-1, 0, 0, 1, 0),
             (-1, 0, 0, 0, 1)]
        """
        return self.negative_even_roots() + self.negative_odd_roots()

    def negative_even_roots(self):
        """
        Return the negative even roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_even_roots()
            [(0, -1, 1, 0, 0), (-1, 0, 1, 0, 0),
             (-1, 1, 0, 0, 0), (0, 0, 0, -1, 1)]
        """
        ct = self.root_system.cartan_type()
        ret = []
        ret += [self.monomial(-i) - self.monomial(-j)
                for i in range(1,ct.m+2)
                for j in range(i+1,ct.m+2)]
        ret += [self.monomial(j) - self.monomial(i)
                for i in range(1,ct.n+2)
                for j in range(i+1,ct.n+2)]
        return ret

    def negative_odd_roots(self):
        """
        Return the negative odd roots of ``self``.

        EXAMPLES::

            sage: e = RootSystem(['A', [2,1]]).ambient_lattice()
            sage: e.negative_odd_roots()
            [(0, 0, -1, 1, 0),
             (0, 0, -1, 0, 1),
             (0, -1, 0, 1, 0),
             (0, -1, 0, 0, 1),
             (-1, 0, 0, 1, 0),
             (-1, 0, 0, 0, 1)]
        """
        ct = self.root_system.cartan_type()
        return [self.monomial(j) - self.monomial(-i)
                for i in range(1,ct.m+2)
                for j in range(1,ct.n+2)]

    def fundamental_weight(self, i):
        """
        Return the fundamental weight `\Lambda_i` of ``self``.

        EXAMPLES::

            sage: L = RootSystem(['A', [3,2]]).ambient_space()
            sage: L.fundamental_weight(-1)
            (1, 1, 1, 0, 0, 0, 0)
            sage: L.fundamental_weight(0)
            (1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2)
            sage: L.fundamental_weight(2)
            (1, 1, 1, 1, -1, -1, 0)
            sage: list(L.fundamental_weights())
            [(1, 0, 0, 0, 0, 0, 0),
             (1, 1, 0, 0, 0, 0, 0),
             (1, 1, 1, 0, 0, 0, 0),
             (1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2),
             (1, 1, 1, 1, -1, 0, 0),
             (1, 1, 1, 1, -1, -1, 0)]
        """
        m = self.root_system.cartan_type().m
        if i < 0:
            return self.sum(self.monomial(j) for j in range(-m-1,i))
        if i == 0:
            return self.sum(self.basis()) / 2
        return (self.sum(self.monomial(j) for j in range(-m-1,1))
                - self.sum(self.monomial(j) for j in range(i+1)))

    def simple_coroot(self, i):
        """
        Return the simple coroot `h_i` of ``self``.

        EXAMPLES::

            sage: L = RootSystem(['A', [3,2]]).ambient_space()
            sage: L.simple_coroot(-2)
            (0, 1, -1, 0, 0, 0, 0)
            sage: L.simple_coroot(0)
            (0, 0, 0, 1, -1, 0, 0)
            sage: L.simple_coroot(2)
            (0, 0, 0, 0, 0, -1, 1)
            sage: list(L.simple_coroots())
            [(1, -1, 0, 0, 0, 0, 0),
             (0, 1, -1, 0, 0, 0, 0),
             (0, 0, 1, -1, 0, 0, 0),
             (0, 0, 0, 1, -1, 0, 0),
             (0, 0, 0, 0, -1, 1, 0),
             (0, 0, 0, 0, 0, -1, 1)]
        """
        if i <= 0:
            return self.simple_root(i)
        return -self.simple_root(i)

    class Element(ambient_space.AmbientSpaceElement):
        def inner_product(self, lambdacheck):
            """
            The scalar product with elements of the coroot lattice
            embedded in the ambient space.

            EXAMPLES::

                sage: e = RootSystem(['A',2]).ambient_space()
                sage: a = e.simple_root(0); a
                (-1, 0, 0)
                sage: a.inner_product(a)
                2
            """
            self_mc = self._monomial_coefficients
            lambdacheck_mc = lambdacheck._monomial_coefficients

            result = self.parent().base_ring().zero()
            for t,c in iteritems(lambdacheck_mc):
                if t not in self_mc:
                    continue
                if t > 0:
                    result -= c*self_mc[t]
                else:
                    result += c*self_mc[t]
            return result

        scalar = inner_product
        dot_product = inner_product

        def associated_coroot(self):
            """
            Return the coroot associated to ``self``, which is
            precisely ``self``.

            EXAMPLES::

                sage: L = RootSystem(['A', [3,2]]).ambient_space()
                sage: al = L.simple_roots()
                sage: al[-1].associated_coroot()
                (0, 0, 1, -1, 0, 0, 0)
                sage: al[0].associated_coroot()
                (0, 0, 0, 1, -1, 0, 0)
                sage: al[1].associated_coroot()
                (0, 0, 0, 0, -1, 1, 0)

                sage: a = al[-1] + al[0] + al[1]; a
                (0, 0, 1, 0, 0, -1, 0)
                sage: a.associated_coroot()
                (0, 0, 1, 0, -2, 1, 0)
                sage: h = L.simple_coroots()
                sage: h[-1] + h[0] + h[1]
                (0, 0, 1, 0, -2, 1, 0)

                sage: (al[-1] + al[0] + al[2]).associated_coroot()
                (0, 0, 1, 0, -1, -1, 1)
            """
            P = self.parent()
            al = P.simple_roots()
            h = P.simple_coroots()
            try:
                return h[al.inverse_family()[self]]
            except KeyError:
                pass
            V = P._dense_free_module()
            dep = V.linear_dependence([self._vector_()] +
                                      [al[i]._vector_() for i in P.index_set()])[0]
            I = P.index_set()
            return P.sum((-c/dep[0]) * h[I[i]] for i,c in dep[1:].iteritems())

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
            sage: ct._repr_(compact=True)
            'A4|2'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.affine() # Not tested -- to be implemented
            ['A', [4, 2], 1]
            sage: ct.dual()
            ['A', [4, 2]]

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

    # A lot of these methods should be implemented by the ABCs of CartanType

    def is_affine(self):
        return False

    def is_finite(self):
        return True

    def dual(self):
        return self

    def type(self):
        return 'A'

    def root_system(self):
        from sage.combinat.root_system.root_system import RootSystem
        return RootSystem(self)

    @cached_method
    def symmetrizer(self):
        from sage.sets.family import Family
        def ell(i): return ZZ.one() if i <= 0 else -ZZ.one()
        return Family(self.index_set(), ell)

    def dynkin_diagram(self): # FIXME
        """
        Returns the Dynkin diagram of super type A.

        EXAMPLES::

            sage: a = CartanType(['A', [4,2]]).dynkin_diagram()
            sage: a
            O---O---O---O---X---O---O
            -4  -3  -2  -1  0   1   2   
            A4|2
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

    def _latex_dynkin_diagram(self, label=lambda i: i, node=None, node_dist=2): # FIXME
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

