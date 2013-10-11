r"""
Library of Hyperplane Arrangements

See :mod:`sage.geometry.hyperplane_arrangement` for details about how
to construct hyperplane arrangements.
"""
#*****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy, deepcopy
from sage.calculus.functional import expand
from sage.calculus.var import var
from sage.combinat.combinat import stirling_number2
from sage.combinat.posets.posets import Poset
from sage.functions.generalized import sign
from sage.functions.other import sqrt
from sage.geometry.polyhedron.all import Polyhedron
from sage.graphs.all import graphs
from sage.matrix.constructor import matrix, random_matrix, zero_matrix
from sage.misc.flatten import flatten
from sage.misc.prandom import random
from sage.misc.misc import powerset
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.plot.line import line
from sage.plot.colors import Color
from sage.plot.graphics import Graphics
from sage.plot.plot import plot, parametric_plot
from sage.plot.point import point
from sage.plot.text import text
from sage.plot.plot3d.parametric_plot3d import parametric_plot3d
from sage.plot.plot3d.shapes2 import text3d
from sage.rings.arith import lcm, binomial
from sage.rings.finite_rings.constructor import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR, var
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangement


class HyperplaneArrangementGenerators():

    def braid(self, n, K=QQ):
        r"""
        The braid arrangement

        This is the set of `n(n-1)/2` hyperplanes: `\{ x_i - x_j = 0 :
        1\leq i \leq j\leq n\\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default: ``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.braid(4)
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        x = polygen(QQ, 'x')
        A = self.graphical(graphs.CompleteGraph(n), K)
        A._characteristic_polynomial = prod(x-i for i in range(n))
        return A

    def bigraphical(self, G, A=None, K=QQ):
        r"""
        The hyperplane arrangement with hyperplanes `x_i - x_j = A[i,j]` and
        `x_j - x_i = A[j,i]` for each edge `v_i, v_j` of ``G``.  The indices
        `i,j` are the indices of elements of ``G.vertices()``.

        INPUT:

        - ``G`` -- Graph
        - ``A`` -- list, matrix, or dictionary (default: None gives semiorder), 'generic'
        - ``K`` -- field (default: `QQ`)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CycleGraph(4)
            sage: G.edges()
            [(0, 1, None), (0, 3, None), (1, 2, None), (2, 3, None)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 3), (1, 2), (2, 3)]
            sage: A = {0:{1:1, 3:2}, 1:{0:3, 2:0}, 2:{1:2, 3:1}, 3:{2:0, 0:2}}
            sage: HA = hyperplane_arrangements.bigraphical(G,A)
            sage: HA.num_regions()
            63
            sage: hyperplane_arrangements.bigraphical(G,'generic').num_regions()
            65
            sage: hyperplane_arrangements.bigraphical(G).num_regions()
            59

        REFERENCES::

        .. [HP] S. Hopkins, D. Perkinson
           "Bigraphical Arrangements"
           :arxiv:`1212.4398`
        """
        n = G.num_verts()
        if A is None:  # default to G-semiorder arrangement
            A = matrix(K, n, lambda i, j: 1)
        elif A == 'generic':
            A = random_matrix(ZZ, n, x=10000)
            A = matrix(K, A)
        hyperplanes = []
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            new[-1] = A[i][j]
            hyperplanes.append(new)
            new = [0]*(n+1)
            new[j] = 1
            new[i] = -1
            new[-1] = A[j][i]
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def Catalan(self, n, K=QQ):
        r"""
        The Catalan arrangement

        This is the set of `3n(n-1)/2` hyperplanes: `\{ x_i - x_j =
        -1,0,1 : 1\leq i \leq j\leq n\\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default: ``QQ`)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.Catalan(5)
            Hyperplane arrangement of 30 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [-1, 0, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        Cn = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        Cn._characteristic_polynomial = x*prod([x-n-i for i in range(1, n)])
        return Cn

    def coordinate(self, n, K=QQ):
        r"""
        The coordinate hyperplane arrangement is the central hyperplane
        arrangement consisting of the coordinate hyperplanes `x_i=0`.

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.coordinate(5)
            Hyperplane arrangement of 5 hyperplanes over Rational Field of dimension 5, rank 5.
        """
        hyperplanes = []
        for i in range(n):
            new = [0]*(n+1)
            new[i] = 1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def G_semiorder(self, G, K=QQ):
        r"""
        The semiorder hyperplane arrangement of a graph G is the arrangement `\{
        x_i - x_j = -1,1\}` where `ij` is an edge of ``G``.

        INPUT:

        - ``G`` -- graph
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_semiorder(G)
            Hyperplane arrangement of 20 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_semiorder(g)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            new[-1] = -1
            hyperplanes.append(new)
            new = deepcopy(new)
            new[-1]=1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def G_Shi(self, G, K=QQ):
        r"""
        Return the Shi hyperplane arrangement of a graph ``G``.

        INPUT:

        - ``G`` -- graph

        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_Shi(G)
            Hyperplane arrangement of 20 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_Shi(g)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: a = hyperplane_arrangements.G_Shi(graphs.WheelGraph(4))
            sage: a.show(frame=false,hyperplane_legend=false,hyperplane_opacities=0.8)
            Displaying the essentialization.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            hyperplanes.append(new)
            new = deepcopy(new)
            new[-1]=1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def graphical(self, G, K=QQ):
        r"""
        The graphical hyperplane arrangement of a graph G is the arrangement `\{
        x_i - x_j = 0\}` where `ij` is an edge of ``G``.

        INPUT:

        - ``G`` -- graph
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.graphical(G)
            Hyperplane arrangement of 10 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.graphical(g)
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            hyperplanes.append(new)
        A = HyperplaneArrangement(hyperplanes, K)
        A._characteristic_polynomial = G.chromatic_polynomial()
        return A

    def Ish(self, n, K=QQ):
        r"""
        The Ish arrangement is the set of `n(n-1)` hyperplanes: ``\{ x_i - x_j =
        0 : 1\leq i \leq j\leq n\} \cup \{x_1 - x_j = i : 1\leq i \leq j\leq n\}:.``

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: a = hyperplane_arrangements.Ish(3)
            sage: a
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 2.
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: b = hyperplane_arrangements.Shi(3)
            sage: b.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x

        REFERENCES::

        .. [AR] D. Armstrong, B. Rhoades
           "The Shi arrangement and the Ish arrangement"
           :arxiv:`1009.1655`
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*(n+1)
                new[i] = 1
                new[j] = -1
                # x_i - x_j = 0
                hyperplanes.append(new)
                # x_1 - x_j = i
                new = [0]*(n+1)
                new[0] = 1
                new[j] = -1
                new[-1] = i + 1
                hyperplanes.append(new)
                A = HyperplaneArrangement(hyperplanes, K)
                x = polygen(QQ, 'x')
                cp = sum([(-1)**k*stirling_number2(n,n-k)*prod([(x-1-j) for j in range(k,n-1)]) for k in range(0,n)])
                cp = x*cp
                cp = expand(cp)
                A._characteristic_polynomial = cp
        return A

    def linial(self, n, K=QQ):
        r"""
        The linial hyperplane arrangement is the set of hyperplanes
        ``\{x_i - x_j = 1 : 1\leq i < j \leq n\}`` 

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: a = hyperplane_arrangements.linial(4)
            sage: a.characteristic_polynomial()
            x^4 - 6*x^3 + 15*x^2 - 14*x
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1,n):
                new = [0]*(n+1)
                new[i] = 1
                new[j] = -1
                new[-1] = 1
                hyperplanes.append(new)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = expand(x*sum(binomial(n,k)*(x-k)**(n-1) for k in range(n+1))/2**n)
        A._characteristic_polynomial = cp
        return A

    def semiorder(self, n, K=QQ):
        r"""
        The semiorder arrangement is the set of `n(n-1)` hyperplanes: `\{ x_i -
        x_j = -1,1 : 1\leq i \leq j\leq n\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.semiorder(4)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [-1, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = x*sum([stirling_number2(n,k)*prod([x-k-i for i in range(1, k)]) for k in range(1,n+1)])
        cp = expand(cp)
        A._characteristic_polynomial = cp
        return A

    def Shi(self, n, K=QQ):
        r"""
        The Shi arrangement is the set of `n(n-1)` hyperplanes: ``\{ x_i - x_j =
        0,1 : 1\leq i \leq j\leq n\}.``

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.Shi(4)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [0, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = sum([(-1)**k*stirling_number2(n,n-k)*prod([(x-1-j) for j in range(k,n-1)]) for k in range(0,n)])
        cp = x*cp
        cp = expand(cp)
        A._characteristic_polynomial = cp
        return A

hyperplane_arrangements = HyperplaneArrangementGenerators()

