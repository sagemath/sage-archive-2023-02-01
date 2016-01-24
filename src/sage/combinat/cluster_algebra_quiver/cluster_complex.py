r"""
Cluster complex (or generalized dual associahedron)

AUTHORS:

- Christian Stump
"""
#*****************************************************************************
#       Copyright (C) 2011      Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.homology.simplicial_complex import Simplex
from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.subword_complex import SubwordComplex, SubwordComplexFacet
from sage.rings.semirings.non_negative_integer_semiring import NN


class ClusterComplex(SubwordComplex):
    r"""
    Cluster complex (or generalized dual associahedron)

    The cluster complex (or generalized dual associahedron) is a
    simplicial complex constructed from a cluster algebra.  Its
    vertices are the cluster variables and its facets are the
    clusters, i.e., maximal subsets of compatible cluster variables.

    The cluster complex of type `A_n` is the simplicial complex with
    vertices being (proper) diagonals in a convex `(n+3)`-gon and with
    facets being triangulations.

    The implementation of the cluster complex depends on its
    connection to subword complexes, see [CLS]_. Let `c` be a Coxeter
    element with reduced word `{\bf c}` in a finite Coxeter group `W`,
    and let `{\bf w}_\circ` be the `c`-sorting word for the longest
    element `w_\circ \in W`.

    The ``multi-cluster complex`` `\Delta(W,k)` has vertices in
    one-to-one correspondence with letters in the word `Q = {\bf c^k
    w}_\circ` and with facets being complements in `Q` of reduced
    expressions for `w_\circ`.

    For `k = 1`, the multi-cluster complex is isomorphic to the
    cluster complex as defined above.

    REFERENCES:

    .. [CLS] C. Ceballos, J.-P. Labbe, C. Stump, ``Subword complexes,
       cluster complexes, and generalized multi-associahedra``,
       :arxiv:`1108.1776`.

    EXAMPLES:

    A first example of a cluster complex::

        sage: C = ClusterComplex(['A', 2]); C
        Cluster complex of type ['A', 2] with 5 vertices and 5 facets

    Its vertices, facets, and minimal non-faces::

        sage: C.vertices()
        (0, 1, 2, 3, 4)

        sage: C.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

        sage: C.minimal_nonfaces()
        [[0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]

    We can do everything we can do on simplicial complexes,
    e.g. computing its homology::

        sage: C.homology()
        {0: 0, 1: Z}
    """
    def __init__(self, W, k=1, coxeter_element=None, algorithm="inductive"):
        """
        EXAMPLES::

            sage: ClusterComplex(['A', 2])
            Cluster complex of type ['A', 2] with 5 vertices and 5 facets
            sage: ClusterComplex(['A', 2], k=2)
            Multi-cluster complex of type ['A', 2] with 7 vertices and 14 facets
        """
        if not k in NN:
            raise ValueError("The additional parameter must be a "
                             "nonnegative integer.")

        if W not in CoxeterGroups:
            W = CoxeterGroup(W)
        if not W.is_finite():
            raise ValueError("The Coxeter group must be finite.")

        if coxeter_element is None:
            c = W.from_reduced_word(W.index_set())
        else:
            c = W.from_reduced_word(coxeter_element)
        w = W.w0
        Q = c.reduced_word() * k + w.coxeter_sorting_word(c)
        SubwordComplex.__init__(self, Q, w, algorithm=algorithm)
        self._W = W
        self._w0 = w
        self._k = k
        if k == 1:
            self.__custom_name = 'Cluster complex'
        else:
            self.__custom_name = 'Multi-cluster complex'

    def __call__(self, F, facet_test=True):
        """
        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: C((0, 1))
            (0, 1)
        """
        return ClusterComplexFacet(self, F, facet_test=facet_test)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2])._repr_()
            "Cluster complex of type ['A', 2] with 5 vertices and 5 facets"
        """
        name = self.__custom_name
        name += ' of type %s with %s vertices and %s facets' \
                % (self.cartan_type(), self.vertices().dimension() + 1,
                   len(self._facets))
        return name

    def k(self):
        r"""
        Return the index `k` of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2]).k()
            1
        """
        return self._k

    def minimal_nonfaces(self):
        """
        Return the minimal non-faces of ``self``.

        EXAMPLES::

            sage: ClusterComplex(['A', 2]).minimal_nonfaces()
            [[0, 2], [0, 3], [1, 3], [1, 4], [2, 4]]
        """
        from sage.combinat.combination import Combinations
        return [X for X in Combinations(self.vertices(), self.k() + 1)
                if not any(set(X).issubset(F) for F in self.facets())]

    def cyclic_action(self):
        """
        EXAMPLES::

            sage: ClusterComplex(['A', 2]).cyclic_action()
            <function act at ...>
        """
        W = self._W
        w = self._w0
        Q = self._Q
        l = len(Q)
        S = W.simple_reflections()
        S_inv = {S[j]: j for j in W.index_set()}
        Q += [S_inv[w * S[k] * w] for k in Q]
        D = {i: (Q[i + 1:].index(Q[i]) + i + 1) % l for i in range(l)}

        def act(F):
            return Simplex(sorted([D[i] for i in F]))
        return act


class ClusterComplexFacet(SubwordComplexFacet):
    def cluster(self):
        """
        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.cluster()
            [(-1, 0), (0, -1)]
        """
        if self.parent().k() != 1:
            raise NotImplementedError("only working for k=1")
        F = self.parent().greedy_facet(side="positive")
        R = F.extended_root_configuration()
        return [-R[i] if i < len(list(F)) else R[i] for i in self]

    def upper_cluster(self):
        """
        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.upper_cluster()
            []
        """
        conf = self._root_configuration_indices()
        W = self.parent().group()
        N = len(W.roots()) / 2
        C = self.cluster()
        return [C[i] for i in xrange(len(conf)) if conf[i] >= N]

    def modified_upper_cluster(self):
        """
        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.modified_upper_cluster()
            []
        """
        W = self.parent().group()
        Upp = self.upper_cluster()
        return [W.prod(W.root_to_reflection(beta)
                       for beta in reversed(Upp[m + 1:])).act_on_root(Upp[m])
                for m in range(len(Upp))]

    def product_of_upper_cluster(self):
        """
        EXAMPLES::

            sage: C = ClusterComplex(['A', 2])
            sage: F = C((0, 1))
            sage: F.product_of_upper_cluster()
            [1 0]
            [0 1]
        """
        W = self.parent().group()
        return W.prod(W.root_to_reflection(beta)
                      for beta in reversed(self.upper_cluster()))
