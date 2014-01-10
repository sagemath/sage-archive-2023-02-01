r"""
Right-Angled Artin Groups

AUTHORS:

- Travis Scrimshaw (2013-09-01): Initial version
"""

##############################################################################
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.misc.cachefunc import cached_method
from sage.structure.list_clone import ClonableArray
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.group import Group
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.graphs.graph import Graph

class RightAngledArtinGroup(Group, UniqueRepresentation):
    r"""
    The right-angled Artin group defined by a graph `G`.

    Let `\Gamma = \{V(\Gamma), E(\Gamma)\}` be a simple graph.
    A *right-angled Artin group* (commonly abbriated as RAAG) is the group

    .. MATH::

        A_{\Gamma} = \langle g_v : v \in V(\Gamma)
        \mid [g_u, g_v] \text{ if } \{u, v\} \notin E(\Gamma) \rangle.

    These are sometimes known as graph groups or partitally commutative groups.
    This RAAG's contains both free groups, given by the complete graphs,
    and free abelian groups, given by disjoint vertices.

    .. NOTE::

        This is the opposite convention of some papers.

    EXAMPLES::

        sage: Gamma = Graph(4)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: a,b,c,d = G.gens()
        sage: a*c*d^4*a^-3*b
        v0^-2*v1*v2*v3^4

        sage: Gamma = graphs.CompleteGraph(4)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: a,b,c,d = G.gens()
        sage: a*c*d^4*a^-3*b
        v0*v2*v3^4*v0^-3*v1

        sage: Gamma = graphs.CycleGraph(5)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: G
        Right-angled Artin group on 5 generators
        sage: a,b,c,d,e = G.gens()
        sage: e^-1*c*b*e*b^-1*c^-4
        v2^-3
    """
    @staticmethod
    def __classcall_private__(cls, G):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: G1 = RightAngledArtinGroup(graphs.CycleGraph(5))
            sage: Gamma = Graph([(0,1),(1,2),(2,3),(3,4),(4,0)])
            sage: G2 = RightAngledArtinGroup(Gamma)
            sage: G3 = RightAngledArtinGroup([(0,1),(1,2),(2,3),(3,4),(4,0)])
            sage: G1 is G2 and G2 is G3
            True
        """
        if not isinstance(G, Graph):
            G = Graph(G, immutable=True)
        else:
            G = G.copy(immutable=True)
        return super(RightAngledArtinGroup, cls).__classcall__(cls, G)

    def __init__(self, G):
        """
        Initialize ``self``.

        INPUT:

        - ``G`` -- a graph

        TESTS::

            sage: G = RightAngledArtinGroup(graphs.CycleGraph(5))
            sage: TestSuite(G).run()
        """
        self._graph = G
        Group.__init__(self)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: RightAngledArtinGroup(graphs.CycleGraph(5))
            Right-angled Artin group on 5 generators
        """
        return "Right-angled Artin group on {} generators".format(self._graph.num_verts())

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gen(2)
            v2
        """
        return self.element_class(self, [(i, 1)])

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gens()
            (v0, v1, v2, v3, v4)
            sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gens()
            (vx, vy, vzeta)
        """
        return tuple(self.gen(i) for i in range(self._graph.num_verts()))

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.ngens()
            5
        """
        return self._graph.num_verts()

    def cardinality(self):
        """
        Return the number of group elements.

        OUTPUT:

        Infinity.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.cardinality()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    order = cardinality
    
    def as_permutation_group(self):
        """
        Raises a ``ValueError`` error since right-angled Artin groups
        are infinite, so they have no isomorphic permutation group.
        
        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.as_permutation_group()
            Traceback (most recent call last):
            ...
            ValueError: the group is infinite
        """
        raise ValueError("the group is infinite")

    def graph(self):
        """
        Return the defining graph of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.graph()
            Graph on 5 vertices
        """
        return self._graph

    @cached_method
    def one(self):
        """
        Return the identity element `1`.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.one()
            1
        """
        return self.element_class(self, [])

    one_element = one

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        TESTS::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: elt = G([[0,3], [3,1], [2,1], [1,1], [3,1]]); elt
            v0^3*v3*v2*v1*v3
            sage: G(elt)
            v0^3*v3*v2*v1*v3
            sage: G(1)
            1
        """
        if isinstance(x, RightAngledArtinGroup.Element):
            if x.parent() is self:
                return x
            raise ValueError("there is no coercion from {} into {}".format(x.parent(), self))
        if x == 1:
            return self.one()
        verts = self._graph.vertices()
        x = map(lambda s: [verts.index(s[0]), s[1]], x)
        return self.element_class(self, self._normal_form(x))

    def _normal_form(self, word):
        """
        Return the normal form of the word ``word``. Helper function for
        creaing elements.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G._normal_form([[0,2], [3,1], [2,1], [0,1], [1,1], [3,1]])
            [[0, 3], [3, 1], [2, 1], [1, 1], [3, 1]]
            sage: a,b,c,d,e = G.gens()
            sage: a^2 * d * c * a * b * d
            v0^3*v3*v2*v1*v3
            sage: a*b*d == d*a*b and a*b*d == a*d*b
            True
            sage: a*c*a^-1*c^-1
            1
            sage: (a*b*c*d*e)^2 * (a*b*c*d*e)^-2
            1
        """
        pos = 0
        G = self._graph
        v = G.vertices()
        w = map(list, word) # Make a (2 level) deep copy
        while pos < len(w):
            comm_set = [w[pos][0]] # The current set of totally commuting elements
            i = pos + 1

            while i < len(w):
                letter = w[i][0] # The current letter
                # Check if this could fit in the commuting set
                if letter in comm_set:
                    # Try to move it in
                    if any(G.has_edge(v[w[j][0]], v[letter]) for j in range(pos + len(comm_set), i)):
                        # We can't, so go onto the next letter
                        i += 1
                        continue
                    j = comm_set.index(letter)
                    w[pos+j][1] += w[i][1]
                    w.pop(i)
                    i -= 1 # Since we removed a syllable
                    # Check cancellations
                    if w[pos+j][1] == 0:
                        w.pop(pos+j)
                        comm_set.pop(j)
                        i -= 1
                        if len(comm_set) == 0:
                            pos = 0 # Start again since cancellation can be pronounced effects
                            break
                elif all( not G.has_edge(v[w[j][0]], v[letter]) for j in range(pos, i)):
                    j = 0
                    for x in comm_set:
                        if x > letter:
                            break
                        j += 1
                    w.insert(pos+j, w.pop(i))
                    comm_set.insert(j, letter)

                i += 1
            pos += len(comm_set)
        return w

    class Element(ClonableArray):
        """
        An element of a right-angled Artin group (RAAG).

        Elements of RAAGs are modeled as lists of pairs ``[i, p]`` where
        ``i`` is the index of a vertex in the defining graph (with some
        fixed order of the vertices) and ``p`` is the power.
        """
        def check(self):
            """
            Check if ``self`` is a valid element. Nothing to check.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: elt = G.gen(2)
                sage: elt.check()
            """
            pass

        def _repr_(self):
            """
            Return a string representation of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: a * b^2 * e^-3
                v0*v1^2*v4^-3
                sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
                sage: G = RightAngledArtinGroup(Gamma)
                sage: x,y,z = G.gens()
                sage: z * y^-2 * x^3
                vzeta*vy^-2*vx^3
            """
            if len(self) == 0:
                return '1'
            v = self.parent()._graph.vertices()
            to_str = lambda i,p: "v{}".format(i) if p == 1 else "v{}^{}".format(i, p)
            return '*'.join(to_str(v[i], p) for i,p in self)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: latex(a*b*e^-4*d^3)
                \sigma_{0}\sigma_{1}\sigma_{4}^{-4}\sigma_{3}^{3}
                sage: latex(G.one())
                1
                sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
                sage: G = RightAngledArtinGroup(Gamma)
                sage: x,y,z = G.gens()
                sage: latex(x^-5*y*z^3)
                \sigma_{\text{\texttt{x}}}^{-5}\sigma_{\text{\texttt{y}}}\sigma_{\text{\texttt{zeta}}}^{3}
            """
            if len(self) == 0:
                return '1'

            from sage.misc.latex import latex
            latexrepr = ''
            v = self.parent()._graph.vertices()
            for i,p in self:
                latexrepr += "\\sigma_{{{}}}".format(latex(v[i]))
                if p != 1:
                    latexrepr += "^{{{}}}".format(p)
            return latexrepr

        def _mul_(self, y):
            """
            Return ``self`` multiplied by ``y``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: a * b
                v0*v1
                sage: b * a
                v1*v0
                sage: a*b*c*d*e
                v0*v1*v2*v3*v4
                sage: a^2*d*c*a*b*d
                v0^3*v3*v2*v1*v3
                sage: e^-1*a*b*d*c*a^-2*e*d*b^2*e*b^-3
                v4^-1*v0*v3*v1*v0^-2*v2*v1^-1*v4*v3*v4
            """
            P = self.parent()
            lst = list(self) + list(y)
            return self.__class__(self.parent(), P._normal_form(lst))

        def __invert__(self):
            """
            Return the inverse of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: (a * b)^-2
                v1^-1*v0^-1*v1^-1*v0^-1
            """
            return self.__class__(self.parent(), map(lambda x: [x[0], -x[1]], reversed(self)))

