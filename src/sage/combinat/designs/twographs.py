r"""
Two-graphs

A two-graph on `n` points is a family `T \subset \binom {[n]}{3}`
of `3`-sets, such that any `4`-set `S\subset [n]` of size four
contains an even number of elements of `T`. Any graph `([n],E)` 
gives rise to a two-graph 
`T(E)=\{t \in \binom {[n]}{3} : \left| \binom {t}{2} \cap E \right|\ odd \}`,
and any two graphs with the same two-graph can be obtained one
from the other by :meth:`Seidel switching <Graph.seidel_switching>`.
This defines an equivalence relation on the graphs on `[n]`, 
called Seidel switching equivalence.
Conversely, given a two-graph `T`, one can construct a graph
`\Gamma` in the corresponding Seidel switching class with an 
isolated vertex `w`. The graph `\Gamma \setminus w` is called
the descendant of `T` w.r.t. `v`.

`T` is called regular if each two-subset of `[n]` is contained
in the same number alpha of triples of `T`.

This module implements a direct construction of a two-graph from a list of
triples, constrution of descendant graphs, regularity checking, and other
things such as constructing the complement two-graph, cf. [BH12]_. 

AUTHORS:

- Dima Pasechnik (Aug 2015)

Index
-----

This module's methods are the following :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~TwoGraph.is_regular_twograph` | returns True if the inc. system is regular twograph
    :meth:`~TwoGraph.complement` | returns the complement of ``self``
    :meth:`~TwoGraph.descendant` | returns the descendant graph at `w`

This module's functions are the following :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`~is_twograph`         | returns True if the incidence system is a two-graph

Methods
---------
"""
from sage.combinat.designs.incidence_structures import IncidenceStructure
from itertools import combinations
from sage.misc.functional import is_odd, is_even

class TwoGraph(IncidenceStructure):
    r"""
    Two-graphs class.

    A two-graph on `n` points is a 3-uniform hypergraph, i.e.  a family
    `T \subset \binom {[n]}{3}` of `3`-sets, such that any
    `4`-set `S\subset [n]` of size four contains an even number of elements of `T`.

    """
    def is_regular_twograph(self, alpha=False, check=False):
        """
        returns True if ``self`` is a regular twograph, i.e. a 2-design:
        each pair of elements of ``self.ground_set()`` is contained in
        exactly ``alpha`` triples.

        INPUT:

            - ``alpha`` -- (optional, default is False) return the value of ``alpha``, if possible.
            - ``check`` -- (optional, default is False), check that we actually have a two-graph.

        EXAMPLES::

            sage: p=graphs.PetersenGraph().twograph()
            sage: p.is_regular_twograph(alpha=True)
            (True, 4)
            sage: p.is_regular_twograph()
            True
            sage: p=graphs.PathGraph(5).twograph()
            sage: p.is_regular_twograph(alpha=True)
            (False, 0)
            sage: p.is_regular_twograph()
            False
        """
        if check:
           from sage.combinat.designs.twographs import is_twograph
           if not is_twograph(self):
               if alpha:
                   return False, 0
               return False
        r, (_,_,_,a) = self.is_t_design(t=2, k=3, return_parameters=True)
        if alpha:
            return r, a
        return r

    def descendant(self, v):
        """
        the descendant graph at ``v``

        The :mod:`switching class of graphs <sage.combinat.designs.twographs>`
        corresponding to ``self`` contains a graph ``D`` with ``v`` its own connected
        component; removing ``v`` from ``D``, one obtains the descendant graph of
        ``self`` at ``v``, which is constructed by this method.

        INPUT:

            - ``v`` -- an element of ``self.ground_set()`` 

        OUTPUT:

            - the descendant :class:`graph <sage.graphs.graph.Graph>` at ``v``
 
        EXAMPLES::

            sage: p=graphs.PetersenGraph().twograph().descendant(0)
            sage: p.is_strongly_regular(parameters=True)
            (9, 4, 1, 2)
        """
        from sage.graphs.graph import Graph
        return Graph(map(lambda y: filter(lambda z: z != v, y),
                            filter(lambda x: v in x, self.blocks())))

    def complement(self):
        """
        the complement of ``self``

        The two-graph constisting exactly of triples not in ``self``.

        EXAMPLES::

            sage: p=graphs.CompleteGraph(8).line_graph().twograph()
            sage: pc = p.complement(); pc
            Incidence structure with 28 points and 1260 blocks

        TESTS::

            sage: from sage.combinat.designs.twographs import is_twograph
            sage: is_twograph(pc)
            True
        """
        return super(TwoGraph, self).complement(uniform=True)


def is_twograph(T):
    """
    True if the incidence system is a two-graph

    EXAMPLES::

        sage: from sage.combinat.designs.twographs import is_twograph
        sage: p=graphs.PetersenGraph().twograph()
        sage: is_twograph(p)
        True
        sage: is_twograph(designs.projective_plane(3))
        False
        sage: is_twograph(designs.projective_plane(2))
        False
    """
    B = map(frozenset, T.blocks())
    return T.is_t_design(k=3) and \
        all(map(lambda f: is_even(sum(map(lambda x: frozenset(x) in B,  combinations(f, 3)))),
                    combinations(T.ground_set(), 4)))
