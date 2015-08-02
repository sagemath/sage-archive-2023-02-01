r"""
Two-graphs

A two-graph on `n` points is a family `T \subset \binom {[n]}{3}`
of `3`-sets, such that any `4`-set `S\subset [n]` of size four
contains an even number of elements of `T`. Any graph `([n],E)` 
gives rise to a two-graph 
`T(E)=\{t \in \binom {[n]}{3} : | \binom {t}{2} \cap E | odd \}`, 
and any two graphs with the same two-graph can be obtained one
from the other by :meth:`Seidel switching <sage.graphs.Graph.seidel_switching>`.
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
things such as constructing the complement two-graph.

REFERENCES:

.. [BH12] A. E. Brouwer, W. H. Haemers, 
  Spectra of Graphs,
  Springer, 2012
  http://dx.doi.org/10.1007/978-1-4614-1939-6

AUTHORS:

- Dima Pasechnik (Aug 2015)

Index
-----

This module's functions are the following :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`~is_regular_twograph` | returns True if the inc. system is regular twograph
    :func:`~is_twograph`         | returns True if the inc.system is a two-graph
    :func:`~twograph_complement` | returns the complement of self 
    :func:`~twograph_descendant` | returns the descendant graph at `w` 

Functions
---------

"""
from sage.combinat.designs.incidence_structures import IncidenceStructure
from itertools import combinations
from sage.misc.functional import is_odd, is_even

def is_regular_twograph(T, alpha=False, check=False):
    """
    returns True if the inc. system is regular twograph
    """
    if check:
       if not is_twograph(T):
           if alpha:
               return False, 0
           return False
    r, (_,_,_,alpha) = T.is_t_design(t=2, return_parameters=True)
    if alpha:
        return r, alpha
    return r

def is_twograph(T):
    """
    True if the inc.system is a two-graph
    """
    return all(map(lambda f: is_even(sum(map(lambda x: x in T.blocks(),  combinations(f, 3)))),
                    combinations(T.ground_set(), 4)))

def twograph_descendant(T,v):
    """
    the descendant graph at `v`
    """
    from sage.graphs.graph import Graph
    edges = map(lambda y: frozenset(filter(lambda z: z != v, y)), filter(lambda x: v in x, T1.blocks()))
    V = T.ground_set()
    V.remove(v)
    return Graph([V, lambda i, j: frozenset((i,j)) in edges])
  
def twograph_complement(T):
    """
    the complement
    """
    Tc = filter(lambda x: not list(x) in T.blocks(), combinations(T.ground_set(), 3))
    return IncidenceStructure(T.ground_set(), Tc)
