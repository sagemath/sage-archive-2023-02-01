r"""
Vertex separation

This module implements several algorithms to compute the vertex separation of a
digraph and the corresponding ordering of the vertices. It also implements tests
functions for evaluation the width of a linear ordering.

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |N^+_G(S)\backslash S|

The *vertex separation* of a digraph `G` is equal to the minimum cost of an
ordering of its vertices.

**Vertex separation and pathwidth**

The vertex separation is defined on a digraph, but one can obtain from a graph
`G` a digraph `D` with the same vertex set, and in which each edge `uv` of `G`
is replaced by two edges `uv` and `vu` in `D`. The vertex separation of `D` is
equal to the pathwidth of `G`, and the corresponding ordering of the vertices of
`D`, also called a *layout*, encodes an optimal path-decomposition of `G`.
This is a result of Kinnersley [Kin92]_ and Bodlaender [Bod98]_.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`path_decomposition` | Returns the pathwidth of the given graph and the ordering of the vertices resulting in a corresponding path decomposition
    :meth:`vertex_separation` | Returns an optimal ordering of the vertices and its cost for vertex-separation
    :meth:`vertex_separation_MILP` | Computes the vertex separation of `G` and the optimal ordering of its vertices using an MILP formulation
    :meth:`lower_bound` | Returns a lower bound on the vertex separation of `G`
    :meth:`is_valid_ordering` | Test if the linear vertex ordering `L` is valid for (di)graph `G`
    :meth:`width_of_path_decomposition` | Returns the width of the path decomposition induced by the linear ordering `L` of the vertices of `G`


Exponential algorithm for vertex separation
-------------------------------------------

In order to find an optimal ordering of the vertices for the vertex separation,
this algorithm tries to save time by computing the function `c'(S)` **at most
once** once for each of the sets `S\subseteq V(G)`. These values are stored in
an array of size `2^n` where reading the value of `c'(S)` or updating it can be
done in constant (and small) time.

Assuming that we can compute the cost of a set `S` and remember it, finding an
optimal ordering is an easy task. Indeed, we can think of the sequence `v_1,
..., v_n` of vertices as a sequence of *sets* `\{v_1\}, \{v_1,v_2\}, ...,
\{v_1,...,v_n\}`, whose cost is precisely `\max c'(\{v_1\}), c'(\{v_1,v_2\}),
... , c'(\{v_1,...,v_n\})`. Hence, when considering the digraph on the `2^n`
sets `S\subseteq V(G)` where there is an arc from `S` to `S'` if `S'=S\cap
\{v\}` for some `v` (that is, if the sets `S` and `S'` can be consecutive in a
sequence), an ordering of the vertices of `G` corresponds to a *path* from
`\emptyset` to `\{v_1,...,v_n\}`. In this setting, checking whether there exists
a ordering of cost less than `k` can be achieved by checking whether there
exists a directed path `\emptyset` to `\{v_1,...,v_n\}` using only sets of cost
less than `k`. This is just a depth-first-search, for each `k`.

**Lazy evaluation of** `c'`

In the previous algorithm, most of the time is actually spent on the computation
of `c'(S)` for each set `S\subseteq V(G)` -- i.e. `2^n` computations of
neighborhoods. This can be seen as a huge waste of time when noticing that it is
useless to know that the value `c'(S)` for a set `S` is less than `k` if all the
paths leading to `S` have a cost greater than `k`. For this reason, the value of
`c'(S)` is computed lazily during the depth-first search. Explanation :

When the depth-first search discovers a set of size less than `k`, the costs of
its out-neighbors (the potential sets that could follow it in the optimal
ordering) are evaluated. When an out-neighbor is found that has a cost smaller
than `k`, the depth-first search continues with this set, which is explored with
the hope that it could lead to a path toward `\{v_1,...,v_n\}`. On the other
hand, if an out-neighbour has a cost larger than `k` it is useless to attempt to
build a cheap sequence going though this set, and the exploration stops
there. This way, a large number of sets will never be evaluated and *a lot* of
computational time is saved this way.

Besides, some improvement is also made by "improving" the values found by
`c'`. Indeed, `c'(S)` is a lower bound on the cost of a sequence containing the
set `S`, but if all out-neighbors of `S` have a cost of `c'(S) + 5` then one
knows that having `S` in a sequence means a total cost of at least `c'(S) +
5`. For this reason, for each set `S` we store the value of `c'(S)`, and replace
it by `\max (c'(S), \min_{\text{next}})` (where `\min_{\text{next}}` is the
minimum of the costs of the out-neighbors of `S`) once the costs of these
out-neighbors have been evaluated by the algrithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on less than 32 vertices. This can be changed to 64 if necessary, but 32
    vertices already require 4GB of memory. Running it on 64 bits is not
    expected to be doable by the computers of the next decade `:-D`

**Lower bound on the vertex separation**

One can obtain a lower bound on the vertex separation of a graph in exponential
time but *small* memory by computing once the cost of each set `S`. Indeed, the
cost of a sequence `v_1, ..., v_n` corresponding to sets `\{v_1\}, \{v_1,v_2\},
..., \{v_1,...,v_n\}` is

.. MATH::

    \max c'(\{v_1\}),c'(\{v_1,v_2\}),...,c'(\{v_1,...,v_n\})\geq\max c'_1,...,c'_n

where `c_i` is the minimum cost of a set `S` on `i` vertices. Evaluating the
`c_i` can take time (and in particular more than the previous exact algorithm),
but it does not need much memory to run.


MILP formulation for the vertex separation
------------------------------------------

We describe below a mixed integer linear program (MILP) for determining an
optimal layout for the vertex separation of `G`, which is an improved version of
the formulation proposed in [SP10]_. It aims at building a sequence `S_t` of
sets such that an ordering `v_1, ..., v_n` of the vertices correspond to
`S_0=\{v_1\}, S_2=\{v_1,v_2\}, ..., S_{n-1}=\{v_1,...,v_n\}`.

**Variables:**


- `y_v^t` -- Variable set to 1 if `v\in S_t`, and 0 otherwise. The order of
  `v` in the layout is the smallest `t` such that `y_v^t==1`.

- `u_v^t` -- Variable set to 1 if `v\not \in S_t` and `v` has an in-neighbor in
  `S_t`. It is set to 0 otherwise.

- `x_v^t` -- Variable set to 1 if either `v\in S_t` or if `v` has an in-neighbor
  in `S_t`. It is set to 0 otherwise.

- `z` -- Objective value to minimize. It is equal to the maximum over all step
  `t` of the number of vertices such that `u_v^t==1`.

**MILP formulation:**

.. MATH::
    :nowrap:

    \begin{alignat}{2}
    \intertext{Minimize:}
    &z&\\
    \intertext{Such that:}
    x_v^t &\leq x_v^{t+1}& \forall v\in V,\ 0\leq t\leq n-2\\
    y_v^t &\leq y_v^{t+1}& \forall v\in V,\ 0\leq t\leq n-2\\
    y_v^t &\leq x_w^t& \forall v\in V,\ \forall w\in N^+(v),\ 0\leq t\leq n-1\\
    \sum_{v \in V} y_v^{t} &= t+1& 0\leq t\leq n-1\\
    x_v^t-y_v^t&\leq u_v^t & \forall v \in V,\ 0\leq t\leq n-1\\
    \sum_{v \in V} u_v^t &\leq z& 0\leq t\leq n-1\\
    0 \leq x_v^t &\leq 1& \forall v\in V,\ 0\leq t\leq n-1\\
    0 \leq u_v^t &\leq 1& \forall v\in V,\ 0\leq t\leq n-1\\
    y_v^t &\in \{0,1\}& \forall v\in V,\ 0\leq t\leq n-1\\
    0 \leq z &\leq n&
    \end{alignat}

The vertex separation of `G` is given by the value of `z`, and the order of
vertex `v` in the optimal layout is given by the smallest `t` for which
`y_v^t==1`.

REFERENCES
----------

.. [Bod98] *A partial k-arboretum of graphs with bounded treewidth*, Hans
  L. Bodlaender, Theoretical Computer Science 209(1-2):1-45, 1998.

.. [Kin92] *The vertex separation number of a graph equals its path-width*,
  Nancy G. Kinnersley, Information Processing Letters 42(6):345-350, 1992.

.. [SP10] *Lightpath Reconfiguration in WDM networks*, Fernando Solano and
  Michal Pioro, IEEE/OSA Journal of Optical Communication and Networking
  2(12):1010-1021, 2010.


AUTHORS
-------

- Nathann Cohen (2011-10): Initial version and exact exponential algorithm

- David Coudert (2012-04): MILP formulation and tests functions



METHODS
-------
"""

include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'
include 'sage/ext/interrupt.pxi'
include 'fast_digraph.pyx'
from libc.stdint cimport uint8_t, int8_t

#*****************************************************************************
#          Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

###############
# Lower Bound #
###############

def lower_bound(G):
    r"""
    Returns a lower bound on the vertex separation of `G`

    INPUT:

    - ``G`` -- a digraph

    OUTPUT:

    A lower bound on the vertex separation of `D` (see the module's
    documentation).

    .. NOTE::

        This method runs in exponential time but has no memory constraint.


    EXAMPLE:

    On a circuit::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = digraphs.Circuit(6)
        sage: lower_bound(g)
        1

    TEST:

    Given anything else than a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = graphs.CycleGraph(5)
        sage: lower_bound(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a DiGraph.
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The graph can have at most 31 vertices.")

    cdef FastDigraph FD = FastDigraph(G)
    cdef int * g = FD.graph
    cdef int n = FD.n

    # minimums[i] is means to store the value of c'_{i+1}
    minimums = <uint8_t *> sage_malloc(sizeof(uint8_t)* n)
    cdef unsigned int i

    # They are initialized to n
    for 0<= i< n:
        minimums[i] = n

    cdef uint8_t tmp, tmp_count

    # We go through all sets
    for 1<= i< <unsigned int> (1<<n):
        tmp_count = <uint8_t> popcount32(i)
        tmp = <uint8_t> compute_out_neighborhood_cardinality(FD, i)

        # And update the costs
        minimums[tmp_count-1] = minimum(minimums[tmp_count-1], tmp)

    # We compute the maximum of all those values
    for 1<= i< n:
        minimums[0] = maximum(minimums[0], minimums[i])

    cdef int min = minimums[0]

    sage_free(minimums)

    return min

################################
# Exact exponential algorithms #
################################

def path_decomposition(G, algorithm = "exponential", verbose = False):
    r"""
    Returns the pathwidth of the given graph and the ordering of the vertices
    resulting in a corresponding path decomposition.

    INPUT:

    - ``G`` -- a digraph

    - ``algorithm`` -- (default: ``"exponential"``) Specify the algorithm to use
      among

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this exponential algorithm only
        works on graphs on less than 32 vertices. This can be changed to 54 if
        necessary, but 32 vertices already require 4GB of memory.

    EXAMPLE:

    The vertex separation of a cycle is equal to 2::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = graphs.CycleGraph(6)
        sage: pw, L = path_decomposition(g); pw
        2
        sage: pwm, Lm = path_decomposition(g, algorithm = "MILP"); pwm
        2

    TEST:

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = digraphs.Circuit(6)
        sage: path_decomposition(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The parameter must be a Graph.")

    from sage.graphs.digraph import DiGraph
    if algorithm == "exponential":
        return vertex_separation(DiGraph(G), verbose = verbose)
    else:
        return vertex_separation_MILP(DiGraph(G), verbosity = (1 if verbose else 0))


def vertex_separation(G, verbose = False):
    r"""
    Returns an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a digraph

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on less than 32 vertices. This can be changed to 54 if necessary,
        but 32 vertices already require 4GB of memory.

    EXAMPLE:

    The vertex separation of a circuit is equal to 1::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: g = digraphs.Circuit(6)
        sage: vertex_separation(g)
        (1, [0, 1, 2, 3, 4, 5])

    TEST:

    Given anything else than a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = graphs.CycleGraph(5)
        sage: lower_bound(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a DiGraph.

    Graphs with non-integer vertices::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: D=digraphs.DeBruijn(2,3)
        sage: vertex_separation(D)
        (2, ['000', '001', '100', '010', '101', '011', '110', '111'])
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The graph should have at most 31 vertices !")

    cdef FastDigraph g = FastDigraph(G)

    if verbose:
        print "Memory allocation"
        g.print_adjacency_matrix()

    sig_on()

    cdef unsigned int mem = 1 << g.n
    cdef int8_t * neighborhoods = <int8_t *> sage_malloc(mem)

    if neighborhoods == NULL:
        sig_off()
        raise MemoryError("Error allocating memory. I just tried to allocate "+str(mem>>10)+"MB, could that be too much ?")

    memset(neighborhoods, <int8_t> -1, mem)

    cdef int i,j , k
    for 0 <= k <g.n:
        if verbose:
            print "Looking for a strategy of cost", str(k)

        if exists(g, neighborhoods, 0, k) <= k:
            break

    if verbose:
        print "... Found !"
        print "Now computing the ordering"

    cdef list order = find_order(g, neighborhoods, k)

    # Relabelling the vertices
    cdef list vertices = G.vertices()
    for i, j in enumerate(order):
        order[i] = vertices[j]

    sage_free(neighborhoods)
    sig_off()

    return k, order

##############################################################################
# Actual algorithm, breadh-first search and updates of the costs of the sets #
##############################################################################

# Check whether an ordering with the given cost exists, and updates data in the
# neighborhoods array at the same time. See the module's documentation

cdef inline int exists(FastDigraph g, int8_t * neighborhoods, int current, int cost):

    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] < 0:
        neighborhoods[current] = compute_out_neighborhood_cardinality(g, current)

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > cost or (current == (1<<g.n)-1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors
    cdef int mini = g.n

    cdef int i
    cdef int next_set


    for 0<= i<g.n:
        if (current >> i)&1:
            continue

        # For each of the out-neighbors next_set of current

        next_set = current | 1<<i

        # Check whether there exists a cheap path toward {1..n}, and updated the
        # cost.
        mini = minimum(mini, exists(g, neighborhoods, next_set, cost))

        # We have found a path !
        if mini <= cost:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]

# Returns the ordering once we are sure it exists
cdef list find_order(FastDigraph g, int8_t * neighborhoods, int cost):
    cdef list ordering = []
    cdef int current = 0
    cdef int n = g.n
    cdef int i

    while n:
        # We look for n vertices

        for 0<= i<g.n:
            if (current >> i)&1:
                continue

            # Find the next set with small cost (we know it exists)
            next_set = current | 1<<i
            if neighborhoods[next_set] <= cost:
                ordering.append(i)
                current = next_set
                break

        # One less to find
        n -= 1

    return ordering

# Min/Max functions

cdef inline int minimum(int a, int b):
    if a<b:
        return a
    else:
        return b

cdef inline int maximum(int a, int b):
    if a>b:
        return a
    else:
        return b


#################################################################
# Function for testing the validity of a linear vertex ordering #
#################################################################

def is_valid_ordering(G, L):
    r"""
    Test if the linear vertex ordering `L` is valid for (di)graph `G`.

    A linear ordering `L` of the vertices of a (di)graph `G` is valid if all
    vertices of `G` are in `L`, and if `L` contains no other vertex and no
    duplicated vertices.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``L`` -- an ordered list of the vertices of ``G``.


    OUTPUT:

    Returns ``True`` if `L` is a valid vertex ordering for `G`, and ``False``
    oterwise.


    EXAMPLE:

    Path decomposition of a cycle::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.is_valid_ordering(G, L)
        True
        sage: vertex_separation.is_valid_ordering(G, [1,2])
        False

    TEST:

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.is_valid_ordering(2, [])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph or a DiGraph.

    Giving anything else than a list::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: vertex_separation.is_valid_ordering(G, {})
        Traceback (most recent call last):
        ...
        ValueError: The second parameter must be of type 'list'.
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("The input parameter must be a Graph or a DiGraph.")
    if not isinstance(L, list):
        raise ValueError("The second parameter must be of type 'list'.")

    return set(L) == set(G.vertices())


####################################################################
# Measurement functions of the widths of some graph decompositions #
####################################################################

def width_of_path_decomposition(G, L):
    r"""
    Returns the width of the path decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `pw_L(G)` of the path decomposition induced by the linear
    ordering `L` of the vertices of `G`. If `G` is a :mod:`DiGraph
    <sage.graphs.digraph>`, it returns instead the width `vs_L(G)` of the
    directed path decomposition induced by the linear ordering `L` of the
    vertices of `G`, where

    .. MATH::

        vs_L(G) & =  \max_{0\leq i< |V|-1} | N^+(L[:i])\setminus L[:i] |\\
        pw_L(G) & =  \max_{0\leq i< |V|-1} | N(L[:i])\setminus L[:i] |\\

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Path decomposition of a cycle::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.width_of_path_decomposition(G, L)
        2

    Directed path decomposition of a circuit::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.width_of_path_decomposition(G, L)
        1

    TESTS:

    Path decomposition of a BalancedTree::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.BalancedTree(3,2)
        sage: pw, L = vertex_separation.path_decomposition(G)
        sage: pw == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: L.reverse()
        sage: pw == vertex_separation.width_of_path_decomposition(G, L)
        False

    Directed path decomposition of a circuit::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(8)
        sage: vs, L = vertex_separation.vertex_separation(G)
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: L = [0,4,6,3,1,5,2,7]
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        False

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = Graph()
        sage: vertex_separation.width_of_path_decomposition(G, ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: The input linear vertex ordering L is not valid for G.
    """
    if not is_valid_ordering(G, L):
        raise ValueError("The input linear vertex ordering L is not valid for G.")

    vsL = 0
    S = set()
    neighbors_of_S_in_V_minus_S = set()

    for u in L:

        # We remove u from the neighbors of S
        neighbors_of_S_in_V_minus_S.discard(u)

        # We add vertex u to the set S
        S.add(u)

        if G._directed:
            Nu = G.neighbors_out(u)
        else:
            Nu = G.neighbors(u)

        # We add the (out-)neighbors of u to the neighbors of S
        for v in Nu:
            if (not v in S):
                neighbors_of_S_in_V_minus_S.add(v)

        # We update the cost of the vertex separation
        vsL = max( vsL, len(neighbors_of_S_in_V_minus_S) )

    return vsL


##########################################
# MILP formulation for vertex separation #
##########################################

def vertex_separation_MILP(G, integrality = False, solver = None, verbosity = 0):
    r"""
    Computes the vertex separation of `G` and the optimal ordering of its
    vertices using an MILP formulation.

    This function uses a mixed integer linear program (MILP) for determining an
    optimal layout for the vertex separation of `G`. This MILP is an improved
    version of the formulation proposed in [SP10]_. See the :mod:`module's
    documentation <sage.graphs.graph_decompositions.vertex_separation>` for more
    details on this MILP formulation.

    INPUTS:

    - ``G`` -- a DiGraph

    - ``integrality`` -- (default: ``False``) Specify if variables `x_v^t` and
      `u_v^t` must be integral or if they can be relaxed. This has no impact on
      the validity of the solution, but it is sometimes faster to solve the
      problem using binary variables only.

    - ``solver`` -- (default: ``None``) Specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve<sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class
      :class:`MixedIntegerLinearProgram<sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``). Sets the level of verbosity. Set
      to 0 by default, which means quiet.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLE:

    Vertex separation of a De Bruijn digraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.DeBruijn(2,3)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs
        2
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: vse, Le = vertex_separation.vertex_separation(G); vse
        2

    The vertex separation of a circuit is 1::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs
        1

    TESTS:

    Comparison with exponential algorithm::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):
        ...       G = digraphs.RandomDirectedGNP(10, 0.2)
        ...       ve, le = vertex_separation.vertex_separation(G)
        ...       vm, lm = vertex_separation.vertex_separation_MILP(G)
        ...       if ve != vm:
        ...          print "The solution is not optimal!"

    Comparison with Different values of the integrality parameter::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):  # long time (11s on sage.math, 2012)
        ...       G = digraphs.RandomDirectedGNP(10, 0.2)
        ...       va, la = vertex_separation.vertex_separation_MILP(G, integrality = False)
        ...       vb, lb = vertex_separation.vertex_separation_MILP(G, integrality = True)
        ...       if va != vb:
        ...          print "The integrality parameter change the result!"

    Giving anything else than a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.vertex_separation_MILP([])
        Traceback (most recent call last):
        ...
        ValueError: The first input parameter must be a DiGraph.
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("The first input parameter must be a DiGraph.")

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
    p = MixedIntegerLinearProgram( maximization = False, solver = solver )

    # Declaration of variables.
    x = p.new_variable( binary = integrality)
    u = p.new_variable( binary = integrality)
    y = p.new_variable( binary = True)
    z = p.new_variable( integer = True, dim = 1 )

    N = G.num_verts()
    V = G.vertices()

    # (2) x[v,t] <= x[v,t+1]   for all v in V, and for t:=0..N-2
    # (3) y[v,t] <= y[v,t+1]   for all v in V, and for t:=0..N-2
    for v in V:
        for t in xrange(N-1):
            p.add_constraint( x[v,t] - x[v,t+1] <= 0 )
            p.add_constraint( y[v,t] - y[v,t+1] <= 0 )

    # (4) y[v,t] <= x[w,t]  for all v in V, for all w in N^+(v), and for all t:=0..N-1
    for v in V:
        for w in G.neighbors_out(v):
            for t in xrange(N):
                p.add_constraint( y[v,t] - x[w,t] <= 0 )

    # (5) sum_{v in V} y[v,t] == t+1 for t:=0..N-1
    for t in xrange(N):
        p.add_constraint( p.sum([ y[v,t] for v in V ]) == t+1 )

    # (6) u[v,t] >= x[v,t]-y[v,t]    for all v in V, and for all t:=0..N-1
    for v in V:
        for t in xrange(N):
            p.add_constraint( x[v,t] - y[v,t] - u[v,t] <= 0 )

    # (7) z >= sum_{v in V} u[v,t]   for all t:=0..N-1
    for t in xrange(N):
        p.add_constraint( p.sum([ u[v,t] for v in V ]) - z['z'] <= 0 )

    # (8)(9) 0 <= x[v,t] and u[v,t] <= 1
    if not integrality:
        for v in V:
            for t in xrange(N):
                p.add_constraint( 0 <= x[v,t] <= 1 )
                p.add_constraint( 0 <= u[v,t] <= 1 )

    # (10) y[v,t] in {0,1}
    p.set_binary( y )

    # (11) 0 <= z <= |V|
    p.add_constraint( z['z'] <= N )

    #  (1) Minimize z
    p.set_objective( z['z'] )

    try:
        obj = p.solve( log=verbosity )
    except MIPSolverException:
        if integrality:
            raise ValueError("Unbounded or unexpected error")
        else:
            raise ValueError("Unbounded or unexpected error. Try with 'integrality = True'.")

    taby = p.get_values( y )
    tabz = p.get_values( z )
    # since exactly one vertex is processed per step, we can reconstruct the sequence
    seq = []
    for t in xrange(N):
        for v in V:
            if (taby[v,t] > 0) and (not v in seq):
                seq.append(v)
                break
    vs = int(round( tabz['z'] ))

    return vs, seq
