r"""
Some examples of posets and lattices.
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from random import random
from sage.combinat.permutation import Permutations, Permutation
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.graph import DiGraph

def BooleanLattice(n):
    """
    Returns the Boolean lattice containing 2n elements.

    EXAMPLES::

        sage: BooleanLattice(5)
        Finite lattice containing 32 elements
    """
    return LatticePoset([[x|(1<<y) for y in range(0,n) if x&(1<<y)==0] for
        x in range(0,2**n)])

def ChainPoset(n):
    """
    Returns a chain (a totally ordered poset) containing n elements.

    EXAMPLES::

        sage: C = ChainPoset(6); C
        Finite lattice containing 6 elements
        sage: C.linear_extension()
        [0, 1, 2, 3, 4, 5]
        sage: for i in range(5):
        ...       for j in range(5):
        ...           if C.covers(C(i),C(j)) and j != i+1:
        ...              print "TEST FAILED"
    """
    c = [[x+1] for x in range(n)]
    c[n-1] = []
    return LatticePoset(c)

def AntichainPoset(n):
    """
    Returns an antichain (a poset with no comparable elements)
    containing n elements.

    EXAMPLES::

        sage: A = AntichainPoset(6); A
        Finite poset containing 6 elements
        sage: for i in range(5):
        ...       for j in range(5):
        ...           if A.covers(A(i),A(j)):
        ...              print "TEST FAILED"
    """
    c = [[] for x in range(n)]
    return Poset(c)

def PentagonPoset():
    """
    Return the "pentagon".

    EXAMPLES::

        sage: PentagonPoset()
        Finite lattice containing 5 elements
    """
    p = LatticePoset([[1,2],[4],[3],[4],[]])
    p.hasse_diagram()._pos = {0:[2,0],1:[0,2],2:[3,1],3:[3,3],4:[2,4]}
    return p

def DiamondPoset(n):
    """
    Returns the lattice of rank two containing n elements.

    EXAMPLES::

        sage: DiamondPoset(7)
        Finite lattice containing 7 elements
    """
    c = [[n-1] for x in range(n)]
    c[0] = [x for x in range(1,n-1)]
    c[n-1] = []
    if n > 2:
        return LatticePoset(c)
    else:
        return Poset(c)

def PosetOfIntegerCompositions(n):
    """
    Returns the poset of integer compositions on the integer n.

    A composition of a positive integer n is a list of positive
    integers that some to n. The order is reverse refinement:
    [p_1,p_2,...,p_l] [q_1,q_2,...,q_m] if q consists of a
    integer composition of p_1, followed by an integer composition of
    p_2, and so on.

    EXAMPLES::

        sage: P = PosetOfIntegerCompositions(7); P
        Finite poset containing 64 elements
        sage: len(P.cover_relations())
        192
    """
    from sage.combinat.composition import Compositions
    C = Compositions(n)
    return Poset((C, [[c,d] for c in C for d in C if d.is_finer(c)]), cover_relations=False)

def PosetOfIntegerPartitions(n):
    """
    Returns the poset of integer partitions on the integer n.

    A partition of a positive integer n is a non-increasing list of
    positive integers that some to n. If p and q are integer partitions
    of n, then p covers q if and only if q is obtained from p by
    joining two parts of p (and sorting, if necessary).

    EXAMPLES::

        sage: P = PosetOfIntegerPartitions(7); P
        Finite poset containing 15 elements
        sage: len(P.cover_relations())
        28
    """
    def lower_covers(partition):
        lc = []
        for i in range(0,len(partition)-1):
            for j in range(i+1,len(partition)):
                new_partition = partition[:]
                del new_partition[j]
                del new_partition[i]
                new_partition.append(partition[i]+partition[j])
                new_partition.sort(reverse=True)
                tup = tuple(new_partition)
                if tup not in lc:
                    lc.append(tup)
        return lc
    from sage.combinat.partition import partitions_list
    H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in
        partitions_list(n)]))
    return Poset(H.reverse())

def PosetOfRestrictedIntegerPartitions(n):
    """
    Returns the poset of integer partitions on the integer n ordered by
    restricted refinement. That is, if p and q are integer partitions
    of n, then p covers q if and only if q is obtained from p by
    joining two distinct parts of p (and sorting, if necessary).

    EXAMPLES::

        sage: P = PosetOfRestrictedIntegerPartitions(7); P
        Finite poset containing 15 elements
        sage: len(P.cover_relations())
        17
    """
    def lower_covers(partition):
        lc = []
        for i in range(0,len(partition)-1):
            for j in range(i+1,len(partition)):
                if partition[i] != partition[j]:
                    new_partition = partition[:]
                    del new_partition[j]
                    del new_partition[i]
                    new_partition.append(partition[i]+partition[j])
                    new_partition.sort(reverse=True)
                    tup = tuple(new_partition)
                    if tup not in lc:
                        lc.append(tup)
        return lc
    from sage.combinat.partition import Partitions
    H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in
        Partitions(n)]))
    return Poset(H.reverse())

def RandomPoset(n,p):
    """
    Generate a random poset on n vertices according to a probability
    distribution p.

    EXAMPLES::

        sage: RandomPoset(17,.15)
        Finite poset containing 17 elements
    """
    p = float(p)
    D = DiGraph(loops=False,multiedges=False)
    D.add_vertices(range(n))
    for i in range(n):
        for j in range(n):
            if random() < p:
                D.add_edge(i,j)
                if not D.is_directed_acyclic():
                    D.delete_edge(i,j)
    return Poset(D,cover_relations=False)

def SymmetricGroupBruhatOrderPoset(n):
    """
    The poset of permutations with respect to Bruhat order.

    EXAMPLES::

        sage: SymmetricGroupBruhatOrderPoset(4)
        Finite poset containing 24 elements
    """
    if n < 10:
        element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
    return Poset(dict([[s,s.bruhat_succ()] for s in Permutations(n)]),element_labels)

def SymmetricGroupBruhatIntervalPoset(start, end):
    """
    The poset of permutations with respect to Bruhat order.

    INPUT:


    -  ``start`` - list permutation

    -  ``end`` - list permutation (same n, of course)


    .. note::

       Must have start = end.

    EXAMPLES::

        sage: from sage.combinat.posets.poset_examples import SymmetricGroupBruhatIntervalPoset

    Any interval is rank symmetric if and only if it avoids these
    permutations::

        sage: P1 = SymmetricGroupBruhatIntervalPoset([0,1,2,3], [2,3,0,1])
        sage: P2 = SymmetricGroupBruhatIntervalPoset([0,1,2,3], [3,1,2,0])
        sage: ranks1 = [P1.rank(v) for v in P1]
        sage: ranks2 = [P2.rank(v) for v in P2]
        sage: [ranks1.count(i) for i in uniq(ranks1)]
        [1, 3, 5, 4, 1]
        sage: [ranks2.count(i) for i in uniq(ranks2)]
        [1, 3, 5, 6, 4, 1]
    """
    start = Permutation(start)
    end = Permutation(end)
    if len(start) != len(end):
        raise TypeError("Start (%s) and end (%s) must have same length."%(start, end))
    if not start.bruhat_lequal(end):
        raise TypeError("Must have start (%s) <= end (%s) in Bruhat order."%(start, end))
    unseen = [start]
    nodes = {}
    while len(unseen) > 0:
        perm = unseen.pop(0)
        nodes[perm] = [succ_perm for succ_perm in perm.bruhat_succ_iterator() if succ_perm.bruhat_lequal(end)]
        for succ_perm in nodes[perm]:
            if not nodes.has_key(succ_perm):
                unseen.append(succ_perm)
    return Poset(nodes)

def SymmetricGroupWeakOrderPoset(n,labels="permutations"):
    """
    The poset of permutations with respect to weak order.

    EXAMPLES::

        sage: SymmetricGroupWeakOrderPoset(4)
        Finite poset containing 24 elements
    """
    if n < 10 and labels == "permutations":
        element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
    if n < 10 and labels == "reduced_words":
        element_labels = dict([[s,"".join(map(str,s.reduced_word_lexmin()))] for s in Permutations(n)])
    def weak_covers(s):
        return [v for v in s.bruhat_succ() if
        	s.length() + (s.inverse()*v).length() == v.length()]
    return Poset(dict([[s,weak_covers(s)] for s in Permutations(n)]),element_labels)

#def PosetExamples(n,uc=None):
#    default_uc = [[[2,3], [], [1], [1], [1], [3,4]],
#        [[1,2,3,4,5],[6],[6],[6],[6],[6],[]]]
#    if uc is None:
#        uc = default_uc[n]
#    Q = Poset(uc)
#    elms = list("abcdefghijklmnopqrstuvwxyz"[0:len(uc)])
#    dag = DiGraph(dict([[i,uc[i]] for i in range(len(uc))]))
#    dag.relabel(elms)
#    P = Poset(dag)
#    return P,Q
