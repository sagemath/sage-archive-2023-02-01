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
from sage.combinat.posets.posets import Poset, Posets_all, FinitePosets_n
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.graph import DiGraph
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject

class PosetsGenerator(object):
    r"""
    A collection of examples of posets.

    EXAMPLES::
        sage: P = Posets()
        sage: P == loads(dumps(P))
        True
    """
    def BooleanLattice(self, n):
        """
        Returns the Boolean lattice containing `2^n` elements.

        EXAMPLES::

            sage: Posets.BooleanLattice(5)
            Finite lattice containing 32 elements
        """
        return LatticePoset([[Integer(x|(1<<y)) for y in range(0,n) if x&(1<<y)==0] for
            x in range(0,2**n)])

    def ChainPoset(self, n):
        """
        Returns a chain (a totally ordered poset) containing ``n`` elements.

        EXAMPLES::

            sage: C = Posets.ChainPoset(6); C
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

    def AntichainPoset(self, n):
        """
        Returns an antichain (a poset with no comparable elements)
        containing ``n`` elements.

        EXAMPLES::

            sage: A = Posets.AntichainPoset(6); A
            Finite poset containing 6 elements
            sage: for i in range(5):
            ...       for j in range(5):
            ...           if A.covers(A(i),A(j)):
            ...              print "TEST FAILED"
        """
        c = [[] for x in range(n)]
        return Poset(c)

    def PentagonPoset(self):
        """
        Return the "pentagon".

        EXAMPLES::

            sage: Posets.PentagonPoset()
            Finite lattice containing 5 elements
        """
        p = LatticePoset([[1,2],[4],[3],[4],[]])
        p.hasse_diagram()._pos = {0:[2,0],1:[0,2],2:[3,1],3:[3,3],4:[2,4]}
        return p

    def DiamondPoset(self, n):
        """
        Returns the lattice of rank two containing ``n`` elements.

        EXAMPLES::

            sage: Posets.DiamondPoset(7)
            Finite lattice containing 7 elements
        """
        c = [[n-1] for x in range(n)]
        c[0] = [x for x in range(1,n-1)]
        c[n-1] = []
        if n > 2:
            return LatticePoset(c)
        else:
            return Poset(c)

    def IntegerCompositions(self, n):
        """
        Returns the poset of integer compositions of the integer ``n``.

        A composition of a positive integer `n` is a list of positive
        integers that sum to `n`. The order is reverse refinement:
        ``[p_1,p_2,...,p_l]`` < ``[q_1,q_2,...,q_m]`` if ``q`` consists
        of an integer composition of ``p_1``, followed by an integer
        composition of ``p_2``, and so on.

        EXAMPLES::

            sage: P = Posets.IntegerCompositions(7); P
            Finite poset containing 64 elements
            sage: len(P.cover_relations())
            192
        """
        from sage.combinat.composition import Compositions
        C = Compositions(n)
        return Poset((C, [[c,d] for c in C for d in C if d.is_finer(c)]), cover_relations=False)

    def IntegerPartitions(self, n):
        """
        Returns the poset of integer partitions on the integer ``n``.

        A partition of a positive integer `n` is a non-increasing list
        of positive integers that sum to `n`. If ``p`` and ``q`` are
        integer partitions of `n`, then ``p`` covers ``q`` if and only
        if ``q`` is obtained from ``p`` by joining two parts of ``p``
        (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.IntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            28
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers
            of elements in the poset of integer paritions.
            """
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

    def RestrictedIntegerPartitions(self, n):
        """
        Returns the poset of integer partitions on the integer ``n``
        ordered by restricted refinement. That is, if ``p`` and ``q``
        are integer partitions of ``n``, then ``p`` covers ``q`` if and
        only if ``q`` is obtained from ``p`` by joining two distinct
        parts of ``p`` (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.RestrictedIntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            17
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers of elements in the
            restricted poset of integer paritions.
            """
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

    def RandomPoset(self, n,p):
        """
        Generate a random poset on ``n`` vertices according to a
        probability distribution ``p``.

        EXAMPLES::

            sage: Posets.RandomPoset(17,.15)
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

    def SymmetricGroupBruhatOrderPoset(self, n):
        """
        The poset of permutations with respect to Bruhat order.

        EXAMPLES::

            sage: Posets.SymmetricGroupBruhatOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10:
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        return Poset(dict([[s,s.bruhat_succ()]
                for s in Permutations(n)]),element_labels)

    def SymmetricGroupBruhatIntervalPoset(self, start, end):
        """
        The poset of permutations with respect to Bruhat order.

        INPUT:

        - ``start`` - list permutation

        - ``end`` - list permutation (same n, of course)

        .. note::

           Must have ``start`` <= ``end``.

        EXAMPLES:

        Any interval is rank symmetric if and only if it avoids these
        permutations::

            sage: P1 = Posets.SymmetricGroupBruhatIntervalPoset([0,1,2,3], [2,3,0,1])
            sage: P2 = Posets.SymmetricGroupBruhatIntervalPoset([0,1,2,3], [3,1,2,0])
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
            nodes[perm] = [succ_perm for succ_perm in perm.bruhat_succ()
                                if succ_perm.bruhat_lequal(end)]
            for succ_perm in nodes[perm]:
                if not nodes.has_key(succ_perm):
                    unseen.append(succ_perm)
        return Poset(nodes)

    def SymmetricGroupWeakOrderPoset(self, n,labels="permutations"):
        """
        The poset of permutations with respect to weak order.

        EXAMPLES::

            sage: Posets.SymmetricGroupWeakOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10 and labels == "permutations":
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        if n < 10 and labels == "reduced_words":
            element_labels = dict([[s,"".join(map(str,s.reduced_word_lexmin()))] for s in Permutations(n)])
        def weak_covers(s):
            r"""
            Nested function for computing the covers of elements in the
            poset of weak order for the symmetric group.
            """
            return [v for v in s.bruhat_succ() if
                s.length() + (s.inverse()*v).length() == v.length()]
        return Poset(dict([[s,weak_covers(s)] for s in Permutations(n)]),element_labels)

    def __call__(self, n=None):
        r"""
        Return either the CombinatorialClass of all posets or the
        CombinatorialClass of all finite posets on ``n`` elements.

        EXAMPLES::
            sage: Posets()
            Posets
            sage: Posets(4)
            Posets containing 4 vertices
        """
        if n is None:
            return Posets_all()
        return FinitePosets_n(n)

posets = PosetsGenerator()
Posets = posets


###########################################################################
##### DEPRECATION WARNINGS ################################################
##### Added 28 April 2009 #################################################
###########################################################################

from sage.misc.misc import deprecation

def BooleanLattice(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.BooleanLattice`` instead.

    TESTS::
        sage: BooleanLattice(3)
        doctest:1: DeprecationWarning: BooleanLattice is deprecated, use Posets.BooleanLattice instead!
        Finite lattice containing 8 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("BooleanLattice", "BooleanLattice"))
    return Posets.BooleanLattice(*args, **kwds)

def ChainPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.ChainPoset`` instead.

    TESTS::
        sage: ChainPoset(3)
        doctest:1: DeprecationWarning: ChainPoset is deprecated, use Posets.ChainPoset instead!
        Finite lattice containing 3 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("ChainPoset","ChainPoset"))
    return Posets.ChainPoset(*args, **kwds)

def AntichainPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.AntichainPoset`` instead.

    TESTS::
        sage: AntichainPoset(3)
        doctest:1: DeprecationWarning: AntichainPoset is deprecated, use Posets.AntichainPoset instead!
        Finite poset containing 3 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("AntichainPoset","AntichainPoset"))
    return Posets.AntichainPoset(*args, **kwds)

def PentagonPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.PentagonPoset`` instead.

    TESTS::
        sage: PentagonPoset()
        doctest:1: DeprecationWarning: PentagonPoset is deprecated, use Posets.PentagonPoset instead!
        Finite lattice containing 5 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("PentagonPoset","PentagonPoset"))
    return Posets.PentagonPoset(*args, **kwds)

def DiamondPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.DiamondPoset`` instead.

    TESTS::
        sage: DiamondPoset(3)
        doctest:1: DeprecationWarning: DiamondPoset is deprecated, use Posets.DiamondPoset instead!
        Finite lattice containing 3 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("DiamondPoset","DiamondPoset"))
    return Posets.DiamondPoset(*args, **kwds)

def PosetOfIntegerCompositions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.IntegerCompositions`` instead.

    TESTS::
        sage: PosetOfIntegerCompositions(3)
        doctest:1: DeprecationWarning: PosetOfIntegerCompositions is deprecated, use Posets.IntegerCompositions instead!
        Finite poset containing 4 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfIntegerCompositions","IntegerCompositions"))
    return Posets.IntegerCompositions(*args, **kwds)

def PosetOfIntegerPartitions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.IntegerPartitions`` instead.

    TESTS::
        sage: PosetOfIntegerPartitions(3)
        doctest:1: DeprecationWarning: PosetOfIntegerPartitions is deprecated, use Posets.IntegerPartitions instead!
        Finite poset containing 3 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfIntegerPartitions","IntegerPartitions"))
    return Posets.IntegerPartitions(*args, **kwds)

def PosetOfRestrictedIntegerPartitions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.RestrictedIntegerPartitions`` instead.

    TESTS::
        sage: PosetOfRestrictedIntegerPartitions(3)
        doctest:1: DeprecationWarning: PosetOfRestrictedIntegerPartitions is deprecated, use Posets.RestrictedIntegerPartitions instead!
        Finite poset containing 3 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfRestrictedIntegerPartitions","RestrictedIntegerPartitions"))
    return Posets.RestrictedIntegerPartitions(*args, **kwds)

def RandomPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.RandomPoset`` instead.

    TESTS::
        sage: RandomPoset(17,.15)
        doctest:1: DeprecationWarning: RandomPoset is deprecated, use Posets.RandomPoset instead!
        Finite poset containing 17 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("RandomPoset","RandomPoset"))
    return Posets.RandomPoset(*args, **kwds)


def SymmetricGroupBruhatOrderPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.SymmetricGroupBruhatOrderPoset`` instead.

    TESTS::
        sage: SymmetricGroupBruhatOrderPoset(3)
        doctest:1: DeprecationWarning: SymmetricGroupBruhatOrderPoset is deprecated, use Posets.SymmetricGroupBruhatOrderPoset instead!
        Finite poset containing 6 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("SymmetricGroupBruhatOrderPoset","SymmetricGroupBruhatOrderPoset"))
    return Posets.SymmetricGroupBruhatOrderPoset(*args, **kwds)

def SymmetricGroupWeakOrderPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.SymmetricGroupWeakOrderPoset`` instead.

    TESTS::
        sage: SymmetricGroupWeakOrderPoset(3)
        doctest:1: DeprecationWarning: SymmetricGroupWeakOrderPoset is deprecated, use Posets.SymmetricGroupWeakOrderPoset instead!
        Finite poset containing 6 elements
    """
    deprecation("%s is deprecated, use Posets.%s instead!" % \
           ("SymmetricGroupWeakOrderPoset","SymmetricGroupWeakOrderPoset"))
    return Posets.SymmetricGroupWeakOrderPoset(*args, **kwds)
