r"""
Crystals of alcove paths
"""

#*****************************************************************************
#       Copyright (C) 2008 Brant Jones <brant at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.crystals.letters import Letter
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.integer import Integer
from sage.misc.misc import math,ellipsis_range
from sage.combinat.partition import Partitions
from sage.combinat.crystals.all import CrystalOfTableaux
from sage.combinat.root_system.root_system import RootSystem
from sage.graphs.graph import DiGraph
import copy

#####################################################################
#  Crystal (parent) class.
#####################################################################

class ClassicalCrystalOfAlcovePaths(UniqueRepresentation, Parent):
    r"""
    Implementation of crystal of alcove paths of the given classical type with
    given highest weight, based on the Lenart--Postnikov model [LP2008].

    These are highest weight crystals for classical types `A_n`, `B_n`, `C_n`,
    `D_n` and the exceptional types `F_4`, `G_2`, `E_6`, `E_7`, `E_8`.

    INPUT:

        - ``cartan_type`` is the Cartan type of a classical Dynkin diagram
        - ``highest_weight`` is a dominant weight as a list of coefficients of
          the fundamental weights `Lambda_i`

    In this model, a chain of roots is associated to the given highest_weight,
    and then the elements of the crystal are indexed by "admissible subsets"
    which indicate "folding positions" along the chain of roots.  See [LP2008]
    for details.

    TODO:

    - Resolve speed issues;  `E_6(\Lambda_1)` takes just under 4 minutes to list().
      To construct the highest-weight node takes 15 sec for `E_6(\Lambda_4)`.
      The initial chain has 42 roots.

    TESTS:

    The following example appears in Figure 2 of [LP2008]::

        sage: C = ClassicalCrystalOfAlcovePaths(['G',2],[0,1]);
        sage: G = C.digraph()
        sage: GG= DiGraph({
        ...       ()        : {(0)         : 2 },
        ...       (0)       : {(0,8)       : 1 },
        ...       (0,1)     : {(0,1,7)     : 2 },
        ...       (0,1,2)   : {(0,1,2,9)   : 1 },
        ...       (0,1,2,3) : {(0,1,2,3,4) : 2 },
        ...       (0,1,2,6) : {(0,1,2,3)   : 1 },
        ...       (0,1,2,9) : {(0,1,2,6)   : 1 },
        ...       (0,1,7)   : {(0,1,2)     : 2 },
        ...       (0,1,7,9) : {(0,1,2,9)   : 2 },
        ...       (0,5)     : {(0,1)       : 1, (0,5,7) : 2 },
        ...       (0,5,7)   : {(0,5,7,9)   : 1 },
        ...       (0,5,7,9) : {(0,1,7,9)   : 1 },
        ...       (0,8)     : {(0,5)       : 1 },
        ...       })
        sage: G.is_isomorphic(GG)
        True

        sage: for edge in sorted([(u.value, v.value, i) for (u,v,i) in G.edges()]):
        ...       print edge
        ([], [0], 2)
        ([0], [0, 8], 1)
        ([0, 1], [0, 1, 7], 2)
        ([0, 1, 2], [0, 1, 2, 9], 1)
        ([0, 1, 2, 3], [0, 1, 2, 3, 4], 2)
        ([0, 1, 2, 6], [0, 1, 2, 3], 1)
        ([0, 1, 2, 9], [0, 1, 2, 6], 1)
        ([0, 1, 7], [0, 1, 2], 2)
        ([0, 1, 7, 9], [0, 1, 2, 9], 2)
        ([0, 5], [0, 1], 1)
        ([0, 5], [0, 5, 7], 2)
        ([0, 5, 7], [0, 5, 7, 9], 1)
        ([0, 5, 7, 9], [0, 1, 7, 9], 1)
        ([0, 8], [0, 5], 1)

    REFERENCES:

        .. [LP2008]  C. Lenart and A. Postnikov. A combinatorial model for crystals of Kac-Moody algebras. Trans. Amer. Math. Soc.  360  (2008), 4349-4381.
    """

    @staticmethod
    def __classcall__(cls, cartan_type, highest_weight):
        """
        cartan_type and heighest_weight are lists, which are mutable, this
        causes a problem for class UniqueRepresentation, the following code
        fixes this problem.

        EXAMPLES::
            sage: ClassicalCrystalOfAlcovePaths.__classcall__(ClassicalCrystalOfAlcovePaths,['A',3],[0,1,0])
            <class 'sage.combinat.crystals.alcove_path.ClassicalCrystalOfAlcovePaths_with_category'>
        """
        cartan_type = CartanType(cartan_type)
        highest_weight = tuple(highest_weight)
        return super(ClassicalCrystalOfAlcovePaths, cls).__classcall__(cls, cartan_type, highest_weight)

    def __init__(self, cartan_type, highest_weight):
        """
        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[1,0,0])
            sage: C.list()
            [[], [0], [0, 1], [0, 1, 2]]
            sage: TestSuite(C).run()
        """
        Parent.__init__(self, category = ClassicalCrystals())
        self._cartan_type = CartanType(cartan_type)
        self._name = "The crystal of alcove paths for type %s"%cartan_type
        self.chain_cache = {}
        self.endweight_cache = {}

        self.R = RootSystem(cartan_type)
        alpha = self.R.root_space().simple_roots()
        Lambda = self.R.weight_space().basis()

        self.positive_roots = sorted(self.R.root_space().positive_roots());

        self.weight = Lambda[Integer(1)] - Lambda[Integer(1)]
        offset = self.R.index_set()[Integer(0)]
        for j in self.R.index_set():
            self.weight = self.weight + highest_weight[j-offset]*Lambda[j]

        self.initial_element = self([])

        self.initial_element.chain = self.get_initial_chain(self.weight)
        rho = (Integer(1)/Integer(2))*sum(self.positive_roots)
        self.initial_element.endweight = rho

        self.chain_cache[ str([]) ] = self.initial_element.chain
        self.endweight_cache[ str([]) ] = self.initial_element.endweight

        self.module_generators = [self.initial_element]

        self._list = super(ClassicalCrystalOfAlcovePaths, self).list()
        self._digraph = super(ClassicalCrystalOfAlcovePaths, self).digraph()
        self._digraph_closure = self.digraph().transitive_closure()

    def get_initial_chain(self, highest_weight):
        """
        Called internally by __init__() to construct the chain of roots
        associated to the highest weight element.

        EXAMPLES::
            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C.get_initial_chain(RootSystem(['A',3]).weight_space().basis()[1])
            [[alpha[1], alpha[1]], [alpha[1] + alpha[2], alpha[1] + alpha[2]], [alpha[1] + alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]]
        """
        pos_roots = self.positive_roots
        tis = self.R.index_set()
        tis.reverse()
        cv_to_pos_root = {}
        to_sort = []
        for r in pos_roots:
            j = highest_weight.scalar( r.associated_coroot() )
            if (int(math.floor(j)) - j == Integer(0)):
                j = j - Integer(1)
            j = int(math.floor(j))
            for k in (ellipsis_range(Integer(0),Ellipsis,j)):
                cv = []
                cv.append((Integer(1)/highest_weight.scalar(r.associated_coroot()))*k)
                for i in tis:
                    cv.append((Integer(1)/highest_weight.scalar(r.associated_coroot()))*r.associated_coroot().coefficient(i))
                cv_to_pos_root[ str(cv) ] = r
                to_sort.append( cv )

        to_sort.sort() # Note:  Python sorts nested lists lexicographically by default.
        lambda_chain = []
        for t in to_sort:
            lambda_chain.append( [ cv_to_pos_root[str(t)], cv_to_pos_root[str(t)] ] )
        return lambda_chain

    def _element_constructor_(self, value):
        """
        Coerces value into self.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[1,0,0])
            sage: C.module_generators
            [[]]
            sage: C([]).e(1)
            sage: C([]).f(1)
            [0]
        """
        return self.element_class(self, value)

    def list(self):
        """
        Returns a list of the elements of self.

        .. warning::

           This can be slow!

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C.list()
            [[], [0], [0, 1], [0, 2], [0, 1, 2], [0, 1, 2, 3]]
        """
        return self._list

    def digraph(self):
        """
        Returns the directed graph associated to self.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C.digraph().vertices()
            [[], [0], [0, 1], [0, 2], [0, 1, 2], [0, 1, 2, 3]]
            sage: C.digraph().edges()
            [([], [0], 2), ([0], [0, 1], 1), ([0], [0, 2], 3), ([0, 1], [0, 1, 2], 3), ([0, 2], [0, 1, 2], 1), ([0, 1, 2], [0, 1, 2, 3], 2)]
        """
        return self._digraph

    def lt_elements(self, x, y):
        r"""
        Returns True if and only if there is a path
        from x to y in the crystal graph.

        Because the crystal graph is classical, it is a directed
        acyclic graph which can be interpreted as a poset. This
        function implements the comparison function of this poset.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: x = C([])
            sage: y = C([0,2])
            sage: C.lt_elements(x,y)
            True
            sage: C.lt_elements(y,x)
            False
            sage: C.lt_elements(x,x)
            False
        """
        assert x.parent() == self and y.parent() == self
        if self._digraph_closure.has_edge(x,y):
            return True
        return False

#####################################################################
#  Element class.
#####################################################################

class ClassicalCrystalOfAlcovePathsElement(Letter):
    r"""
    Crystal of alcove paths element

    These elements are indexed by "admissible subsets" which indicate "folding positions"
    along the chain of roots.  Not every subset corresponds to an element.

    TESTS::

        sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
        sage: C.list()
        [[], [0], [0, 1], [0, 2], [0, 1, 2], [0, 1, 2, 3]]
        sage: [ [ x < y for y in C ] for x in C ]
        [[False, True, True, True, True, True], [False, False, True, True, True, True], [False, False, False, False, True, True], [False, False, False, False, True, True], [False, False, False, False, False, True], [False, False, False, False, False, False]]
        sage: C([]).parent()
        <class 'sage.combinat.crystals.alcove_path.ClassicalCrystalOfAlcovePaths_with_category'>
        sage: C([])
        []
        sage: C([]) == C([0])
        False
        sage: C([]) == C([])
        True
        sage: C([]) < C([0])
        True
        sage: C([0,1]) < C([0,2])
        False
        sage: C([0,1]) < C([0,1,2,3])
        True
        sage: TestSuite(C).run()
    """

    def __init__(self, parent, value):
        """
        INPUT:

        - ``value`` is an admissible subset, stored as a list.

        chain is a chain of roots.  This is computed from value, using
        get_chain_from_subset().  We do this lazily because it can be
        expensive.

        endweight is the weight that lies at the end of the chain of roots.

        chain and endweight are cached (in the parent), and they are only
        computed when required by E(), F() or weight().

        The weight of the element can be theoretically be computed from
        endweight, although this is not yet implemented.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C.list()
            [[], [0], [0, 1], [0, 2], [0, 1, 2], [0, 1, 2, 3]]
            sage: c = C([0,1])
            sage: TestSuite(c).run()
        """
        Letter.__init__(self, parent, value)
        if (str(value) in parent.chain_cache.keys()):
            self.chain = parent.chain_cache[ str(value) ]
            self.endweight = parent.endweight_cache[ str(value) ]
        else:
            self.chain = None
            self.endweight = None

    def get_chain_from_subset(self, verbose = False):
        r"""
        Returns a list of pairs of roots, from an admissible subset of
        an initial chain of pairs of roots.  Called internally by
        __init__().

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C([]).get_chain_from_subset()
            ([[alpha[2], alpha[2]], [alpha[1] + alpha[2], alpha[1] + alpha[2]], [alpha[2] + alpha[3], alpha[2] + alpha[3]], [alpha[1] + alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]], 3/2*alpha[1] + 2*alpha[2] + 3/2*alpha[3])
        """
        rchain = []
        # NT: Is a deep copy required?
        rchain = copy.deepcopy(self.parent().initial_element.chain)
        rend = copy.deepcopy(self.parent().initial_element.endweight)
        J = sorted(copy.deepcopy(self.value))
        J.reverse()

        if verbose == True:
            print "**** get_chain_from_subset():  initial weight at infinity ", rend
        for i in J:
            (rchain, rend) = self.fold(i, rchain, rend)
            if verbose == True:
                print "**** get_chain_from_subset() after fold ", i, ":  weight at infinity ", rend
        return (rchain, rend)

    def fold(self, i, chain, endweight):
        r"""
        Returns a new chain that results from folding at position i.
        Called internally by __init__().

        TESTS::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C([0]).get_chain_from_subset()
            ([[alpha[2], -alpha[2]], [alpha[1], alpha[1]], [alpha[3], alpha[3]], [alpha[1] + alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]], 3/2*alpha[1] + alpha[2] + 3/2*alpha[3])
        """
        s = self.parent().R.root_space().reflection(chain[i][Integer(0)])
        chain[i][Integer(1)] = s(chain[i][Integer(1)])
        for j in (ellipsis_range((i+Integer(1)),Ellipsis,len(chain)-Integer(1))):
            pair = chain[j]
            pair[Integer(0)] = s(pair[Integer(0)])
            pair[Integer(1)] = s(pair[Integer(1)])
        return (chain, s(endweight))

    def __hash__(self):
        """
        Return hash value.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0])
            sage: C([0]).__hash__() == C([]).f(2).__hash__()
            True
        """
        return hash(tuple(self.value))

    def e(self, i, verbose=False):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['E',6],[1,0,0,0,0,0])
            sage: C.module_generators
            [[]]
            sage: C([]).f(1)
            [0]
            sage: C([]).f(1).e(1)
            []
        """
        assert i in self.index_set()

        # Construct the chain associated to the subset self.value
        if self.chain is None:
            if (str(self.value) in self.parent().chain_cache.keys()):
                self.chain = self.parent().chain_cache[ str(self.value) ]
                self.endweight = self.parent().endweight_cache[ str(self.value) ]
            else:
                (self.chain, self.endweight) = self.get_chain_from_subset()
                self.parent().chain_cache[str(self.value)] = copy.deepcopy(self.chain)
                self.parent().endweight_cache[str(self.value)] = copy.deepcopy(self.endweight)

        ai = self.parent().R.root_space().simple_root(i)

        if verbose == True:
            print "( self.value, i = ", self.value, i, " ) "
            print "( chain = ", self.chain, " ) "
        IJp = filter( lambda j: self.chain[j][Integer(0)] == ai or self.chain[j][Integer(0)] == (-Integer(1))*ai, (ellipsis_range(Integer(0) ,Ellipsis, len(self.chain)-Integer(1))) )
        if verbose == True:
            print "( IJp = ", IJp, ") "
        SJp = []
        for j in IJp:
            pair = self.chain[j]
            if ( pair[Integer(0)].coefficients()[Integer(0)] < Integer(0) ):
                SJp.append( -Integer(1) )
            else:
                SJp.append( Integer(1) )
            if ( pair[Integer(1)].coefficients()[Integer(0)] < Integer(0) ):
                SJp.append( -Integer(1) )
            else:
                SJp.append( Integer(1) )

        if (self.endweight.scalar( ai.associated_coroot() ) < Integer(0)):
            SJp.append( -Integer(1) )
        else:
            SJp.append( Integer(1) )

        if verbose == True:
            print "( SJp = ", SJp, ") "

        LJp = []
        t = Integer(0)-(Integer(1)/Integer(2))
        for j in (ellipsis_range(Integer(0),Ellipsis,len(SJp)-Integer(1))):
            if SJp[j] == Integer(1):
                t = t+(Integer(1)/Integer(2))
            elif SJp[j] == -Integer(1):
                t = t-(Integer(1)/Integer(2))
            if ((j % Integer(2)) == Integer(0)):
                LJp.append(t)

        if verbose == True:
            print "( LJp = ", LJp, " ) "

        MJp = max(LJp)
        if verbose == True:
            print "( MJp = ", MJp, " ) "

        if (MJp <= LJp[ len(LJp)-1 ]):
            return None

        rLJp = copy.deepcopy(LJp)
        rLJp.reverse()
        k = rLJp.index(MJp)
        k = (len(LJp)-1) - k
        m = k + 1

        rsubseq = []
        rsubseq.extend(self.value)
        rsubseq.remove(IJp[k])
        if (m < len(IJp)):
            rsubseq.append(IJp[m])
        rsubseq.sort()

        return self.parent()(rsubseq)

    def f(self, i, verbose=False):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = ClassicalCrystalOfAlcovePaths(['E',6],[1,0,0,0,0,0])
            sage: C.module_generators
            [[]]
            sage: C([]).f(1)
            [0]
            sage: C([]).f(1).f(3)
            [0, 1]
            sage: C([]).f(1).f(3).f(4)
            [0, 1, 2]
            sage: C([]).f(1).f(3).f(4).f(5)
            [0, 1, 2, 4]
            sage: C([]).f(1).f(3).f(4).f(2)
            [0, 1, 2, 3]
        """
        assert i in self.index_set()

        # Construct the chain associated to the subset self.value
        if self.chain is None:
            if (str(self.value) in self.parent().chain_cache.keys()):
                self.chain = self.parent().chain_cache[ str(self.value) ]
                self.endweight = self.parent().endweight_cache[ str(self.value) ]
            else:
                (self.chain, self.endweight) = self.get_chain_from_subset()
                self.parent().chain_cache[str(self.value)] = copy.deepcopy(self.chain)
                self.parent().endweight_cache[str(self.value)] = copy.deepcopy(self.endweight)

        ai = self.parent().R.root_space().simple_root(i)

        if verbose == True:
            print "( self.value, i = ", self.value, i, " ) "
            print "( chain = ", self.chain, " ) "
        IJp = filter( lambda j: self.chain[j][Integer(0)] == ai or self.chain[j][Integer(0)] == (-Integer(1))*ai, (ellipsis_range(Integer(0) ,Ellipsis, len(self.chain)-Integer(1))) )
        if verbose == True:
            print "( IJp = ", IJp, ") "
        SJp = []
        for j in IJp:
            pair = self.chain[j]
            if ( pair[Integer(0)].coefficients()[Integer(0)] < Integer(0) ):
                SJp.append( -Integer(1) )
            else:
                SJp.append( Integer(1) )
            if ( pair[Integer(1)].coefficients()[Integer(0)] < Integer(0) ):
                SJp.append( -Integer(1) )
            else:
                SJp.append( Integer(1) )

        # Incidentally, the vector at infinity in LJp is <mu(J), a_p^check>.
        # This gives an alternative way to calculate the last SJp/LJp entries.
        if (self.endweight.scalar( ai.associated_coroot() ) < Integer(0)):
            SJp.append( -Integer(1) )
        else:
            SJp.append( Integer(1) )

        if verbose == True:
            print "( SJp = ", SJp, ") "

        LJp = []
        t = Integer(0)-(Integer(1)/Integer(2))
        for j in (ellipsis_range(Integer(0),Ellipsis,len(SJp)-Integer(1))):
            if SJp[j] == Integer(1):
                t = t+(Integer(1)/Integer(2))
            elif SJp[j] == -Integer(1):
                t = t-(Integer(1)/Integer(2))
            if ((j % Integer(2)) == Integer(0)):
                LJp.append(t)

        if verbose == True:
            print "( LJp = ", LJp, " ) "

        MJp = max(LJp)
        if verbose == True:
            print "( MJp = ", MJp, " ) "

        if (MJp <= Integer(0)):
            return None

        m = LJp.index(MJp)
        k = m - Integer(1)

        rsubseq = []
        rsubseq.extend(self.value)
        if (m < len(IJp)):
            rsubseq.remove(IJp[m])
        rsubseq.append(IJp[k])
        rsubseq.sort()

        return self.parent()(rsubseq)

ClassicalCrystalOfAlcovePaths.Element = ClassicalCrystalOfAlcovePathsElement

#####################################################################
#  Code to test by comparing with existing crystal implementations.
#####################################################################

def compare_graphs(g1, g2, node1, node2):
    r"""
    Returns two edge-labeled DiGraphs obtained from Crystal.digraph(), starting
    from the root nodes of each graph.

    EXAMPLES:
        sage: G1 = sage.combinat.crystals.all.CrystalOfTableaux(['A',3], shape = [1,1]).digraph()
        sage: G2 = ClassicalCrystalOfAlcovePaths(['A',3],[0,1,0]).digraph()
        sage: sage.combinat.crystals.alcove_path.compare_graphs(G1, G2, G1.vertices()[0], G2.vertices()[0])
        True
    """
    for out_edge in g1.outgoing_edges( node1 ):
        matched = False
        for o2 in g2.outgoing_edges( node2 ):
            if o2[2] == out_edge[2]:
                if matched == True:
                    print "ERROR:  Two edges with the same label for ", out_edge, " exist."
                    return False
                matched = True
                result = compare_graphs(g1, g2, out_edge[1], o2[1])
                if result == False:
                    return False
        if matched == False:
            print "ERROR:  No matching edge for ", out_edge, "."
            return False
    return True

def test_against_tableaux(R, N, k):
    r"""
    Tests our ClassicalCrystalOfAlcovePaths against all of the tableaux crystals
    of type R in rank N with highest weight given by a partition of k.

    EXAMPLES::

        sage: sage.combinat.crystals.alcove_path.test_against_tableaux(['A',3], 3, 2)
        ** Shape  [2]
          T has  10  nodes.
          C weight  [2, 0, 0]
          C has  10  nodes.
          Compare graphs:  True
        ** Shape  [1, 1]
          T has  6  nodes.
          C weight  [0, 1, 0]
          C has  6  nodes.
          Compare graphs:  True
    """
    shapes = Partitions(k).list()
    for shape in shapes:
        print "** Shape ", shape
        T = CrystalOfTableaux(R, shape = shape)
        ct = len(T.list())
        print "  T has ", ct, " nodes."
        #T.digraph().show(edge_labels=True)
        H = T.digraph()
        weight = T.module_generators[0].weight()
        w = [ weight.scalar(RootSystem(R).ambient_space().simple_coroot(i)) for i in range(1,N+1) ]
        print "  C weight ", w

        C = ClassicalCrystalOfAlcovePaths(R , w)

        cc = len(C.list())
        #C.digraph().show(edge_labels=True)
        G = C.digraph()
        print "  C has ", cc, " nodes."
        if cc != ct:
            print "FAIL: number of nodes differ.", cc, ct
            return
        print "  Compare graphs: ", compare_graphs(G, H, G.vertices()[0], H.vertices()[0])



def test_some_specific_examples():
    r"""
    Tests our ClassicalCrystalOfAlcovePaths against some specific examples.

    EXAMPLES::

        sage: sage.combinat.crystals.alcove_path.test_some_specific_examples()
        G2 example passed.
        C3 example passed.
        B3 example 1 passed.
        B3 example 2 passed.
        True
    """

    # This appears in Lenart.
    C = ClassicalCrystalOfAlcovePaths(['G',2],[0,1]);
    G = C.digraph();

    GT = DiGraph({
        ()        : {(0)         : 2 },
        (0)       : {(0,8)       : 1 },
        (0,1)     : {(0,1,7)     : 2 },
        (0,1,2)   : {(0,1,2,9)   : 1 },
        (0,1,2,3) : {(0,1,2,3,4) : 2 },
        (0,1,2,6) : {(0,1,2,3)   : 1 },
        (0,1,2,9) : {(0,1,2,6)   : 1 },
        (0,1,7)   : {(0,1,2)     : 2 },
        (0,1,7,9) : {(0,1,2,9)   : 2 },
        (0,5)     : {(0,1)       : 1, (0,5,7) : 2 },
        (0,5,7)   : {(0,5,7,9)   : 1 },
        (0,5,7,9) : {(0,1,7,9)   : 1 },
        (0,8)     : {(0,5)       : 1 }
        })


    if (G.is_isomorphic(GT) != True):
        return False;
    else:
        print "G2 example passed."

    # Some examples from Hong--Kang:

    # type C, ex. 8.3.5, pg. 189
    C = ClassicalCrystalOfAlcovePaths(['C',3],[0,0,1]);
    G = C.digraph();
    GT = DiGraph({
        ():{ (0): 3},
        (0):{ (0, 6): 2},
        (0, 1):{ (0, 1, 3): 3, (0, 1, 7): 1},
        (0, 1, 2):{ (0, 1, 2, 3): 3},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 8): 2},
        (0, 1, 2, 3, 4):{ (0, 1, 2, 3, 4, 5): 3},
        (0, 1, 2, 3, 8):{ (0, 1, 2, 3, 4): 2},
        (0, 1, 3):{ (0, 1, 3, 7): 1},
        (0, 1, 3, 7):{ (0, 1, 2, 3): 1, (0, 1, 3, 7, 8): 2},
        (0, 1, 3, 7, 8):{ (0, 1, 2, 3, 8): 1},
        (0, 1, 7):{ (0, 1, 2): 1, (0, 1, 3, 7): 3},
        (0, 6):{ (0, 1): 2, (0, 6, 7): 1},
        (0, 6, 7):{ (0, 1, 7): 2}
        })

    if (G.is_isomorphic(GT) != True):
        return False;
    else:
        print "C3 example passed."

    # type B, fig. 8.1 pg. 172
    C = ClassicalCrystalOfAlcovePaths(['B',3],[2,0,0]);
    G = C.digraph();

    GT = DiGraph({
        ():{ (6): 1},
        (0):{ (0, 7): 2},
        (0, 1):{ (0, 1, 11): 3},
        (0, 1, 2):{ (0, 1, 2, 9): 2},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 10): 1},
        (0, 1, 2, 3, 10):{ (0, 1, 2, 3, 4): 1},
        (0, 1, 2, 9):{ (0, 1, 2, 3): 2, (0, 1, 2, 9, 10): 1},
        (0, 1, 2, 9, 10):{ (0, 1, 2, 3, 10): 2},
        (0, 1, 5):{ (0, 1, 2): 3, (0, 1, 5, 9): 2},
        (0, 1, 5, 9):{ (0, 1, 2, 9): 3, (0, 1, 5, 9, 10): 1},
        (0, 1, 5, 9, 10):{ (0, 1, 2, 9, 10): 3},
        (0, 1, 8):{ (0, 1, 5): 3},
        (0, 1, 8, 9):{ (0, 1, 5, 9): 3, (0, 1, 8, 9, 10): 1},
        (0, 1, 8, 9, 10):{ (0, 1, 5, 9, 10): 3},
        (0, 1, 11):{ (0, 1, 8): 3},
        (0, 7):{ (0, 1): 2, (0, 7, 11): 3},
        (0, 7, 8):{ (0, 7, 8, 9): 2},
        (0, 7, 8, 9):{ (0, 1, 8, 9): 2},
        (0, 7, 8, 9, 10):{ (0, 1, 8, 9, 10): 2},
        (0, 7, 11):{ (0, 1, 11): 2, (0, 7, 8): 3},
        (6):{ (0): 1, (6, 7): 2},
        (6, 7):{ (0, 7): 1, (6, 7, 11): 3},
        (6, 7, 8):{ (0, 7, 8): 1, (6, 7, 8, 9): 2},
        (6, 7, 8, 9):{ (6, 7, 8, 9, 10): 1},
        (6, 7, 8, 9, 10):{ (0, 7, 8, 9, 10): 1},
        (6, 7, 11):{ (0, 7, 11): 1, (6, 7, 8): 3}
        })

    if (G.is_isomorphic(GT) != True):
        return False;
    else:
        print "B3 example 1 passed."

    C = ClassicalCrystalOfAlcovePaths(['B',3],[0,1,0]);
    G = C.digraph();

    GT = DiGraph({
        ():{ (0): 2},
        (0):{ (0, 1): 1, (0, 7): 3},
        (0, 1):{ (0, 1, 7): 3},
        (0, 1, 2):{ (0, 1, 2, 8): 2},
        (0, 1, 2, 3):{ (0, 1, 2, 3, 5): 1, (0, 1, 2, 3, 9): 3},
        (0, 1, 2, 3, 4):{ (0, 1, 2, 3, 4, 5): 1},
        (0, 1, 2, 3, 4, 5):{ (0, 1, 2, 3, 4, 5, 6): 2},
        (0, 1, 2, 3, 5):{ (0, 1, 2, 3, 5, 9): 3},
        (0, 1, 2, 3, 5, 9):{ (0, 1, 2, 3, 4, 5): 3},
        (0, 1, 2, 3, 9):{ (0, 1, 2, 3, 4): 3, (0, 1, 2, 3, 5, 9): 1},
        (0, 1, 2, 5):{ (0, 1, 2, 3, 5): 2},
        (0, 1, 2, 8):{ (0, 1, 2, 3): 2},
        (0, 1, 2, 8, 9):{ (0, 1, 2, 3, 9): 2},
        (0, 1, 7):{ (0, 1, 2): 3, (0, 1, 7, 8): 2},
        (0, 1, 7, 8):{ (0, 1, 7, 8, 9): 3},
        (0, 1, 7, 8, 9):{ (0, 1, 2, 8, 9): 3},
        (0, 2):{ (0, 1, 2): 1, (0, 2, 5): 2},
        (0, 2, 5):{ (0, 2, 5, 8): 1},
        (0, 2, 5, 8):{ (0, 1, 2, 5): 1},
        (0, 7):{ (0, 1, 7): 1, (0, 2): 3}
        })

    if (G.is_isomorphic(GT) != True):
        return False;
    else:
        print "B3 example 2 passed."

    # type B, fig. 8.3 pg. 174

    return True;

