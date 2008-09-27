"""
Incidence structures.

An incidence structure is specified by a list of points, blocks,
and an incidence matrix ([1], [2]).

Classes:

    IncidenceStructure

This software is released under the terms of the GNU General Public License,
version 2 or above (your choice). For details on licencing, see the
accompanying documentation.

REFERENCES:
  [1] Block designs and incidence structures from wikipedia,
      http://en.wikipedia.org/wiki/Block_design
      http://en.wikipedia.org/wiki/Incidence_structure
  [2] E. Assmus, J. Key, Designs and their codes, CUP, 1992.

This is a significantly modified form of part of the module
block_design.py (version 0.6) written by Peter Dobcsanyi
<peter@designtheory.org>.

Copyright 2007-2008 by David Joyner <wdjoyner@gmail.com>,
Peter Dobcsanyi <peter@designtheory.org>.

"""

import types
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial, integer_floor

###  utility functions  -------------------------------------------------------

def block_cmp(b1, b2):
    """
    Auxiliary comparison function to sort blocks which are sorted lists.
    It compares two blocks first by length then by lexicographic order.
    It returns -1 if b1 is "smaller" than b2; 0 if they are equal; and
    +1 otherwise.

    EXAMPLES:
        sage: from sage.combinat.designs.incidence_structures import block_cmp
        sage: L1 = [1,2,3,4]
        sage: L2 = [5,6,7]
        sage: L1<L2
        True
        sage: block_cmp(L1, L2)
        1
        sage: block_cmp(L2, L1)
        -1

    """
    if len(b1) < len(b2):
        return -1
    elif len(b1) == len(b2):
        if b1 < b2:
            return -1
        elif b1 == b2:
            return 0
        else:
            return 1
    else:
        return 1

def IncidenceStructureFromMatrix(M, name=None):
    """
    M must be a (0,1)-matrix. Creates a set of "points" from the rows
    and a set of "blocks" from the columns.

    EXAMPLES:
        sage: from sage.combinat.designs.block_design import BlockDesign
        sage: BD1 = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: M = BD1.incidence_matrix()
        sage: BD2 = IncidenceStructureFromMatrix(M)
        sage: BD1 == BD2
        True

    """
    nm = name
    v = len(M.rows())
    b = len(M.columns())
    #points = range(v)
    blocks = []
    for i in range(b):
        B = []
        for j in range(v):
            if M[i,j]!=0:
                B.append(j)
        blocks.append(B)
    return IncidenceStructure(range(v), blocks, name=nm)

class IncidenceStructure(object):
    """
    This the base class for block designs.

    """
    def __init__(self, pts, blks, inc_mat=None, name=None):
        """
        The parameters are a pair pts, blks, both of which are a list
        (blks is a list of lists). If each B in blks is contained in pts
        then the incidence matrix inc_mat need not (and should not) be
        given. Otherwise, inc_mat should be the |pts|x|blks| (0,1)-matrix A
        for which A_{i,j}=1 iff blks[j] is incident with pts[i].

        Optional keywords are:
            "inc_mat" (for giving the (0,1)-incidence matrix), and
            "name" (a string, such as "Fano plane").

        EXAMPLES:
            sage: IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            IncidenceStructure<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

        REFERENCES:
            E. Assmus, J. Key, Designs and their codes, CUP, 1992.

        """
        from sage.combinat.designs.incidence_structures import block_cmp
        bs = []
        self.pnts = pts
        v, blocks = len(pts), blks
        for block in blocks:
            for x in block:
                if v <= x or x < 0:
                    raise ValueError('Point %s is not in the base set.'%x)
            y = block[:]
            y.sort()
            bs.append(y)
        bs.sort(block_cmp)
        self.v = v
        self.blcks = bs
        self.name = name
        self._incidence_matrix = inc_mat

    def __repr__(self):
        """
        A print method.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            BlockDesign<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

        """
        if self.name:
            repr = '%s<v=%d, blocks=%s>'%(self.name, self.v, self.blcks)
        else:
            repr = 'IncidenceStructure<v=%d, blocks=%s>'%( self.v, self.blcks)
        return repr

    def automorphism_group(self):
        """
        Returns the subgroup of the automorphism group of the incidence graph
        which respects the P\cup B partition. This is (isomorphic to)
        the automorphism group of the block design, although the degrees differ.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: G = BD.automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: BD = BlockDesign(4,[[0],[0,1],[1,2],[3,3]])
            sage: G = BD.automorphism_group(); G
            Permutation Group with generators []
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: G = BD.automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]

        """
        from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        M1 = self.incidence_matrix()
        M2 =  MatrixStruct(M1)
        M2.run()
        gens = M2.automorphism_group()[0]
        v = len(self.points())
        G = SymmetricGroup(v)
        gns = []
        for g in gens:
            L = [j+1 for j in g]
            gns.append(G(L))
        return PermutationGroup(gns)

    def blocks(self):
        """
        Return the list of blocks.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]

        """
        B = self.blcks
        B.sort()
        return B

    def __eq__(self, other):
        """
        Returns true if their points and blocks are equal (resp.).

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD1 = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()
            sage: BD2 = IncidenceStructureFromMatrix(M)
            sage: BD1 == BD2
            True

        """
        bool1 = self.points() == other.points()
        bool2 = self.blocks() == other.blocks()
        return (bool1 and bool2)

    def block_sizes(self):
        """
        Return a list of block's sizes.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]

        """
        bs = []
        for b in self.blocks():
            bs.append(len(b))
        self._block_sizes = bs
        return bs

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD._gap_()
            'BlockDesign(7,[[1, 2, 3], [1, 4, 5], [1, 6, 7], [2, 4, 6], [2, 5, 7], [3, 4, 7], [3, 5, 6]])'

        """
        from sage.sets.set import Set
        B = self.blocks()
        v = len(self.points())
        gB = []
        for b in B:
           gB.append([x+1 for x in b])
        return "BlockDesign("+str(v)+","+str(gB)+")"

    def dual_block_design(self):
        """
        Wraps GAP Design's DualBlockDesign (see [1])

        REQUIRES: GAP's Design package.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
           sage: D = BlockDesign(4, [[0,2],[1,2,3],[2,3]])
           sage: D
           BlockDesign<v=4, blocks=[[0, 2], [2, 3], [1, 2, 3]]>
           sage: D.dual_block_design()          # requires optional gap package
           BlockDesign<v=3, blocks=[[0], [1], [1, 2], [0, 1, 2]]>
           sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="FanoPlane")
           sage: BD
           FanoPlane<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
           sage: BD.dual_block_design()         # requires optional gap package
           BlockDesign<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

        REFERENCE:
          Soicher, Leonard, Design package manual, available at
          http://www.gap-system.org/Manuals/pkg/design/htm/CHAP003.htm
        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        from sage.misc.flatten import flatten
        from sage.combinat.designs.block_design import BlockDesign
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        gap.eval("DD:=DualBlockDesign("+gD+")")
        v = eval(gap.eval("DD.v"))
        gblcks = eval(gap.eval("DD.blocks"))
        gB = []
        for b in gblcks:
           gB.append([x-1 for x in b])
        return BlockDesign(v, gB)

    def incidence_matrix(self):
        '''
        Return the incidence matrix A of the design.
        A is a (v x b) matrix defined by:
            A[i,j] = 1   if i is in block B_j
                     0   otherwise

        EXAMPLES:
            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
            sage: BD.incidence_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [1 0 0 0 0 1 1]
            [0 1 0 1 0 1 0]
            [0 1 0 0 1 0 1]
            [0 0 1 1 0 0 1]
            [0 0 1 0 1 1 0]

        '''
        if self._incidence_matrix!=None:
            return self._incidence_matrix
        else:
            v = len(self.points())
            blks = self.blocks()
            b = len(blks)
            MS = MatrixSpace(ZZ,v,b)
            A = MS(0)
            #A = NUM.zeros((v,b), NUM.Int)
            for i in range(v):
                for j, b in enumerate(blks):
                    if i in b:
                        A[i,j] = 1
            self._incidence_matrix = A
            return A

    def incidence_graph(self):
        """
        Returns the incidence graph of the design, where the incidence matrix
        of the design is the adjacency matrix of the graph.

        EXAMPLE:
            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.incidence_graph()
            Bipartite graph on 14 vertices
            sage: A = BD.incidence_matrix()
            sage: Graph(block_matrix([A*0,A,A.transpose(),A*0])) == BD.incidence_graph()
            True

        REFERENCE:
           Sage Reference Manual on Graphs,
           http://www.sagemath.org/doc/ref/node44.html
        """
        from sage.graphs.bipartite_graph import BipartiteGraph
        A = self.incidence_matrix()
        return BipartiteGraph(A)
        #same as return Graph(block_matrix([A*0,A,A.transpose(),A*0]))

    def is_block_design(self, t, v, k, lmbda, type=None):
        """
        This is *not* just a wrapper for GAP Design's IsBlockDesign.
        The GAP Design function IsBlockDesign
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm#SSEC001.1
        apparently simply checks the record structure and no mathematical
        properties. Instead, the function below checks the definition.

        OPTIONS:
            type can be "simple or "binary" or "connected"
            Depending on the option, this wraps IsBinaryBlockDesign,
            IsSimpleBlockDesign, or IsConnectedBlockDesign.

        Binary: no block has a repeated element.
        Simple: no block is repeated.
        Connected: its incidence graph is a connected graph.

        REQUIRES: GAP's Design package.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.parameters()
            (2, 7, 3, 1)
            sage: BD.is_block_design(2, 7, 3, 1)
            True
            sage: BD.is_block_design(2, 7, 3, 1,"binary")
            True
            sage: BD.is_block_design(2, 7, 3, 1,"connected")
            True
            sage: BD.is_block_design(2, 7, 3, 1,"simple")
            True

        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        if not(v == len(self.points())):
            return False
        b = lmbda*binomial(v,t)/binomial(k,t)
        r = int(b*k/v)
        if not(b == len(self.blocks())):
            return False
        if not(v.divides(b*k)):
            return False
        A = self.incidence_matrix()
        #k = sum(A.columns()[0])
        #r = sum(A.rows()[0])
        for i in range(b):
            if not(sum(A.columns()[i]) == k):
                return False
        for i in range(v):
            if not(sum(A.rows()[i]) == r):
                return False
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        if type==None:
            return True
        if type=="binary":
            for b in self.blocks():
                if len(b)!=len(Set(b)):
                     return False
            return True
        if type=="simple":
            B = self.blocks()
            for b in B:
                 if B.count(b)>1:
                     return False
            return True
        if type=="connected":
            Gamma = self.incidence_graph()
            if Gamma.is_connected():
                return True
            else:
                return False

    def pairwise_balanced_lambda(self):
        """
        Returns 0 if D is not (binary or) pairwise balanced, and otherwise the positive
        constant lambda such that every pair of distinct points of D is in exactly
        lambda blocks.

        A binary block design D is pairwise balanced if D has at least two
        points and every pair of distinct points is contained in exactly lambda
        blocks, for some positive constant lambda.

        This does not check if the input is a binary design.

        Wraps GAP Design's PairwiseBalancedLambda
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm

        REQUIRES: GAP's Design package.

        EXAMPLES:
            sage: blks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = IncidenceStructure(range(7), blks)
            sage: BD.pairwise_balanced_lambda()      # requires optional gap package
            1

        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        out = gap.eval("PairwiseBalancedLambda("+gD+")")
        if out=="fail":
            return 0
        return eval(out)

    def parameters(self, t=2):
        """
        Returns (t,v,k,lambda). Does not check if the input is a
        block design. Uses t=2 by default.


        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="FanoPlane")
            sage: BD.parameters()
            (2, 7, 3, 1)
            sage: BD.parameters(t=3)
            (3, 7, 3, 0)

        """
        v = len(self.points())
        blks = self.blocks()
        k = len(blks[int(0)])
        b = len(blks)
        #A = self.incidence_matrix()
        #r = sum(A.rows()[0])
        lmbda = int(b/(binomial(v,t)/binomial(k,t)))
        return (t,v,k,lmbda)

    def points(self):
        """
        Returns the list of points.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.points()
            [0, 1, 2, 3, 4, 5, 6]

        """
        return self.pnts

    def points_from_gap(self):
        """
        Literally pushes this block design over to GAP and
        returns the points of that. Other than debugging, usefulness
        is unclear.

        REQUIRES: GAP's Design package.

        EXAMPLES:
            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.points_from_gap()      # requires optional gap package
            [1, 2, 3, 4, 5, 6, 7]

        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        gP = gap.eval("BlockDesignPoints("+gD+")").replace("..",",")
        return range(eval(gP)[0],eval(gP)[1]+1)



