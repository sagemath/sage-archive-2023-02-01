"""
Incidence structures.

An incidence structure is specified by a list of points, blocks, and
an incidence matrix ([1]_, [2]_).

REFERENCES:

.. [1] Block designs and incidence structures from wikipedia,
  :wikipedia:`Block_design`
  :wikipedia:`Incidence_structure`

.. [2] E. Assmus, J. Key, Designs and their codes, CUP, 1992.

AUTHORS:

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of part of the module block_design.py
  (version 0.6) written by Peter Dobcsanyi peter@designtheory.org.


Classes and methods
-------------------
"""
#***************************************************************************
#                              Copyright (C) 2007                          #
#                                                                          #
#                Peter Dobcsanyi       and         David Joyner            #
#           <peter@designtheory.org>          <wdjoyner@gmail.com>         #
#                                                                          #
#                                                                          #
#    Distributed under the terms of the GNU General Public License (GPL)   #
#    as published by the Free Software Foundation; either version 2 of     #
#    the License, or (at your option) any later version.                   #
#                    http://www.gnu.org/licenses/                          #
#***************************************************************************

from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial
from sage.misc.decorators import rename_keyword

###  utility functions  -------------------------------------------------------

def coordinatewise_product(L):
    """
    Returns the coordinatewise product of a list of vectors.

    INPUT:

    - ``L`` is a list of `n`-vectors or lists all of length `n` with a common
      parent. This returns the vector whose `i`-th coordinate is the product of
      the `i`-th coordinates of the vectors.

    EXAMPLES::

        sage: from sage.combinat.designs.incidence_structures import coordinatewise_product
        sage: L = [[1,2,3],[-1,-1,-1],[5,7,11]]
        sage: coordinatewise_product(L)
        [-5, -14, -33]
    """
    n = len(L[0])
    ans = [1]*n
    for x in L:
        ans = [ans[i]*x[i] for i in range(n)]
    return ans

def IncidenceStructureFromMatrix(M, name=None):
    """
    Builds and incidence structure from a matrix.

    INPUT:

    - ``M`` -- a binary matrix. Creates a set of "points" from the rows and a
      set of "blocks" from the columns.

    EXAMPLES::

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
    def __init__(self, pts, blks, inc_mat=None, name=None, test=True):
        """
        INPUT:

        - ``pts, blks`` -- a list of points, and a list of lists (list of blocks).

          If each `B` in ``blks`` is contained in ``pts`` then the incidence
          matrix ` inc_mat`` need not (and should not) be given.  Otherwise,
          ``inc_mat`` should be the ``ptsxblks`` `(0,1)`-matrix `A` for which
          `A_i,j=1` iff ``blks[j]`` is incident with ``pts[i]``.

        - ``inc_mat`` (for giving the `(0,1)`-incidence matrix)

        - ``name`` (a string, such as "Fano plane").

        - ``test`` (boolean) - if ``True``, then each block must be a list of pts.

        EXAMPLES::

            sage: IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            Incidence structure with 7 points and 7 blocks

        Points are sorted  ::

            sage: BD1 = IncidenceStructure([4,6,0,3,2,5,1],[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD1.points()
            [0, 1, 2, 3, 4, 5, 6]

        TESTS:

        The following shows that Trac Ticket #11333  is fixed.  ::

            sage: A = IncidenceStructure([0,1],[[0]])
            sage: B = IncidenceStructure([1,0],[[0]])
            sage: B==A
            True

        REFERENCES:

        - E. Assmus, J. Key, Designs and their codes, CUP, 1992.
        """
        bs = []
        self.pnts = pts
        self.pnts.sort()
        v, blocks = len(pts), blks
        for block in blocks:
            if test:
                for x in block:
                    if not(x in self.pnts):
                        raise ValueError('Point %s is not in the base set.'%x)
            try:
                y = block[:]
                y.sort()
                bs.append(y)
            except StandardError:
                bs.append(block)
        bs.sort(cmp)
        self.v = v
        self.blcks = bs
        self.name = name
        self._incidence_matrix = inc_mat

    def __iter__(self):
        """
        Iterator over the blocks.

        EXAMPLE::

            sage: sts = designs.steiner_triple_system(9)
            sage: list(sts)
            [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7], [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]
        """

        return iter(self.blcks)


    def __repr__(self):
        """
        A print method.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD
            Incidence structure with 7 points and 7 blocks
        """
        repr = 'Incidence structure with %s points and %s blocks'%(len(self.pnts),len(self.blcks))
        return repr

    def __str__(self):
        """
        A print method.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            BlockDesign<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
        """
        if self.name:
            repr = '%s<points=%s, blocks=%s>'%(self.name, self.pnts, self.blcks)
        else:
            repr = 'IncidenceStructure<points=%s, blocks=%s>'%( self.pnts, self.blcks)
        return repr

    def automorphism_group(self):
        """
        Returns the subgroup of the automorphism group of the incidence graph
        which respects the P B partition. It is (isomorphic to) the automorphism
        group of the block design, although the degrees differ.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: G = BD.automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: BD = BlockDesign(4,[[0],[0,1],[1,2],[3,3]],test=False)
            sage: G = BD.automorphism_group(); G
            Permutation Group with generators [()]
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

    def block_design_checker(self, t, v, k, lmbda, type=None):
        """
        This is *not* a wrapper for GAP Design's IsBlockDesign. The GAP
        Design function IsBlockDesign
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm
        apparently simply checks the record structure and no mathematical
        properties. Instead, the function below checks some necessary (but
        not sufficient) "easy" identities arising from the identity.

        INPUT:

        - ``t`` - the t as in "t-design"

        - ``v`` - the number of points

        - ``k`` - the number of blocks incident to a point

        - ``lmbda`` - each t-tuple of points should be incident with
          lmbda blocks

        - ``type`` - can be 'simple' or 'binary' or 'connected'
          Depending on the option, this wraps IsBinaryBlockDesign,
          IsSimpleBlockDesign, or IsConnectedBlockDesign.

          - Binary: no block has a repeated element.

          - Simple: no block is repeated.

          - Connected: its incidence graph is a connected graph.

        WARNING: This is very fast but can return false positives.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.parameters()
            (2, 7, 3, 1)
            sage: BD.block_design_checker(2, 7, 3, 1)
            True
            sage: BD.block_design_checker(2, 7, 3, 1,"binary")
            True
            sage: BD.block_design_checker(2, 7, 3, 1,"connected")
            True
            sage: BD.block_design_checker(2, 7, 3, 1,"simple")
            True
        """
        from sage.sets.set import Set
        if not(v == len(self.points())):
            return False
        b = lmbda*binomial(v,t)/binomial(k,t)
        r = int(b*k/v)
        if not(b == len(self.blocks())):
            return False
        if not(ZZ(v).divides(b*k)):
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

    def blocks(self):
        """
        Return the list of blocks.

        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
        """
        self._block_sizes = map(len,self.blocks())
        return self._block_sizes

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES::

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

    def dual_incidence_structure(self, algorithm=None):
        """
        Wraps GAP Design's DualBlockDesign (see [1]). The dual of a block
        design may not be a block design.

        Also can be called with ``dual_design``.

        .. NOTE:

        The algorithm="gap" option requires GAP's Design package (included in
        the gap_packages Sage spkg).

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: D = BlockDesign(4, [[0,2],[1,2,3],[2,3]], test=False)
            sage: D
            Incidence structure with 4 points and 3 blocks
            sage: D.dual_design()
            Incidence structure with 3 points and 4 blocks
            sage: print D.dual_design(algorithm="gap")       # optional - gap_design
            IncidenceStructure<points=[0, 1, 2], blocks=[[0], [0, 1, 2], [1], [1, 2]]>
            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="FanoPlane")
            sage: BD
            Incidence structure with 7 points and 7 blocks
            sage: print BD.dual_design(algorithm="gap")         # optional - gap_design
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD.dual_incidence_structure()
            Incidence structure with 7 points and 7 blocks

        REFERENCE:

        - Soicher, Leonard, Design package manual, available at
          http://www.gap-system.org/Manuals/pkg/design/htm/CHAP003.htm
        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        from sage.misc.flatten import flatten
        from sage.combinat.designs.block_design import BlockDesign
        from sage.misc.functional import transpose
        if algorithm=="gap":
            gap.load_package("design")
            gD = self._gap_()
            gap.eval("DD:=DualBlockDesign("+gD+")")
            v = eval(gap.eval("DD.v"))
            gblcks = eval(gap.eval("DD.blocks"))
            gB = []
            for b in gblcks:
                gB.append([x-1 for x in b])
            return IncidenceStructure(range(v), gB, name=None, test=False)
        pts = self.blocks()
        M = transpose(self.incidence_matrix())
        blks = self.points()
        return IncidenceStructure(pts, blks, M, name=None, test=False)

    dual_design = dual_incidence_structure  # to preserve standard terminology

    def incidence_matrix(self):
        """
        Return the incidence matrix `A` of the design. A is a `(v \times b)`
        matrix defined by: ``A[i,j] = 1`` if ``i`` is in block ``B_j`` and 0
        otherwise.

        EXAMPLES::

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
        """
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
        Returns the incidence graph of the design, where the incidence
        matrix of the design is the adjacency matrix of the graph.

        EXAMPLE::

            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.incidence_graph()
            Bipartite graph on 14 vertices
            sage: A = BD.incidence_matrix()
            sage: Graph(block_matrix([[A*0,A],[A.transpose(),A*0]])) == BD.incidence_graph()
            True

        REFERENCE:

        - Sage Reference Manual on Graphs
        """
        from sage.graphs.bipartite_graph import BipartiteGraph
        A = self.incidence_matrix()
        return BipartiteGraph(A)
        #same as return Graph(block_matrix([[A*0,A],[A.transpose(),A*0]]))

    def is_block_design(self):
        """
        Returns a pair ``True, pars`` if the incidence structure is a
        `t`-design, for some `t`, where ``pars`` is the list of parameters `(t,
        v, k, lmbda)`.  The largest possible `t` is returned, provided `t=10`.

        EXAMPLES::

            sage: BD = IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.is_block_design()
            (True, [2, 7, 3, 1])
            sage: BD.block_design_checker(2, 7, 3, 1)
            True
            sage: BD = designs.WittDesign(9)   # optional - gap_packages (design package)
            sage: BD.is_block_design()         # optional - gap_packages (design package)
            (True, [2, 9, 3, 1])
            sage: BD = designs.WittDesign(12)  # optional - gap_packages (design package)
            sage: BD.is_block_design()         # optional - gap_packages (design package)
            (True, [5, 12, 6, 1])
            sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: BD.is_block_design()
            (True, [2, 8, 2, 2])
        """
        from sage.combinat.designs.incidence_structures import coordinatewise_product
        from sage.combinat.combinat import unordered_tuples
        from sage.combinat.combination import Combinations
        A = self.incidence_matrix()
        v = len(self.points())
        b = len(self.blocks())
        k = sum(A.columns()[0])
        rowsA = A.rows()
        VS = rowsA[0].parent()
        r = sum(rowsA[0])
        for i in range(b):
            if not(sum(A.columns()[i]) == k):
                return False
        for i in range(v):
            if not(sum(A.rows()[i]) == r):
                return False
        t_found_yet = False
        lambdas = []
        for t in range(2,min(v,11)):
            #print t
            L1 = Combinations(range(v),t)
            L2 = [[rowsA[i] for i in L] for L in L1]
            #print t,len(L2)
            lmbda = VS(coordinatewise_product(L2[0])).hamming_weight()
            lambdas.append(lmbda)
            pars = [t,v,k,lmbda]
            #print pars
            for ell in L2:
                a = VS(coordinatewise_product(ell)).hamming_weight()
                if not(a == lmbda) or a==0:
                    if not(t_found_yet):
                        pars = [t-1,v,k,lambdas[t-3]]
                        return False, pars
                    else:
                        #print pars, lambdas
                        pars = [t-1,v,k,lambdas[t-3]]
                        return True, pars
                t_found_yet = True
        pars = [t-1,v,k,lambdas[t-3]]
        return True, pars

    def parameters(self, t=2):
        """
        Returns `(t,v,k,lambda)`. Does not check if the input is a block
        design. Uses `t=2` by default.

        EXAMPLES::

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

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.points()
            [0, 1, 2, 3, 4, 5, 6]
        """
        return self.pnts

    def points_from_gap(self):
        """
        Literally pushes this block design over to GAP and returns the
        points of that. Other than debugging, usefulness is unclear.

        REQUIRES: GAP's Design package.

        EXAMPLES::

            sage: from sage.combinat.designs.block_design import BlockDesign
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.points_from_gap()      # optional - gap_packages (design package)
            doctest:1: DeprecationWarning: Unless somebody protests this method will be removed, as nobody seems to know why it is there.
            See http://trac.sagemath.org/14499 for details.
            [1, 2, 3, 4, 5, 6, 7]
        """
        from sage.misc.superseded import deprecation
        deprecation(14499, ('Unless somebody protests this method will be '
                            'removed, as nobody seems to know why it is there.'))
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        gap.load_package("design")
        gD = self._gap_()
        gP = gap.eval("BlockDesignPoints("+gD+")").replace("..",",")
        return range(eval(gP)[0],eval(gP)[1]+1)



