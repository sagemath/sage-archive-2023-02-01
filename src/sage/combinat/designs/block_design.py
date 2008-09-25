"""
A module to help with constructions and computations of block designs and
other incidence structures.

A block design is an incidence structure consisting of a set of points P
and a set of blocks B, where each block is considered as a subset of P.
More precisely, a *block design* B is a class of k-element subsets of P such that
the number r of blocks that contain any point x in P is independent of x,
and the number lambda of blocks that contain any given t-element
subset T is independent of the choice of T (see [1] for more). Such a block
design is also called a t-(v,k,lambda)-design, and v (the number of points),
b (the number of blocks), k, r, and lambda are the parameters of the design.
(In Python, lambda is reserved, so we sometimes use lmbda or L instead.)

In Sage, sets are replaced by (ordered) lists and the standard
representation of a block design uses P = [0,1,..., v-1], so a block
design is specified by (v,B). However, an incidence structure can be
specified with points, blocks, and an incidence matrix.

Classes:

    IncidenceStructure

This software is released under the terms of the GNU General Public License,
version 2 or above (your choice). For details on licencing, see the
accompanying documentation.

REFERENCES:
  [1] Block design from wikipedia, http://en.wikipedia.org/wiki/Block_design
  [2] What is a block design?, http://designtheory.org/library/extrep/html/node4.html
      (in "The External Representation of Block Designs" by Peter J. Cameron,
      Peter Dobcsanyi, John P. Morgan, Leonard H. Soicher)

This is a significantly modified form of the module block_design.py (version 0.6)
written by Peter Dobcsanyi <peter@designtheory.org>.

Copyright 2007-2008 by Peter Dobcsanyi <peter@designtheory.org>,
and David Joyner <wdjoyner@gmail.com>.

TODO: Implement DerivedDesign, ComplementaryDesign, Hadamard3Design

"""

import types
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial
from sage.rings.finite_field import FiniteField

###  utility functions  -------------------------------------------------------

def block_cmp(b1, b2):
    """
    Auxiliary comparison function to sort blocks which are sorted lists.
    It compares two blocks first by length then by lexicographic order.
    It returns -1 if b1 is "smaller" than b2; 0 if they are equal; and
    +1 otherwise.

    EXAMPLES:
        sage: from sage.combinat.designs.block_design import block_cmp
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

def tdesign_params(t, v, k, L):
    '''
    Return the design's parameters: (t, v, b, r , k, L)
    '''
    x = binomial(v, t)
    y = binomial(k, t)
    b = divmod(L * x, y)[0]
    x = binomial(v-1, t-1)
    y = binomial(k-1, t-1)
    r = int_div(L * x, y)
    return (t, v, b, r, k, L)

def ProjectiveGeometryDesign(n, d, F, method=None):
    """
    Input: n is the projective dimension, so the number of points is
             v = |PP^n(GF(q))|
           d is the dimension of the subspaces of P = PP^n(GF(q))
             which make up the blocks, so b is the number of d-dimensional
             subspaces of P

    Wraps GAP Design's PGPointFlatBlockDesign. Does *not* require GAP's Design.

    EXAMPLES:
        sage: ProjectiveGeometryDesign(2, 1, GF(2))
        ProjectiveGeometryDesign<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
        sage: ProjectiveGeometryDesign(2, 1, GF(2), method="gap")      # requires optional gap package
        ProjectiveGeometryDesign<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>


    """
    q = F.order()
    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    if method == None:
        gap.eval("V:=GaloisField(%s)^%s"%(q,n+1))
        gap.eval("points:=AsSet(Subspaces(V,1))")
        gap.eval("flats:=AsSet(Subspaces(V,%s));"%(d+1))
        gblcks = eval(gap.eval("AsSortedList(List(flats,f->Filtered([1..Length(points)],i->IsSubset(f,points[i]))));"))
        v = (q**(n+1)-1)/(q-1)
        gB = []
        for b in gblcks:
            gB.append([x-1 for x in b])
        return BlockDesign(v, gB, name="ProjectiveGeometryDesign")
    if method == "gap":   # Requires GAP's Design
        gap.eval('LoadPackage("design")')
        gap.eval("D := PGPointFlatBlockDesign( %s, %s, %d )"%(n,q,d))
        v = eval(gap.eval("D.v"))
        gblcks = eval(gap.eval("D.blocks"))
        gB = []
        for b in gblcks:
            gB.append([x-1 for x in b])
        return BlockDesign(v, gB, name="ProjectiveGeometryDesign")

def AffineGeometryDesign(n, d, F):
    """
    Input: n is the Euclidian dimension, so the number of points is
             v = |GF(q)^n|
           d is the dimension of the (affine) subspaces of P = GF(q)^n
             which make up the blocks.

    Wraps some functions used in GAP Design's PGPointFlatBlockDesign.
    Does *not* require GAP's Design.

    EXAMPLES:
        sage: BD = AffineGeometryDesign(3, 1, GF(2))
        sage: BD.parameters()
        (2, 8, 2, 2)
        sage: BD.is_block_design(2,8,2,2)
        True
        sage: BD = AffineGeometryDesign(3, 2, GF(2))
        sage: BD.parameters()
        (2, 8, 4, 12)
        sage: BD.is_block_design(2,8,4,12)
        True

    """
    q = F.order()
    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    gap.eval("V:=GaloisField(%s)^%s"%(q,n))
    gap.eval("points:=AsSet(V)")
    gap.eval("Subs:=AsSet(Subspaces(V,%s));"%d)
    gap.eval("CP:=Cartesian(points,Subs)")
    flats = gap.eval("flats:=List(CP,x->Sum(x))") # affine spaces
    gblcks = eval(gap.eval("AsSortedList(List(flats,f->Filtered([1..Length(points)],i->points[i] in f)));"))
    v = q**n
    gB = []
    for b in gblcks:
       gB.append([x-1 for x in b])
    return BlockDesign(v, gB, name="AffineGeometryDesign")

def WittDesign(n):
    """
    Input: n is in {9,10,11,12,21,22,23,24}.

    Wraps GAP Design's WittDesign. If n=24 then this function returns
    the large Witt design W24, the unique (up to isomorphism) 5-(24,8,1) design.
    If n=12 then this function returns the small Witt design W12, the
    unique (up to isomorphism) 5-(12,6,1) design. The other values of n return
    a block design derived from these.

    REQUIRES: GAP's Design package.

    EXAMPLES:
        sage: BD = WittDesign(9)   # requires optional gap package
        sage: BD.parameters()      # requires optional gap package
        (2, 9, 3, 1)
        sage: BD                   # requires optional gap package
        WittDesign<v=9, blocks=[[0, 1, 7], [0, 2, 5], [0, 3, 4], [0, 6, 8], [1, 2, 6], [1, 3, 5], [1, 4, 8], [2, 3, 8], [2, 4, 7], [3, 6, 7], [4, 5, 6], [5, 7, 8]]>
        sage: BD = WittDesign(12)  # requires optional gap package
        sage: BD.parameters(t=5)   # requires optional gap package
        (5, 12, 6, 1)

    """
    from sage.interfaces.gap import gap, GapElement
    gap.eval('LoadPackage("design")')
    gap.eval("B:=WittDesign(%s)"%n)
    v = eval(gap.eval("B.v"))
    gblcks = eval(gap.eval("B.blocks"))
    gB = []
    for b in gblcks:
       gB.append([x-1 for x in b])
    return BlockDesign(v, gB, name="WittDesign")

def HadamardDesign(n):
    """
    As described in \S 1, p. 10, in [CvL]. The input n must
    have the property that there is a Hadamard matrix of order
    n+1 (and that a construction of that Hadamard matrix has been
    implemented...).


    EXAMPLES:
        sage: HadamardDesign(7)
        HadamardDesign<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    REFERENCES:
        [CvL] P. Cameron, J. H. van Lint, Designs, graphs, codes
        and their links, London Math. Soc., 1991.
    """
    from sage.combinat.combinat import hadamard_matrix
    from sage.matrix.constructor import matrix
    H = hadamard_matrix(n+1)
    H1 = H.matrix_from_columns(range(1,n+1))
    H2 = H1.matrix_from_rows(range(1,n+1))
    J = matrix(ZZ,n,n,[1]*n*n)
    MS = J.parent()
    A = MS((H2+J)/2) # convert -1's to 0's; coerce entries to ZZ
    # A is the incidence matrix of the block design
    return IncidenceStructureFromMatrix(A,name="HadamardDesign")

def IncidenceStructureFromMatrix(M, name=None):
    """
    M must be a (0,1)-matrix. Creates a set of "points" from the rows
    and a set of "blocks" from the columns.

    EXAMPLES:
        sage: BD1 = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: M = BD1.incidence_matrix()
        sage: BD2 = IncidenceStructureFromMatrix(M)
        sage: BD1 == BD2
        True

    """
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
    return BlockDesign(v, blocks, name)

def BlockDesign(max_pt, blks, name=None):
    """
    Returns an instance of the IncidenceStructure class. Requires each
    B in blks to be contained in range(max_pt). Does not test if
    the result is a block design.

    EXAMPLES:
        sage: BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="Fano plane")
        Fano plane<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    """
    nm = name
    if nm == None:
        nm = "BlockDesign"
    BD = IncidenceStructure( range(max_pt), blks, inc_mat=None, name=nm )
    #if BD.is_block_design():
    #    return BD
    #else:
    #    raise TypeError("parameters are not those of a block design.")
    return BD


###  classes  -----------------------------------------------------------------

class IncidenceStructure(object):
    '''
    This the base class for block designs.

    EXAMPLES:
    sage: IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
    IncidenceStructure<v=7, blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    '''
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

        REFERENCES:
            E. Assmus, J. Key, Designs and their codes, CUP, 1992.

        """
        bs = []
        self.pnts = pts
        v, blocks = len(pts), blks
        for block in blocks:
            for x in block:
                if v <= x or x < 0:
                    raise sage.combinat.designs.Error(
                        'Point "%s" is not in the base set.' % x)
            y = block[:]
            y.sort()
            bs.append(y)
        bs.sort(block_cmp)
        self.v = v
        self.blcks = bs
        self.name = name
        self._incidence_matrix = inc_mat

    def __repr__(self):
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
        from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct
        #from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        M1 = self.incidence_matrix()
        #M2 =  MatrixStruct(M1)
        M2 = NonlinearBinaryCodeStruct(M1)
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
        """
        B = self.blcks
        B.sort()
        return B

    def blocks_from_gap(self):
        """
        Literally pushes this block design over to GAP and
        returns the blocks of that. Usefulness is unclear...
        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        return eval(gap.eval("BlockDesignBlocks("+gD+")"))

    def __eq__(self, other):
        """
        Returns true if their points and blocks are equal (resp.).

        EXAMPLES:
            sage: BD1 = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()
            sage: BD2 = IncidenceStructureFromMatrix(M)
            sage: BD1 == BD2
            True

        """
        bool1 = self.points() == other.points()
        bool2 = self.blocks() == other.blocks()
        return (bool1 and bool2)

    def _get_block_sizes(self):
        '''
        Return a list of block's sizes.
        '''
        try:
            return self._block_sizes
        except AttributeError:
            bs = []
            for b in self.blocks:
                bs.append(len(b))
            self._block_sizes = bs
        return bs

    block_sizes = _get_block_sizes

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES:
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
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        gap.eval("DD:=DualBlockDesign("+gD+")")
        v = eval(gap.eval("DD.v"))
        gblcks = eval(gap.eval("DD.blocks"))
        gB = []
        for b in gblcks:
           gB.append([x-1 for x in b])
        return BlockDesign(v, gB)

    def _get_incidence_matrix(self):
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

    incidence_matrix = _get_incidence_matrix

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
        The GAP Design function IsBlockDesign
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm#SSEC001.1
        apparently simply checks the record structure and no mathematical
        properties, so it is not wrapped. Instead, the function below
        checks the definition.

        OPTIONS:
            type can be "simple or "binary" or "connected"
            Depending on the option, this wraps IsBinaryBlockDesign,
            IsSimpleBlockDesign, or IsConnectedBlockDesign.

        Binary: no block has a repeated element.
        Simple: no block is repeated.
        Connected: its incidence graph is a connected graph.

        REQUIRES: GAP's Design package.

        EXAMPLES:
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.parameters()
            (2, 7, 3, 1)
            sage: BD.is_block_design(2, 7, 3, 1)
            True
            sage: BD.is_block_design(2, 7, 3, 1,"binary")      # requires optional gap package
            True
            sage: BD.is_block_design(2, 7, 3, 1,"connected")   # requires optional gap package
            True
            sage: BD.is_block_design(2, 7, 3, 1,"simple")      # requires optional gap package
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
            #if gap.eval("IsBlockDesign("+gD+")")=="true":
            #    return True
            #else:
            #    return False
        if type=="binary":
            if gap.eval("IsBinaryBlockDesign("+gD+")")=="true":
                return True
            else:
                return False
        if type=="simple":
            if gap.eval("IsSimpleBlockDesign("+gD+")")=="true":
                return True
            else:
                return False
        if type=="connected":
            if gap.eval("IsConnectedBlockDesign("+gD+")")=="true":
                return True
            else:
                return False

    def linear_code(self, F = FiniteField(2)):
        """
        Returns the linear code obtained from the F-linear span
        of the rows of the incidence matrix (section 2.4 in [1]).

        EXAMPLES:
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.linear_code()
            Linear code of length 7, dimension 4 over Finite Field of size 2

        Indeed, this is the Hamming [7,4,3] code (Example 2.4.1 in [1]).

        REFERENCES:
           [1] Assmus and Key, Designs and their codes, CUP, 1992.
        """
        from sage.coding.linear_code import LinearCode
        M = self.incidence_matrix().change_ring(F)
        return LinearCode(M)

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
        block design.


        EXAMPLES:
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

        return self.pnts

    def points_from_gap(self):
        """
        Literally pushes this block design over to GAP and
        returns the points of that. Usefulness is unclear
        but could be used for checking Sage/Python constructions
        against GAP Design.

        EXAMPLES:
            sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.points_from_gap()      # requires optional gap package
            [1, 2, 3, 4, 5, 6, 7]
            sage: BD.blocks_from_gap()      # requires optional gap package
            [[1, 2, 3], [1, 4, 5], [1, 6, 7], [2, 4, 6], [2, 5, 7], [3, 4, 7], [3, 5, 6]]

        """
        from sage.interfaces.gap import gap, GapElement
        from sage.sets.set import Set
        gap.eval('LoadPackage("design")')
        gD = self._gap_()
        gP = gap.eval("BlockDesignPoints("+gD+")").replace("..",",")
        return range(eval(gP)[0],eval(gP)[1]+1)



