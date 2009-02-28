"""
Block designs.

A module to help with constructions and computations of block
designs and other incidence structures.

A block design is an incidence structure consisting of a set of
points P and a set of blocks B, where each block is considered as a
subset of P. More precisely, a *block design* B is a class of
k-element subsets of P such that the number r of blocks that
contain any point x in P is independent of x, and the number lambda
of blocks that contain any given t-element subset T is independent
of the choice of T (see [1] for more). Such a block design is also
called a t-(v,k,lambda)-design, and v (the number of points), b
(the number of blocks), k, r, and lambda are the parameters of the
design. (In Python, lambda is reserved, so we sometimes use lmbda
or L instead.)

In Sage, sets are replaced by (ordered) lists and the standard
representation of a block design uses P = [0,1,..., v-1], so a
block design is specified by (v,B).

This software is released under the terms of the GNU General Public
License, version 2 or above (your choice). For details on
licencing, see the accompanying documentation.

REFERENCES:

- [1] Block design from wikipedia,
  http://en.wikipedia.org/wiki/Block_design

- [2] What is a block design?,
  http://designtheory.org/library/extrep/html/node4.html (in 'The
  External Representation of Block Designs' by Peter J. Cameron, Peter
  Dobcsanyi, John P. Morgan, Leonard H. Soicher)

This is a significantly modified form of the module
block_design.py (version 0.6) written by Peter Dobcsanyi
peter@designtheory.org. Thanks go to Robert Miller for lots of good
design suggestions.

Copyright 2007-2008 by Peter Dobcsanyi peter@designtheory.org, and
David Joyner wdjoyner@gmail.com.

TODO: Implement DerivedDesign, ComplementaryDesign,
Hadamard3Design
"""

import types
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial, integer_floor
from sage.rings.finite_field import FiniteField
from sage.combinat.designs.incidence_structures import IncidenceStructure, IncidenceStructureFromMatrix

###  utility functions  -------------------------------------------------------

def tdesign_params(t, v, k, L):
    """
    Return the design's parameters: (t, v, b, r , k, L). Note t must be
    given.

    EXAMPLES::

        sage: BD = BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: from sage.combinat.designs.block_design import tdesign_params
        sage: tdesign_params(2,7,3,1)
        (2, 7, 7, 3, 3, 1)
    """
    x = binomial(v, t)
    y = binomial(k, t)
    b = divmod(L * x, y)[0]
    x = binomial(v-1, t-1)
    y = binomial(k-1, t-1)
    r = integer_floor(L * x/y)
    return (t, v, b, r, k, L)

def ProjectiveGeometryDesign(n, d, F, method=None):
    """
    Input: n is the projective dimension, so the number of points is v
    = PPn(GF(q)) d is the dimension of the subspaces of P = PPn(GF(q))
    which make up the blocks, so b is the number of d-dimensional
    subspaces of P

    Wraps GAP Design's PGPointFlatBlockDesign. Does *not* require
    GAP's Design.

    EXAMPLES::

        sage: ProjectiveGeometryDesign(2, 1, GF(2))
        Incidence structure with 7 points and 7 blocks
        sage: BD = ProjectiveGeometryDesign(2, 1, GF(2), method="gap")      # requires optional gap package
        sage: BD.is_block_design()                                     # requires optional gap package
        (True, [2, 7, 3, 1])
    """
    q = F.order()
    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    if method == None:
        V = VectorSpace(F, n+1)
        points = list(V.subspaces(1))
        flats = list(V.subspaces(d+1))
        blcks = []
        for p in points:
            b = []
            for i in range(len(flats)):
                if p.is_subspace(flats[i]):
                    b.append(i)
            blcks.append(b)
        v = (q**(n+1)-1)/(q-1)
        return BlockDesign(v, blcks, name="ProjectiveGeometryDesign")
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
    r"""
    Input: n is the Euclidean dimension, so the number of points is
    `v = |F^n|` (F = GF(q), some q) d is the dimension of the
    (affine) subspaces of `P = GF(q)^n` which make up the
    blocks.

    `AG_{n,d} (F)`, as it is sometimes denoted, is a
    `2` - `(v, k, \lambda)` design of points and
    `d`- flats (cosets of dimension n) in the affine geometry
    `AG_n (F)`, where

    .. math::

             v = q^n,\  k = q^d ,
             \lambda =\frac{(q^{n-1}-1) \cdots (q^{n+1-d}-1)}{(q^{n-1}-1) \cdots (q-1)}.



    Wraps some functions used in GAP Design's PGPointFlatBlockDesign.
    Does *not* require GAP's Design.

    EXAMPLES::

        sage: BD = AffineGeometryDesign(3, 1, GF(2))
        sage: BD.parameters()
        (2, 8, 2, 2)
        sage: BD.is_block_design()
        (True, [2, 8, 2, 2])
        sage: BD = AffineGeometryDesign(3, 2, GF(2))
        sage: BD.parameters()
        (2, 8, 4, 12)
        sage: BD.is_block_design()
        (True, [3, 8, 4, 4])
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
    Input: n is in 9,10,11,12,21,22,23,24.

    Wraps GAP Design's WittDesign. If n=24 then this function returns
    the large Witt design W24, the unique (up to isomorphism)
    5-(24,8,1) design. If n=12 then this function returns the small
    Witt design W12, the unique (up to isomorphism) 5-(12,6,1) design.
    The other values of n return a block design derived from these.

    REQUIRES: GAP's Design package.

    EXAMPLES::

        sage: BD = WittDesign(9)   # requires optional gap package
        sage: BD.parameters()      # requires optional gap package
        (2, 9, 3, 1)
        sage: BD                   # requires optional gap package
        Incidence structure with 9 points and 12 blocks
        sage: print BD             # requires optional gap package
        WittDesign<points=[0, 1, 2, 3, 4, 5, 6, 7, 8], blocks=[[0, 1, 7], [0, 2, 5], [0, 3, 4], [0, 6, 8], [1, 2, 6], [1, 3, 5], [1, 4, 8], [2, 3, 8], [2, 4, 7], [3, 6, 7], [4, 5, 6], [5, 7, 8]]>
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
    return BlockDesign(v, gB, name="WittDesign", test=True)

def HadamardDesign(n):
    """
    As described in Section 1, p. 10, in [CvL]. The input n must have the
    property that there is a Hadamard matrix of order n+1 (and that a
    construction of that Hadamard matrix has been implemented...).

    EXAMPLES::

        sage: HadamardDesign(7)
        Incidence structure with 7 points and 7 blocks
        sage: print HadamardDesign(7)
        HadamardDesign<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    REFERENCES:

    - [CvL] P. Cameron, J. H. van Lint, Designs, graphs, codes and
      their links, London Math. Soc., 1991.
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

def BlockDesign(max_pt, blks, name=None, test=True):
    """
    Returns an instance of the IncidenceStructure class. Requires each
    B in blks to be contained in range(max_pt). Does not test if the
    result is a block design.

    EXAMPLES::

        sage: BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="Fano plane")
        Incidence structure with 7 points and 7 blocks
        sage: print BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]], name="Fano plane")
        Fano plane<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
    """
    nm = name
    tst = test
    if nm == None and test:
        nm = "BlockDesign"
    BD = BlockDesign_generic( range(max_pt), blks, name=nm, test=tst )
    if not(test):
        return BD
    else:
        pars = BD.parameters()
        if BD.block_design_checker(pars[0],pars[1],pars[2],pars[3]):
            return BD
        else:
            raise TypeError("parameters are not those of a block design.")

BlockDesign_generic = IncidenceStructure
"""
    Possibly in the future there will be methods which apply to
    block designs and not incidence structures. None have been
    implemented yet though. The class name BlockDesign_generic
    is reserved in the name space in case more specialized
    methods are implemented later. In that case, BlockDesign_generic
    should inherit from IncidenceStructure.
"""
