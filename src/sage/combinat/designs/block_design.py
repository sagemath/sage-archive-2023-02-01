"""
Block designs.

A module to help with constructions and computations of block
designs and other incidence structures.

A block design is an incidence structure consisting of a set of points `P` and a
set of blocks `B`, where each block is considered as a subset of `P`. More
precisely, a *block design* `B` is a class of `k`-element subsets of `P` such
that the number `r` of blocks that contain any point `x` in `P` is independent
of `x`, and the number `\lambda` of blocks that contain any given `t`-element
subset `T` is independent of the choice of `T` (see [1]_ for more). Such a block
design is also called a `t-(v,k,\lambda)`-design, and `v` (the number of
points), `b` (the number of blocks), `k`, `r`, and `\lambda` are the parameters
of the design. (In Python, ``lambda`` is reserved, so we sometimes use ``lmbda``
or ``L`` instead.)

In Sage, sets are replaced by (ordered) lists and the standard representation of
a block design uses `P = [0,1,..., v-1]`, so a block design is specified by
`(v,B)`.

REFERENCES:

.. [1] Block design from wikipedia,
  :wikipedia:`Block_design`

.. [2] What is a block design?,
  http://designtheory.org/library/extrep/extrep-1.1-html/node4.html (in 'The
  External Representation of Block Designs' by Peter J. Cameron, Peter
  Dobcsanyi, John P. Morgan, Leonard H. Soicher)

AUTHORS:

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of the module block_design.py (version
  0.6) written by Peter Dobcsanyi peter@designtheory.org. Thanks go to Robert
  Miller for lots of good design suggestions.


Functions and methods
---------------------
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

from sage.modules.free_module import VectorSpace
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial, integer_floor
from sage.combinat.designs.incidence_structures import IncidenceStructure, IncidenceStructureFromMatrix
from sage.misc.decorators import rename_keyword
from sage.rings.finite_rings.constructor import FiniteField

###  utility functions  -------------------------------------------------------

def tdesign_params(t, v, k, L):
    """
    Return the design's parameters: `(t, v, b, r , k, L)`. Note that `t` must be
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

def ProjectiveGeometryDesign(n, d, F, algorithm=None):
    """
    Returns a projective geometry design.

    A projective geometry design of parameters `n,d,F` has for points the lines
    of `F^{n+1}`, and for blocks the `d+1`-dimensional subspaces of `F^{n+1}`,
    each of which contains `\\frac {|F|^{d+1}-1} {|F|-1}` lines.

    INPUT:

    - ``n`` is the projective dimension

    - ``d`` is the dimension of the subspaces of `P = PPn(F)` which
      make up the blocks.

    - ``F`` is a finite field.

    - ``algorithm`` -- set to ``None`` by default, which results in using Sage's
      own implementation. In order to use GAP's implementation instead (i.e. its
      ``PGPointFlatBlockDesign`` function) set ``algorithm="gap"``. Note that
      GAP's "design" package must be available in this case.

    EXAMPLES:

    The points of the following design are the `\\frac {2^{2+1}-1} {2-1}=7`
    lines of `\mathbb{Z}_2^{2+1}`. It has `7` blocks, corresponding to each
    2-dimensional subspace of `\mathbb{Z}_2^{2+1}`::

        sage: designs.ProjectiveGeometryDesign(2, 1, GF(2))
        Incidence structure with 7 points and 7 blocks
        sage: BD = designs.ProjectiveGeometryDesign(2, 1, GF(2), algorithm="gap") # optional - gap_packages (design package)
        sage: BD.is_block_design()                                     # optional - gap_packages (design package)
        (True, [2, 7, 3, 1])
    """
    q = F.order()
    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    if algorithm is None:
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
    if algorithm == "gap":   # Requires GAP's Design
        gap.load_package("design")
        gap.eval("D := PGPointFlatBlockDesign( %s, %s, %d )"%(n,q,d))
        v = eval(gap.eval("D.v"))
        gblcks = eval(gap.eval("D.blocks"))
        gB = []
        for b in gblcks:
            gB.append([x-1 for x in b])
        return BlockDesign(v, gB, name="ProjectiveGeometryDesign")

def ProjectivePlaneDesign(n, type="Desarguesian"):
    r"""
    Returns a projective plane of order `n`.

    A finite projective plane is a 2-design with `n^2+n+1` lines (or blocks) and
    `n^2+n+1` points. For more information on finite projective planes, see the
    :wikipedia:`Projective_plane#Finite_projective_planes`.

    INPUT:

    - ``n`` -- the finite projective plane's order

    - ``type`` -- When set to ``"Desarguesian"``, the method returns
      Desarguesian projective planes, i.e. a finite projective plane obtained by
      considering the 1- and 2- dimensional spaces of `F_n^3`.

      For the moment, no other value is available for this parameter.

    .. SEEALSO::

        :meth:`ProjectiveGeometryDesign`

    EXAMPLES::

        sage: designs.ProjectivePlaneDesign(2)
        Incidence structure with 7 points and 7 blocks

    Non-existent ones::

        sage: designs.ProjectivePlaneDesign(10)
        Traceback (most recent call last):
        ...
        ValueError: No projective plane design of order 10 exists.
        sage: designs.ProjectivePlaneDesign(14)
        Traceback (most recent call last):
        ...
        ValueError: By the Bruck-Ryser-Chowla theorem, no projective plane of order 14 exists.

    An unknown one::

        sage: designs.ProjectivePlaneDesign(12)
        Traceback (most recent call last):
        ...
        ValueError: If such a projective plane exists, we do not know how to build it.

    TESTS::

        sage: designs.ProjectivePlaneDesign(10, type="AnyThingElse")
        Traceback (most recent call last):
        ...
        ValueError: The value of 'type' must be 'Desarguesian'.
    """
    from sage.rings.arith import two_squares

    if type != "Desarguesian":
        raise ValueError("The value of 'type' must be 'Desarguesian'.")

    try:
        F = FiniteField(n, 'x')
    except ValueError:
        if n == 10:
            raise ValueError("No projective plane design of order 10 exists.")
        try:
            if (n%4) in [1,2]:
                two_squares(n)
        except ValueError:
            raise ValueError("By the Bruck-Ryser-Chowla theorem, no projective"
                             " plane of order "+str(n)+" exists.")

        raise ValueError("If such a projective plane exists, "
                         "we do not know how to build it.")

    return ProjectiveGeometryDesign(2,1,F)

def AffineGeometryDesign(n, d, F):
    r"""
    Returns an Affine Geometry Design.

    INPUT:

    - `n` (integer) -- the Euclidean dimension. The number of points is
      `v=|F^n|`.

    - `d` (integer) -- the dimension of the (affine) subspaces of `P = GF(q)^n`
      which make up the blocks.

    - `F` -- a Finite Field (i.e. ``FiniteField(17)``), or a prime power
      (i.e. an integer)

    `AG_{n,d} (F)`, as it is sometimes denoted, is a `2` - `(v, k, \lambda)`
    design of points and `d`- flats (cosets of dimension `n`) in the affine
    geometry `AG_n (F)`, where

    .. math::

             v = q^n,\  k = q^d ,
             \lambda =\frac{(q^{n-1}-1) \cdots (q^{n+1-d}-1)}{(q^{n-1}-1) \cdots (q-1)}.

    Wraps some functions used in GAP Design's ``PGPointFlatBlockDesign``.  Does
    *not* require GAP's Design package.

    EXAMPLES::

        sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
        sage: BD.parameters(t=2)
        (2, 8, 2, 1)
        sage: BD.is_block_design()
        (True, [2, 8, 2, 1])
        sage: BD = designs.AffineGeometryDesign(3, 2, GF(2))
        sage: BD.parameters(t=3)
        (3, 8, 4, 1)
        sage: BD.is_block_design()
        (True, [3, 8, 4, 1])

    A 3-design::

        sage: D = IncidenceStructure(range(32),designs.steiner_quadruple_system(32))
        sage: D.is_block_design()
        (True, [3, 32, 4, 1])

    With an integer instead of a Finite Field::

        sage: BD = designs.AffineGeometryDesign(3, 2, 4)
        sage: BD.parameters(t=2)
        (2, 64, 16, 5)
    """
    try:
        q = int(F)
    except TypeError:
        q = F.order()

    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    gap.eval("V:=GaloisField(%s)^%s"%(q,n))
    gap.eval("points:=AsSet(V)")
    gap.eval("Subs:=AsSet(Subspaces(V,%s));"%d)
    gap.eval("CP:=Cartesian(points,Subs)")
    flats = gap.eval("flats:=List(CP,x->Sum(x))") # affine spaces
    gblcks = eval(gap.eval("Set(List(flats,f->Filtered([1..Length(points)],i->points[i] in f)));"))
    v = q**n
    gB = []
    for b in gblcks:
       gB.append([x-1 for x in b])
    return BlockDesign(v, gB, name="AffineGeometryDesign")

def WittDesign(n):
    """
    INPUT:

    - ``n`` is in `9,10,11,12,21,22,23,24`.

    Wraps GAP Design's WittDesign. If ``n=24`` then this function returns the
    large Witt design `W_{24}`, the unique (up to isomorphism) `5-(24,8,1)`
    design. If ``n=12`` then this function returns the small Witt design
    `W_{12}`, the unique (up to isomorphism) `5-(12,6,1)` design.  The other
    values of `n` return a block design derived from these.

    .. NOTE:

        Requires GAP's Design package (included in the gap_packages Sage spkg).

    EXAMPLES::

        sage: BD = designs.WittDesign(9)   # optional - gap_packages (design package)
        sage: BD.is_block_design()      # optional - gap_packages (design package)
        (2, 9, 3, 1)
        sage: BD                   # optional - gap_packages (design package)
        Incidence structure with 9 points and 12 blocks
        sage: print BD             # optional - gap_packages (design package)
        WittDesign<points=[0, 1, 2, 3, 4, 5, 6, 7, 8], blocks=[[0, 1, 7], [0, 2, 5], [0, 3, 4], [0, 6, 8], [1, 2, 6], [1, 3, 5], [1, 4, 8], [2, 3, 8], [2, 4, 7], [3, 6, 7], [4, 5, 6], [5, 7, 8]]>
        sage: BD = designs.WittDesign(12)  # optional - gap_packages (design package)
        sage: BD.is_block_design()   # optional - gap_packages (design package)
        (5, 12, 6, 1)
    """
    from sage.interfaces.gap import gap, GapElement
    gap.load_package("design")
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
    property that there is a Hadamard matrix of order `n+1` (and that a
    construction of that Hadamard matrix has been implemented...).

    EXAMPLES::

        sage: designs.HadamardDesign(7)
        Incidence structure with 7 points and 7 blocks
        sage: print designs.HadamardDesign(7)
        HadamardDesign<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    REFERENCES:

    - [CvL] P. Cameron, J. H. van Lint, Designs, graphs, codes and
      their links, London Math. Soc., 1991.
    """
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix
    from sage.matrix.constructor import matrix
    H = hadamard_matrix(n+1)
    H1 = H.matrix_from_columns(range(1,n+1))
    H2 = H1.matrix_from_rows(range(1,n+1))
    J = matrix(ZZ,n,n,[1]*n*n)
    MS = J.parent()
    A = MS((H2+J)/2) # convert -1's to 0's; coerce entries to ZZ
    # A is the incidence matrix of the block design
    return IncidenceStructureFromMatrix(A,name="HadamardDesign")

def steiner_triple_system(n):
    r"""
    Returns a Steiner Triple System

    A Steiner Triple System (STS) of a set `\{0,...,n-1\}`
    is a family `S` of 3-sets such that for any `i \not = j`
    there exists exactly one set of `S` in which they are
    both contained.

    It can alternatively be thought of as a factorization of
    the complete graph `K_n` with triangles.

    A Steiner Triple System of a `n`-set exists if and only if
    `n \equiv 1 \pmod 6` or `n \equiv 3 \pmod 6`, in which case
    one can be found through Bose's and Skolem's constructions,
    respectively [AndHon97]_.

    INPUT:

    - ``n`` returns a Steiner Triple System of `\{0,...,n-1\}`

    EXAMPLE:

    A Steiner Triple System on `9` elements ::

        sage: sts = designs.steiner_triple_system(9)
        sage: sts
        Incidence structure with 9 points and 12 blocks
        sage: list(sts)
        [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7], [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]

    As any pair of vertices is covered once, its parameters are ::

        sage: sts.parameters(t=2)
        (2, 9, 3, 1)

    An exception is raised for invalid values of ``n`` ::

        sage: designs.steiner_triple_system(10)
        Traceback (most recent call last):
        ...
        ValueError: Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6

    REFERENCE:

    .. [AndHon97] A short course in Combinatorial Designs,
      Ian Anderson, Iiro Honkala,
      Internet Editions, Spring 1997,
      http://www.utu.fi/~honkala/designs.ps
    """

    name = "Steiner Triple System on "+str(n)+" elements"

    if n%6 == 3:
        t = (n-3)/6
        Z = range(2*t+1)

        T = lambda (x,y) : x + (2*t+1)*y

        sts = [[(i,0),(i,1),(i,2)] for i in Z] + \
            [[(i,k),(j,k),(((t+1)*(i+j)) % (2*t+1),(k+1)%3)] for k in range(3) for i in Z for j in Z if i != j]

    elif n%6 == 1:

        t = (n-1)/6
        N = range(2*t)
        T = lambda (x,y) : x+y*t*2 if (x,y) != (-1,-1) else n-1

        L1 = lambda i,j : (i+j) % (int((n-1)/3))
        L = lambda i,j : L1(i,j)/2 if L1(i,j)%2 == 0 else t+(L1(i,j)-1)/2

        sts = [[(i,0),(i,1),(i,2)] for i in range(t)] + \
            [[(-1,-1),(i,k),(i-t,(k+1) % 3)] for i in range(t,2*t) for k in [0,1,2]] + \
            [[(i,k),(j,k),(L(i,j),(k+1) % 3)] for k in [0,1,2] for i in N for j in N if i < j]

    else:
        raise ValueError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")

    from sage.sets.set import Set
    sts = Set(map(lambda x: Set(map(T,x)),sts))

    return BlockDesign(n, sts, name=name)

def BlockDesign(max_pt, blks, name=None, test=True):
    """
    Returns an instance of the :class:`IncidenceStructure` class.

    Requires each B in blks to be contained in range(max_pt). Does not test if
    the result is a block design.

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
        pars = BD.parameters(t=2)
        if BD.block_design_checker(pars[0],pars[1],pars[2],pars[3]):
            return BD
        else:
            raise TypeError("parameters are not those of a block design.")

# Possibly in the future there will be methods which apply to block designs and
# not incidence structures. None have been implemented yet though. The class
# name BlockDesign_generic is reserved in the name space in case more
# specialized methods are implemented later. In that case, BlockDesign_generic
# should inherit from IncidenceStructure.
BlockDesign_generic = IncidenceStructure
