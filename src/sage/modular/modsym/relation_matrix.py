"""
Relation matrices for ambient modular symbols spaces

This file contains functions that are used by the various ambient modular
symbols classes to compute presentations of spaces in terms of generators and
relations, using the standard methods based on Manin symbols.
"""

#*****************************************************************************
#       Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

SPARSE=True

import sage.matrix.matrix_space as matrix_space
import sage.matrix.all
import sage.rings.all as rings
from   sage.misc.search import search
from sage.rings.rational_field import is_RationalField


import sage.misc.misc as misc

import manin_symbols


# S = [0,-1; 1,0]
# T = [0,-1; 1,-1],
# T^2 = [-1, 1, -1, 0]
# I = [-1,0; 0,1]


######################################################################
# The following four functions are used to compute the quotient
# modulo the S, I, and T relations more efficiently that the generic
# code in the relation_matrix file:
#    modS_relations -- compute the S relations.
#    modI_quotient --  compute the I relations.
#    T_relation_matrix -- matrix whose echelon form gives
#                                 the quotient by 3-term T relations.
#    gens_to_basis_matrix -- compute echelon form of 3-term
#                                    relation matrix, and read off each
#                                    generator in terms of basis.
# These four functions are orchestrated in the function
#    compute_presentation
# which is defined below.  See the comment at the beginning
# of that function for an overall description of the algorithm.
######################################################################
def modS_relations(syms):
    """
    Compute quotient of Manin symbols by the S relations.

    Here S is the 2x2 matrix [0, -1; 1, 0].

    INPUT:


    -  ``syms`` - manin_symbols.ManinSymbols


    OUTPUT:


    -  ``rels`` - set of pairs of pairs (j, s), where if
       mod[i] = (j,s), then x_i = s\*x_j (mod S relations)


    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma0
        sage: from sage.modular.modsym.relation_matrix import modS_relations

    ::

        sage: syms = ManinSymbolList_gamma0(2, 4); syms
        Manin Symbol List of weight 4 for Gamma0(2)
        sage: modS_relations(syms)
        set([((3, -1), (4, 1)), ((5, -1), (5, 1)), ((1, 1), (6, 1)), ((0, 1), (7, 1)), ((3, 1), (4, -1)), ((2, 1), (8, 1))])

    ::

        sage: syms = ManinSymbolList_gamma0(7, 2); syms
        Manin Symbol List of weight 2 for Gamma0(7)
        sage: modS_relations(syms)
        set([((3, 1), (4, 1)), ((2, 1), (7, 1)), ((5, 1), (6, 1)), ((0, 1), (1, 1))])

    Next we do an example with Gamma1::

        sage: from sage.modular.modsym.manin_symbols import ManinSymbolList_gamma1
        sage: syms = ManinSymbolList_gamma1(3,2); syms
        Manin Symbol List of weight 2 for Gamma1(3)
        sage: modS_relations(syms)
        set([((3, 1), (6, 1)), ((0, 1), (5, 1)), ((0, 1), (2, 1)), ((3, 1), (4, 1)), ((6, 1), (7, 1)), ((1, 1), (2, 1)), ((1, 1), (5, 1)), ((4, 1), (7, 1))])
    """
    if not isinstance(syms, manin_symbols.ManinSymbolList):
        raise TypeError("syms must be a ManinSymbolList")
    tm = misc.verbose()
    # We will fill in this set with the relations x_i + s*x_j = 0,
    # where the notation is as in _sparse_2term_quotient.
    rels = set()
    for i in xrange(len(syms)):
        j, s = syms.apply_S(i)
        assert j != -1
        if i < j:
            rels.add( ((i,1),(j,s)) )
        else:
            rels.add( ((j,s),(i,1)) )
    misc.verbose("finished creating S relations",tm)
    return rels

def modI_relations(syms, sign):
    """
    Compute quotient of Manin symbols by the I relations.

    INPUT:

    -  ``syms`` - ManinSymbols

    -  ``sign`` - int (either -1, 0, or 1)

    OUTPUT:

    -  ``rels`` - set of pairs of pairs (j, s), where if
       mod[i] = (j,s), then x_i = s\*x_j (mod S relations)

    EXAMPLE::

        sage: L = sage.modular.modsym.manin_symbols.ManinSymbolList_gamma1(4, 3)
        sage: sage.modular.modsym.relation_matrix.modI_relations(L, 1)
        set([((14, 1), (20, 1)), ((0, 1), (0, -1)), ((7, 1), (7, -1)), ((9, 1), (3, -1)), ((3, 1), (9, -1)), ((16, 1), (22, 1)), ((10, 1), (4, -1)), ((1, 1), (1, -1)), ((19, 1), (19, 1)), ((8, 1), (2, -1)), ((12, 1), (12, 1)), ((20, 1), (14, 1)), ((21, 1), (15, 1)), ((5, 1), (11, -1)), ((15, 1), (21, 1)), ((22, 1), (16, 1)), ((6, 1), (6, -1)), ((2, 1), (8, -1)), ((17, 1), (23, 1)), ((4, 1), (10, -1)), ((18, 1), (18, 1)), ((11, 1), (5, -1)), ((23, 1), (17, 1)), ((13, 1), (13, 1))])

    .. warning::

       We quotient by the involution eta((u,v)) = (-u,v), which has
       the opposite sign as the involution in Merel's Springer LNM
       1585 paper! Thus our +1 eigenspace is his -1 eigenspace,
       etc. We do this for consistency with MAGMA.
    """
    tm = misc.verbose()
    # We will fill in this set with the relations x_i - sign*s*x_j = 0,
    # where the notation is as in _sparse_2term_quotient.
    rels = set()
    for i in xrange(len(syms)):
        j, s = syms.apply_I(i)
        assert j != -1
        rels.add( ((i,1),(j,-sign*s)) )
    misc.verbose("finished creating I relations",tm)
    return rels

def T_relation_matrix_wtk_g0(syms, mod, field, sparse):
    r"""
    Compute a matrix whose echelon form gives the quotient by 3-term T
    relations. Despite the name, this is used for all modular symbols spaces
    (including those with character and those for `\Gamma_1` and `\Gamma_H`
    groups), not just `\Gamma_0`.

    INPUT:

    -  ``syms`` - ManinSymbols

    -  ``mod`` - list that gives quotient modulo some two-term relations, i.e.,
       the S relations, and if sign is nonzero, the I relations.

    -  ``field`` - base_ring

    -  ``sparse`` - (True or False) whether to use sparse rather than dense
       linear algebra

    OUTPUT: A sparse matrix whose rows correspond to the reduction of
    the T relations modulo the S and I relations.

    EXAMPLE::

        sage: from sage.modular.modsym.relation_matrix import *
        sage: L = sage.modular.modsym.manin_symbols.ManinSymbolList_gamma_h(GammaH(36, [17,19]), 2)
        sage: modS = sparse_2term_quotient(modS_relations(L), 216, QQ)
        sage: T_relation_matrix_wtk_g0(L, modS, QQ, False)
        72 x 216 dense matrix over Rational Field (use the '.str()' method to see the entries)
        sage: T_relation_matrix_wtk_g0(L, modS, GF(17), True)
        72 x 216 sparse matrix over Finite Field of size 17 (use the '.str()' method to see the entries)
    """
    tm = misc.verbose()
    row = 0
    entries = {}
    already_seen = set()
    w = syms.weight()
    for i in xrange(len(syms)):
        if i in already_seen:
            continue
        iT_plus_iTT = syms.apply_T(i) + syms.apply_TT(i)
        j0, s0 = mod[i]
        v = {j0:s0}
        for j, s in iT_plus_iTT:
            if w==2: already_seen.add(j)
            j0, s0 = mod[j]
            s0 = s*s0
            if j0 in v:
                v[j0] += s0
            else:
                v[j0] = s0
        for j0 in v.keys():
            entries[(row, j0)] = v[j0]
        row += 1

    MAT = matrix_space.MatrixSpace(field, row,
                                len(syms), sparse=True)
    R = MAT(entries)
    if not sparse:
        R = R.dense_matrix()
    misc.verbose("finished (number of rows=%s)"%row, tm)
    return R

def gens_to_basis_matrix(syms, relation_matrix, mod, field, sparse):
    """
    Compute echelon form of 3-term relation matrix, and read off each
    generator in terms of basis.

    INPUT:

    -  ``syms`` - a ManinSymbols object

    -  ``relation_matrix`` - as output by
       ``__compute_T_relation_matrix(self, mod)``

    -  ``mod`` - quotient of modular symbols modulo the
       2-term S (and possibly I) relations

    -  ``field`` - base field

    -  ``sparse`` - (bool): whether or not matrix should be
       sparse

    OUTPUT:

    -  ``matrix`` - a matrix whose ith row expresses the
       Manin symbol generators in terms of a basis of Manin symbols
       (modulo the S, (possibly I,) and T rels) Note that the entries of
       the matrix need not be integers.

    -  ``list`` - integers i, such that the Manin symbols `x_i` are a basis.

    EXAMPLE::

        sage: from sage.modular.modsym.relation_matrix import *
        sage: L = sage.modular.modsym.manin_symbols.ManinSymbolList_gamma1(4, 3)
        sage: modS = sparse_2term_quotient(modS_relations(L), 24, GF(3))
        sage: gens_to_basis_matrix(L, T_relation_matrix_wtk_g0(L, modS, GF(3), 24), modS, GF(3), True)
        (24 x 2 sparse matrix over Finite Field of size 3, [13, 23])
    """
    from sage.matrix.matrix import is_Matrix
    if not is_Matrix(relation_matrix):
        raise TypeError("relation_matrix must be a matrix")
    if not isinstance(mod, list):
        raise TypeError("mod must be a list")

    misc.verbose(str(relation_matrix.parent()))

    try:
        h = relation_matrix.height()
    except AttributeError:
        h = 9999999
    tm = misc.verbose("putting relation matrix in echelon form (height = %s)"%h)
    if h < 10:
        A = relation_matrix.echelon_form(algorithm='multimodular', height_guess=1)
    else:
        A = relation_matrix.echelon_form()
    A.set_immutable()
    tm = misc.verbose('finished echelon', tm)

    tm = misc.verbose("Now creating gens --> basis mapping")

    basis_set = set(A.nonpivots())
    pivots = A.pivots()

    basis_mod2 = set([j for j,c in mod if c != 0])

    basis_set = basis_set.intersection(basis_mod2)
    basis = sorted(basis_set)

    ONE = field(1)

    misc.verbose("done doing setup",tm)


    tm = misc.verbose("now forming quotient matrix")
    M = matrix_space.MatrixSpace(field, len(syms), len(basis), sparse=sparse)

    B = M(0)
    cols_index = dict([(basis[i], i) for i in range(len(basis))])

    for i in basis_mod2:
        t, l = search(basis, i)
        if t:
            B[i,l] = ONE
        else:
            _, r = search(pivots, i)    # so pivots[r] = i
            # Set row i to -(row r of A), but where we only take
            # the non-pivot columns of A:
            B._set_row_to_negative_of_row_of_A_using_subset_of_columns(i, A, r, basis, cols_index)

    misc.verbose("done making quotient matrix",tm)

    # The following is very fast (over Q at least).
    tm = misc.verbose('now filling in the rest of the matrix')
    k = 0
    for i in range(len(mod)):
        j, s = mod[i]
        if j != i and s != 0:   # ignored in the above matrix
            k += 1
            B.set_row_to_multiple_of_row(i, j, s)
    misc.verbose("set %s rows"%k)
    tm = misc.verbose("time to fill in rest of matrix", tm)

    return B, basis

def compute_presentation(syms, sign, field, sparse=None):
    r"""
    Compute the presentation for self, as a quotient of Manin symbols
    modulo relations.

    INPUT:

    -  ``syms`` - manin_symbols.ManinSymbols

    -  ``sign`` - integer (-1, 0, 1)

    -  ``field`` - a field


    OUTPUT:

    -  sparse matrix whose rows give each generator
       in terms of a basis for the quotient

    -  list of integers that give the basis for the
       quotient

    -  mod: list where mod[i]=(j,s) means that x_i
       = s\*x_j modulo the 2-term S (and possibly I) relations.


    ALGORITHM:

    #. Let `S = [0,-1; 1,0], T = [0,-1; 1,-1]`, and
       `I = [-1,0; 0,1]`.

    #. Let `x_0,\ldots, x_{n-1}` by a list of all
       non-equivalent Manin symbols.

    #. Form quotient by 2-term S and (possibly) I relations.

    #. Create a sparse matrix `A` with `m` columns,
       whose rows encode the relations

       .. math::

                          [x_i] + [x_i T] + [x_i T^2] = 0.


       There are about n such rows. The number of nonzero entries per row
       is at most 3\*(k-1). Note that we must include rows for *all* i,
       since even if `[x_i] = [x_j]`, it need not be the case
       that `[x_i T] = [x_j T]`, since `S` and
       `T` do not commute. However, in many cases we have an a
       priori formula for the dimension of the quotient by all these
       relations, so we can omit many relations and just check that there
       are enough at the end--if there aren't, we add in more.

    #. Compute the reduced row echelon form of `A` using sparse
       Gaussian elimination.

    #. Use what we've done above to read off a sparse matrix R that
       uniquely expresses each of the n Manin symbols in terms of a subset
       of Manin symbols, modulo the relations. This subset of Manin
       symbols is a basis for the quotient by the relations.


    EXAMPLE::

        sage: L = sage.modular.modsym.manin_symbols.ManinSymbolList_gamma0(8,2)
        sage: sage.modular.modsym.relation_matrix.compute_presentation(L, 1, GF(9,'a'), True)
        (
        [2 0 0]
        [1 0 0]
        [0 0 0]
        [0 2 0]
        [0 0 0]
        [0 0 2]
        [0 0 0]
        [0 2 0]
        [0 0 0]
        [0 1 0]
        [0 1 0]
        [0 0 1], [1, 9, 11], [(1, 2), (1, 1), (0, 0), (9, 2), (0, 0), (11, 2), (0, 0), (9, 2), (0, 0), (9, 1), (9, 1), (11, 1)]
        )
    """
    if sparse is None:
        if syms.weight() >= 6:
            sparse = False
        else:
            sparse = True
    R, mod = relation_matrix_wtk_g0(syms, sign, field, sparse)
    B, basis = gens_to_basis_matrix(syms, R, mod, field, sparse)
    return B, basis, mod

def relation_matrix_wtk_g0(syms, sign, field, sparse):
    r"""
    Compute the matrix of relations. Despite the name, this is used for all
    spaces (not just for Gamma0). For a description of the algorithm, see the
    docstring for ``compute_presentation``.

    INPUT:

    - ``syms``: sage.modular.modsym.manin_symbols.ManinSymbolList object

    - ``sign``: integer (0, 1 or -1)

    - ``field``: the base field (non-field base rings not supported at present)

    - ``sparse``: (True or False) whether to use sparse arithmetic.

    Note that ManinSymbolList objects already have a specific weight, so there
    is no need for an extra ``weight`` parameter.

    OUTPUT: a pair (R, mod) where

    - R is a matrix as output by ``T_relation_matrix_wtk_g0``

    - mod is a set of 2-term relations as output by ``sparse_2term_quotient``

    EXAMPLE::

        sage: L =  sage.modular.modsym.manin_symbols.ManinSymbolList_gamma0(8,2)
        sage: A = sage.modular.modsym.relation_matrix.relation_matrix_wtk_g0(L, 0, GF(2), True); A
        (
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 1 1 1 0]
        [0 0 0 0 0 0 1 0 0 1 1 0]
        [0 0 0 0 0 0 1 0 0 0 0 0], [(1, 1), (1, 1), (8, 1), (10, 1), (6, 1), (11, 1), (6, 1), (9, 1), (8, 1), (9, 1), (10, 1), (11, 1)]
        )
        sage: A[0].is_sparse()
        True
    """
    rels = modS_relations(syms)
    if sign != 0:
        # Let rels = rels union I relations.
        rels.update(modI_relations(syms,sign))

    if syms._apply_S_only_0pm1() and is_RationalField(field):
        import relation_matrix_pyx
        mod = relation_matrix_pyx.sparse_2term_quotient_only_pm1(rels, len(syms))
    else:
        mod = sparse_2term_quotient(rels, len(syms), field)

    R = T_relation_matrix_wtk_g0(syms, mod, field, sparse)
    return R, mod

def sparse_2term_quotient(rels, n, F):
    r"""
    Performs Sparse Gauss elimination on a matrix all of whose columns
    have at most 2 nonzero entries. We use an obvious algorithm, which
    runs fast enough. (Typically making the list of relations takes
    more time than computing this quotient.) This algorithm is more
    subtle than just "identify symbols in pairs", since complicated
    relations can cause generators to surprisingly equal 0.

    INPUT:


    -  ``rels`` - set of pairs ((i,s), (j,t)). The pair
       represents the relation s\*x_i + t\*x_j = 0, where the i, j must
       be Python int's.

    -  ``n`` - int, the x_i are x_0, ..., x_n-1.

    -  ``F`` - base field


    OUTPUT:


    -  ``mod`` - list such that mod[i] = (j,s), which means
       that x_i is equivalent to s\*x_j, where the x_j are a basis for
       the quotient.


    EXAMPLE: We quotient out by the relations

    .. math::

                    3*x0 - x1 = 0,\qquad  x1 + x3 = 0,\qquad   x2 + x3 = 0,\qquad  x4 - x5 = 0


    to get

    ::

        sage: v = [((int(0),3), (int(1),-1)), ((int(1),1), (int(3),1)), ((int(2),1),(int(3),1)), ((int(4),1),(int(5),-1))]
        sage: rels = set(v)
        sage: n = 6
        sage: from sage.modular.modsym.relation_matrix import sparse_2term_quotient
        sage: sparse_2term_quotient(rels, n, QQ)
        [(3, -1/3), (3, -1), (3, -1), (3, 1), (5, 1), (5, 1)]
    """
    if not isinstance(rels, set):
        raise TypeError("rels must be a set")
    n = int(n)
    if not isinstance(F, rings.Ring):
        raise TypeError("F must be a ring.")

    tm = misc.verbose("Starting sparse 2-term quotient...")
    free = range(n)
    ONE = F(1)
    ZERO = F(0)
    coef = [ONE for i in xrange(n)]
    related_to_me = [[] for i in xrange(n)]
    for v0, v1 in rels:
        c0 = coef[v0[0]] * F(v0[1])
        c1 = coef[v1[0]] * F(v1[1])

        # Mod out by the following relation:
        #
        #    c0*free[v0[0]] + c1*free[v1[0]] = 0.
        #
        die = None
        if c0 == ZERO and c1 == ZERO:
            pass
        elif c0 == ZERO and c1 != ZERO:  # free[v1[0]] --> 0
            die = free[v1[0]]
        elif c1 == ZERO and c0 != ZERO:
            die = free[v0[0]]
        elif free[v0[0]] == free[v1[0]]:
            if c0 + c1 != 0:
                # all xi equal to free[v0[0]] must now equal to zero.
                die = free[v0[0]]
        else:  # x1 = -c1/c0 * x2.
            x = free[v0[0]]
            free[x] = free[v1[0]]
            coef[x] = -c1/c0
            for i in related_to_me[x]:
                free[i] = free[x]
                coef[i] *= coef[x]
                related_to_me[free[v1[0]]].append(i)
            related_to_me[free[v1[0]]].append(x)
        if die is not None:
            for i in related_to_me[die]:
                free[i] = 0
                coef[i] = ZERO
            free[die] = 0
            coef[die] = ZERO

    mod = [(free[i], coef[i]) for i in xrange(len(free))]
    misc.verbose("finished",tm)
    return mod



#############################################################
## The following two sparse_relation_matrix are not
## used by any modular symbols code.  They're here for
## historical reasons, and can probably be safely deleted.
#############################################################

## def sparse_relation_matrix_wt2_g0n(list, field, sign=0):
##     r"""
##     Create the sparse relation matrix over $\Q$ for Manin symbols of
##     weight 2 on $\Gamma_0(N)$, with given sign.

##     INPUT:
##         list -- sage.modular.modsym.p1list.List
##     OUTPUT:
##         A -- a sparse matrix that gives the 2-term and 3-term
##              relations between Manin symbols.

##     MORE DETAILS:
##     \begin{enumerate}
##       \item Create an empty sparse matrix.

##       \item Let $S = [0,-1; 1,0]$, $T = [0,-1; 1,-1]$, $I = [-1,0; 0,1]$.

##       \item Enter the T relations:
##            $$
##                    x + x T = 0.
##            $$
##            Remove x and x*T from reps to consider.

##       \item If sign $\neq 0$, enter the I relations:
##            $$
##                    x - sign\cdot x\cdot I = 0.
##            $$

##       \item Enter the S relations in the matrix:
##            $$
##                    x + x S + x S^2 = 0
##            $$
##            by putting 1s at cols corresponding to $x$, $x S$, and $x S^2$.
##            Remove $x$, $x S$, and $x S^2$ from list of reps to consider.
##     \end{enumerate}
##     """
##     ZERO = field(0)
##     ONE =  field(1)
##     TWO =  field(2)

##     # This will be a dict of the entries of the sparse matrix, where
##     # the notation is entries[(i,j)]=x.
##     entries = {}

##     # The current row
##     row = 0

##     ##  The S relations
##     already_seen= set([])
##     for i in range(len(list)):
##         if i in already_seen:
##             continue
##         u,v = list[i]
##         j = list.index(v,-u)
##         already_seen.add(j)
##         if i != j:
##             entries[(row,i)] = ONE
##             entries[(row,j)] = ONE
##         else:
##             entries[(row,i)] = TWO
##         row += 1
##     number_of_S_relations = row
##     misc.verbose("There were %s S relations"%(number_of_S_relations))

##     ##  The eta relations:
##     ##    eta((u,v)) = -(-u,v)
##     if sign != 0:
##         SIGN = field(sign)
##         already_seen= set([])
##         for i in range(len(list)):
##             if i in already_seen:
##                 continue
##             u, v = list[i]
##             j = list.index(-u,v)
##             already_seen.add(j)
##             if i != j:
##                 entries[(row,i)] = ONE
##                 entries[(row,j)] = SIGN*ONE
##             else:
##                 entries[(row,i)] = ONE + SIGN
##             row += 1
##     number_of_I_relations = row - number_of_S_relations
##     misc.verbose("There were %s I relations"%(number_of_I_relations))

##     ## The three-term T relations
##     already_seen = set([])
##     for i in range(len(list)):
##         if i in already_seen:
##             continue
##         u,v = list[i]
##         j1 = list.index(v,-u-v)
##         already_seen.add(j1)
##         j2 = list.index(-u-v,u)
##         already_seen.add(j2)
##         v = {i:ZERO, j1:ZERO, j2:ZERO}
##         v[i] = ONE
##         v[j1] += ONE
##         v[j2] += ONE
##         for x in v.keys():
##             entries[(row,x)] = v[x]
##         row += 1

##     number_of_T_relations = row - number_of_I_relations - number_of_S_relations
##     misc.verbose("There were %s T relations"%(number_of_T_relations))

##     M = matrix_space.MatrixSpace(RationalField(), row,
##                     len(list), sparse=True)
##     if not sparse:
##         M = M.dense_matrix()

##     return M(entries)

## def sparse_relation_matrix_wtk_g0n(M, field, sign=0):
##     r"""
##     Create the sparse relation matrix over $\Q$ for Manin symbols of
##     given weight on $\Gamma_0(N)$, with given sign.

##     INPUT:
##         M -- manin_symbols.ManinSymbolList
##         field -- base field
##         weight -- the weight, an integer > 2
##         sign -- element of [-1,0,1]

##     OUTPUT:
##         A -- a SparseMatrix that gives the 2-term and 3-term relations
##              between Manin symbols.

##     MORE DETAILS:
##     \begin{enumerate}
##        \item Create an empty sparse matrix.

##         \item Let $S = [0,-1; 1,0]$, $T = [0,-1; 1,-1]$, $I = [-1,0; 0,1]$.

##         \item Enter the $T$ relations:
##                  $$  x + x*T = 0  $$
##            Remove $x$ and $x T$ from reps to consider.

##         \item If sign $\neq 0$, enter the I relations:
##         $$
##                    x + sign x I = 0.
##         $$

##         \item Enter the $S$ relations in the matrix:
##            $$
##                    x + x S + x S^2 = 0
##            $$
##            by putting 1's at cols corresponding to $x$, $x S$, and $x S^2$.
##            Remove x from list of reps to consider.
##     \end{enumerate}
##     """
##     weight = M.weight()
##     if not (isinstance(weight, int) and weight > 2):
##         raise TypeError, "weight must be an int > 2"

##     ZERO = field(0)
##     ONE =  field(1)
##     TWO =  field(2)

##     # This will be a dict of the entries of the sparse matrix, where
##     # the notation is entries[(i,j)]=x.
##     entries = {}

##     # The current row
##     row = 0

##     # The list of Manin symbol triples (i,u,v)
##     n = len(M)

##     ##  The S relations
##     already_seen= set([])
##     for i in xrange(n):
##         if i in already_seen:
##             continue
##         j, s = M.apply_S(i)
##         already_seen.add(j)
##         if i != j:
##             entries[(row,i)] = ONE
##             entries[(row,j)] = field(s)
##         else:
##             entries[(row,i)] = ONE+field(s)
##         row += 1
##     number_of_S_relations = row
##     misc.verbose("There were %s S relations"%(number_of_S_relations))
##     cnt = row
##     ##  The I relations
##     if sign != 0:
##         SIGN = field(sign)
##         already_seen= set([])
##         for i in xrange(n):
##             if i in already_seen:
##                 continue
##             j, s = M.apply_I(i)
##             already_seen.add(j)
##             if i != j:
##                 entries[(row,i)] = ONE
##                 entries[(row,j)] = -SIGN*field(s)
##             else:
##                 entries[(row,i)] = ONE-SIGN*field(s)
##             row += 1
##     number_of_I_relations = row - number_of_S_relations
##     misc.verbose("There were %s I relations"%(number_of_I_relations))
##     cnt = row

##     ## The T relations
##     already_seen = set([])
##     for i in xrange(n):
##         if i in already_seen:
##             continue
##         iT_plus_iTT = M.apply_T(i) + M.apply_TT(i)
##         v = {i:ONE}
##         for j, s in iT_plus_iTT:
##             if j in v:
##                 v[j] += field(s)
##             else:
##                 v[j] = field(s)
##         for j in v.keys():
##             entries[(row, j)] = v[j]
##         row += 1

