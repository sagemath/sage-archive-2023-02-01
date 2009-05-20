r"""
Guava error-correcting code constructions.

This module only contains Guava wrappers (Guava is an optional GAP package).

AUTHOR:
    -- David Joyner (2005-11-22, 2006-12-03): initial version
    -- Nick Alexander (2006-12-10): factor GUAVA code to guava.py
    -- David Joyner (2007-05): removed Golay codes, toric and trivial
                               codes and placed them in code_constructions;
                               renamed RandomLinearCode->RandomLinearCodeGuava
    -- David Joyner (2008-03): removed QR, XQR, cyclic and ReedSolomon codes
    --  " (2009-05): added "optional package" comments, fixed some docstrings to
                     to be sphinx compatible

"""

#*****************************************************************************
#       Copyright (C) 2007 David Joyner <wdj@usna.edu>
#                     2006 Nick Alexander <ncalexan@math.uci.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.interfaces.all import gap
from sage.misc.randstate import current_randstate
from sage.misc.preparser import *
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_field import FiniteField as GF
from sage.interfaces.gap import gfq_gap_to_sage
from sage.groups.perm_gps.permgroup import *
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod, add
from sage.misc.functional import log
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from linear_code import *

def BinaryReedMullerCode(r,k):
    r"""
    The binary 'Reed-Muller code' with dimension k and
    order r is a code with length `2^k` and minimum distance `2^k-r`
    (see for example, section 1.10 in [HP]_). By definition, the
    `r^{th}` order binary Reed-Muller code of length `n=2^m`, for
    `0 \leq r \leq m`, is the set of all vectors `(f(p)\ |\ p \\in GF(2)^m)`,
    where `f` is a multivariate polynomial of degree at most `r`
    in `m` variables.

    INPUT:
        r, k -- positive integers with `2^k>r`.

    OUTPUT:
        Returns the binary 'Reed-Muller code' with dimension k and order r.

    EXAMPLE::
        sage: C = BinaryReedMullerCode(2,4); C  # requires optional package
        Linear code of length 16, dimension 11 over Finite Field of size 2
        sage: C.minimum_distance()              # requires optional package
        4
        sage: C.gen_mat()                       # requires optional package
        [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
        [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]
        [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]
        [0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1]
        [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1]
        [0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1]
        [0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1]
        [0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1]
        [0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1]
        [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1]

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=ReedMullerCode("+str(r)+", "+str(k)+")")
    gap.eval("G:=GeneratorMat(C)")
    k = int(gap.eval("Length(G)"))
    n = int(gap.eval("Length(G[1])"))
    G = [[gfq_gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def QuasiQuadraticResidueCode(p):
    r"""
    A (binary) quasi-quadratic residue code (or QQR code), as defined by
    Proposition 2.2 in [BM]_, has a generator matrix in the block form `G=(Q,N)`.
    Here `Q` is a `p \times p` circulant matrix whose top row
    is `(0,x_1,...,x_{p-1})`, where `x_i=1` if and only if `i`
    is a quadratic residue `\mod p`, and `N` is a `p \times p` circulant
    matrix whose top row is `(0,y_1,...,y_{p-1})`, where `x_i+y_i=1` for all `i`.

    INPUT:
        p -- a prime >2.

    OUTPUT:
        Returns a QQR code of length 2p.

    EXAMPLES::
        sage: C = QuasiQuadraticResidueCode(11); C   # requires optional package
        Linear code of length 22, dimension 11 over Finite Field of size 2

    REFERENCES:
    ..[BM] Bazzi and Mitter, {\it Some constructions of codes from group actions}, (preprint
             March 2003, available on Mitter's MIT website).
    ..[J]  D. Joyner, {\it On quadratic residue codes and hyperelliptic curves},
             (preprint 2006)

    These are self-orthogonal in general and self-dual when $p \\equiv 3 \\pmod 4$.

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=QQRCode("+str(p)+")")
    gap.eval("G:=GeneratorMat(C)")
    k = int(gap.eval("Length(G)"))
    n = int(gap.eval("Length(G[1])"))
    G = [[gfq_gap_to_sage(gap.eval("G[%s][%s]" % (i,j)),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))



def RandomLinearCodeGuava(n,k,F):
    r"""
    The method used is to first construct a `k \times n` matrix of the block
    form `(I,A)`, where `I` is a `k \times k` identity matrix and `A` is a
    `k \times (n-k)` matrix constructed using random elements of `F`. Then
    the columns are permuted using a randomly selected element of the symmetric
    group `S_n`.

    INPUT:
        Integers `n,k`, with `n>k>1`.

    OUTPUT:
        Returns a "random" linear code with length n, dimension k over field F.

    EXAMPLES::
        sage: C = RandomLinearCodeGuava(30,15,GF(2)); C      # requires optional package
        Linear code of length 30, dimension 15 over Finite Field of size 2
        sage: C = RandomLinearCodeGuava(10,5,GF(4,'a')); C      # requires optional package
        Linear code of length 10, dimension 5 over Finite Field in a of size 2^2

    AUTHOR: David Joyner (11-2005)
    """
    current_randstate().set_seed_gap()

    q = F.order()
    gap.eval("C:=RandomLinearCode("+str(n)+","+str(k)+", GF("+str(q)+"))")
    gap.eval("G:=GeneratorMat(C)")
    k = int(gap.eval("Length(G)"))
    n = int(gap.eval("Length(G[1])"))
    G = [[gfq_gap_to_sage(gap.eval("G[%s][%s]" % (i,j)),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))


