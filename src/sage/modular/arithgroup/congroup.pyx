"""
Cython helper functions for congruence subgroups

This file contains optimized Cython implementations of a few functions related
to the standard congruence subgroups `\Gamma_0, \Gamma_1, \Gamma_H`.  These
functions are for internal use by routines elsewhere in the Sage library.
"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

include "cysignals/memory.pxi"

import random
from congroup_gamma1 import Gamma1_constructor as Gamma1
from congroup_gamma0 import Gamma0_constructor as Gamma0

cimport sage.rings.fast_arith
import sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int arith_int
arith_int = sage.rings.fast_arith.arith_int()
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ
Mat2Z = MatrixSpace(ZZ,2)

cdef Matrix_integer_dense genS, genT, genI
genS = Matrix_integer_dense(Mat2Z, [0,-1, 1, 0], True, True)
genT = Matrix_integer_dense(Mat2Z, [1, 1, 0, 1], True, True)
genI = Matrix_integer_dense(Mat2Z, [1, 0, 0, 1], True, True)


# This is the C version of a function formerly implemented in python in
# sage.modular.congroup.  It is orders of magnitude faster (e.g., 30
# times).  The key speedup is in replacing looping through the
# elements of the Python list R with looping through the elements of a
# C-array.

def degeneracy_coset_representatives_gamma0(int N, int M, int t):
    r"""
    Let `N` be a positive integer and `M` a divisor of `N`.  Let `t` be a
    divisor of `N/M`, and let `T` be the `2 \times 2` matrix `(1, 0; 0, t)`.
    This function returns representatives for the orbit set `\Gamma_0(N)
    \backslash T \Gamma_0(M)`, where `\Gamma_0(N)` acts on the left on `T
    \Gamma_0(M)`.

    INPUT:

    - ``N`` -- int
    - ``M`` -- int (divisor of `N`)
    - ``t`` -- int (divisor of `N/M`)

    OUTPUT:

    list -- list of lists ``[a,b,c,d]``, where ``[a,b,c,d]`` should be viewed
    as a 2x2 matrix.

    This function is used for computation of degeneracy maps between
    spaces of modular symbols, hence its name.

    We use that `T^{-1} \cdot (a,b;c,d) \cdot T = (a,bt; c/t,d)`, that the
    group `T^{-1} \Gamma_0(N) T` is contained in `\Gamma_0(M)`, and that
    `\Gamma_0(N) T` is contained in `T \Gamma_0(M)`.

    ALGORITHM:

    1. Compute representatives for $\Gamma_0(N/t,t)$ inside of $\Gamma_0(M)$:

      + COSET EQUIVALENCE: Two right cosets represented by `[a,b;c,d]` and
        `[a',b';c',d']` of `\Gamma_0(N/t,t)` in `{\rm SL}_2(\ZZ)` are equivalent if
        and only if `(a,b)=(a',b')` as points of `\mathbf{P}^1(\ZZ/t\ZZ)`,
        i.e., `ab' \cong ba' \pmod{t}`, and `(c,d) = (c',d')` as points of
        `\mathbf{P}^1(\ZZ/(N/t)\ZZ)`.

      + ALGORITHM to list all cosets:

        a) Compute the number of cosets.
        b) Compute a random element `x` of `\Gamma_0(M)`.
        c) Check if x is equivalent to anything generated so far; if not, add x
           to the list.
        d) Continue until the list is as long as the bound
           computed in step (a).

    2. There is a bijection between `\Gamma_0(N)\backslash T \Gamma_0(M)` and
       `\Gamma_0(N/t,t) \backslash \Gamma_0(M)` given by `T r \leftrightarrow
       r`. Consequently we obtain coset representatives for
       `\Gamma_0(N)\backslash T \Gamma_0(M)` by left multiplying by `T` each
       coset representative of `\Gamma_0(N/t,t) \backslash \Gamma_0(M)` found
       in step 1.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import degeneracy_coset_representatives_gamma0
        sage: len(degeneracy_coset_representatives_gamma0(13, 1, 1))
        14
        sage: len(degeneracy_coset_representatives_gamma0(13, 13, 1))
        1
        sage: len(degeneracy_coset_representatives_gamma0(13, 1, 13))
        14
    """
    if N % M != 0:
        raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)

    if (N/M) % t != 0:
        raise ArithmeticError, "t (=%s) must be a divisor of N/M (=%s)"%(t,N/M)

    cdef int n, i, j, k, aa, bb, cc, dd, g, Ndivt, halfmax, is_new
    cdef int* R

    # total number of coset representatives that we'll find
    n = Gamma0(N).index() / Gamma0(M).index()
    k = 0   # number found so far
    Ndivt = N / t
    R = <int*> sig_malloc(sizeof(int) * (4*n))
    if R == <int*>0:
        raise MemoryError
    halfmax = 2*(n+10)
    while k < n:
        # try to find another coset representative.
        cc = M*random.randrange(-halfmax, halfmax+1)
        dd =   random.randrange(-halfmax, halfmax+1)
        g = arith_int.c_xgcd_int(-cc,dd,&bb,&aa)
        if g == 0: continue
        cc = cc / g
        if cc % M != 0: continue
        dd = dd / g
        # Test if we've found a new coset representative.
        is_new = 1
        for i from 0 <= i < k:
            j = 4*i
            if (R[j+1]*aa - R[j]*bb)%t == 0 and \
               (R[j+3]*cc - R[j+2]*dd)%Ndivt == 0:
                is_new = 0
                break
        # If our matrix is new add it to the list.
        if is_new:
            R[4*k] = aa
            R[4*k+1] = bb
            R[4*k+2] = cc
            R[4*k+3] = dd
            k = k + 1

    # Return the list left multiplied by T.
    S = []
    for i from 0 <= i < k:
        j = 4*i
        S.append([R[j], R[j+1], R[j+2]*t, R[j+3]*t])
    sig_free(R)
    return S

def degeneracy_coset_representatives_gamma1(int N, int M, int t):
    r"""
    Let `N` be a positive integer and `M` a divisor of `N`.  Let `t` be a
    divisor of `N/M`, and let `T` be the `2 \times 2` matrix `(1,0; 0,t)`.
    This function returns representatives for the orbit set `\Gamma_1(N)
    \backslash T \Gamma_1(M)`, where `\Gamma_1(N)` acts on the left on `T
    \Gamma_1(M)`.

    INPUT:

    - ``N`` -- int
    - ``M`` -- int (divisor of `N`)
    - ``t`` -- int (divisor of `N/M`)

    OUTPUT:

    list -- list of lists ``[a,b,c,d]``, where ``[a,b,c,d]`` should be viewed
    as a 2x2 matrix.

    This function is used for computation of degeneracy maps between
    spaces of modular symbols, hence its name.

    ALGORITHM:

    Everything is the same as for
    :func:`~degeneracy_coset_representatives_gamma0`, except for coset
    equivalence.   Here `\Gamma_1(N/t,t)` consists of matrices that are of the
    form `(1,*; 0,1) \bmod N/t` and `(1,0; *,1) \bmod t`.

    COSET EQUIVALENCE: Two right cosets represented by `[a,b;c,d]` and
    `[a',b';c',d']` of `\Gamma_1(N/t,t)` in `{\rm SL}_2(\ZZ)` are equivalent if
    and only if

    .. math::

        a \cong a' \pmod{t},
        b \cong b' \pmod{t},
        c \cong c' \pmod{N/t},
        d \cong d' \pmod{N/t}.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import degeneracy_coset_representatives_gamma1
        sage: len(degeneracy_coset_representatives_gamma1(13, 1, 1))
        168
        sage: len(degeneracy_coset_representatives_gamma1(13, 13, 1))
        1
        sage: len(degeneracy_coset_representatives_gamma1(13, 1, 13))
        168
    """

    if N % M != 0:
        raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)

    if (N/M) % t != 0:
        raise ArithmeticError, "t (=%s) must be a divisor of N/M (=%s)"%(t,N/M)

    cdef int d, g, i, j, k, n, aa, bb, cc, dd, Ndivt, halfmax, is_new
    cdef int* R


    # total number of coset representatives that we'll find
    n = Gamma1(N).index() / Gamma1(M).index()
    d = arith_int.c_gcd_int(t, N/t)
    n = n / d
    k = 0   # number found so far
    Ndivt = N / t
    R = <int*> sig_malloc(sizeof(int) * (4*n))
    if R == <int*>0:
        raise MemoryError
    halfmax = 2*(n+10)
    while k < n:
        # try to find another coset representative.
        cc =     M*random.randrange(-halfmax, halfmax+1)
        dd = 1 + M*random.randrange(-halfmax, halfmax+1)
        g = arith_int.c_xgcd_int(-cc,dd,&bb,&aa)
        if g == 0: continue
        cc = cc / g
        if cc % M != 0: continue
        dd = dd / g
        if M != 1 and dd % M != 1: continue
        # Test if we've found a new coset representative.
        is_new = 1
        for i from 0 <= i < k:
            j = 4*i
            if (R[j] - aa)%t == 0 and \
               (R[j+1] - bb)%t == 0 and \
               (R[j+2] - cc)%(Ndivt) == 0 and \
               (R[j+3] - dd)%(Ndivt) == 0:
                is_new = 0
                break
        # If our matrix is new add it to the list.
        if is_new:
            if k > n:
                sig_free(R)
                raise RuntimeError, "bug!!"
            R[4*k] = aa
            R[4*k+1] = bb
            R[4*k+2] = cc
            R[4*k+3] = dd
            k = k + 1

    # Return the list left multiplied by T.
    S = []
    for i from 0 <= i < k:
        j = 4*i
        S.append([R[j], R[j+1], R[j+2]*t, R[j+3]*t])
    sig_free(R)
    return S

def generators_helper(coset_reps, level):
    r"""
    Helper function for generators of Gamma0, Gamma1 and GammaH.

    These are computed using coset representatives, via an "inverse
    Todd-Coxeter" algorithm, and generators for `{\rm SL}_2(\ZZ)`.

    ALGORITHM: Given coset representatives for a finite index
    subgroup `G` of `{\rm SL}_2(\ZZ)` we compute generators for `G` as follows.
    Let `R` be a set of coset representatives for `G`.  Let `S, T \in {\rm
    SL}_2(\ZZ)` be defined by `(0,-1; 1,0)` and `(1,1,0,1)`, respectively.
    Define maps `s, t: R \to G` as follows. If `r \in R`, then there exists a
    unique `r' \in R` such that `GrS = Gr'`. Let `s(r) = rSr'^{-1}`. Likewise,
    there is a unique `r'` such that `GrT = Gr'` and we let `t(r) = rTr'^{-1}`.
    Note that `s(r)` and `t(r)` are in `G` for all `r`.  Then `G` is generated
    by `s(R)\cup t(R)`.

    There are more sophisticated algorithms using group actions on trees (and
    Farey symbols) that give smaller generating sets -- this code is now
    deprecated in favour of the newer implementation based on Farey symbols.

    EXAMPLES::

        sage: Gamma0(7).generators(algorithm="todd-coxeter") # indirect doctest
        [
        [1 1]  [-1  0]  [ 1 -1]  [1 0]  [1 1]  [-3 -1]  [-2 -1]  [-5 -1]
        [0 1], [ 0 -1], [ 0  1], [7 1], [0 1], [ 7  2], [ 7  3], [21  4],
        <BLANKLINE>
        [-4 -1]  [-1  0]  [ 1  0]
        [21  5], [ 7 -1], [-7  1]
        ]
    """
    cdef Matrix_integer_dense x,y,z,v,vSmod,vTmod

    crs = coset_reps.list()
    try:
        reps = [Matrix_integer_dense(Mat2Z,lift_to_sl2z(c, d, level),False,True) for c,d in crs]
    except Exception:
        raise ArithmeticError("Error lifting to SL2Z: level=%s crs=%s" % (level, crs))
    ans = []
    cdef Py_ssize_t i
    for i in range(len(crs)):
        x = reps[i]
        v = Matrix_integer_dense(Mat2Z,[crs[i][0],crs[i][1],0,0],False,True)
        vSmod = (v*genS)
        vTmod = (v*genT)
        y_index = coset_reps.normalize(vSmod[0,0],vSmod[0,1])
        z_index = coset_reps.normalize(vTmod[0,0],vTmod[0,1])
        y_index = crs.index(y_index)
        z_index = crs.index(z_index)
        y = reps[y_index]
        z = reps[z_index]
        y = y._invert_unit()
        z = z._invert_unit()
        ans.append(x*genS*y)
        ans.append(x*genT*z)
    return [x for x in ans if x != genI]
