"""
Helper functions for congruence group class
"""

import random

import sage.rings.arith

cimport arith
import arith
cdef arith.arith_int arith_int
arith_int  = arith.arith_int()

include "cdefs.pxi"


# This is the C version of a function also implemented in python in
# sage.modular.congroup.  It is orders of magnitude faster (e.g., 30
# times).  The key speedup is in replacing looping through the
# elements of the Python list R with looping through the elements of a
# C-array.

def degeneracy_coset_representatives_gamma0(int N, int M, int t):
    r"""
    Let $N$ be a positive integer and $M$ a multiple of $N$.  Let $t$ be a
    divisor of $N/M$, and let $T$ be the $2x2$ matrix $T=[0,0; 0,t]$.
    This function returns representatives for the orbit set
    $\Gamma_0(N) \backslash T \Gamma_0(M)$, where $\Gamma_0(N)$
    acts on the left on $T \Gamma_0(M)$.

    INPUT:
        N -- int
        M -- int (divisor of N)
        t -- int (divisors of N/M)

    OUTPUT:
        list -- list of lists [a,b,c,d], where [a,b,c,d] should
                be viewed as a 2x2 matrix.

    This function is used for computation of degeneracy maps between
    spaces of modular symbols, hence its name.

    We use that $T^{-1}*(a,b;c,d)*T = (a,bt,c/t,d),$
    that the group $T^{-1}Gamma_0(N) T$ is contained in $\Gamma_0(M)$,
    and that $\Gamma_0(N) T$ is contained in $T \Gamma_0(M)$.

    ALGORITHM:
    \begin{enumerate}
    \item Compute representatives for $\Gamma_0(N/t,t)$ inside of $\Gamma_0(M)$:
          COSET EQUIVALENCE:
           Two right cosets represented by $[a,b;c,d]$ and
           $[a',b';c',d']$ of $\Gamma_0(N/t,t)$ in $\SL_2(\Z)$
           are equivalent if and only if
           $(a,b)=(a',b')$ as points of $\P^1(\Z/t\Z)$,
           i.e., $ab' \con ba' \pmod{t}$,
           and $(c,d) = (c',d')$ as points of $\P^1(\Z/(N/t)\Z)$.

        ALGORITHM to list all cosets:
        \begin{enumerate}
           \item Compute the number of cosets.
           \item Compute a random element x of Gamma_0(M).
           \item Check if x is equivalent to anything generated so
               far; if not, add x to the list.
           \item Continue until the list is as long as the bound
               computed in A.
        \end{enumerate}

    \item There is a bijection between $\Gamma_0(N)\backslash T
          \Gamma_0(M)$ and $\Gamma_0(N/t,t) \backslash \Gamma_0(M)$
          given by $T r$ corresponds to $r$.  Consequently we obtain
          coset representatives for $\Gamma_0(N)\backslash T
          \Gamma_0(M)$ by left multiplying by $T$ each coset
          representative of $\Gamma_0(N/t,t) \Gamma_0(M)$ found in
          step 1.

    \end{enumerate}
    """
    import sage.modular.dims

    if N % M != 0:
        raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)

    if (N/M) % t != 0:
        raise ArithmeticError, "t (=%s) must be a divisor of N/M (=%s)"%(t,N/M)

    cdef int n, i, j, k, aa, bb, cc, dd, g, Ndivt, halfmax, is_new
    cdef int* R

    # total number of coset representatives that we'll find
    n = sage.modular.dims.idxG0(N) / sage.modular.dims.idxG0(M)
    k = 0   # number found so far
    Ndivt = N / t
    R = <int*> PyMem_Malloc(sizeof(int) * (4*n))
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
    PyMem_Free(R)
    return S




def degeneracy_coset_representatives_gamma1(int N, int M, int t):
    r"""
    Let $N$ be a positive integer and $M$ a multiple of $N$.  Let $t$ be a
    divisor of $N/M$, and let $T$ be the $2x2$ matrix $T=[0,0; 0,t]$.
    This function returns representatives for the orbit set
    $\Gamma_1(N) \backslash T \Gamma_1(M)$, where $\Gamma_1(N)$
    acts on the left on $T \Gamma_1(M)$.

    INPUT:
        N -- int
        M -- int (divisor of N)
        t -- int (divisors of N/M)

    OUTPUT:
        list -- list of lists [a,b,c,d], where [a,b,c,d] should
                be viewed as a 2x2 matrix.

    This function is used for computation of degeneracy maps between
    spaces of modular symbols, hence its name.

    ALGORITHM:
    Everything is the same as for degeneracy_coset_representatives_gamma1,
    except for coset equivalent.   Here $\Gamma_1(N/t,t)$ consists
    of matrices that are of the form $[1,*;0,1]$ module $N/t$ and
    $[1,0;*,1]$ modulo $t$.

           COSET EQUIVALENCE:
           Two right cosets represented by $[a,b;c,d]$ and
           $[a',b';c',d']$ of $\Gamma_1(N/t,t)$ in $\SL_2(\Z)$
           are equivalent if and only if
           $$
            a \con a' \pmod{t},
            b \con b' \pmod{t},
            c \con c' \pmod{N/t},
            d \con d' \pmod{N/t}.
           $$

    """
    import sage.modular.dims

    if N % M != 0:
        raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)

    if (N/M) % t != 0:
        raise ArithmeticError, "t (=%s) must be a divisor of N/M (=%s)"%(t,N/M)

    cdef int d, g, i, j, k, n, aa, bb, cc, dd, Ndivt, halfmax, is_new
    cdef int* R


    # total number of coset representatives that we'll find
    n = sage.modular.dims.idxG1(N) / sage.modular.dims.idxG1(M)
    d = arith_int.c_gcd_int(t, N/t)
    n = n / d
    k = 0   # number found so far
    Ndivt = N / t
    R = <int*> PyMem_Malloc(sizeof(int) * (4*n))
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
        if dd % M != 1: continue
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
            if k > n: raise RuntimeError, "bug!!"
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
    PyMem_Free(R)
    return S


