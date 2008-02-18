"""
Saturation over ZZ
"""

from sage.rings.all import ZZ, gcd, binomial, GF
from sage.matrix.constructor import identity_matrix, random_matrix
from sage.misc.misc import verbose
import matrix_integer_dense_hnf
from copy import copy


def p_saturation(A, p, proof=True):
    """
    INPUT:
        A -- a matrix over ZZ
        p -- a prime
        proof -- bool (default: True)

    OUTPUT:
        The p-saturation of the matrix A, i.e., a new matrix in Hermite form
        whose row span a ZZ-module that is p-saturated.

    ISSUES:
         * What if p is large? -- Sage linear algebra is maybe
           currently only implemented for small p.  This could be a
           major problem.  Maybe Linbox will help.
    """
    tm = verbose("%s-saturating a %sx%s matrix"%(p, A.nrows(), A.ncols()))
    H = A.hermite_form(include_zero_rows=False, proof=proof)
    while True:
        if p == 2:
            A = H.change_ring(GF(p))
        else:
            try:
                # Faster than change_ring
                A = H._reduce(p)
            except OverflowError:
                # fall back to generic GF(p) matrices
                A = H.change_ring(GF(p))
        assert A.nrows() <= A.ncols()
        K = A.kernel()
        if K.dimension() == 0:
            verbose("done saturating", tm)
            return H
        B = K.basis_matrix().lift()
        C = ((B * H) / p).change_ring(ZZ)
        H = H.stack(C).hermite_form(include_zero_rows=False, proof=proof)
    verbose("done saturating", tm)

import random
def random_sublist_of_size(k, n):
    """
    INPUT:
        k -- an integer
        n -- an integer

    OUTPUT:
        a randomly chosen sublist of range(k) of size n.
    """
    if n > k:
        raise ValueError, "n must be <= len(v)"
    if n == k:
        return range(k)
    if n >= k//2+5:
        # use complement
        w = random_sublist_of_size(k, k-n)
        m = set(w)
        w = [z for z in xrange(k) if z not in m]
        assert(len(w)) == n
        return w

    w = set([])
    while len(w) < n:
        z = random.randrange(k)
        if not z in w:
            w.add(z)
    w = list(w)
    w.sort()
    return w

def solve_system_with_difficult_last_row(B, A):
    """
    Solve the matrix equation B*Z = A when the last row of $B$
    contains huge entries.

    INPUT:
        B -- a square n x n nonsingular matrix with painful big bottom row.
        A -- an n x k matrix.
    OUTPUT:
        the unique solution to B*Z = A.

    EXAMPLES:
    """
    # See the comments in the function of the same name in matrix_integer_dense_hnf.py.
    # This function is just a generalization of that one to A a matrix.
    C = copy(B)
    while True:
        C[C.nrows()-1] = random_matrix(ZZ,1,C.ncols()).row(0)
        try:
            X = C.solve_right(A)
        except ValueError:
            verbose("Try difficult solve again with different random vector")
        else:
            break
    D = B.matrix_from_rows(range(C.nrows()-1))
    N = D._rational_kernel_iml()
    if N.ncols() != 1:
        verbose("Difficult solve quickly failed.  Using direct approach.")
        return B.solve_right(A)

    tm = verbose("Recover correct linear combinations")
    k = N.matrix_from_columns([0])

    # The sought for solution Z to B*Z = A is some linear combination
    #       Z = X + alpha*k
    #  Let w be the last row of B; then Z satisfies
    #       w * Z = A'
    # where A' is the last row of A.  Thus
    #       w * (X + alpha*k) = A'
    # so    w * X + alpha*w*k = A'
    # so    alpha*w*k  = A' - w*X.
    w = B[-1]  # last row of B
    A_prime = A[-1]  # last row of A
    lhs = w*k
    rhs = A_prime - w * X

    if lhs[0] == 0:
        verbose("Difficult solve quickly failed.  Using direct approach.")
        return B.solve_right(A)

    for i in range(X.ncols()):
        alpha = rhs[i] / lhs[0]
        X.set_column(i, (X.matrix_from_columns([i]) + alpha*k).list())
    verbose("Done getting linear combinations.", tm)
    return X

def saturation(A, proof=True, p=0, max_dets=5):
    """
    INPUT:
        A     -- a matrix over ZZ
        proof -- bool (default: True)
        p     -- int (default: 0); if not 0
                 only guarantees that output is p-saturated
        max_dets -- int (default: 4) max number of dets of
                 submatrices to compute.

    OUTPUT:
        the saturation of A
    """
    # Find a submatrix of full rank and instead saturate that matrix.
    r = A.rank()
    if A.is_square() and r == A.nrows():
        return identity_matrix(ZZ, r)
    if A.nrows() > r:
        P = []
        while len(P) < r:
            P = matrix_integer_dense_hnf.probable_pivot_rows(A)
        A = A.matrix_from_rows(P)

    # Factor out all common factors from all rows, just in case.
    A._factor_out_common_factors_from_each_row()

    if A.nrows() <= 1:
        return A

    if max_dets > 0:
        # Take the GCD of at most num_dets randomly chosen determinants.
        nr = A.nrows(); nc = A.ncols()
        d = 0
        trials = min(binomial(nc, nr), max_dets)
        already_tried = []
        while len(already_tried) < trials:
            v = random_sublist_of_size(nc, nr)
            tm = verbose('saturation -- checking det condition on submatrix')
            d = gcd(d, A.matrix_from_columns(v).determinant(proof=proof))
            verbose('saturation -- got det down to %s'%d, tm)
            if gcd(d, p) == 1:
                return A
            already_tried.append(v)

        if gcd(d, p) == 1:
            # already p-saturated
            return A

        # Factor and p-saturate at each p.
        # This is not a good algorithm, because all the HNF's in it are really slow!
        #
        #tm = verbose('factoring gcd %s of determinants'%d)
        #limit = 2**31-1
        #F = d.factor(limit = limit)
        #D = [p for p, e in F if p <= limit]
        #B = [n for n, e in F if n > limit]  # all big factors -- there will only be at most one
        #assert len(B) <= 1
        #C = B[0]
        #for p in D:
        #    A = p_saturation(A, p=p, proof=proof)

    # This is a really simple but powerful algorithm.
    # FACT: If A is a matrix of full rank, then hnf(transpose(A))^(-1)*A is a saturation of A.
    # To make this practical we use solve_system_with_difficult_last_row, since the
    # last column of HNF's are typically the only really big ones.
    B = A.transpose().hermite_form(include_zero_rows=False, proof=proof)
    if B.nrows() > r:
        B = B.matrix_from_rows(range(r))
    B = B.transpose()

    # Now compute B^(-1) * A
    C = solve_system_with_difficult_last_row(B, A)
    return C.change_ring(ZZ)



################################################################
# Saturation
# David Kohel sent me the following a couple of years ago.
# It's probably the algorithm to use.
## function pAdicSaturation(B,p)
##     if Type(B[1]) eq SeqEnum then
##         V := RSpace(Rationals(),#B[1]);
## 	B := [ V | v : v in B];
##     end if;
##     V := Universe(B);
##     n := Degree(V);
##     for i in [1..#B] do
## 	B[i] *:= LCM([ Denominator(c) : c in Eltseq(B[i]) ]);
##     end for;
##     ZZ := Integers();
##     FF := FiniteField(p);
##     B := RMatrixSpace(ZZ,#B,n)!Matrix(B);
##     m := Rank(B);
##     B := Submatrix(HermiteForm(B),1,1,m,n);
##     N := RMatrixSpace(FF,m,n)!B;
##     while Rank(N) lt m do
## 	K := Kernel(N);
## 	vprintf pAdicSaturation :
## 	    "Rank(N) + Rank(K) = %o + %o = %o\n", Rank(N), Rank(K), m;
## 	C := RMatrixSpace(ZZ,#Basis(K),n)!
## 	Matrix([ (1/p)*V!&+[ ZZ!u[i]*B[i] : i in [1..m] ] : u in Basis(K) ]);
## 	vtime pAdicSaturation, 2 :
## 	    B := Submatrix(HermiteForm(VerticalJoin(B,C)),1,1,m,n);
## 	N := RMatrixSpace(FF,m,n)!B;
##     end while;
##     vprintf pAdicSaturation : "Rank(N) = %o \n", Rank(N), m;
##     return [ B[i] : i in [1..m] ];
## end function;
#################################################################


##########################
# Allan also says:
## > How does the MAGMA command PureLattice work?  What is the
## > algorithm, etc.?
## > Do you do this:
## >
## > 1. Find echelon form of basis of lattice.
## >
## > 2. Write down matrix over $\Z$ that has saturation of lattice
## >    as kernel.
## >
## > 3. Find the kernel using algorithm 2.7.2 of Cohen's book (Kernel
##    over Z using LLL).
## More complicated than this.  That would work, but requires 2 kernels
## and the 2nd one can't done by a modular algorithm: I don't want to
## compute kernels, because I do that by modular methods and they only do
## it over Q and then you need this very saturation alg to get the kernel
## over Z!!!
## Here is basic form of one "standard" saturation algorithm, which
## I used to do:
##     Given basis B.
##     H = HermiteForm(B);
##     for (;;)
##     {
## 	Get Smith form S of H and P so that S = P*H*?;
## 	    [right transformation mat ? not needed]
## 	If diag of S is all ones, then return H;
## 	H = P*H;
## 	Remove content from all rows of H;
## 	    [if entry (i,i) of S has val d > 1, then d will divide all
## 	     entries of row i of H]
##     }
## Then H is basis at the end.
## I have a new modular algorithm for this done about a year ago which is
## complicated -- I may publish this if I get a chance.  It is now
## used by all funcs in Magma which need saturation.  It finds the largest
## elem divisor D of B by modular method, then partially factors D and
## for small divisors, uses modular method to get rid of those primes, and
## then does the big primes another way I think.
