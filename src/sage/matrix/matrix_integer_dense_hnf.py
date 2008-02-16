"""
Modular algorithm to compute Hermite normal forms of integer matrices.

AUTHORS:
    -- Clement Pernet and William Stein (2008-02-07): initial version

TODO:

   [ ] provably correct det
   [ ] provably correct hnf (since we use rank profile -- if wrong need to retry!)
   [ ] transformation matrix
   [ ] fix memory leaks:
sage..: a = random_matrix(ZZ, 300,400)
sage..: get_memory_usage()
'626M+'
sage..: time h = hnf(a
sage..: get_memory_usage()
'634M+'

   [ ] Make this the default in Sage and test free modules, etc.
   [ ] automated testing for correctness
"""

from copy import copy

from sage.misc.misc import verbose, prod
from sage.matrix.constructor import random_matrix, matrix, matrix

from sage.rings.all import ZZ, QQ, previous_prime, CRT_list
import math

#MAX_DET_PRIME = 67108879   # next prime after 2^26 -- biggest for linbox (?)
MAX_DET_PRIME=16777259      # but this prime much faster for linbox

def det_from_modp_and_divisor(A, d, p, z_mod, moduli):
    tm = verbose("Multimodular stage of det calculation -- using p = %s"%p, level=1)
    z = A._linbox_modn_det(p) / d
    z = z.lift()
    z_mod.append(z)
    moduli.append(p)
    z = CRT_list(z_mod, moduli)
    N = prod(moduli)
    if z > N//2:
        z = z - N
    verbose("finished multimodular det for p = %s"%p, tm, level=1)
    return d * z

def det_given_divisor(A, d, proof=True, stabilize=3):
    """
    Given a divisor d of the determinant of A, compute the
    determinant of A.

    INPUT:
        A -- a square integer matrix
        d -- an integer that is assumed to divide the determinant of A
        proof -- bool (default True) compute det modulo enough primes
                 so that the determinant is computed provably correctly
                 (via the Hadamard bound).
        stabilize -- int (default: 3) if proof = False, then compute det
                 mod p until stabilize successive modulo det computations
                 stabilize.
    """
    p = MAX_DET_PRIME
    z_mod = []
    moduli = []
    if proof:
        N = 1
        B = (10**A.hadamard_bound()) // d + 1
        dd = d
        verbose("Multimodular det -- need to use about %s primes."%(int(math.log(B)/math.log(p))))
        while N < B:
            dd = det_from_modp_and_divisor(A, d, p, z_mod, moduli)
            N *= p
            p = previous_prime(p)
        return dd
    else:
        val = []
        while True:
            dd = det_from_modp_and_divisor(A, d, p, z_mod, moduli)
            val.append(dd)
            if len(val) >= stabilize and len(set(val[-stabilize:])) == 1:
                return val[-1]
            p = previous_prime(p)

def det_padic(A, proof=True, stabilize=3):
    """
    Return the determinant of A, computed using a p-adic/multimodular
    algorithm.

    INPUTS:
        A -- a square matrix
        proof -- boolean
        stabilize -- if proof False, number of successive primes so that
                     CRT det must stabilize.

    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as h
        sage: a = matrix(ZZ, 3, [1..9])
        sage: h.det_padic(a)
        0
        sage: a = matrix(ZZ, 3, [1,2,5,-7,8,10,192,5,18])
        sage: h.det_padic(a)
        -3669
        sage: a.determinant(algorithm='ntl')
        -3669
    """
    if not A.is_square():
        raise ValueError, "A must be a square matrix"
    r = A.rank()
    if r < A.nrows():
        return ZZ(0)
    v = random_matrix(ZZ, A.nrows(), 1)
    d = A.solve_right(v, check_rank=False).denominator()
    return det_given_divisor(A, d, proof=proof, stabilize=stabilize)

def double_det (A, b, c, proof):
    """
    Compute the determinants of the stacked integer matrices
    A.stack(b) and A.stack(c).

    INPUT:
        A -- an (n-1) x n matrix
        b -- an 1 x n matrix
        c -- an 1 x n matrix
        proof -- whether or not to compute the det modulo enough times
                 to provably compute the determinant.

    OUTPUT:
        a pair of two integers.

    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: A = matrix(ZZ, 2, 3, [1,2,3, 4,-2,5])
        sage: b = matrix(ZZ, 1, 3, [1,-2,5])
        sage: c = matrix(ZZ, 1, 3, [8,2,10])
        sage: A.stack(b).det()
        -48
        sage: A.stack(c).det()
        42
        sage: hnf.double_det(A, b, c, False)
        (-48, 42)
    """
    # We use the "two for the price of one" algorithm, which I made up. (William Stein)

    # This is a clever trick!  First we transpose everything.  Then
    # we use that if [A|b]*v = c then [A|c]*w = b with w easy to write down!
    # In fact w is got from v by dividing all entries by -v[n], where n is the
    # number of rows of v, and *also* dividing the last entry of w by v[n] again.
    # See this as an algebra exercise where you have to think of matrix vector
    # multiply as "linear combination of columns".
    A = A.transpose()
    b = b.transpose()
    c = c.transpose()
    t = verbose('starting double det')
    B = A.augment(b)
    v = B.solve_right(-c, check_rank=True)  # infinite loop if not full rank and don't do this.

    # we use stabilize=2, since the det is typically so small with this construction
    db = det_given_divisor(B, v.denominator(), proof=proof, stabilize=2)

    n = v.nrows()
    vn = v[n-1,0]
    w = (-1/vn)*v
    w[n-1] = w[n-1]/vn
    dc = det_given_divisor(A.augment(c), w.denominator(), proof=proof, stabilize=2)

    verbose('finished double det', t)

    return (db, dc)


def add_column(B, H_B, a):
    """
    The add column procedure.

    INPUT:
        B   -- a non-singular square matrix
        H_B -- the Hermite normal form of B
        a   -- a column vector

    OUTPUT:
        x   -- a vector such that H' = H_B.augment(x) is the HNF of A = B.augment(a).

    EXAMPLES:
        sage: B = matrix(ZZ, 3, 3, [1,2,5, 0,-5,3, 1,1,2])
        sage: H_B = B.echelon_form()
        sage: a = matrix(ZZ, 3, 1, [1,8,-2])
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: x = hnf.add_column(B, H_B, a); x
        [18]
        [ 3]
        [23]
        sage: H_B.augment(x)
        [ 1  0 17 18]
        [ 0  1  3  3]
        [ 0  0 18 23]
        sage: B.augment(a).echelon_form()
        [ 1  0 17 18]
        [ 0  1  3  3]
        [ 0  0 18 23]
    """
    t0 = verbose('starting add_column')

    # We use a direct solve method without inverse.  This
    # is more clever than what is in Allan Steel's talk and
    # what is in that paper, in 2 ways -- (1) no inverse need
    # to be computed, and (2) we cleverly solve a vastly easier
    # system and recover the solution to the original system.

    # Here's how:
    # 1. We make a copy of B but with the last *nasty* row of B replaced
    #    by a random very nice row.
    C = copy(B)
    C[C.nrows()-1] = [1]*C.ncols()

    # 2. Then we find the unique solution to C * x = a
    #    (todo -- recover from bad case.)
    x = C.solve_right(a)

    # 3. We next delete the last row of B and find a basis vector k
    #    for the 1-dimensional kernel.
    D = B.matrix_from_rows(range(C.nrows()-1))
    N = D._rational_kernel_iml()
    if N.ncols() != 1:
        raise NotImplementedError, "need to recover gracefully from rank issues with matrix."
    k = N.matrix_from_columns([0])

    # 4. The sought for solution z to B*z = a is some linear combination
    #       z = x + alpha*k
    # and setting w to be the last row of B, this column vector z satisfies
    #       w * z = a'
    # where a' is the last entry of a.  Thus
    #       w * (x + alpha*k) = a'
    # so    w * x + alpha*w*k = a'
    # so    alpha*w*k  = a' - w*x.

    w = B[-1]  # last row of B
    a_prime = a[-1]
    lhs = w*k
    rhs = a_prime - w * x
    alpha = rhs[0] / lhs[0]
    z = x + alpha*k

    zd, d = z._clear_denom()
    x = H_B * zd
    if d != 1:
        for i in range(x.nrows()):
            x[i,0] = x[i,0]/d

    return x

def add_row(A, b, pivots, include_zero_rows):
    """
    The add row procedure.

    INPUT:
        A -- a matrix in Hermite normal form with n column
        b -- an n x 1 row matrix
        pivots -- sorted list of integers; the pivot positions of A.

    OUTPUT:
        H -- the Hermite normal form of A.stack(b).
        new_pivots -- the pivot columns of H.

    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: A = matrix(ZZ, 2, 3, [-21, -7, 5, 1,20,-7])
        sage: b = matrix(ZZ, 1,3, [-1,1,-1])
        sage: hnf.add_row(A, b, A.pivots(), True)
        ([ 1  6 29]
        [ 0  7 28]
        [ 0  0 46], [0, 1, 2])
        sage: A.stack(b).echelon_form()
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
    """
    t = verbose('add hnf row')
    v = b.row(0)
    H, pivs = A._add_row_and_maintain_echelon_form(b.row(0), pivots)
    if include_zero_rows and H.nrows() != A.nrows()+1:
        H = H.matrix_from_rows(range(A.nrows()+1))
    verbose('finished add hnf row', t)
    return H, pivs

def hnf_square(A, proof):
    """
    INPUT:
        a nonsingular n x n matrix A over the integers.
    OUTPUT:
        the Hermite normal form of A.

    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: A = matrix(ZZ, 3, [-21, -7, 5, 1,20,-7, -1,1,-1])
        sage: hnf.hnf_square(A, False)
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
        sage: A.echelon_form()
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
    """
    n = A.nrows()
    m = A.ncols()
    if n != m:
        raise NotImplementedError, "A must be square."

    # Small cases -- don't use this algorithm
    if n <= 3:
        return A.echelon_form(algorithm="pari")


    t = verbose("starting slicings")
    B = A.matrix_from_rows(range(m-2)).matrix_from_columns(range(n-1))
    c = A.matrix_from_rows([m-2]).matrix_from_columns (range(n-1))
    d = A.matrix_from_rows([m-1]).matrix_from_columns (range(n-1))
    b = A.matrix_from_columns([n-1]).matrix_from_rows(range(m-2))
    verbose("done slicing", t)

    try:
        (d1,d2) = double_det (B,c,d, proof=proof)
    except (ValueError, ZeroDivisionError):
        verbose("det computation failed -- we compute hnf of submatrix directly.")
        d1 = d2 = ZZ(0)
    (g,k,l) = d1._xgcd (d2, minimal=True)
    W = B.stack (k*c + l*d)
    verbose("submatrix det: g=%s"%g)
    if g == 0 or g > 2**30:  # inconceivable that g will be at all large on even slightly random input
        H = hnf(W)
    else:
        H = W._hnf_mod(2*g)
    x = add_column(B.stack(k*c+l*d), H, b.stack(matrix(1,1,[k*A[m-2,m-1] + l*A[m-1,m-1]])))
    Hprime = H.augment(x)
    pivots = range(Hprime.nrows())
    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m-2]), pivots, include_zero_rows=False)
    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m-1]), pivots, include_zero_rows=False)
    H = Hprime.matrix_from_rows(range(m))
    return H

def interleave_matrices(A, B, cols1, cols2):
    """
    INPUT:
        A, B -- matrices with the same number of rows
        cols1, cols2 -- disjoint lists of integers
    OUTPUT:
        construct a new matrix C by sticking the columns
        of A at the positions specified by cols1 and the
        columns of B at the positions specified by cols2.

    EXAMPLES:
        sage: A = matrix(ZZ, 2, [1,2,3,4]); B = matrix(ZZ, 2, [-1,5,2,3])
        sage: A
        [1 2]
        [3 4]
        sage: B
        [-1  5]
        [ 2  3]
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.interleave_matrices(A, B, [1,3], [0,2])
        [-1  1  5  2]
        [ 2  3  3  4]
    """
    D = A.augment(B)
    w = cols1 + cols2
    v = [w.index(i) for i in range(len(cols1) + len(cols2))]
    return D.matrix_from_columns(v)

def probable_pivot_rows(A):
    """
    INPUT:
        A -- a matrix
    OUTPUT:
        a list of integers
    """
    return probable_pivot_columns(A.transpose())

def probable_pivot_columns(A):
    """
    INPUT:
        A -- a matrix
    OUTPUT:
        a list of integers
    """
    p = ZZ.random_element(10007, 46000).next_prime()
    return A._reduce(p).pivots()

def ones(H, pivots):
    # Find the "onecol" pivot columns of H, i.e., the columns
    # that contain exactly one "1" entry and all other entries 0.
    onecol = []
    onerow = []
    i = 0
    for c in pivots:
        if H[i,c] == 1:
            onecol.append(c)
            onerow.append(i)
        i += 1
    onecol_set = set(onecol)
    non_onecol = [i for i in range(H.ncols()) if i not in onecol_set]
    non_onerow = [i for i in range(len(pivots)) if i not in onerow]
    return onecol, onerow, non_onecol, non_onerow

def extract_ones_data(H, pivots):
        onecol, onerow, non_onecol, non_onerow = ones(H, pivots)
        verbose('extract_ones -- got submatrix of size %s'%len(non_onecol))
        if len(non_onecol) in [1, 2]:
            # Extract submatrix of all non-onecol columns and onecol rows
            C = H.matrix_from_rows_and_columns(onerow, non_onecol)
            # Extract submatrix of all non-onecol columns and other rows
            D = H.matrix_from_rows_and_columns(non_onerow, non_onecol).transpose()
            tt = verbose("extract ones -- INVERT %s x %s"%(len(non_onecol), len(non_onecol)), level=1)
            E = D**(-1)
            verbose("done inverting", tt, level=1)
            return C, D, E, onecol, onerow, non_onecol, non_onerow
        else:
            return None, None, None, onecol, onerow, non_onecol, non_onerow

def is_in_hnf_form(H):
    """
    Return True precisely if the matrix H is in Hermite normal form.
    """
    tt = verbose('testing if matrix is in HNF')
    for i in xrange(H.nrows()):
        for j in xrange(H.ncols()):
            if i > j:
                if H[i,j]:
                    verbose('done verifying in HNF - not in HNF', tt)
                    return False
            elif i < j and j < H.nrows():
                if H[i,j] < 0 or H[i,j] >= H[j,j]:
                    verbose('done verifying in HNF - not in HNF', tt)
                    return False
            elif i == j:
                if H[i,i] < 0:
                    verbose('done verifying in HNF - not in HNF', tt)
                    return False
    verbose('done verifying in HNF', tt)
    return True

def hnf(A, include_zero_rows=True, proof=False):
    """
    INPUT:
        A -- an n x m matrix A over the integers.
        include_zero_rows -- bool (default: True) whether or not to
                             include zero rows in the output matrix
        proof -- whether or not to prove the result correct.

    OUTPUT:
        the Hermite normal form of A.
    """
    if proof == False:
        return probable_pivots_hnf(A, include_zero_rows = include_zero_rows, proof=False)

    while True:
        try:
            H = probable_pivots_hnf(A, include_zero_rows = include_zero_rows, proof=True)
        except (AssertionError, ZeroDivisionError):
            verbose("Assertion occured when computing HNF; guessed pivot columns likely wrong.")
            pass
        else:
            if is_in_hnf_form(H):
                return H
        verbose("After attempt the return matrix is not in HNF form since pivots must have been wrong.  We try again.")


def probable_pivots_hnf(A, include_zero_rows, proof):
    """
    OUTPUT:
        the Hermite normal form of A.
    """
    # Find left-most full rank submatrix by working modulo a prime
    rows = probable_pivot_rows(A)
    B    = A.matrix_from_rows(rows)
    cols = probable_pivot_columns(B)
    C   = B.matrix_from_columns(cols)
    # Now C is a submatrix of A that has full rank and is square.

    # We compute the HNF of C.
    H = hnf_square(C, proof=proof)

    # The transformation matrix to HNF is the unique
    # matrix U such that U * C = H, i.e., U = H*C^(-1).

    if len(cols) < B.ncols():
        # We compute the HNF of B by finding the
        # HNF of B by multiplying the matrix D
        # got from the columns not in C by U:
        # We want to compute X = U*D.  But U = H*C^(-1),
        # so X = U*D = H*C^(-1)*D.
        # So C*H^(-1)*X = D

        # find y s.t C*y = D
        #   H^(-1)*X = y ===> X = H*y
        #
        cols_set = set(cols)
        cols2 = [i for i in range(B.ncols()) if not i in cols_set]
        D = B.matrix_from_columns(cols2)
        Y = C.solve_right(D)
        H2 = H*Y
        H2 = H2.change_ring(ZZ)
        # The HNF of B is got by assembling together
        # the matrices H and H2.
        H = interleave_matrices(H, H2, cols, cols2)

    # Now H is the HNF of the matrix B
    # Finally we add all remaining rows of A to H using
    # the add_row function.
    C, D, E, onecol, onerow, non_onecol, non_onerow = extract_ones_data(H, cols)
    if len(non_onecol) == 0:
        # Identity matrix -- done
        verbose("hnf -- got identity matrix -- early abort")
        return pad_zeros(H, A.nrows() - H.nrows())

    rows_set = set(rows)
    nz = 0
    for i in range(A.nrows()):
        if not i in rows_set:
            if E is None:
                H, cols = add_row(H, A.matrix_from_rows([i]), cols, include_zero_rows)
                C, D, E, onecol, onerow, non_onecol, non_onerow = extract_ones_data(H, cols)
                if len(non_onecol) == 0:
                    # Identity matrix -- done
                    verbose("hnf -- got identity matrix -- early abort")
                    return pad_zeros(H, A.nrows() - H.nrows())
            else:
                z = A.matrix_from_rows_and_columns([i], non_onecol)
                w = A.matrix_from_rows_and_columns([i], onecol)
                tt = verbose("checking denom (%s x %s)"%(D.nrows(), D.ncols()))
                Y = (z - w*C).transpose()
                k = E*Y
                verbose("done checking denom",tt)
                if k.denominator() != 1:
                    H, cols = add_row(H, A.matrix_from_rows([i]), cols, include_zero_rows)
                    D = H.matrix_from_rows_and_columns(non_onerow, non_onecol).transpose()
                else:
                    nz += 1
                nn = ones(H, cols)
                if len(nn[2]) == 0:
                    verbose("hnf -- got identity matrix -- early abort")
                    return pad_zeros(H, A.nrows() - H.nrows())


    if include_zero_rows:
        return pad_zeros(H, nz)
    return H

def pad_zeros(H, nz):
    if nz == 0:
        return H
    return H.stack(matrix(ZZ, nz, H.ncols()))


def benchmark_hnf(nrange, bits=4):
    """
    Run benchmark program.

    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.benchmark_hnf([50,100],32)
        ('sage', 50, 32, ...),
        ('sage', 100, 32, ...),
    """
    from sage.misc.misc import cputime
    b = 2**bits
    for n in nrange:
        a = random_matrix(ZZ, n, x=-b,y=b)
        t = cputime()
        h = hnf(a)
        tm = cputime(t)
        print '%s,'%(('sage', n, bits, tm),)

def benchmark_magma_hnf(nrange, bits=4):
    """
    EXAMPLES:
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.benchmark_magma_hnf([50,100],32)     # optional -- requires magma
        ('magma', 50, 32, ...),
        ('magma', 100, 32, ...),
    """
    from sage.misc.misc import cputime
    from sage.interfaces.all import magma
    b = 2**bits
    for n in nrange:
        a = magma('MatrixAlgebra(IntegerRing(),%s)![Random(%s,%s) : i in [1..%s]]'%(n,-b,b,n**2))
        t = magma.cputime()
        h = a.EchelonForm()
        tm = magma.cputime(t)
        print '%s,'%(('magma', n, bits, tm),)


################################################################
# Integer Kernel
#################################################################

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


