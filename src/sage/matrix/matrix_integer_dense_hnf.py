"""
Modular algorithm to compute Hermite normal forms of integer matrices.

AUTHORS:
    -- Clement Pernet and William Stein (2008-02-07): initial version

TODO:

   [ ] optimize double det
   [ ] more rows than columns
   [ ] more columns than rows
   [ ] degenerate cases -- fail nicely

"""

from copy import copy

from sage.misc.misc import verbose
from sage.matrix.constructor import random_matrix, matrix

from sage.rings.all import ZZ, QQ, previous_prime

def fastdet(A, times=1):
    """
    Likely to be correct fast det (no guarantee).
    """
    v = random_matrix(ZZ,A.nrows(), 1, x=-1,y=1)
    w = A.solve_right(v, check_rank=False)
    d = w.denominator()
    p = 46337
    for i in range(times):
        Amod = A._reduce(p)
        det = Amod.determinant() / d
        m = ZZ(det)
        if m >= p//2:
            m -= p
        d *= m
        p = previous_prime(p)
    return d

def xxx_doubleDet (A, b, c):
    """
    Compute the determinants of the stacked matrices A.stack(b) and A.stack(c).

    INPUT:
        A -- an (n-1) x n matrix
        b -- an 1 x n matrix
        c -- an 1 x n matrix

    OUTPUT:
        a pair of two integers.
    """
    t = verbose('starting double det')
    d1 = fastdet(A.stack(b))
    d2 = fastdet(A.stack(c))
    verbose('finished double det', t)
    return (d1,d2)

def doubleDet (A, b, c):
    """
    Compute the determinants of the stacked matrices A.stack(b) and A.stack(c).

    INPUT:
        A -- an (n-1) x n matrix
        b -- an 1 x n matrix
        c -- an 1 x n matrix

    OUTPUT:
        a pair of two integers.
    """
    # We use the "two for the price of one" algorithm.

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
    v = B.solve_right(-c, check_rank=False)  # infinite loop if not full rank!
    db = v.denominator()
    p = 46337
    Bmod = B._reduce(p)
    z = Bmod.determinant() / db
    z = z.lift()
    if z > p//2:
        z = z - p
    db *= z

    n = v.nrows()
    vn = v[n-1,0]
    w = (-1/vn)*v
    w[n-1] = w[n-1]/vn
    dc = w.denominator()

    Cmod = A.augment(c)._reduce(p)
    z = Cmod.determinant() / dc
    z = z.lift()
    if z > p//2:
        z = z - p
    dc *= z

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
    if d > 1:
        for i in range(x.ncols()):
            x[i,0] = x[i,0]/d

    verbose('finished add column', t0)

    return x

def add_row(A, b, pivots):
    """
    The add row procedure.

    INPUT:
        A -- an n x (n-1) matrix in Hermite normal form
        b -- an n x 1 matrix
        pivots -- sorted list of integers; the pivot positions of A.

    OUTPUT:
        H -- the Hermite normal form of A.stack(b)
    """
    t = verbose('first add row')
    global X, v
    X = A
    v = b.row(0)
    z = A._add_row_and_maintain_echelon_form(b.row(0), pivots)
    verbose('finished add row', t)
    return z

def hnf(A):
    """
    INPUT:
        an n x m matrix A over the integers.
    OUTPUT:
        the Hermite normal form of A.
    """
    n = A.nrows()
    m = A.ncols()
    if n != m:
        raise NotImplementedError, "A must be square."

    t = verbose("starting slicings")
    B = A.matrix_from_rows(range(m-2)).matrix_from_columns(range(n-1))
    c = A.matrix_from_rows([m-2]).matrix_from_columns (range(n-1))
    d = A.matrix_from_rows([m-1]).matrix_from_columns (range(n-1))
    b = A.matrix_from_columns([n-1]).matrix_from_rows(range(m-2))
    verbose("done slicing", t)

    (d1,d2) = doubleDet (B,c,d)

    (g,k,l) = d1._xgcd (d2, minimal=True)

    W = B.stack (k*c + l*d)

    H = W._hnf_mod(2*g)

    x = add_column(B.stack(k*c+l*d), H, b.stack(matrix(1,1,[k*A[m-2,m-1] + l*A[m-1,m-1]])))

    Hprime = H.augment(x)

    pivots = range(Hprime.nrows())

    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m-2]), pivots)
    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m-1]), pivots)

    H = Hprime.matrix_from_rows(range(m))
    return H




def benchmark(nrange, bits=4):
    from sage.misc.misc import cputime
    b = 2**bits
    for n in nrange:
        a = random_matrix(ZZ, n, x=-b,y=b)
        t = cputime()
        h = hnf(a)
        tm = cputime(t)
        print '%s,'%(('sage', n, bits, tm),)

def benchmark_magma(nrange, bits=4):
    from sage.misc.misc import cputime
    b = 2**bits
    for n in nrange:
        a = magma('MatrixAlgebra(IntegerRing(),%s)![Random(%s,%s) : i in [1..%s]]'%(n,-b,b,n**2))
        t = magma.cputime()
        h = a.EchelonForm()
        tm = magma.cputime(t)
        print '%s,'%(('magma', n, bits, tm),)
