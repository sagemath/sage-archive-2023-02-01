from sage.matrix.constructor import matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

@cached_function
def compute_M(p, t, A):
    r"""
    Find generators of set `M_t(A)=\{f\in\mathbb{Z}[X]^d\mid Af
    \equiv 0\pmod{p^t}\}` as a `\mathbb{Z}[X]`-module.

    INPUT:

    - `p` -- a prime number

    - `t` -- a non-negative integer

    - `A` -- an immutable matrix over `\mathbb{Z}[X]`

    OUTPUT:

    A matrix `F`. The columns of `(p^tI F)` are generators of `M_t(A)`.

    EXAMPLES::

        sage: X = polygen(ZZ)
        sage: A = matrix([[X-1, X-2, X], [X-1, X-3, X+2]])
        sage: A.set_immutable()
        sage: for t in range(3):
        ....:     print compute_M(3, t, A)
    """

    if t == 0:
        return matrix(A.parent().base(), A.ncols(), 0)

    P = A.parent()
    ZZX = P.base()
    (X,) = ZZX.gens()
    d = A.ncols()

    G = compute_M(p, t-1, A)
    R = A*G/p**(t-1)
    R.change_ring(ZZX)

    AR = matrix.block([[A, R]])
    Fp = GF(p)
    FpX = PolynomialRing(Fp, name=X)

    ARb = AR.change_ring(FpX)
    (Db, Sb, Tb) = ARb.smith_form()
    assert Sb * ARb * Tb == Db
    assert all(i == j or Db[i, j].is_zero()
               for i in range(Db.nrows())
               for j in range(Db.ncols()))

    r = Db.rank()
    T = Tb.change_ring(ZZX)

    F1 = matrix.block([[p**(t-1) * matrix.identity(d), G]])*T
    F = F1.matrix_from_columns(range(r, F1.ncols()))
    assert all((entry % (p**t))== 0
               for row in A*F
               for entry in row), "A*F=%s" % str(A*F)
    return matrix.block([[F, p*G]])
