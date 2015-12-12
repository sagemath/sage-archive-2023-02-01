"""
Calculate nu
"""

from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject

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

    A matrix `F`. The columns of `\begin{pmatrix}p^tI& F\end{pmatrix}` are generators of `M_t(A)`.

    EXAMPLES::

        sage: from calculate_nu import compute_M # not tested
        sage: X = polygen(ZZ)
        sage: A = matrix([[X-1, X-2, X], [X-1, X-3, X+2]])
        sage: A.set_immutable()
        sage: for t in range(3):
        ....:     print compute_M(3, t, A)
        []
        [      2|]
        [  x + 2|]
        [2*x + 1|]
        [              6         6*x + 2|              6]
        [        3*x + 6 6*x^2 + 4*x + 8|        3*x + 6]
        [        6*x + 3 3*x^2 + 2*x + 4|        6*x + 3]
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
    assert (A*F % (p**t)).is_zero(), "A*F=%s" % str(A*F)
    return matrix.block([[F, p*G]])


class Compute_nu(SageObject):
    r"""
    Compute nu for a given matrix.

    INPUT:

    - ``B`` -- an immutable matrix over `\mathbb{Z}`.

    OUTPUT:

    An object which can be called with integer arguments yielding the
    actual nu values.

    EXAMPLES::

        sage: from calculate_nu import compute_nu # not tested
        sage: B = matrix([[1, 2], [3, 4]])
        sage: C = Compute_nu(B)
        sage: for t in range(3): # not implemented
        ....:     print C(t)
    """
    def __init__(self, B):
        super(Compute_nu, self).__init__()
        if not B.is_square():
            raise TypeError("square matrix required.")

        self._B = B
        X = polygen(B.base_ring())
        adjunct = (X-B).adjoint()
        d = B.nrows()**2
        b = matrix(d, 1, adjunct.list())
        chi_B = B.charpoly(X)
        self._A = matrix.block([[b, chi_B*matrix.identity(d)]])
        self._A.set_immutable()
        self._ZX = X.parent()
        # we also need the multivariate polynomial ring in one variable
        # because Groebner bases are only implemented there.
        self._ZXm = PolynomialRing(X.base_ring(), X.variable_name(), 1)

    def ideal(self, p, t):
        r"""
        Return the ideal `N_t(B)=\{ f\in \mathbb{Z}[X] \mid \exists M\in\mathbb{Z}^{n\times n}\colon f \operatorname{adj}(X-B) \equiv
        \chi_B M \pmod{p^t}\}`.

        INPUT:

        - ``t`` -- a nonnegative integer.

        OUTPUT:

        An ideal in `\mathbb{Z}[X]`.

        EXAMPLES::

            sage: from calculate_nu import compute_nu # not tested
            sage: B = matrix([[1, 2], [3, 4]])
            sage: C = Compute_nu(B)
            sage: for t in range(3):
            ....:     print C.ideal(3, t)
            Principal ideal (1) of Univariate Polynomial Ring in x over Integer Ring
            Ideal (x^2 + x + 1, 3) of Univariate Polynomial Ring in x over Integer Ring
            Ideal (x^2 + 4*x + 7, 9) of Univariate Polynomial Ring in x over Integer Ring
        """
        M = compute_M(p, t, self._A)
        generators = [self._ZX(f) for f in [p**t] + list(M.row(0))]
        assert all((g(self._B) % p**t).is_zero()
                   for g in generators)
        # switch to multivariate polynomial ring in one variable
        # for Groebner basis computation
        I = self._ZXm.ideal(generators)
        GB = I.groebner_basis()
        return self._ZX.ideal(*GB)
