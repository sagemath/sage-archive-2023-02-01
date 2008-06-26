include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"

include "../../libs/flint/flint.pxi"
include "../../libs/flint/fmpz_poly.pxi"

from sage.rings.integer cimport Integer

## from sage.rings.all import QQ
## R = QQ['X']
## def apply_to_monomial_0(int i, int j, int a, int b, int c, int d):
##     f = R([b,a])
##     g = R([d,c])
##     if i < 0 or j-i < 0:
##         raise ValueError, "i (=%s) and j-i (=%s) must both be nonnegative."%(i,j-i)
##     h = (f**i)*(g**(j-i))
##     return [int(a) for a in h.padded_list(j+1)]

def apply_to_monomial(int i, int j, int a, int b, int c, int d):
    """
    Returns a list of the coefficients of
    $$
             (aX + bY)^i (cX + dY)^{j-i},
    $$
    where $0 \leq i \leq j$, and $a,b,c,d$ are integers.

    One should think of $j$ as being $k-2$ for the application to
    modular symbols.

    INPUT:
        i, j, a, b, c, d -- all ints

    OUTPUT:
        list of ints, which are the coefficients
        of $Y^j$, $Y^{j-1}*X$, \ldots, $X^j$, respectively.

    EXAMPLE:
    We compute that $(X+Y)^2(X-Y) = X^3 + X^2Y - XY^2 - Y^3$.
        sage: from sage.modular.modsym.manin_symbols import apply_to_monomial
        sage: apply_to_monomial(2, 3, 1,1,1,-1)
        [-1, -1, 1, 1]
        sage: apply_to_monomial(5, 8, 1,2,3,4)
        [2048, 9728, 20096, 23584, 17200, 7984, 2304, 378, 27]
        sage: apply_to_monomial(6,12, 1,1,1,-1)
        [1, 0, -6, 0, 15, 0, -20, 0, 15, 0, -6, 0, 1]
    """
    cdef Py_ssize_t k
    if i < 0 or j-i < 0:
        raise ValueError, "i (=%s) and j-i (=%s) must both be nonnegative."%(i,j-i)

    # initialize flint polys
    cdef fmpz_poly_t f, g, ff, gg, pr
    fmpz_poly_init(f)
    fmpz_poly_init(g)
    fmpz_poly_init(ff)
    fmpz_poly_init(gg)
    fmpz_poly_init(pr)

    # f = b+a*x, g = d+c*x
    fmpz_poly_set_coeff_si(f, 0, b)
    fmpz_poly_set_coeff_si(f, 1, a)
    fmpz_poly_set_coeff_si(g, 0, d)
    fmpz_poly_set_coeff_si(g, 1, c)

    # h = (f**i)*(g**(j-i))
    fmpz_poly_power(ff, f, i)
    fmpz_poly_power(gg, g, j-i)
    fmpz_poly_mul(pr, ff, gg)

    cdef Integer res
    v = []
    for k from 0 <= k <= j:
        res = <Integer>PY_NEW(Integer)
        fmpz_poly_get_coeff_mpz(res.value, pr, k)
        v.append(int(res))

    # clear flint polys
    fmpz_poly_clear(f)
    fmpz_poly_clear(g)
    fmpz_poly_clear(ff)
    fmpz_poly_clear(gg)
    fmpz_poly_clear(pr)
    return v
