from sage.functions.other import binomial, floor, sqrt
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer

IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

def poly2list(p, len):
    r"""
    Returns ``p`` as a list of its coefficients of length ``len``.

    INPUT:

    - ``p`` -- a polynomial

    - ``len`` -- an integer

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import poly2list
        sage: F.<x> = GF(41)[]
        sage: p = 9*x^2 + 8*x + 37
        sage: poly2list(p, 4)
        [37, 8, 9, 0]
    """
    return list(p) + [0]*max(0, len-p.degree()-1)

def list_decoding_range(n, d, q=None):
    r"""
    Returns the minimal and maximal number of errors correctable by a
    Johnson-distance list decoder beyond half the minimal distance.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code
    - ``q`` -- (default: ``None``) ????????

    EXAMPLES::

        sage: from sage.coding.grs import list_decoding_range
        sage: list_decoding_range(250, 181)
        (91, 118)
    """
    if q is None:
        return (ligt((d-1)/2), gilt(n - sqrt(n*(n-d))))
    else:
        f = (q-1.)/q
        return (ligt((d-1)/2), gilt(f*(n-sqrt(n*(n-d/f)))))

def ligt(x):
    r"""
    Returns the least integer greater than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import ligt
        sage: ligt(41)
        42

    It works with any type of numbers (not only integers)::

        sage: ligt(41.041)
        42
    """
    return floor(x+1)

def gilt(x):
    r"""
    Returns the greatest integer smaller than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import gilt
        sage: gilt(43)
        42

    It works with any type of numbers (not only integers)::

        sage: gilt(43.041)
        43
    """
    if x in ZZ:
        return Integer(x-1)
    else:
        return floor(x)

def solve2deg_int(a,b,c):
    r"""
    Returns the greatest integer range `[i1, i2]` such that
    `i1 > x1` and `i2 < x2` where `x1`,`x2` are the two zeroes of the equation in `x`:
    `ax^2+bx+c=0`.

    If there is no real solution to the equation, it returns an empty, range with negative coefficients.

    INPUT:

    - ``a``, ``b`` and ``c`` -- coefficients of a second degree equation, ``a`` being the coefficient of
      the higher degree term.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import solve2deg_int
        sage: solve2deg_int(1, -5, 1)
        (1, 4)

    If there is no real solution::

        sage: solve2deg_int(50, 5, 42)
        (-2, -1)
    """
    D = b**2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx, maxx = (-b-sD)/2.0/a , (-b+sD)/2.0/a
    mini, maxi = (ligt(minx), gilt(maxx))
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def find_minimal_satisfiable(f, startn=1, contiguous=True):
    r"""
    Returns the minimal integral ``n``, `n > 0` such that ``f(n) == True``.

    If the interval for which `f` is true is contiguous and open
    towards infinity, a logarithmic algorithm is used, otherwise linear.
    `startn` can be given as a hint to a value that might be true.

    INPUT:

    - ``f`` -- a function

    - ``startn`` -- (default: ``1``) the starting point of the algorithm.

    - ``contiguous`` -- (default: ``True``) boolean describing the contiguousity of ``f``'s
      interval

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import find_minimal_satisfiable
        sage: def f(x):
        ....:    return None if x > 10 or x == 1 else x + 1

        sage: find_minimal_satisfiable(f)
        2
    """
    if not contiguous:
        n = startn
        if f(n):
            while f(n) and n > 0:
                n = n - 1
            return n + 1
        else:
            while not f(n):
                n = n + 1
            return n
    else:
        maxn = startn
        minn = 1
        # Keep doubling n to find one that works and then binary
        while not f(maxn):
            minn = maxn + 1
            maxn *= 2
        while minn < maxn:
            tryn = minn + floor((maxn - minn) * 0.5)
            if f(tryn):
                maxn = tryn
            else:
                minn = tryn + 1
        return maxn


