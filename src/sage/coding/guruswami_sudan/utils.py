from sage.functions.other import binomial, floor, sqrt
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer

IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

def poly2list(p, len):
    """Convert the polynomial p into a list of coefficients of length len"""
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

def gs_satisfactory(tau, s, l, C = None, n_k = None):
    r"""
    Returns whether input parameters satisfy the governing equation of
    Guruswami-Sudan.

    See [N13]_ page 49, definition 3.3 and proposition 3.4 for details.

    INPUT:

    - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
      to correct
    - ``s`` -- an integer, multiplicity parameter of Guruswami-Sudan algorithm
    - ``l`` -- an integer, list size parameter
    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import gs_satisfactory as gs_sat
        sage: tau, s, l = 97, 1, 2
        sage: n, k = 250, 70
        sage: gs_sat(tau, s, l, n_k = (n, k))
        True

    One can also pass a GRS code::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: gs_sat(tau, s, l, C = C)
        True

    Another example where ``s`` and ``l`` does not satisfy the equation::

        sage: tau, s, l = 118, 47, 80
        sage: gs_sat(tau, s, l, n_k = (n, k))
        False

    If one provides both ``C`` and ``n_k`` an exception is returned::

        sage: tau, s, l = 97, 1, 2
        sage: n, k = 250, 70
        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: gs_sat(tau, s, l, C = C, n_k = (n, k))
        Traceback (most recent call last):
        ...
        ValueError: Please provide only the code or its length and dimension

    Same if one provides none of these::

        sage: gs_sat(tau, s, l)
        Traceback (most recent call last):
        ...
        ValueError: Please provide either the code or its length and dimension
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]
    return l > 0 and s > 0 and n * s * (s+1) < (l+1) * (2*s*(n-tau) - (k-1) * l)

def s_l_from_tau(tau, C = None, n_k = None):
    r"""
    Returns suitable ``s`` and ``l`` according to input parameters.

    See [N13]_ pages 53-54, proposition 3.11 for details.

    INPUT:

    - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
      to correct
    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`

    OUTPUT:

    - ``(s, l)`` -- a couple of integers, where:
        - ``s`` is the multiplicity parameter of Guruswami-Sudan algorithm and
        - ``l`` is the list size parameter

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import s_l_from_tau as s_l
        sage: tau = 97
        sage: n, k = 250, 70
        sage: s_l(tau, n_k = (n, k))
        (2, 3)

    Same one with a GRS code::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: s_l(tau, C = C)
        (2, 3)

    Another one with a bigger ``tau``::

        sage: s_l(118, C = C)
        (47, 89)
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]

    w = k - 1
    atau = n - tau
    smin = tau * w / (atau ** 2 - n * w)
    s = floor(1 + smin)
    D = (s - smin) * (atau ** 2 - n * w) * s + (w**2) /4
    l = floor(atau / w * s + 0.5 - sqrt(D)/w)
    assert gs_satisfactory(tau,s,l, n_k = (n, k)) , IMPOSSIBLE_PARAMS
    return (s, l)

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


