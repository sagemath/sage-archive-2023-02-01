import weakref

cache = {}

def MPolynomialRing(base_ring, names, n=1, order='degrevlex', inject_variables=True):
    r"""
    Create a Multivariate polynomial ring over a commutative base ring.

    INPUT:
        base_ring -- CommutativeRing
        names -- list or tuple or string:
                   - list or tuple of n variable names
                   - if string, names the variables the characters in the string, or
                     if there are commas in the string by the variables between commas.
        n -- int, number of variables  (default: 1)
        order -- string; the term order, or an object of type TermOrder:
                 'degrevlex' (default) -- degree reverse lexicographic
                 'revlex' -- reverse lexicographic
                 'lex'  -- lexicographic
                 'deglex' -- degree lexicographic
                 'wp(w1,...,wn)' -- weight reverse lexicographic
                 'Wp(w1,...,wn)' -- weight lexicographic

    NOTE: For convenience names and n may be switched in the input.

    EXAMPLES:
        sage: R = MPolynomialRing(RationalField(), 'x', 3); R
        Polynomial Ring in x0, x1, x2 over Rational Field
        sage: x0,x1,x2 = R.gens()
        sage: x0.element()
        PolyDict with representation {(1, 0, 0): 1}
        sage: x0 + x1 + x2
        x2 + x1 + x0
        sage: (x0 + x1 + x2)**2
        x2^2 + 2*x1*x2 + x1^2 + 2*x0*x2 + 2*x0*x1 + x0^2

    This example illustrates the quick shorthand for naming several
    variables one-letter names.
        sage: MPolynomialRing(ZZ, 'xyzw', 4)
        Polynomial Ring in x, y, z, w over Integer Ring

    To obtain both the ring and its generators, use the \code{objgens} function.
        sage: R, (x,y,z,w) = MPolynomialRing(ZZ, 4, 'xyzw').objgens()
        sage: (x+y+z+w)^2
        w^2 + 2*z*w + z^2 + 2*y*w + 2*y*z + y^2 + 2*x*w + 2*x*z + 2*x*y + x^2

    We can construct multi-variate polynomials rings over completely
    arbitrary SAGE rings.  In this example, we construct a polynomial
    ring S in 3 variables over a polynomial ring in 2 variables over
    GF(9).  Then we construct a polynomial ring in 20 variables over S!

        sage: R, (n1,n2) = MPolynomialRing(GF(9),2,['n1','n2']).objgens()
        sage: R = GF(9)['n1, n2']
        sage: n1^2 + 2*n2
        2*n2 + n1^2
        sage: S = MPolynomialRing(R,3, names='a'); a0,a1,a2=S.gens()
        sage: S
        Polynomial Ring in a0, a1, a2 over Polynomial Ring in n1, n2 over Finite Field in a of size 3^2
        sage: x = (n1+n2)*a0 + 2*a1**2
        sage: x
        2*a1^2 + (n2 + n1)*a0
        sage: x**3
        2*a1^6 + (n2^3 + n1^3)*a0^3
        sage: T = MPolynomialRing(S, 20)
        sage: T
        Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19 over Polynomial Ring in a0, a1, a2 over Polynomial Ring in n1, n2 over Finite Field in a of size 3^2

    """
    import multi_polynomial_ring as m

    if isinstance(names, (int, long, m.Integer)):
        # swap them.
        names, n = n, names

    T = m.TermOrder(order)
    if isinstance(names, list):
        names = tuple(names)

    elif isinstance(names, str):
        if len(names) > 1:
            names = tuple(names)

    key = (base_ring, n, names, T)
    if cache.has_key(key):
        R = cache[key]()
        if not (R is None):
            if inject_variables:
                R.inject_variables()
            return R

    if not isinstance(base_ring, m.commutative_ring.CommutativeRing):
        raise TypeError, "Base ring must be a commutative ring."

    if m.integral_domain.is_IntegralDomain(base_ring):
        R = m.MPolynomialRing_polydict_domain(base_ring, n, names, T)
    else:
        R = m.MPolynomialRing_polydict(base_ring, n, names, T)

    cache[key] = weakref.ref(R)
    if inject_variables:
        R.inject_variables()
    return R

def is_MPolynomialRing(x):
    import multi_polynomial_ring as m
    return isinstance(x, m.MPolynomialRing_generic)

