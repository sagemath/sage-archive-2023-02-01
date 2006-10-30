import weakref

cache = {}

def PolynomialRing(base_ring, arg1=None, arg2=None,
                   sparse=False, order='degrevlex'):
    r"""
    Return the globally unique univariate or multivariate polynomial
    ring with given properties and variable name or names.

    There are three ways to call the polynomial ring constructor:
          1. PolynomialRing(base_ring, name,    sparse=False)
          2. PolynomialRing(base_ring, names,   order='degrevlex')
          3. PolynomialRing(base_ring, name, n, order='degrevlex')
          4. PolynomialRing(base_ring, n, name, order='degrevlex')

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:
         base_ring -- a commutative ring
         name -- a string
         names -- a list or tuple of names, or a comma separated string
         n -- an integer
         sparse -- bool (default: False), whether or not elements are sparse
         order -- string or TermOrder, e.g.,
                 'degrevlex' (default) -- degree reverse lexicographic
                 'revlex' -- reverse lexicographic
                 'lex'  -- lexicographic
                 'deglex' -- degree lexicographic
                 'wp(w1,...,wn)' -- weight reverse lexicographic
                 'Wp(w1,...,wn)' -- weight lexicographic

    OUTPUT:
        PolynomialRing(base_ring, name, sparse=False) returns a univariate
        polynomial ring; all other input formats return a multivariate
        polynomial ring.

    UNIQUENESS and IMMUTABILITY: In SAGE there is exactly one
    single-variate polynomial ring over each base ring in each choice
    of variable and sparsenes.  There is also exactly one multivariate
    polynomial ring over each base ring for each choice of names of
    variables and term order.  The names of the generators cannot be
    changed after the ring has been created (immutability).

    SQUARE BRACKETS NOTATION: You can alternatively create a single or
    multivariate polynomial ring over a ring $R$ by writing
    \code{R['varname']} or \code{R['var1,var2,var3,...']}.  This
    square brackets notation doesn't allow for setting any of the
    optional arguments.  Also, it must *never* be used in library
    code, since it injects variables into the global scope.

    EXAMPLES:
        sage: PolynomialRing(ZZ, 'y')
        Univariate Polynomial Ring in y over Integer Ring
        sage: PolynomialRing(QQ['z'], 'y')
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in z over Rational Field
        sage: PolynomialRing(QQ, 'abc')
        Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, 'abc', sparse=True)
        Sparse Univariate Polynomial Ring in abc over Rational Field
        sage: PolynomialRing(QQ, 'y', 3, sparse=True)
        Polynomial Ring in y0, y1, y2 over Rational Field
    """
    import polynomial_ring as m

    if isinstance(arg1, (int, long, m.integer.Integer)):
        arg1, arg2 = arg2, arg1

    if not m.ring.is_Ring(base_ring):
        raise TypeError, 'base_ring must be a ring'

    R = None
    if isinstance(arg2, (int, long, m.integer.Integer)):
        # 3. PolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError, "You *must* specify the names of the variables (as a list of strings or a comma-seperated string)."
        n = int(arg2)
        names = arg1
        R = multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str):
        if not ',' in arg1:
            # 1. PolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError, "if second arguments is a string with no commas, then there must be no other non-optional arguments"
            name = arg1
            R = single_variate(base_ring, name, sparse)
        else:
            # 2-4. PolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
            names = arg1.split(',')
            n = len(names)
            R = multi_variate(base_ring, names, n, sparse, order)
    elif isinstance(arg1, (list, tuple)):
            # PolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
            names = arg1
            n = len(names)
            R = multi_variate(base_ring, names, n, sparse, order)

    if R is None:
        raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"

    return R


def _get_from_cache(key):
    try:
        if cache.has_key(key):
            return cache[key]()
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)
    return None

def _save_in_cache(key, R):
    try:
        cache[key] = weakref.ref(R)
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)


cdef single_variate(base_ring, name, sparse):
    import polynomial_ring as m
    key = (base_ring, name, sparse)
    R = _get_from_cache(key)
    if not R is None: return R

    if m.integer_mod_ring.is_IntegerModRing(base_ring) and not sparse:
        n = base_ring.order()
        if n.is_prime():
            R = m.PolynomialRing_dense_mod_p(base_ring, name)
        else:
            R = m.PolynomialRing_dense_mod_n(base_ring, name)

    elif base_ring.is_field():
        R = m.PolynomialRing_field(base_ring, name, sparse)

    elif base_ring.is_integral_domain():
        R = m.PolynomialRing_integral_domain(base_ring, name, sparse)
    else:
        R = m.PolynomialRing_generic(base_ring, name, sparse)

    _save_in_cache(key, R)
    return R

cdef multi_variate(base_ring, names, n, sparse, order):
    import multi_polynomial_ring as m

    order = m.TermOrder(order)

    if isinstance(names, list):
        names = tuple(names)

    elif isinstance(names, str):
        if len(names) > 1:
            names = tuple(names)

    key = (base_ring, names, n, sparse, order)
    R = _get_from_cache(key)
    if not R is None:
        return R

    if m.integral_domain.is_IntegralDomain(base_ring):
        R = m.MPolynomialRing_polydict_domain(base_ring, n, names, order)
    else:
        R = m.MPolynomialRing_polydict(base_ring, n, names, order)

    _save_in_cache(key, R)
    return R


