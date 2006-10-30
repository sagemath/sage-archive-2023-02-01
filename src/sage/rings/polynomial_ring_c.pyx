import weakref

cache = {}

def PolynomialRing(base_ring, arg1=None, arg2=None, sparse=False,
                   order='degrevlex', scope=None):
    r"""
    Return the globally unique univariate or multivariate polynomial
    ring with given properties and variable name or names.

    In SAGE there is exactly one single-variate polynomial ring over
    each base ring in each choice of variable and sparsenes.  There is
    also exactly one multivariate polynomial ring over each base ring
    for each choice of names of variables and term order.

    There are three ways to call the polynomial ring constructor
    function.  The optional arguments sparse, order, and scope *must*
    be explicitly named.

    VARIABLE INJECTION: By default, if this function is called from
    the interactive interpreter or notebook, then the optional scope
    argumetn is set to the global scope of the interpreter.  If this
    function is called from Python code, scope is None, and the
    indeterminates are *not* injected into the module-scope namespace.

    You can alternatively create a single or multivariate polynomial
    ring over a ring $R$ by writing \code{R['varname']} or
    \code{R['var1,var2,var3,...']}.  This square brackets notation
    doesn't allow for setting any of the optional arguments.  Also, it
    must *never* be used in library code, since it injects variables
    into the global scope.

    1. PolynomialRing(base_ring, name, sparse=False):
        INPUT:
            base_ring -- the base ring
            name -- (str) the name of the generator
            sparse -- (bool; default: False) whether or not elements are sparse.
        OUTPUT:
            A univariate polynomial ring.

    2. PolynomialRing(base_ring, names, order='degrevlex', inject_variables=True):
        INPUT:
            base_ring -- the base ring
            names -- list of strings, tuple of strings, or comma separaed string
                     with at least one comma
            order -- string or term order (default: 'degrevlex'); see docstring
                     for MPolynomialRing.
            inject_variables -- (default: True) whether or not to inject variables
                     into current scope.
        OUTPUT:
            A multivariate polynomoial ring

    3. PolynomialRing(base_ring, name, n, order='degrevlex', inject_variables=True):
        INPUT:
            base_ring -- the base ring
            name -- single string
            n -- an integer, the number of variables
            order -- string or term order (default: 'degrevlex'); see docstring
                     for MPolynomialRing.
            inject_variables -- (default: True) whether or not to inject
                                variables into current scope.
        OUTPUT:
            A multivariate polynomoial ring with variables named, e.g., x0, x1,
            x2, etc., if name="x".

        NOTE: You can also call this n and the variable name switched.
           PolynomialRing(base_ring, n, name, order='degrevlex', inject_variables=True):

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
    if isinstance(arg1, str):
        if not ',' in arg1:
            # 1. PolynomialRing(base_ring, name, sparse=False, inject_variables=True):
            if not arg2 is None:
                raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
            name = arg1
            R = single_variate(base_ring, name, sparse, scope)
        else:
            # 2. PolynomialRing(base_ring, names, order='degrevlex', inject_variables=True):
            if not arg2 is None:
                raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
            names = arg1.split(',')
            n = len(names)
            R = multi_variate(base_ring, names, n, sparse, order, scope)

    elif isinstance(arg2, (int, long, m.integer.Integer)):
        # 3. PolynomialRing(base_ring, names, n, order='degrevlex', inject_variables=True):
        if not isinstance(arg1, str):
            raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
        n = int(arg2)
        names = arg1
        R = multi_variate(base_ring, names, n, sparse, order, scope)

    if R is None:
        raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"

    if scope:
        R.inject_variables(scope)

    return R


def _get_from_cache(key, scope):
    if cache.has_key(key):
        return cache[key]()
    return None

def _save_in_cache(key, R):
    cache[key] = R


cdef single_variate(base_ring, name, sparse, scope):
    import polynomial_ring as m
    key = (base_ring, name, sparse)
    R = _get_from_cache(key, scope)
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

cdef multi_variate(base_ring, names, n, sparse, order, scope):
    import multi_polynomial_ring as m

    order = m.TermOrder(order)

    if isinstance(names, list):
        names = tuple(names)

    elif isinstance(names, str):
        if len(names) > 1:
            names = tuple(names)

    key = (base_ring, names, n, sparse, order, scope)
    R = _get_from_cache(key, scope)
    if not R is None:
        return R

    if m.integral_domain.is_IntegralDomain(base_ring):
        R = m.MPolynomialRing_polydict_domain(base_ring, n, names, order)
    else:
        R = m.MPolynomialRing_polydict(base_ring, n, names, order)

    _save_in_cache(key, R)
    return R


