#########################################################################################
# Factory function for making polynomial rings
#########################################################################################

import weakref

_cache = {}

def PolynomialRing(base_ring, arg1=None, arg2=None,
                   sparse=False, order='degrevlex'):
    r"""
    Return the globally unique univariate or multivariate polynomial
    ring with given properties and variable name or names.

    There are four ways to call the polynomial ring constructor:
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
    1. PolynomialRing(base_ring, name,    sparse=False):
        sage: PolynomialRing(QQ, 'w')
        Univariate Polynomial Ring in w over Rational Field

    In the interactive interpreter the variable is immediately
    ready for use after you define the ring:
        sage: (1 + w)^3
        w^3 + 3*w^2 + 3*w + 1

    You must specify a name:
        sage: PolynomialRing(QQ)
        Traceback (most recent call last):
        ...
        TypeError: invalid input to PolynomialRing function; please see the docstring for that function

        sage: PolynomialRing(QQ, 'abc', sparse=True)
        Sparse Univariate Polynomial Ring in abc over Rational Field

        sage: PolynomialRing(PolynomialRing(GF(7),'k'), 'w')
        Univariate Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

    The square bracket notation:
        sage: R = QQ['y']; R
        Univariate Polynomial Ring in y over Rational Field
        sage: y^2 + y
        y^2 + y

    This is exactly the same ring as what PolynomialRing returns:
        sage: R is PolynomialRing(QQ,'y')
        True

    However, rings with different variables are different:
        sage: QQ['x'] == QQ['y']
        False

    2. PolynomialRing(base_ring, names,   order='degrevlex'):
        sage: R = PolynomialRing(QQ, 'a,b,c'); R
        Polynomial Ring in a, b, c over Rational Field

        sage: S = PolynomialRing(QQ, ['a','b','c']); S
        Polynomial Ring in a, b, c over Rational Field

        sage: T = PolynomialRing(QQ, ('a','b','c')); T
        Polynomial Ring in a, b, c over Rational Field

    All three rings are identical.
        sage: (R is S) and  (S is T)
        True

    There is a unique polynomial ring with each term order:
        sage: R = PolynomialRing(QQ, 'x,y,z', order='degrevlex'); R
        Polynomial Ring in x, y, z over Rational Field
        sage: S = PolynomialRing(QQ, 'x,y,z', order='revlex'); S
        Polynomial Ring in x, y, z over Rational Field
        sage: S is PolynomialRing(QQ, 'x,y,z', order='revlex')
        True
        sage: R == S
        False


    3. PolynomialRing(base_ring, name, n, order='degrevlex'):

    If you specify a single name as a string and a number of
    variables, then variables labeled with numbers are created.
        sage: PolynomialRing(QQ, 'x', 10)
        Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 over Rational Field

        sage: PolynomialRing(GF(7), 'y', 5)
        Polynomial Ring in y0, y1, y2, y3, y4 over Finite Field of size 7

        sage: PolynomialRing(QQ, 'y', 3, sparse=True)
        Polynomial Ring in y0, y1, y2 over Rational Field

    It is easy in Python to create fairly aribtrary variable names.
    For example, here is a ring with generators labeled by the first
    100 primes:

        sage: PolynomialRing(ZZ, ['x%s'%p for p in primes(100)])
        Polynomial Ring in x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97 over Integer Ring
        sage: (x2 + x41 + x71)^2
        x71^2 + 2*x41*x71 + x41^2 + 2*x2*x71 + 2*x2*x41 + x2^2
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
        R = _multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str):
        if not ',' in arg1:
            # 1. PolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError, "if second arguments is a string with no commas, then there must be no other non-optional arguments"
            name = arg1
            R = _single_variate(base_ring, name, sparse)
        else:
            # 2-4. PolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError, "invalid input to PolynomialRing function; please see the docstring for that function"
            names = arg1.split(',')
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)
    elif isinstance(arg1, (list, tuple)):
            # PolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
            names = arg1
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)

    if arg1 is None and arg2 is None:
        raise TypeError, "you *must* specify the indeterminates."
    if R is None:
        raise TypeError, "invalid input (%s, %s, %s) to PolynomialRing function; please see the docstring for that function"%(
            base_ring, arg1, arg2)

    return R


def _get_from_cache(key):
    try:
        if _cache.has_key(key):
            return _cache[key]()
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)
    return None

def _save_in_cache(key, R):
    try:
        _cache[key] = weakref.ref(R)
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)


def _single_variate(base_ring, name, sparse):
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

def _multi_variate(base_ring, names, n, sparse, order):
    import multi_polynomial_ring as m

    order = m.TermOrder(order)

    if isinstance(names, list):
        names = tuple(names)

    elif isinstance(names, str):
        if ',' in names:
            names = tuple(names.split(','))

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

#########################################################################################
# END (Factory function for making polynomial rings)
#########################################################################################

