#########################################################################################
# Factory function for making polynomial rings
#########################################################################################

from sage.structure.parent_gens import normalize_names
import ring
import weakref
import padics.padic_ring_generic
import padics.padic_field_generic
import padics.padic_ring_lazy
import padics.padic_field_lazy

_cache = {}

def PolynomialRing(base_ring, arg1=None, arg2=None,
                   sparse=False, order='degrevlex',
                   names=None, name=None):
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
    variables and term order.  The names of the generators can only
    be temporarily changed after the ring has been created.  Do this
    using the localvars context:

        EXAMPLES of VARIABLE NAME CONTEXT:
            sage: R.<x,y> = PolynomialRing(QQ,2); R
            Polynomial Ring in x, y over Rational Field
            sage: f = x^2 - 2*y^2

        You can't just globally change the names of those variables.
        This is because objects all over SAGE could have pointers to
        that polynomial ring.
            sage: R._assign_names(['z','w'])
            Traceback (most recent call last):
            ...
            ValueError: variable names cannot be changed after object creation.

        However, you can very easily change the names within a "with" block:
            sage: with localvars(R, ['z','w']):
            ...     print f
            ...
            -2*w^2 + z^2

        After the with block the names revert to what they were before.
            sage: print f
            -2*y^2 + x^2


    SQUARE BRACKETS NOTATION: You can alternatively create a single or
    multivariate polynomial ring over a ring $R$ by writing
    \code{R['varname']} or \code{R['var1,var2,var3,...']}.  This
    square brackets notation doesn't allow for setting any of the
    optional arguments.

    EXAMPLES:
    1. PolynomialRing(base_ring, name,    sparse=False):
        sage: PolynomialRing(QQ, 'w')
        Univariate Polynomial Ring in w over Rational Field

    Use the diamond brackets notation to make the variable
    ready for use after you define the ring:
        sage: R.<w> = PolynomialRing(QQ)
        sage: (1 + w)^3
        w^3 + 3*w^2 + 3*w + 1

    You must specify a name:
        sage: PolynomialRing(QQ)
        Traceback (most recent call last):
        ...
        TypeError: You must specify the names of the variables.

        sage: R.<abc> = PolynomialRing(QQ, sparse=True); R
        Sparse Univariate Polynomial Ring in abc over Rational Field

        sage: R.<w> = PolynomialRing(PolynomialRing(GF(7),'k')); R
        Univariate Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

    The square bracket notation:
        sage: R.<y> = QQ['y']; R
        Univariate Polynomial Ring in y over Rational Field
        sage: y^2 + y
        y^2 + y

    In fact, since the diamond brackets on the left determine the variable name, you can omit the
    variable from the square brackets:
        sage: R.<zz> = QQ[ ]; R
        Univariate Polynomial Ring in zz over Rational Field
        sage: (zz + 1)^2
        zz^2 + 2*zz + 1

    This is exactly the same ring as what PolynomialRing returns:
        sage: R is PolynomialRing(QQ,'zz')
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

        sage: R = PolynomialRing(ZZ, ['x%s'%p for p in primes(100)]); R
        Polynomial Ring in x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97 over Integer Ring

    By calling the \code{inject_variables()} method all those variable
    names are available for interactive use:
        sage: R.inject_variables()
        Defining x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97
        sage: (x2 + x41 + x71)^2
        x71^2 + 2*x41*x71 + x41^2 + 2*x2*x71 + 2*x2*x41 + x2^2

    You can also call \code{injvar}, which is a convenient shortcut for \code{inject_variables()}.
        sage: R = PolynomialRing(GF(7),15,'w'); R
        Polynomial Ring in w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14 over Finite Field of size 7
        sage: R.injvar()
        Defining w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14
        sage: (w0 + 2*w8 + w13)^2
        w13^2 + 4*w8*w13 + 4*w8^2 + 2*w0*w13 + 4*w0*w8 + w0^2
    """
    import polynomial_ring as m
    if isinstance(arg1, (int, long, m.integer.Integer)):
        arg1, arg2 = arg2, arg1

    if not names is None:
        arg1 = names
    elif not name is None:
        arg1 = name

    if not m.ring.is_Ring(base_ring):
        raise TypeError, 'base_ring must be a ring'

    if arg1 is None:
        raise TypeError, "You must specify the names of the variables."

    R = None
    if isinstance(arg2, (int, long, m.integer.Integer)):
        # 3. PolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError, "You *must* specify the names of the variables."
        n = int(arg2)
        names = arg1
        R = _multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str) or len(arg1) == 1:
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
        raise TypeError, "you *must* specify the indeterminates (as not None)."
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
    name = normalize_names(1, name)
    key = (base_ring, name, sparse)
    R = _get_from_cache(key)
    if not R is None: return R

    if isinstance(base_ring, ring.CommutativeRing):
        if m.integer_mod_ring.is_IntegerModRing(base_ring) and not sparse:
            n = base_ring.order()
            if n.is_prime():
                R = m.PolynomialRing_dense_mod_p(base_ring, name)
            elif n > 1:
                R = m.PolynomialRing_dense_mod_n(base_ring, name)
            else:  # n == 1!
                R = m.PolynomialRing_integral_domain(base_ring, name)   # specialized code breaks in this case.

        elif isinstance(base_ring, padics.padic_ring_lazy.pAdicRingLazy):
            R = m.PolynomialRing_padic_ring_lazy(base_ring, name)

        elif isinstance(base_ring, padics.padic_field_lazy.pAdicFieldLazy):
            R = m.PolynomialRing_padic_field_lazy(base_ring, name)

        elif isinstance(base_ring, padics.padic_field_generic.pAdicFieldGeneric):
            R = m.PolynomialRing_padic_field(base_ring, name)

        elif isinstance(base_ring, padics.padic_ring_generic.pAdicRingGeneric):
            R = m.PolynomialRing_padic_ring(base_ring, name)

        elif base_ring.is_field():
            R = m.PolynomialRing_field(base_ring, name, sparse)

        elif base_ring.is_integral_domain():
            R = m.PolynomialRing_integral_domain(base_ring, name, sparse)
        else:
            R = m.PolynomialRing_commutative(base_ring, name, sparse)
    else:
        R = m.PolynomialRing_general(base_ring, name, sparse)

    _save_in_cache(key, R)
    return R

def _multi_variate(base_ring, names, n, sparse, order):
    names = normalize_names(n, names)

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

