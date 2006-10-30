import weakref

cache = {}

def FiniteField(order, name=None, modulus=None):
    """
    Return the globally unique finite field of given order with generator
    labeled by the given name and possibly with given modulus.

    INPUT:
        order --   int
        name --    string; must be specified in not a prime field
        modulus -- (optional) defining polynomial for field, i.e.,
                   generator of the field will be a root of this
                   polynomial; if not specified the choice of
                   definining polynomials can be arbitrary.

    EXAMPLES:
        sage: k = FiniteField(9, 'a'); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

    You can also use GF instead of FiniteField -- they are identical.
    """
    order = int(order)

    key = (order, name, modulus)
    if cache.has_key(key):
        K = cache[key]()
        if not K is None:
            if inject_variable:
                K.inject_variables()
            return K

    if finite_field.arith.is_prime(order):
        K = finite_field.integer_mod_ring.IntegerModRing(order)
    else:
        if name is None:
            raise TypeError, "you must specify the generator name"
        if False and order < 2**16:   # todo -- re-enable
            K = finite_field.FiniteField_givaro(order, name, modulus)
        else:
            K = finite_field.FiniteField_ext_pari(order, name, modulus)

    cache[key] = weakref.ref(K)
    if inject_variable:
        K.inject_variables()
    return K


def is_FiniteField(x):
    return isinstance(x, finite_field.FiniteField_generic)

def is_PrimeFiniteField(x):
    return isinstance(x, finite_field.FiniteField_prime_modn)

