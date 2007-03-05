import weakref
import sage.rings.padics.padic_field_capped_relative
import sage.rings.padics.padic_field_lazy

pAdicFieldCappedRelative = sage.rings.padics.padic_field_capped_relative.pAdicFieldCappedRelative
pAdicFieldLazy = sage.rings.padics.padic_field_lazy.pAdicFieldLazy


padic_field_cache = {}
def Qp(p, prec = 20, type = 'capped-rel', print_mode = None, halt = 40):
    """
    A creation function for p-adic fields.

    INPUT:
        p -- integer: the p in Q_p
        prec -- integer (default: 20) the precision cap of the field.  Individual elements keep track of their own precision.
        type -- string (default: 'capped-rel') see Notes
        print_mode -- string (default: None) the print mode
    OUTPUT:
        the corresponding p-adic field

    EXAMPLES:
        sage: K = Qp(5); a = K(4); a
        4 + O(5^20)
        sage: L = Qp(5, 10, type = 'lazy'); b = L(2); b
        2 + O(5^10)
        sage: a + b
        1 + 5 + O(5^20)

    NOTES:
         values of type:
         'capped-rel' -> pAdicFieldCappedRelative.  This is the default, considers precision as the precision of the unit part.  Tracks precision of individual elements, but bounds the precision of any element with a precision cap.
        'lazy' -> pAdicFieldLazy.  Uses lazy elements so that additional precision can be requested during a computation.  There is some amount of performance penalty because of this ability.

        values of print_mode:
        'val-unit' -- elements are displayed as p^k*u
        'integer' -- elements are displayed as an integer
        'series' -- elements are displayed as series in p
        'val-unit-p' -- same as val-unit, except that p is written as "p"
        'integer-p' -- same as integer, except that p is written as "p"
        'series-p' -- same as series, except that p is written as "p"
    """
    if not p.is_prime():
        raise ValueError, "p must be prime"
    if type != 'lazy':
        key = (p, prec, type)
    else:
        key = (p, prec, halt)
    if padic_field_cache.has_key(key):
        K = padic_field_cache[key]()
        if K != None:
            if print_mode != None:
                K.set_print_mode(print_mode)
            else:
                K.set_print_mode('series')
            return K
    if print_mode == None:
        print_mode = 'series'
    if (type == 'capped-rel'):
        K = pAdicFieldCappedRelative(p, prec, print_mode)
    elif (type == 'lazy'):
        K = pAdicFieldLazy(p, prec, print_mode, halt)
    else:
        raise ValueError, "type must be either 'capped-rel' or 'lazy'"
    padic_field_cache[key] = weakref.ref(K)
    return K

pAdicField = Qp # for backwards compatibility; and it's not hard.

qadic_field_cache = {}
def Qq(q, name=None, prec=20, type='capped-rel', print_mode=None, halt=40, modulus=None, check=True):
    r"""
    Given a prime power q = p^n, return the unique unramified extension
    of Qp of degree n.

    Currently, there's no code for unramified field extensions, so
    we just return the UnramifiedRingExtension.
    """

    from sage.rings.integer import Integer
    from sage.rings.polynomial_ring import PolynomialRing
    from sage.rings.padics.unramified_ring_extension import UnramifiedRingExtension
    from sage.rings.integer_ring import ZZ

    if name is None:
        raise TypeError, "You must specify the name of the generator."

    q = Integer(q)
    F = q.factor()
    if len(F) != 1:
        raise ValueError, "q must be a prime power"
    if F[0][1] == 1:
        return Qp(q, prec, type, print_mode, halt)

    if type != 'lazy':
        key = (q, name, prec, type)
    else:
        key = (q, name, prec, halt)
    if qadic_field_cache.has_key(key):
        K = qadic_field_cache[key]()
        if not (K is None):
            if not (print_mode is None):
                K.set_print_mode(print_mode)
            return K

    if modulus is None:
        check = False
        from sage.rings.finite_field import GF
        modulus = PolynomialRing(Qp(F[0][0], prec, type, print_mode, halt), name)(GF(q,name).modulus().change_ring(ZZ))
    if print_mode is None:
        print_mode = 'series'
    K = UnramifiedRingExtension(modulus, prec, print_mode, check)
    qadic_field_cache[key] = weakref.ref(K)

    return K

