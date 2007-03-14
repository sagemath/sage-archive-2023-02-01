import weakref
import sage.rings.padics.padic_field_capped_relative
import sage.rings.padics.padic_field_lazy

Integer = sage.rings.integer.Integer
pAdicFieldCappedRelative = sage.rings.padics.padic_field_capped_relative.pAdicFieldCappedRelative
pAdicFieldLazy = sage.rings.padics.padic_field_lazy.pAdicFieldLazy


padic_field_cache = {}
def Qp(p, prec = 20, type = 'capped-rel', print_mode = 'series', halt = 40, check = True):
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
    """
    if check:
        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        if not isinstance(prec, (int, long, Integer)):
            raise TypeError, "prec must be an integer"
        elif isinstance(prec, (int, long)):
            prec = Integer(prec)
        if not isinstance(halt, (int, long, Integer)):
            raise TypeError, "prec must be an integer"
        elif isinstance(halt, (int, long)):
            halt = Integer(halt)
    if type != 'lazy':
        key = (p, prec, type, print_mode)
    else:
        key = (p, prec, halt, print_mode)
    if padic_field_cache.has_key(key):
        K = padic_field_cache[key]()
        if K != None:
            return K
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
def Qq(q, names=None, prec=20, type='capped-rel', print_mode='series', halt=40, modulus=None, check=True):
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

    if check:
        if names is None:
            raise TypeError, "You must specify the name of the generator."
        if isinstance(names, (list, tuple)):
            names = names[0]
        if not isinstance(prec, (int, long, Integer)):
            raise TypeError, "prec must be an integer"
        elif isinstance(prec, (int, long)):
            prec = Integer(prec)
        if not (modulus is None or isinstance(modulus, Polynomial)):
            raise TypeError, "modulus must be a polynomial"
        if not isinstance(names, str):
            raise TypeError, "names must be a string"
        if not isinstance(halt, (int, long, Integer)):
            raise TypeError, "halt must be an integer"
        elif isinstance(halt, (int, long)):
            halt = Integer(halt)

    q = Integer(q)
    F = q.factor()
    if len(F) != 1:
        raise ValueError, "q must be a prime power"
    if F[0][1] == 1:
        return Qp(q, prec, type, print_mode, halt)

    if type != 'lazy':
        key = (q, names, prec, type, print_mode)
    else:
        key = (q, names, prec, halt, print_mode)
    if qadic_field_cache.has_key(key):
        K = qadic_field_cache[key]()
        if not (K is None):
            return K

    if modulus is None:
        check = False
        from sage.rings.finite_field import GF
        modulus = PolynomialRing(Qp(F[0][0], prec, type, print_mode, halt), names)(GF(q,names).modulus().change_ring(ZZ))
    K = UnramifiedRingExtension(modulus, prec, print_mode, check)
    qadic_field_cache[key] = weakref.ref(K)

    return K


######################################################
# Short constructor names for different types
######################################################

def QpCR(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def QpCA(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def QpFM(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def QpL(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')


def QqCR(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def QqCA(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def QqFM(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def QqL(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')
