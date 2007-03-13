import weakref
import sage.rings.padics.padic_ring_capped_relative
import sage.rings.padics.padic_ring_capped_absolute
import sage.rings.padics.padic_ring_fixed_mod
import sage.rings.padics.padic_ring_lazy
import sage.rings.integer
import sage.rings.padics.unramified_ring_extension_capped_relative
import sage.rings.padics.unramified_ring_extension_capped_absolute
import sage.rings.padics.unramified_ring_extension_fixed_mod
import sage.rings.padics.unramified_ring_extension_lazy


from sage.rings.integer_ring import ZZ

from sage.rings.polynomial_ring import PolynomialRing

UnramifiedRingExtensionCappedRelative = sage.rings.padics.unramified_ring_extension_capped_relative.UnramifiedRingExtensionCappedRelative
UnramifiedRingExtensionCappedAbsolute = sage.rings.padics.unramified_ring_extension_capped_absolute.UnramifiedRingExtensionCappedAbsolute
UnramifiedRingExtensionFixedMod = sage.rings.padics.unramified_ring_extension_fixed_mod.UnramifiedRingExtensionFixedMod
UnramifiedRingExtensionLazy = sage.rings.padics.unramified_ring_extension_lazy.UnramifiedRingExtensionLazy
pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy
Integer = sage.rings.integer.Integer

padic_ring_cache = {}
def Zp(p, prec = 20, type = 'capped-rel', print_mode = None, halt = 40, names = None, check=True):
    """
    A creation function for p-adic rings.
    INPUT:
        p -- integer the p in Z_p
        prec -- integer (default: 20) the precision cap of the ring.  Except for the fixed modulus case, individual elements keep track of their own precision.
        type -- string (default: 'capped-rel') see Notes
        print_mode -- string (default: None) the print mode
    OUTPUT:
        the corresponding p-adic ring
    EXAMPLES:
       sage: R = Zp(5); a = R(4); a
            4 + O(5^20)
       sage: S = Zp(5, 10, type = 'capped-abs'); b = S(2); b
            2 + O(5^10)
       sage: a + b
            1 + 5 + O(5^10)

    NOTES:
         values of type:
         'capped-rel' -> pAdicRingCappedRelative.  This is the default, considers precision as the precision of the unit part.  Tracks precision of individual elements, but bounds the precision of any element with a precision cap.
        'fixed-mod'  -> pAdicRingFixedMod.  This is basically a wrapper around $\Z / p^n \Z$, adding functions appropriate to p-adics.  This is the fastest option.
        'capped-abs' -> pAdicRingCappedAbsolute.  The same as pAdicRingFixedMod, but keeps track of precision.
        'lazy' -> pAdicRingLazy.  Uses lazy elements so that additional precision can be requested during a computation.  There is some amount of performance penalty because of this ability.

        values of print_mode:
        'val-unit' -- elements are displayed as p^k*u
        'integer' -- elements are displayed as an integer
        'series' -- elements are displayed as series in p
        'val-unit-p' -- same as val-unit, except that p is written as "p"
        'integer-p' -- same as integer, except that p is written as "p"
        'series-p' -- same as series, except that p is written as "p"
    """
    # if such a ring already exists reset it's print mode and return it
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
    if names is None:
        names = str(p)
    if type != 'lazy':
        key = (p, prec, type, names)
    else:
        key = (p, prec, halt, names)
    if padic_ring_cache.has_key(key):
        K = padic_ring_cache[key]()
        if not (K is None):
            if not (print_mode is None):
                K.set_print_mode(print_mode)
            return K
    if print_mode == None:
        print_mode = 'series'
    if (type == 'capped-rel'):
        K = pAdicRingCappedRelative(p, prec, print_mode, names)
    elif (type == 'fixed-mod'):
        K = pAdicRingFixedMod(p, prec, print_mode, names)
    elif (type == 'capped-abs'):
        K = pAdicRingCappedAbsolute(p, prec, print_mode, names)
    elif (type == 'lazy'):
        K = pAdicRingLazy(p, prec, print_mode, halt, names)
    else:
        raise ValueError, "type must be one of 'capped-rel', 'fixed-mod', 'capped-abs' or 'lazy'"
    padic_ring_cache[key] = weakref.ref(K)
    return K

def Zq(q, prec = 20, type = 'capped-abs', modulus = None, names=None, print_mode = None, halt = 40, zp_name = None, check = True):
    r"""
    The creation function for unramified extensions of $\Z_p$
    """
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
        return Zp(q, prec, type, print_mode, halt, names, check)
    base = Zp(F[0][0], prec, type, print_mode, halt, zp_name, check = False)
    if modulus is None:
        from sage.rings.finite_field import GF
        if zp_name is None:
            zp_name = str(p)
        modulus = PolynomialRing(base, 'x')(GF(q, names).modulus().change_ring(ZZ))
    return ExtensionFactory(base, modulus, prec, names, print_mode, halt, check, unram = True)
