import weakref
import sage.rings.padics.padic_ring_capped_relative
import sage.rings.padics.padic_ring_capped_absolute
import sage.rings.padics.padic_ring_fixed_mod
import sage.rings.padics.padic_ring_lazy

pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy

padic_ring_cache = {}
def Zp(p, prec = 20, type = 'capped-rel', print_mode = None, halt = 40):
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
            6 + O(5^10)

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
    if not p.is_prime():
        raise ValueError, "p must be prime"
    if type != 'lazy':
        key = (p, prec, type)
    else:
        key = (p, prec, halt)
    if padic_ring_cache.has_key(key):
        K = padic_ring_cache[key]()
        if K != None:
            if print_mode != None:
                K.set_print_mode(print_mode)
            return K
    if print_mode == None:
        print_mode = 'val-unit'
    if (type == 'capped-rel'):
        K = pAdicRingCappedRelative(p, prec, print_mode)
    elif (type == 'fixed-mod'):
        K = pAdicRingFixedMod(p, prec, print_mode)
    elif (type == 'capped-abs'):
        K = pAdicRingCappedAbsolute(p, prec, print_mode)
    elif (type == 'lazy'):
        K = pAdicRingLazy(p, prec, print_mode, halt)
    else:
        raise ValueError, "type must be one of 'capped-rel', 'fixed-mod', 'capped-abs' or 'lazy'"
    padic_ring_cache[key] = weakref.ref(K)
    return K
