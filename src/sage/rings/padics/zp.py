import weakref
import sage.rings.padics.padic_ring_capped_relative
import sage.rings.padics.padic_ring_capped_absolute
import sage.rings.padics.padic_ring_fixed_mod
import sage.rings.padics.padic_ring_lazy
import sage.rings.integer

from sage.rings.integer_ring import ZZ

from sage.rings.polynomial_ring import PolynomialRing
import sage.rings.padics.unramified_ring_extension

pAdicRingCappedRelative = sage.rings.padics.padic_ring_capped_relative.pAdicRingCappedRelative
pAdicRingCappedAbsolute = sage.rings.padics.padic_ring_capped_absolute.pAdicRingCappedAbsolute
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingLazy = sage.rings.padics.padic_ring_lazy.pAdicRingLazy
Integer = sage.rings.integer.Integer
UnramifiedRingExtension = sage.rings.padics.unramified_ring_extension.UnramifiedRingExtension

padic_ring_cache = {}
def Zp(p, prec = 20, type = 'capped-rel', print_mode = 'series', halt = 40, check=True):
    """
    Return a model of the $p$-adic integer $\Z_p$.

    INPUT:
        p -- integer the p in Z_p
        prec -- integer (default: 20) the precision cap of the ring.
                Except for the fixed modulus case, individual elements keep
                track of their own precision.
        type -- string (default: 'capped-rel') see notes section below for options.
        print_mode -- string (default: series) the print mode; see notes section
                below for options.
        halt -- integer (default: 40): only applicable for type='lazy'
        check -- bool (default: True): wether to verify that the input is valid.

    OUTPUT:
        the corresponding p-adic ring

    EXAMPLES:
    We create rings with various parameters
       sage: Zp(7)
       7-adic Ring with capped relative precision 20
       sage: Zp(9)
       Traceback (most recent call last):
       ...
       ValueError: p must be prime
       sage: Zp(17, 5)
       17-adic Ring with capped relative precision 5
       sage: Zp(17, 5)(-1)
       16 + 16*17 + 16*17^2 + 16*17^3 + 16*17^4 + O(17^5)
       sage: Zp(next_prime(10^50), 10000)
       100000000000000000000000000000000000000000000000151-adic Ring with capped relative precision 10000

    We create each type of ring:
        sage: Zp(7, 20, 'capped-rel')
        7-adic Ring with capped relative precision 20
        sage: Zp(7, 20, 'fixed-mod')
        7-adic Ring of fixed modulus 7^20
        sage: Zp(7, 20, 'capped-abs')
        7-adic Ring with capped absolute precision 20
        sage: Zp(7, 20, 'lazy')
        Lazy 7-adic Ring

    We create a capped relative ring with each print mode:
        sage: k = Zp(7, 8, print_mode='series'); k
        7-adic Ring with capped relative precision 8
        sage: k(7*(19))
        5*7 + 2*7^2 + O(7^9)
        sage: k(7*(-19))
        2*7 + 4*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + O(7^9)

        sage: k = Zp(7, print_mode='val-unit'); k
        7-adic Ring with capped relative precision 20
        sage: k(7*(19))
        7 * 19 + O(7^21)
        sage: k(7*(-19))
        7 * 79792266297611982 + O(7^21)

        sage: k = Zp(7, print_mode='integer'); k
        7-adic Ring with capped relative precision 20
        sage: k(7*(19))
        133 + O(7^21)
        sage: k(7*(-19))
        558545864083283874 + O(7^21)

    Note that $p$-adic rings are cached (via weak references):
        sage: a = Zp(7); b = Zp(7)
        sage: a is b
        True

    We create some elements in various rings:
        sage: R = Zp(5); a = R(4); a
        4 + O(5^20)
        sage: S = Zp(5, 10, type = 'capped-abs'); b = S(2); b
        2 + O(5^10)
        sage: a + b
        1 + 5 + O(5^10)

    NOTES:
       type -- string (default: 'capped-rel'), the type of p-adic ring.

           'capped-rel' -- pAdicRingCappedRelative.  This is the
                           default, considers precision as the
                           precision of the unit part.  Tracks
                           precision of individual elements, but
                           bounds the precision of any element with a
                           precision cap.
           'fixed-mod' --  pAdicRingFixedMod.  This is basically a
                           wrapper around $\Z / p^n \Z$, adding
                           functions appropriate to p-adics.  This is
                           the fastest option.
           'capped-abs' -- pAdicRingCappedAbsolute.  The same as
                           pAdicRingFixedMod, but keeps track of
                           precision.
           'lazy' -- pAdicRingLazy.  Uses lazy elements so that
                     additional precision can be requested during a
                     computation.  There is some amount of performance
                     penalty because of this ability.

       print_mode -- string (default: 'series', unless it has been
                     previously specified for a cached version of this ring)
           'val-unit' -- elements are displayed as p^k*u
           'integer' -- elements are displayed as an integer
           'series' -- elements are displayed as series in p
    """
    # if such a ring already exists reset it's print mode (unless the input print mode is None) and return it
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
    if padic_ring_cache.has_key(key):
        K = padic_ring_cache[key]()
        if not (K is None):
            return K
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

qadic_ring_cache = {}
def Zq(q, prec = 20, type = 'capped-abs', modulus = None, names=None,
          print_mode='series', halt = 40, check = True):
    r"""
    Return an unramified extension of $\Z_p$.

    INPUT:
        q -- prime power
        prec -- integer (default: 20)
        type -- string (default: 'capped-abs'); see the documentation for Zq
        modulus -- polynomial (default: None)
        names -- tuple
        print_mode -- string (default: 'series'); see the documentation for print_mode
        halt -- integer (default: 40): only applicable for type='lazy'
        check -- bool (default: True): wether to verify that the input is valid.

    OUTPUT:
        -- an unramified extension of Z_p

    EXAMPLES:

    TODO: This printing is all completely backwards -- a and x must be switched.
    We
        sage: k.<a> = Zq(4); k
        Unramified Extension of 2-adic Ring with capped absolute precision 20
        in x defined by (1 + O(2^20))*a^2 + (1 + O(2^20))*a + 1 + O(2^20)
        sage: k.<a> = Zq(3^10); k
        Unramified Extension of 3-adic Ring with capped absolute precision 20 in x
        defined by (1 + O(3^20))*a^10 + O(3^20)*a^9 + O(3^20)*a^8 + O(3^20)*a^7 +
        (2 + O(3^20))*a^6 + (2 + O(3^20))*a^5 + (2 + O(3^20))*a^4 + O(3^20)*a^3 +
        O(3^20)*a^2 + (1 + O(3^20))*a + 2 + O(3^20)
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
        return Zp(q, prec, type, print_mode, halt)
    if type != 'lazy':
        key = (q, names, prec, type, modulus, print_mode)
    else:
        key = (q, names, prec, halt, modulus, print_mode)
    if qadic_ring_cache.has_key(key):
        K = qadic_ring_cache[key]()
        if not (K is None):
            return K
    if modulus is None:
        check = False
        from sage.rings.finite_field import GF
        modulus = PolynomialRing(Zp(F[0][0], prec, type, print_mode, halt), names)(GF(q, names).modulus().change_ring(ZZ))

    K = UnramifiedRingExtension(modulus, prec, print_mode, check)
    qadic_ring_cache[key] = weakref.ref(K)
    return K


######################################################
# Short constructor names for different types
######################################################

def ZpCR(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def ZpCA(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def ZpFM(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def ZpL(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')


def ZqCR(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def ZqCA(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def ZqFM(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def ZqL(p, prec = 20, print_mode = 'series', halt = 40, check=True):
    return Zq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')
