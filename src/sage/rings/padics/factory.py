from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.rings.padics.padic_ring_capped_relative import pAdicRingCappedRelative
from sage.rings.padics.padic_ring_capped_absolute import pAdicRingCappedAbsolute
from sage.rings.padics.padic_ring_fixed_mod import pAdicRingFixedMod
from sage.rings.padics.padic_ring_lazy import pAdicRingLazy
from sage.rings.padics.padic_field_capped_relative import pAdicFieldCappedRelative
from sage.rings.padics.padic_field_lazy import pAdicFieldLazy

######################################################
# ext_table --
# This dictionary controls what class is created by the extension
# factory when it finds a given class in the ground ring of the tower.
######################################################

from padic_extension_leaves import *
#This imports all of the classes used in the ext_table below.

ext_table = {}
ext_table['e', pAdicFieldCappedRelative] = EisensteinExtensionFieldCappedRelative
#ext_table['e', pAdicFieldLazy] = EisensteinExtensionFieldLazy
ext_table['e', pAdicRingCappedAbsolute] = EisensteinExtensionRingCappedAbsolute
ext_table['e', pAdicRingCappedRelative] = EisensteinExtensionRingCappedRelative
ext_table['e', pAdicRingFixedMod] = EisensteinExtensionRingFixedMod
#ext_table['e', pAdicRingLazy] = EisensteinExtensionRingLazy
#ext_table['p', pAdicFieldCappedRelative] = pAdicGeneralExtensionFieldCappedRelative
#ext_table['p', pAdicFieldLazy] = pAdicGeneralExtensionFieldLazy
#ext_table['p', pAdicRingCappedAbsolute] = pAdicGeneralExtensionRingCappedAbsolute
#ext_table['p', pAdicRingCappedRelative] = pAdicGeneralExtensionRingCappedRelative
#ext_table['p', pAdicRingFixedMod] = pAdicGeneralExtensionRingFixedMod
#ext_table['p', pAdicRingLazy] = pAdicGeneralExtensionRingLazy
ext_table['u', pAdicFieldCappedRelative] = UnramifiedExtensionFieldCappedRelative
#ext_table['u', pAdicFieldLazy] = UnramifiedExtensionFieldLazy
ext_table['u', pAdicRingCappedAbsolute] = UnramifiedExtensionRingCappedAbsolute
ext_table['u', pAdicRingCappedRelative] = UnramifiedExtensionRingCappedRelative
ext_table['u', pAdicRingFixedMod] = UnramifiedExtensionRingFixedMod
#ext_table['u', pAdicRingLazy] = UnramifiedExtensionRingLazy

import weakref

#######################################################################################################
#
#  p-Adic Fields
#  Qp -- base field
#  Qq -- unramified extension field of Qp
#  QpCR, QpL, QqCR, QqL -- shortcuts for capped relative and lazy versions of Qp and Qq
#
#######################################################################################################

padic_field_cache = {}
def Qp(p, prec = 20, type = 'capped-rel', print_mode = None, halt = 40, names = None, check = True):
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

        #sage: L = Qp(5, 10, type = 'lazy'); b = L(2); b
        #2 + O(5^10)
        #sage: a + b
        #1 + 5 + O(5^20)

    NOTES:
        values of type:
        'capped-rel' -> pAdicFieldCappedRelative.  This is the default, considers precision as the precision of the unit part.  Tracks precision of individual elements, but bounds the precision of any element with a precision cap.
        'lazy' -> pAdicFieldLazy.  Uses lazy elements so that additional precision can be requested during a computation.  There is some amount of performance penalty because of this ability.

        values of print_mode:
        Leaving print_mode as None uses the global default print mode.  Other allowable values are:
        'val-unit' -- elements are displayed as p^k*u
        'terse' -- elements are displayed as an integer in base 10 or the quotient of an integer by a power of p (still in base 10)
        'series' -- elements are displayed as series in p
        'digits' -- elements are displayed as a string of base p digits
        'bars' -- elements are displayed as a string of base p digits with separators
        For more details and more control, see sage.rings.padics.padic_printing or look at padic_printing.<tab> from the command line.
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
    if names is None:
        name = str(p)
    elif isinstance(names, tuple):
        name = names[0]
    else:
        name = str(names)
    if type != 'lazy':
        key = (p, prec, type, name, print_mode)
    else:
        key = (p, prec, halt, name, print_mode)
    if padic_field_cache.has_key(key):
        K = padic_field_cache[key]()
        if K != None:
            return K
    if (type == 'capped-rel'):
        K = pAdicFieldCappedRelative(p, prec, print_mode, name)
    elif (type == 'lazy'):
        raise NotImplementedError, "lazy p-adics need more work.  Sorry."
        K = pAdicFieldLazy(p, prec, print_mode, halt, name)
    else:
        raise ValueError, "type must be either 'capped-rel' or 'lazy'"
    padic_field_cache[key] = weakref.ref(K)
    return K

######################################################
# Qq -- unramified extensions
######################################################

def Qq(q, prec = 20, type = 'capped-rel', modulus = None, names=None, print_mode=None, halt=40, qp_name = None, check=True):
    r"""
    Given a prime power q = p^n, return the unique unramified extension
    of Qp of degree n.

    Currently, there's no code for unramified field extensions, so
    we just return the UnramifiedRingExtension.

        values of print_mode:
        Leaving print_mode as None uses the global default print mode.  Other allowable values are:
        'val-unit' -- elements are displayed as p^k*u
        'terse' -- elements are displayed as an integer in base 10 or the quotient of an integer by a power of p (still in base 10)
        'series' -- elements are displayed as series in p
        'digits' -- elements are displayed as a string of base p digits
        'bars' -- elements are displayed as a string of base p digits with separators
        For more details and more control, see sage.rings.padics.padic_printing or look at padic_printing.<tab> from the command line.
    """

    from sage.rings.integer import Integer
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
        return Qp(q, prec, type, print_mode, halt, names, check)
    base = Qp(F[0][0], prec, type, print_mode, halt, qp_name, check = False)
    if modulus is None:
        from sage.rings.finite_field import GF
        from sage.rings.integer_ring import ZZ
	from sage.rings.polynomial.polynomial_ring import PolynomialRing
        if qp_name is None:
            qp_name = (str(F[0][0]),)
        modulus = PolynomialRing(base, 'x')(GF(q, names).modulus().change_ring(ZZ))
    return ExtensionFactory(modulus, prec, print_mode, halt, names, check, unram = True)

######################################################
# Short constructor names for different types
######################################################

def QpCR(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def QpL(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Qp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')


def QqCR(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def QqL(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Qq(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')

#######################################################################################################
#
#  p-Adic Rings
#  Zp -- base rings
#  Zq -- unramified extension ring of Zp
#  ZpCR, ZpCA, ZpFM, ZpL, ZqCR, ZqCA, ZqFM, ZqL -- shortcuts for capped relative and lazy versions of Zp and Zq
#
#######################################################################################################

padic_ring_cache = {}
def Zp(p, prec = 20, type = 'capped-rel', print_mode = None, halt = 40, names = None, check=True):
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

    It works even with a fairly huge cap:
       sage: Zp(next_prime(10^50), 100000)
       100000000000000000000000000000000000000000000000151-adic Ring with capped relative precision 100000

    We create each type of ring:
        sage: Zp(7, 20, 'capped-rel')
        7-adic Ring with capped relative precision 20
        sage: Zp(7, 20, 'fixed-mod')
        7-adic Ring of fixed modulus 7^20
        sage: Zp(7, 20, 'capped-abs')
        7-adic Ring with capped absolute precision 20

        #sage: Zp(7, 20, 'lazy')
        #Lazy 7-adic Ring

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

        sage: k = Zp(7, print_mode='terse'); k
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

       print_mode -- string (default: None)
        Leaving print_mode as None uses the global default print mode.  Other allowable values are:
        'val-unit' -- elements are displayed as p^k*u
        'terse' -- elements are displayed as an integer in base 10
        'series' -- elements are displayed as series in p
        'digits' -- elements are displayed as a string of base p digits
        'bars' -- elements are displayed as a string of base p digits with separators
        For more details and more control, see sage.rings.padics.padic_printing or look at padic_printing.<tab> from the command line.
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
    if names is None:
        name = str(p)
    elif isinstance(names, tuple):
        name = names[0]
    else:
        name = str(names)
    if type != 'lazy':
        key = (p, prec, type, name, print_mode)
    else:
        key = (p, prec, halt, name, print_mode)
    if padic_ring_cache.has_key(key):
        K = padic_ring_cache[key]()
        if not (K is None):
            return K
    if (type == 'capped-rel'):
        K = pAdicRingCappedRelative(p, prec, print_mode, name)
    elif (type == 'fixed-mod'):
        K = pAdicRingFixedMod(p, prec, print_mode, name)
    elif (type == 'capped-abs'):
        K = pAdicRingCappedAbsolute(p, prec, print_mode, name)
    elif (type == 'lazy'):
        raise NotImplementedError, "lazy p-adics need more work.  Sorry."
        K = pAdicRingLazy(p, prec, print_mode, halt, name)
    else:
        raise ValueError, "type must be one of 'capped-rel', 'fixed-mod', 'capped-abs' or 'lazy'"
    padic_ring_cache[key] = weakref.ref(K)
    return K

######################################################
# Zq -- unramified extensions
######################################################

def Zq(q, prec = 20, type = 'capped-abs', modulus = None, names=None,
          print_mode=None, halt = 40, zp_name = None, check = True):
    r"""
    Return an unramified extension of $\Z_p$.

    INPUT:
        q -- prime power
        prec -- integer (default: 20)
        type -- string (default: 'capped-abs'); see the documentation for Zp
        modulus -- polynomial (default: None)
        names -- tuple
        print_mode -- string (default: None); see the documentation for Zp
        halt -- integer (default: 40): only applicable for type='lazy'
        zp_name -- string (default: None): a name for the underlying Zp's prime
        check -- bool (default: True): whether to verify that the input is valid.

    OUTPUT:
        -- an unramified extension of Z_p

    EXAMPLES:
        sage: k.<a> = Zq(4); k
        Unramified Extension of 2-adic Ring with capped absolute precision 20
        in a defined by (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20))
        sage: k.<a> = Zq(3^10); k
        Unramified Extension of 3-adic Ring with capped absolute precision 20 in a
        defined by (1 + O(3^20))*x^10 + (2 + O(3^20))*x^6 + (2 + O(3^20))*x^5 +
        (2 + O(3^20))*x^4 + (1 + O(3^20))*x + (2 + O(3^20))
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
            names = str(names)
            #raise TypeError, "names must be a string"
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
            zp_name = (str(F[0][0]),)
        modulus = PolynomialRing(base, 'x')(GF(q, names).modulus().change_ring(ZZ))
    return ExtensionFactory(base, modulus, prec, print_mode, halt, names, check, unram = True)

######################################################
# Short constructor names for different types
######################################################

def ZpCR(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def ZpCA(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def ZpFM(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def ZpL(p, prec = 20, print_mode = None, halt = 40, check=True):
    return Zp(p=p, prec=prec, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')


def ZqCR(q, prec = 20, modulus = None, names = None, print_mode = 'series', halt = 40, zp_name = None, check=True):
    return Zq(q=q, prec=prec, modulus = modulus, names = names, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-rel')

def ZqCA(q, prec = 20, modulus = None, names = None, print_mode = 'series', halt = 40, zp_name = None, check=True):
    return Zq(q=q, prec=prec, modulus = modulus, names = names, print_mode=print_mode, halt=halt, check=check,
              type = 'capped-abs')

def ZqFM(q, prec = 20, modulus = None, names = None, print_mode = 'series', halt = 40, zp_name = None, check=True):
    return Zq(q=q, prec=prec, modulus = modulus, names = names, print_mode=print_mode, halt=halt, check=check,
              type = 'fixed-mod')

def ZqL(q, prec = 20, modulus = None, names = None, print_mode = 'series', halt = 40, zp_name = None, check=True):
    return Zq(q=q, prec=prec, modulus = modulus, names = names, print_mode=print_mode, halt=halt, check=check,
              type = 'lazy')

#######################################################################################################
#
#  The Extension Factory -- creates extensions of p-adic rings and fields
#
#######################################################################################################

extension_cache = {}
def ExtensionFactory(base, premodulus, prec = None, print_mode = None, halt = None, names = None, check = True, unram = False):
    from sage.calculus.all import is_SymbolicExpression
    from sage.rings.polynomial.all import is_Polynomial
    if check:
        if is_SymbolicExpression(premodulus):
            if len(premodulus.variables()) != 1:
                raise ValueError, "symbolic expression must be in only one variable"
            modulus = premodulus.polynomial(base)
        elif is_Polynomial(premodulus):
            if premodulus.parent().ngens() != 1:
                raise ValueError, "must use univariate polynomial"
            modulus = premodulus.change_ring(base)
        else:
            raise ValueError, "modulus must be a polynomial"
        if modulus.degree() <= 1:
            raise NotImplementedError, "degree of modulus must be at least 2"
        # need to add more checking here.
        if not unram: #this is not quite the condition we want for not checking these things; deal with fixed-mod sanely
            if not modulus.is_monic():
                if base.is_field():
                    modulus = modulus / modulus.leading_coefficient()
                elif modulus.leading_coefficient().valuation() <= min(c.valuation() for c in modulus.list()):
                    modulus = modulus.parent()(modulus / modulus.leading_coefficient())
                else:
                    modulus = modulus / modulus.leading_coefficient()
                    base = base.fraction_field()
            #Now modulus is monic
            if not krasner_check(modulus, prec):
                raise ValueError, "polynomial does not determine a unique extension.  Please specify more precision or use parameter check=False."
        if print_mode is None:
            print_mode = base.print_mode()
        if names is None:
            raise ValueError, "must specify name of generator of extension"
        if isinstance(names, tuple):
            names = names[0]
    else:
        modulus = premodulus
    #print type(base)
    # We now decide on the extension class: unramified, eisenstein, two-step or general
    if unram or is_unramified(modulus):
        polytype = 'u'
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter()
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        if prec is None:
            prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()])
        else:
            prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()] + [prec])
        shift_seed = None
    elif is_eisenstein(modulus):
        polytype = 'e'
        e = modulus.degree()
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter() * e
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        # The precision of an eisenstein extension is governed both by the absolute precision of the polynomial,
        # and also by the precision of polynomial with the leading term removed (for shifting).
        # The code below is to determine the correct prec for the extension, and possibly to obtain
        # the information needed to shift right with full precision from the premodulus.
        if is_SymbolicExpression(premodulus):
            # Here we assume that the output of coeffs is sorted in increasing order by exponent:
            coeffs = premodulus.coeffs()
            preseed = premodulus / coeffs[-1][0]
            preseed -= preseed.variables()[0]**coeffs[-1][1]
            preseed /= base.prime() # here we assume that the base is unramified over Qp
            shift_seed = -preseed.polynomial(base)
        else: # a polynomial
            if not premodulus.is_monic():
                preseed = preseed / premodulus.leading_coefficient()
            else:
                preseed = premodulus
            preseed = preseed[:preseed.degree()]
            if base.is_fixed_mod():
                shift_seed = -preseed.change_ring(base)
                shift_seed = shift_seed.parent()([a >> 1 for a in shift_seed.list()])
            else:
                if base.e() == 1:
                    try:
                        preseed *= 1/base.prime()
                        shift_seed = -preseed.change_ring(base)
                    except TypeError:
                        # give up on getting more precision
                        shift_seed = -preseed.change_ring(base)
                        shift_seed /= base.uniformizer()
                else:
                    # give up on getting more precision
                    shift_seed = -preseed.change_ring(base)
                    shift_seed /= base.uniformizer()
        if prec is None:
            prec = min([c.precision_absolute() for c in shift_seed.list() if not c._is_exact_zero()] + [modulus.leading_coefficient().precision_absolute()] + [base.precision_cap()]) * e
        else:
            prec = min([c.precision_absolute() * e for c in shift_seed.list() if not c._is_exact_zero()] + [modulus.leading_coefficient().precision_absolute() * e] + [base.precision_cap() * e] + [prec])
    else:
        polytype = 'p'
    #print "polytype = %s"%polytype
    if polytype == 'u' or polytype == 'e':
        modulus = truncate_to_prec(modulus, prec)
        key = (base, modulus, names, prec, halt, print_mode)
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        #print ext_table[polytype, type(base.ground_ring_of_tower())]
        K = ext_table[polytype, type(base.ground_ring_of_tower())](premodulus, modulus, prec, halt, print_mode, shift_seed, names)
        extension_cache[key] = weakref.ref(K)
        return K
    else:
        upoly, epoly, prec = split(modulus, prec)
        key = (base, upoly, epoly, names, prec, halt, print_mode)
        precmult = epoly.degree()
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        K = ext_table['p', type(base.ground_ring_of_tower())](premodulus, upoly, epoly, prec*precmult, halt, print_mode, names)
        extension_cache[key] = weakref.ref(K)
        return K

######################################################
# Helper functions for the Extension Factory
######################################################

def split(poly, prec):
    raise NotImplementedError, "Extensions by general polynomials not yet supported.  Please use an unramified or eisenstein polynomial."

def truncate_to_prec(poly, absprec):
    return poly.parent()(poly, absprec = absprec) # Is this quite right?  We don't want flat necessarily...

def krasner_check(poly, prec):
    return True #This needs to be implemented

def is_eisenstein(poly):
    if poly[0].valuation() != 1:
        return False
    if reduce(lambda a, b: a or b, [(c.valuation() < 1) for c in poly.list()[1:poly.degree()]]):
        return False
    return True

def is_unramified(poly):
    if poly[0].valuation() > 0:
        return False
    if reduce(lambda a, b: a or b, [(c.valuation() < 0) for c in poly.list()[1:poly.degree()]]):
        return False
    F = poly.parent().change_ring(poly.base_ring().residue_class_field())(poly).factor()
    if len(F) != 1 or F[0][1] != 1:
        return False
    return True
