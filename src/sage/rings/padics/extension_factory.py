from padic_extension_leaves import *
#This imports all of the classes used in the ext_table below.

from sage.rings.padics.padic_ring_capped_relative import pAdicRingCappedRelative
from sage.rings.padics.padic_ring_capped_absolute import pAdicRingCappedAbsolute
from sage.rings.padics.padic_ring_fixed_mod import pAdicRingFixedMod
from sage.rings.padics.padic_ring_lazy import pAdicRingLazy
from sage.rings.padics.padic_field_capped_relative import pAdicFieldCappedRelative
from sage.rings.padics.padic_field_lazy import pAdicFieldLazy

ext_table = {}
ext_table['e', pAdicFieldCappedRelative] = EisensteinExtensionFieldCappedRelative
ext_table['e', pAdicFieldLazy] = EisensteinExtensionFieldLazy
ext_table['e', pAdicRingCappedAbsolute] = EisensteinExtensionRingCappedAbsolute
ext_table['e', pAdicRingCappedRelative] = EisensteinExtensionRingCappedRelative
ext_table['e', pAdicRingFixedMod] = EisensteinExtensionRingFixedMod
ext_table['e', pAdicRingLazy] = EisensteinExtensionRingLazy
ext_table['p', pAdicFieldCappedRelative] = pAdicGeneralExtensionFieldCappedRelative
ext_table['p', pAdicFieldLazy] = pAdicGeneralExtensionFieldLazy
ext_table['p', pAdicRingCappedAbsolute] = pAdicGeneralExtensionRingCappedAbsolute
ext_table['p', pAdicRingCappedRelative] = pAdicGeneralExtensionRingCappedRelative
ext_table['p', pAdicRingFixedMod] = pAdicGeneralExtensionRingFixedMod
ext_table['p', pAdicRingLazy] = pAdicGeneralExtensionRingLazy
ext_table['u', pAdicFieldCappedRelative] = UnramifiedExtensionFieldCappedRelative
ext_table['u', pAdicFieldLazy] = UnramifiedExtensionFieldLazy
ext_table['u', pAdicRingCappedAbsolute] = UnramifiedExtensionRingCappedAbsolute
ext_table['u', pAdicRingCappedRelative] = UnramifiedExtensionRingCappedRelative
ext_table['u', pAdicRingFixedMod] = UnramifiedExtensionRingFixedMod
ext_table['u', pAdicRingLazy] = UnramifiedExtensionRingLazy

import weakref


extension_cache = {}
def ExtensionFactory(modulus, prec = None, print_mode = None, halt = None, names = None, check = True, unram = False):
    if check:
        pass # need to add sanity checking here.  Poly degree > 1.  Irreducible
    base = modulus.base_ring()
    if not unram: #this is not quite the condition we want for not checking these things; deal with fixed-mod sanely
        if not modulus.is_monic():
            if modulus.base_ring().is_field():
                modulus = modulus / modulus.leading_coefficient()
            elif modulus.leading_coefficient().valuation() <= min(c.valuation() for c in modulus.list()):
                modulus = modulus.parent()(modulus / modulus.leading_coefficient())
            else:
                modulus = modulus / modulus.leading_coefficient()
                base = base.fraction_field()
        #Now modulus is monic
        if not krasner_check(modulus, prec):
            raise ValueError, "polynomial does not determine a unique extension.  Please specify more precision."
    if print_mode is None:
        print_mode = base.print_mode()
    if names is None:
        raise ValueError, "must specify name of generator of extension"
    if isinstance(names, tuple):
        names = names[0]
    # We now decide on the extension class: unramified, eisenstein or general (padic)
    polytype = 'p'
    if unram or is_unramified(modulus):
        polytype = 'u'
        precmult = 1
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter()
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        if prec is None:
            prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()])
        else:
            prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()] + [prec])
    elif is_eisenstein(modulus):
        polytype = 'e'
        e = modulus.degree()
        precmult = e
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter() * e
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        if prec is None:
            prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()]) * e
        else:
            prec = min([c.precision_absolute() * e for c in modulus.list()] + [base.precision_cap() * e] + [prec])
    print "polytype = %s"%polytype
    if polytype != 'p':
        modulus = truncate_to_prec(modulus, prec)
        key = (modulus, names, prec, halt, print_mode)
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        K = ext_table[polytype, type(base.ground_ring_of_tower())](modulus, prec*precmult, halt, print_mode, names)
        extension_cache[key] = weakref.ref(K)
        return K
    else:
        upoly, epoly, prec = split(modulus, prec)
        key = (upoly, epoly, names, prec, halt, print_mode)
        precmult = epoly.degree()
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        K = ext_table['p', type(base.ground_ring_of_tower())](upoly, epoly, prec*precmult, halt, print_mode, names)
        extension_cache[key] = weakref.ref(K)
        return K

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
