import sage.rings.padics.eisenstein_field_extension_capped_relative
import sage.rings.padics.eisenstein_field_extension_lazy
import sage.rings.padics.eisenstein_ring_extension_capped_absolute
import sage.rings.padics.eisenstein_ring_extension_capped_relative
import sage.rings.padics.eisenstein_ring_extension_fixed_mod
import sage.rings.padics.eisenstein_ring_extension_lazy
import sage.rings.padics.padic_field_extension_capped_relative
import sage.rings.padics.padic_field_extension_lazy
import sage.rings.padics.padic_ring_extension_capped_absolute
import sage.rings.padics.padic_ring_extension_capped_relative
import sage.rings.padics.padic_ring_extension_fixed_mod
import sage.rings.padics.padic_ring_extension_lazy
import sage.rings.padics.unramified_field_extension_capped_relative
import sage.rings.padics.unramified_field_extension_lazy
import sage.rings.padics.unramified_ring_extension_capped_absolute
import sage.rings.padics.unramified_ring_extension_capped_relative
import sage.rings.padics.unramified_ring_extension_fixed_mod
import sage.rings.padics.unramified_ring_extension_lazy

from sage.rings.padics.padic_ring_capped_relative import pAdicRingCappedRelative
from sage.rings.padics.padic_ring_capped_absolute import pAdicRingCappedAbsolute
from sage.rings.padics.padic_ring_fixe_mod import pAdicRingFixedMod
from sage.rings.padics.padic_ring_lazy import pAdicRingLazy
from sage.rings.padics.padic_field_capped_relative import pAdicFieldCappedRelative
from sage.rings.padics.padic_field_lazy import pAdicFieldLazy

ext_table = {}
ext_table['e', pAdicFieldCappedRelative] = sage.rings.padics.eisenstein_field_extension_capped_relative.EisensteinFieldExtensionCappedRelative
ext_table['e', pAdicFieldLazy] = sage.rings.padics.eisenstein_field_extension_lazy.EisensteinFieldExtensionLazy
ext_table['e', pAdicRingCappedAbsolute] = sage.rings.padics.eisenstein_ring_extension_capped_absolute.EisensteinRingExtensionCappedAbsolute
ext_table['e', pAdicRingCappedRelative] = sage.rings.padics.eisenstein_ring_extension_capped_relative.EisensteinRingExtensionCappedRelative
ext_table['e', pAdicRingFixedMod] = sage.rings.padics.eisenstein_ring_extension_fixed_mod.EisensteinRingExtensionFixedMod
ext_table['e', pAdicRingLazy] = sage.rings.padics.eisenstein_ring_extension_lazy.EisensteinRingExtensionLazy
ext_table['p', pAdicFieldCappedRelative] = sage.rings.padics.padic_field_extension_capped_relative.pAdicFieldExtensionCappedRelative
ext_table['p', pAdicFieldLazy] = sage.rings.padics.padic_field_extension_lazy.pAdicFieldExtensionLazy
ext_table['p', pAdicRingCappedAbsolute] = sage.rings.padics.padic_ring_extension_capped_absolute.pAdicRingExtensionCappedAbsolute
ext_table['p', pAdicRingCappedRelative] = sage.rings.padics.padic_ring_extension_capped_relative.pAdicRingExtensionCappedRelative
ext_table['p', pAdicRingFixedMod] = sage.rings.padics.padic_ring_extension_fixed_mod.pAdicRingExtensionFixedMod
ext_table['p', pAdicRingLazy] = sage.rings.padics.padic_ring_extension_lazy.pAdicRingExtensionLazy
ext_table['u', pAdicFieldCappedRelative] = sage.rings.padics.unramified_field_extension_capped_relative.UnramifiedFieldExtensionCappedRelative
ext_table['u', pAdicFieldLazy] = sage.rings.padics.unramified_field_extension_lazy.UnramifiedFieldExtensionLazy
ext_table['u', pAdicRingCappedAbsolute] = sage.rings.padics.unramified_ring_extension_capped_absolute.UnramifiedRingExtensionCappedAbsolute
ext_table['u', pAdicRingCappedRelative] = sage.rings.padics.unramified_ring_extension_capped_relative.UnramifiedRingExtensionCappedRelative
ext_table['u', pAdicRingFixedMod] = sage.rings.padics.unramified_ring_extension_fixed_mod.UnramifiedRingExtensionFixedMod
ext_table['u', pAdicRingLazy] = sage.rings.padics.unramified_ring_extension_lazy.UnramifiedRingExtensionLazy

import weakref


extension_cache = {}
def ExtensionFactory(base, modulus, prec = None, names = None, print_mode = None, halt = None, check = True, unram = False):
    if check:
        pass # need to add sanity checking here
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
        if not krasner_check(poly, prec):
            raise ValueError, "polynomial does not determine a unique extension.  Please specify more precision."
    if print_mode is None:
        print_mode = base.print_mode()
    if names is None:
        raise ValueError, "must specify name of generator of extension"
    # We now decide on the extension class: unramified, eisenstein or general (padic)
    polytype = 'p'
    if unram or is_unramified(poly):
        polytype = 'u'
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter()
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        if prec is None:
            prec = min([c.precision_absolute() for c in poly.list()].append(base.precision_cap()))
        else:
            prec = min([c.precision_absolute() for c in poly.list()].append(base.precision_cap()).append(prec))
    elif is_eisenstein(poly):
        polytype = 'e'
        e = poly.degree()
        if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = base.halting_paramter() * e
        elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
            halt = None
        if prec is None:
            prec = min([c.precision_absolute() for c in poly.list()].append(base.precision_cap())) * e
        else:
            prec = min([c.precision_absolute() * e for c in poly.list()].append(base.precision_cap() * e).append(prec))
    if polytype != 'p':
        modulus = truncate_to_prec(modulus, prec)
        key = (base, modulus, names, prec, halt, print_mode)
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        K = ext_table[polytype, type(base.ground_ring_of_tower())](base, modulus, names, prec, halt, print_mode)
        extension_cache[key] = weakref.ref(K)
        return K
    else:
        upoly, intermediate_ext, epoly, prec = split(modulus, prec)
        key = (base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode)
        if extension_cache.has_key(key):
            K = extension_cache[key]()
            if not (K is None):
                return K
        K = ext_table['p', type(base.ground_ring_of_tower())](base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode)
        extension_cache[key] = weakref.ref(K)
        return K

def split(poly, prec):
    raise NotImplementedError, "Extensions by general polynomials not yet supported.  Please use an unramified or eisenstein polynomial."

def truncate_to_prec(poly, absprec):
    return modulus.parent()(modulus, absprec = prec) # Is this quite right?  We don't want flat necessarily...

def krasner_check(poly, prec):
    return True #This needs to be implemented

def is_eisnstein(poly):
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
