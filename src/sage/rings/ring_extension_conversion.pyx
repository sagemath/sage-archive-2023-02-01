#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.categories.pushout import construction_tower
from sage.categories.map cimport Map, FormalCompositeMap
from sage.categories.morphism import IdentityMorphism
from sage.rings.ring_extension cimport RingExtension_generic
from sage.rings.ring_extension_element cimport RingExtensionElement
from sage.rings.ring_extension_morphism cimport RingExtensionHomomorphism
from sage.rings.ring_extension_morphism cimport RingExtensionBackendIsomorphism
from sage.rings.ring_extension_morphism cimport RingExtensionBackendReverseIsomorphism


# For parents
#############

cpdef backend_parent(R):
    if isinstance(R, RingExtension_generic):
        return (<RingExtension_generic>R)._backend
    else:
        return R

cpdef from_backend_parent(R, RingExtension_generic E):
    tower = construction_tower(R)
    bases_tower = [ parent for (_, parent) in tower ]
    for base in E.bases():
        backend = backend_parent(base)
        try:
            s = bases_tower.index(backend)
            ans = base
            for i in range(s, 0, -1):
                functor = tower[i][0]
                ans = functor(ans)
            return ans
        except ValueError:
            pass
    return R


# For elements
##############

cpdef backend_element(x):
    if isinstance(x, RingExtensionElement):
        return (<RingExtensionElement>x)._backend
    else:
        return x

cpdef from_backend_element(x, RingExtension_generic E):
    parent = from_backend_parent(x.parent(),E)
    if parent is None:
        return x
    else:
        return parent(x)


# For morphisms
###############

cdef _backend_morphism(f):
    domain = f.domain()
    if not isinstance(f.domain(), RingExtension_generic) and not isinstance(f.codomain(), RingExtension_generic):
        return f
    elif isinstance(f, RingExtensionHomomorphism):
        return (<RingExtensionHomomorphism>f)._backend
    elif isinstance(f, FormalCompositeMap):
        return _backend_morphism(f.then()) * _backend_morphism(f.first())
    elif isinstance(f, IdentityMorphism):
        ring = backend_parent(domain)
        return ring.Hom(ring).identity()
    elif domain is domain.base_ring():
        ring = f.codomain()
        if isinstance(ring, RingExtension_generic):
            ring = ring._backend
        if ring.has_coerce_map_from(domain):
            return ring.coerce_map_from(domain)
    raise NotImplementedError

cpdef backend_morphism(f, forget="all"):
    try:
        g = _backend_morphism(f)
        if forget is None and (isinstance(f.domain(), RingExtension_generic) or isinstance(f.codomain(), RingExtension_generic)):
            g = RingExtensionHomomorphism(f.domain().Hom(f.codomain()), g)
        if forget == "domain" and isinstance(f.codomain(), RingExtension_generic):
            g = RingExtensionHomomorphism(g.domain().Hom(f.codomain()), g)
        if forget == "codomain" and isinstance(f.domain(), RingExtension_generic):
            g = RingExtensionHomomorphism(f.domain().Hom(g.codomain()), g)
    except NotImplementedError:
        g = f
        if (forget == "all" or forget == "domain") and isinstance(f.domain(), RingExtension_generic):
            ring = f.domain()._backend
            g = g * RingExtensionBackendIsomorphism(ring.Hom(f.domain()))
        if (forget == "all" or forget == "codomain") and isinstance(f.codomain(), RingExtension_generic):
            ring = f.codomain()._backend
            g = RingExtensionBackendReverseIsomorphism(f.codomain().Hom(ring)) * g
    return g

cpdef from_backend_morphism(f, RingExtension_generic E):
    cdef domain = from_backend_parent(f.domain(), E)
    cdef codomain = from_backend_parent(f.codomain(), E)
    return domain.Hom(codomain)(f)


# Generic
#########

def to_backend(arg):
    if isinstance(arg, list):
        return [ to_backend(x) for x in arg ]
    elif isinstance(arg, tuple):
        return tuple([ to_backend(x) for x in arg ])
    elif isinstance(arg, dict):
        return { to_backend(key): to_backend(value) for (key, value) in arg.items() }
    elif isinstance(arg, RingExtension_generic):
        return (<RingExtension_generic>arg)._backend
    elif isinstance(arg, RingExtensionElement):
        return (<RingExtensionElement>arg)._backend
    return arg

def from_backend(arg, E):
    ans = None
    if isinstance(arg, list):
        ans = [ from_backend(x,E) for x in arg ]
    elif isinstance(arg, tuple):
        ans = tuple(from_backend(x,E) for x in arg)
    elif isinstance(arg, dict):
        ans = { from_backend(key,E): from_backend(value,E) for (key, value) in arg.items() }
    elif isinstance(arg, Parent):
        ans = from_backend_parent(arg,E)
    elif isinstance(arg, Map):
        ans = from_backend_morphism(arg,E)
    elif isinstance(arg, Element):
        ans = from_backend_element(arg,E)
    if ans is None:
        return arg
    else:
        return ans
