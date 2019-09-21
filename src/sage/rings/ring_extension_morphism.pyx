#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.element cimport Element
from sage.categories.map import Map
from sage.rings.morphism cimport RingHomomorphism
from sage.rings.ring_extension import RingExtension_class, RingExtensionWithBasis
from sage.rings.ring_extension import _common_base
from sage.rings.ring_extension_element cimport RingExtensionElement


# Helper functions

cpdef _backend_morphism(f):
    from sage.categories.map import FormalCompositeMap
    from sage.categories.morphism import IdentityMorphism
    domain = f.domain()
    if not isinstance(f.domain(), RingExtension_class) and not isinstance(f.codomain(), RingExtension_class):
        return f
    elif isinstance(f, RingExtensionHomomorphism):
        return (<RingExtensionHomomorphism>f)._backend
    elif isinstance(f, FormalCompositeMap):
        return _backend_morphism(f.then()) * _backend_morphism(f.first())
    elif isinstance(f, IdentityMorphism):
        ring = domain._backend
        return ring.Hom(ring).identity()
    elif domain is domain.base_ring():
        ring = f.codomain()
        if isinstance(ring, RingExtension_class):
            ring = ring._backend
        if ring.has_coerce_map_from(domain):
            return ring.coerce_map_from(domain)
    raise NotImplementedError

def backend_morphism(f, forget="all"):
    try:
        g = _backend_morphism(f)
        if forget is None and (isinstance(f.domain(), RingExtension_class) or isinstance(f.codomain(), RingExtension_class)):
            g = RingExtensionHomomorphism(f.domain().Hom(f.codomain()), g)
        if forget == "domain" and isinstance(f.codomain(), RingExtension_class):
            g = RingExtensionHomomorphism(g.domain().Hom(f.codomain()), g)
        if forget == "codomain" and isinstance(f.domain(), RingExtension_class):
            g = RingExtensionHomomorphism(f.domain().Hom(g.codomain()), g)
    except NotImplementedError:
        g = f
        if (forget == "all" or forget == "domain") and isinstance(f.domain(), RingExtension_class):
            ring = f.domain()._backend
            g = g * RingExtensionBackendIsomorphism(ring.Hom(f.domain()))
        if (forget == "all" or forget == "codomain") and isinstance(f.codomain(), RingExtension_class):
            ring = f.codomain()._backend
            g = RingExtensionBackendReverseIsomorphism(f.codomain().Hom(ring)) * g
    return g

# I don't trust the operator ==
def are_equal_morphisms(f, g):
    if f is None and g is None:
        return True
    if f is None:
        f, g = g, f
    if g is None:
        for x in f.domain().gens():
            if f(x) != x: return False
    else:
        for x in f.domain().gens():
            if f(x) != g(x): return False
    return True


# Classes

cdef class RingExtensionHomomorphism(RingHomomorphism):
    r"""
    Homomorphisms between extensions

    EXAMPLES:

        sage: F = GF(5^2)
        sage: K = GF(5^4)
        sage: L = GF(5^8)
        sage: E1 = RingExtension(K,F)
        sage: E2 = RingExtension(L,K)

    """
    def __init__(self, parent, defn, base_map=None, check=True):
        RingHomomorphism.__init__(self, parent)
        backend_domain = domain = self.domain()
        if isinstance(backend_domain, RingExtension_class):
            backend_domain = backend_domain._backend
        backend_codomain = codomain = self.codomain()
        if isinstance(backend_codomain, RingExtension_class):
            backend_codomain = backend_codomain._backend
        # We construct the backend morphism
        if isinstance(defn, Map):
            if base_map is not None:
                raise ValueError("base_map must not be set when providing the backend morphism")
            backend = backend_morphism(defn)
            if backend.domain() is not backend_domain:
                raise TypeError("the domain of the backend morphism is not correct")
            if backend.codomain() is not backend_codomain:
                raise TypeError("the codomain of the backend morphism is not correct")
            self._backend = backend
            self._im_gens = None
            self._base_map_construction = False
        elif isinstance(defn, (list, tuple)):
            # We figure out what is the base
            if base_map is not None:
                base = base_map.domain()
                gens = domain.gens(base)
            else:
                base = domain
                gens = tuple([])
                while True:
                    if len(gens) == len(defn): 
                        break  
                    if len(gens) > len(defn) or base is base.base_ring():
                        raise ValueError("the number of images does not match the number of generators")
                    gens += base.gens()
                    base = base.base_ring()
            # We construct the backend morphism
            im_gens = [ codomain(x) for x in defn ]
            backend_bases = [ backend_domain ]
            b = backend_domain.base_ring()
            while b is not b.base_ring():
                backend_bases.append(b)
                b = b.base_ring()
            backend_bases.reverse()
            current_morphism = None
            for current_domain in backend_bases:
                current_im_gens = [ ]
                for x in current_domain.gens():
                    pol = domain(backend_domain(x)).polynomial(base)
                    if base_map is not None:
                        pol = pol.map_coefficients(base_map)
                    y = pol(im_gens)
                    if isinstance(codomain, RingExtension_class):
                        y = (<RingExtensionElement>y)._backend
                    current_im_gens.append(y)
                current_morphism = current_domain.hom(current_im_gens, base_map=current_morphism, check=check)
            # We check that everything went well
            if check:
                for i in range(len(gens)):
                    x = domain(gens[i])
                    if isinstance(domain, RingExtension_class):
                        x = (<RingExtensionElement>x)._backend
                    y = im_gens[i]
                    if isinstance(codomain, RingExtension_class):
                        y = (<RingExtensionElement>y)._backend
                    if current_morphism(x) != y:
                        raise ValueError("images do not define a valid homomorphism")
                coercion_morphism = backend_morphism(domain.defining_morphism(base))
                if base_map is None:
                    backend_base_map = coercion_morphism
                else:
                    backend_base_map = backend_morphism(base_map)
                restriction_current_morphism = current_morphism * coercion_morphism
                if not are_equal_morphisms(restriction_current_morphism, backend_base_map):
                    raise ValueError("images do not define a valid homomorphism")
            self._backend = current_morphism
            self._im_gens = im_gens[:domain.ngens()]
            if base is domain.base_ring():
                self._base_map_construction = base_map
            else:
                self._base_map_construction = { 
                    'im_gens': defn[domain.ngens():], 
                    'base_map': base_map, 
                    'check': False
                }
        else:
            raise TypeError

    cpdef Element _call_(self, x):
        if isinstance(self.domain(), RingExtension_class):
            x = (<RingExtensionElement>x)._backend
        y = self._backend(x)
        if isinstance(self.codomain(), RingExtension_class):
            y = self._codomain(y)
        return y

    @cached_method
    def base_map(self):
        domain = self.domain()
        codomain = self.codomain()
        base = domain.base_ring()
        if base is base.base_ring():
            return None
        base_map = self._base_map_construction
        if base_map is False:
            if base is domain:
                base_map = None
            else:
                base_map = self * domain.coerce_map_from(base)
        elif isinstance(base_map, dict):
            base_map = base.hom(**self._base_map_construction)
        if base_map is None:
            return None
        if (codomain.has_coerce_map_from(base) 
        and are_equal_morphisms(backend_morphism(base_map), backend_morphism(codomain.coerce_map_from(base)))):
            return None
        if base_map.codomain() is not self.codomain():
            base_map = base_map.extend_codomain(self.codomain())
        return base_map

    cpdef _richcmp_(self, other, int op):
        return self._backend._richcmp_(backend_morphism(other), op)

    def is_identity(self):
        if self.domain() is not self.codomain():
            return False
        return self._backend.is_identity()

    def is_injective(self):
        return self._backend.is_injective()

    def is_surjective(self):
        return self._backend.is_surjective()

    def _repr_defn(self):
        import re
        s = ""
        gens = self.domain().gens()
        if self._im_gens is None:
            self._im_gens = [ self(x) for x in gens ]
        for i in range(len(gens)):
            s += "%s |--> %s\n" % (gens[i], self._im_gens[i])
        if self.base_map() is not None:
            s += "with map on base ring"
            ss = self.base_map()._repr_defn()
            ss = re.sub('\nwith map on base ring:?$', '', ss, 0, re.MULTILINE)
            if ss != "": s += ":\n" + ss
        if s != "" and s[-1] == "\n":
            s = s[:-1]
        return s

    def _composition(self, right):
        domain = right.domain()
        codomain = self.codomain()
        backend_right = backend_morphism(right)
        if self._backend is None:
            print(self)
        backend = self._backend * backend_right
        if isinstance(domain, RingExtension_class) or isinstance(codomain, RingExtension_class):
            return RingExtensionHomomorphism(domain.Hom(codomain), backend)
        else:
            return backend

    cdef _update_slots(self, dict _slots):
        self._backend = _slots['_backend']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        slots = RingHomomorphism._extra_slots(self)
        slots['_backend'] = self._backend
        return slots


cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    def __init__(self, parent):
        RingHomomorphism.__init__(self, parent)
        domain = self.domain()
        self._backend = domain.Hom(domain).identity()

    def _repr_type(self):
        return "Coercion"

    def _repr_defn(self):
        return ""

    cpdef Element _call_(self, x):
        codomain = self.codomain()
        return codomain.element_class(codomain, x)


cdef class RingExtensionBackendReverseIsomorphism(RingExtensionHomomorphism):
    def __init__(self, parent):
        RingHomomorphism.__init__(self, parent)
        codomain = self.codomain()
        self._backend = codomain.Hom(codomain).identity()

    def _repr_type(self):
        return "Canonical"

    def _repr_defn(self):
        return ""

    cpdef Element _call_(self, x):
        return (<RingExtensionElement>x)._backend


class MapVectorSpaceToRelativeField(Map):
    def __init__(self, E, K):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.degree(K)
        self._basis = [ (<RingExtensionElement>x)._backend for x in E.basis_over(K) ]
        self._f = backend_morphism(E.defining_morphism(K), forget="codomain")
        domain = K ** self._degree
        parent = domain.Hom(E)
        Map.__init__(self, parent)

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

    #cpdef Element _call_(self, v):
    def _call_(self, v):
        elt = sum(self._f(v[i]) * self._basis[i] for i in range(self._degree))
        return self.codomain()(elt)


class MapRelativeFieldToVectorSpace(Map):
    def __init__(self, E, K):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.degree(K)
        self._basis = [ (<RingExtensionElement>x)._backend for x in E.basis_over(K) ]
        f = backend_morphism(E.defining_morphism(K), forget="codomain")
        codomain = K ** self._degree
        parent = E.Hom(codomain)
        Map.__init__(self, parent)

        if isinstance(K, RingExtension_class):
            K = K._backend
        L = E._backend

        # We compute the matrix of our isomorphism (over base)
        base = _common_base(K,L)
        EK, iK, jK = K.free_module(base, map=True)
        EL, iL, jL = L.free_module(base, map=True)

        self._dimK = EK.dimension()
        self._iK = iK
        self._jL = jL

        M = [ ]
        for x in self._basis:
            for v in EK.basis():
                y = x * f(iK(v))
                M.append(jL(y))
        from sage.matrix.matrix_space import MatrixSpace
        self._matrix = MatrixSpace(base,len(M))(M).inverse_of_unit()

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

    #cpdef Element _call_(self, x):
    def _call_(self, x):
        coeffs = self.backend_coefficients(x)
        return self.codomain()(coeffs)

    def backend_coefficients(self, x):
        dK = self._dimK
        w = (self._jL((<RingExtensionElement>x)._backend) * self._matrix).list()
        coeffs = [ ]
        for i in range(self._degree):
            coeff = self._iK(w[i*dK:(i+1)*dK])
            coeffs.append(coeff)
        return coeffs
