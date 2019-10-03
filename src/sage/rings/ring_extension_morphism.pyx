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
from sage.structure.richcmp import op_EQ, op_NE

from sage.structure.element cimport Element
from sage.categories.map import Map
from sage.rings.ring cimport CommutativeRing
from sage.rings.morphism cimport RingMap
from sage.rings.ring_extension cimport RingExtension_generic, RingExtensionWithBasis
from sage.rings.ring_extension_element cimport RingExtensionElement
from sage.rings.ring_extension_conversion cimport backend_parent, backend_element, backend_morphism


# I don't trust the operator ==
cdef are_equal_morphisms(f, g):
    r"""
    Return ``True`` if ``f`` and ``g`` coincide on the
    generators of the domain, ``False`` otherwise.

    INPUT:

    - ``f`` - a ring homomorphism or ``None``; if ``None``,
      we consider that ``f`` is a coercion map

    - ``g`` - a ring homomorphism or ``None``

    TESTS::

        sage: S.<x> = QQ[]
        sage: T.<y> = S[]
        sage: TT = T.over(QQ)
        sage: H = End(TT)

        sage: cc = S.hom([-x])
        sage: f = T.hom([x^2 + y^2], base_map=cc)
        sage: g = T.hom([x^2 + y^2])

        sage: H(f) == H(g)      # indirect doctest
        False
        sage: H(f^2) == H(g^2)  # indirect doctest
        True
    """
    cdef CommutativeRing b
    cdef tuple gens
    if f is None and g is None:
        return True
    if f is None:
        f, g = g, f
    gens = tuple()
    b = f.domain()
    while b is not b._base:
        gens += b.gens()
        b = b._base
    if g is None:
        for x in gens:
            if f(x) != x: return False
    else:
        for x in gens:
            if f(x) != g(x): return False
    return True


cdef class RingExtensionHomomorphism(RingMap):
    r"""
    A class for ring homomorphisms between extensions.

    EXAMPLES:

        sage: K.<a> = GF(5^2).over()
        sage: L.<b> = GF(5^4).over(K)
        sage: phi = L.hom([b^5, a^5])
        sage: phi
        Ring endomorphism of Field in b with defining polynomial x^2 + (3 - a)*x + a over its base
          Defn: b |--> (2 + a) + 2*b
                with map on base ring:
                a |--> 1 - a

        sage: TestSuite(phi).run()

    """
    def __init__(self, parent, defn, base_map=None, check=True):
        r"""
        Initialize this morphism.
        """
        RingMap.__init__(self, parent)
        domain = self.domain()
        backend_domain = backend_parent(domain)
        codomain = self.codomain()
        backend_codomain = backend_parent(codomain)
        # We construct the backend morphism
        if isinstance(defn, Map):
            if base_map is not None:
                raise ValueError("base_map cannot be set when passing in the backend morphism")
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
                    current_im_gens.append(backend_element(y))
                current_morphism = current_domain.hom(current_im_gens, base_map=current_morphism, check=check)
            # We check that everything went well
            if check:
                for i in range(len(gens)):
                    x = backend_element(domain(gens[i]))
                    y = backend_element(im_gens[i])
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

    def _repr_type(self):
        return "Ring"

    cpdef Element _call_(self, x):
        y = self._backend(backend_element(x))
        if isinstance(self.codomain(), RingExtension_generic):
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
        eq = are_equal_morphisms(self._backend, backend_morphism(other))
        if op == op_EQ:
            return eq
        if op == op_NE:
            return not eq
        raise NotImplemented

    def is_identity(self):
        return are_equal_morphisms(self._backend, None)

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
        backend = self._backend * backend_right
        if isinstance(domain, RingExtension_generic) or isinstance(codomain, RingExtension_generic):
            return RingExtensionHomomorphism(domain.Hom(codomain), backend)
        else:
            return backend

    cdef _update_slots(self, dict _slots):
        self._backend = _slots['_backend']
        RingMap._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        slots = RingMap._extra_slots(self)
        slots['_backend'] = self._backend
        return slots


cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    def __init__(self, parent):
        RingMap.__init__(self, parent)
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
        RingMap.__init__(self, parent)
        codomain = self.codomain()
        self._backend = codomain.Hom(codomain).identity()

    def _repr_type(self):
        return "Canonical"

    def _repr_defn(self):
        return ""

    cpdef Element _call_(self, x):
        return (<RingExtensionElement>x)._backend


cdef class MapVectorSpaceToRelativeField(Map):
    def __init__(self, E, K):
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

    cpdef Element _call_(self, v):
        cdef Element elt
        elt = self._f(v[0]) * self._basis[0]
        for i in range(1, self._degree):
            elt += self._f(v[i]) * self._basis[i]
        return self.codomain()(elt)


cdef class MapRelativeFieldToVectorSpace(Map):
    def __init__(self, E, K):
        cdef CommutativeRing L, base

        self._degree = (<RingExtensionWithBasis>E)._degree_over(K)
        self._basis = [ (<RingExtensionElement>x)._backend for x in E.basis_over(K) ]
        f = backend_morphism(E.defining_morphism(K), forget="codomain")
        codomain = K ** self._degree
        Map.__init__(self, E.Hom(codomain))

        K = backend_parent(K)
        L = (<RingExtensionWithBasis>E)._backend

        # We compute the matrix of our isomorphism (over base)
        from sage.rings.ring_extension import common_base
        base = common_base(K, L, False)
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

    cpdef Element _call_(self, x):
        coeffs = self.backend_coefficients(x)
        return self.codomain()(coeffs)

    cdef list backend_coefficients(self, RingExtensionElement x):
        dK = self._dimK
        w = (self._jL(x._backend) * self._matrix).list()
        coeffs = [ ]
        for i in range(self._degree):
            coeff = self._iK(w[i*dK:(i+1)*dK])
            coeffs.append(coeff)
        return coeffs
