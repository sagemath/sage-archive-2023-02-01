include "../../libs/ntl/decl.pxi"

cimport sage.structure.element
from sage.structure.element cimport Element, RingElement
from sage.rings.morphism cimport RingHomomorphism_coercion
from sage.categories.morphism cimport Section

cdef class Section_ZpFM_ZZ(Section):
    pass
cdef class Morphism_ZZ_ZpFM(RingHomomorphism_coercion):
    pass
cdef class Section_ZpCA_ZZ(Section):
    pass
cdef class Morphism_ZZ_ZpCA(RingHomomorphism_coercion):
    pass
cdef class Section_ZpCR_ZZ(Section):
    pass
cdef class Morphism_ZZ_ZpCR(RingHomomorphism_coercion):
    pass
cdef class Section_QpCR_ZZ(Section):
    pass
cdef class Morphism_ZZ_QpCR(RingHomomorphism_coercion):
    pass
cdef class Section_QpCR_QQ(Section):
    pass
cdef class Morphism_QQ_QpCR(RingHomomorphism_coercion):
    pass
cdef class Section_QpCR_ZpCR(Section):
    pass
cdef class Morphism_ZpCR_QpCR(RingHomomorphism_coercion):
    pass
cdef class Section_ZpCA_ZpCR(Section):
    pass
cdef class Morphism_ZpCR_ZpCA(RingHomomorphism_coercion):
    pass
cdef class Section_QpCR_ZpCA(Section):
    pass
cdef class Morphism_ZpCA_QpCR(RingHomomorphism_coercion):
    pass
cdef class Section_ZpFM_ZpCA(Section):
    pass
cdef class Morphism_ZpCA_ZpFM(RingHomomorphism_coercion):
    pass
cdef class Morphism_ZpCR_ZpFM(RingHomomorphism_coercion):
    pass
## Extension types
cdef class Morphism_ZZ_EisFM(RingHomomorphism_coercion):
    pass
cdef class Morphism_ZZ_UnrFM(RingHomomorphism_coercion):
    pass
cdef class Morphism_ZpFM_EisFM(RingHomomorphism_coercion):
    pass
cdef class Morphism_ZpFM_UnrFM(RingHomomorphism_coercion):
    pass
