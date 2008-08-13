

from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.rings.ring cimport Ring

cdef class RingMap(Morphism):
    pass

cdef class RingMap_lift(RingMap):
    cdef Ring S

cdef class RingHomomorphism(RingMap):
    cdef RingMap _lift

cdef class RingHomomorphism_coercion(RingHomomorphism):
    pass

cdef class RingHomomorphism_im_gens(RingHomomorphism):
    cdef __im_gens

cdef class RingHomomorphism_cover(RingHomomorphism):
    pass

cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    cdef phi
