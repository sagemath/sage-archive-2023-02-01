

from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.rings.ring cimport Ring

cdef class RingMap(Morphism):
    pass

cdef class RingMap_lift(RingMap):
    cdef Ring S
    cdef Element _call_c_impl(self, Element x)
    cdef _update_slots(self, _slots)
    cdef _extra_slots(self, _slots)

cdef class RingHomomorphism(RingMap):
    cdef RingMap _lift
    cdef _update_slots(self, _slots)
    cdef _extra_slots(self, _slots)

cdef class RingHomomorphism_coercion(RingHomomorphism):
    cdef Element _call_c_impl(self, Element x)

cdef class RingHomomorphism_im_gens(RingHomomorphism):
    cdef __im_gens
    cdef Element _call_c_impl(self, Element x)
    cdef _update_slots(self, _slots)
    cdef _extra_slots(self, _slots)

cdef class RingHomomorphism_cover(RingHomomorphism):
    cdef Element _call_c_impl(self, Element x)

cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    cdef phi
    cdef Element _call_c_impl(self, Element x)
    cdef _update_slots(self, _slots)
    cdef _extra_slots(self, _slots)