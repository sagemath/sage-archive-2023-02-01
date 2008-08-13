from sage.rings.morphism cimport RingHomomorphism
from sage.structure.element cimport Element
from sage.rings.integer_mod cimport NativeIntStruct

cdef class ResidueFieldHomomorphism(RingHomomorphism):
    pass

cdef class NFResidueFieldHomomorphism(ResidueFieldHomomorphism):
    cdef object im_gen
    cdef object p
