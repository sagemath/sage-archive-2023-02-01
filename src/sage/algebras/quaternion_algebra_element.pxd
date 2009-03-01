include "../ext/cdefs.pxi"


import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element
from sage.categories.morphism cimport Morphism

cdef class QuaternionAlgebraElement_abstract(AlgebraElement):
    pass

cdef class QuaternionAlgebraElement_generic(QuaternionAlgebraElement_abstract):
    cdef object x, y, z, w, d, a, b
        # we will assume that our element has the representation
        # x + yi + zj + wk, where i^2 = a, j^2 = b

    pass

cdef class QuaternionAlgebraElement_number_field(QuaternionAlgebraElement_abstract):
    pass

cdef class QuaternionAlgebraElement_rational_field(QuaternionAlgebraElement_abstract):
    cdef mpz_t x, y, z, w, a, b, d
