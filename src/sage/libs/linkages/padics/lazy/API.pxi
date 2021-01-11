# Operations on digits (intended to be small elements in the exact subring)

cdef inline void digit_init(cdigit a):
    pass

cdef inline void digit_clear(cdigit a):
    pass

cdef inline void digit_set(cdigit a, cdigit b):
    pass

cdef inline void digit_set_ui(cdigit a, slong b):
    pass

cdef inline bint digit_equal(cdigit a, cdigit b):
    pass

cdef inline bint digit_equal_ui(cdigit a, slong b):
    pass

cdef inline void digit_set_sage(cdigit a, elt):
    pass

cdef inline digit_get_sage(cdigit a):
    pass

cdef inline void digit_add(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_sub(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_mul(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_mod(cdigit res, cdigit a, cdigit modulus):
    pass

cdef inline void digit_quorem(cdigit quo, cdigit rem, cdigit a, cdigit modulus):
    pass


# Operations on elements (represented as series of digits)

cdef inline void element_init(cdigit a):
    pass

cdef inline void element_clear(cdigit a):
    pass

cdef inline void element_set_coeff(celement x, cdigit a, slong i):
    pass

cdef inline void cdigit_ptr element_get_coeff(celement x, slong i):
    pass

cdef inline void cdigit_ptr element_get_slice(celement x, slong start, slong length):
    pass

cdef inline void element_iadd_coeff(celement x, cdigit a, slong i):
    pass

cdef inline void element_isub_coeff(celement x, cdigit a, slong i):
    pass

cdef inline void element_iadd_slice(celement x, celement slice, slong start):
    pass

cdef inline void element_isub_slice(celement x, celement slice, slong start):
    pass

cdef inline void element_reduce_coeff(celement x, slong i):
    pass
