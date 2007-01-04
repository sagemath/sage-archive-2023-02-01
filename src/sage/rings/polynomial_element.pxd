import sage.rings.commutative_algebra_element as commutative_algebra_element
cimport sage.rings.commutative_algebra_element as commutative_algebra_element


cdef class Polynomial(commutative_algebra_element.CommutativeAlgebraElement):
    cdef ssize_t degree
    cdef char _is_gen

cdef class Polynomial_generic_dense(Polynomial):
    cdef object __coeffs # a python list

cdef class Polynomial_generic_sparse(Polynomial):
    cdef object __coeffs # a python dict (for now)
