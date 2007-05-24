from sage.misc.sagex_ds import BinaryTree
from sage.misc.sagex_ds cimport BinaryTree
include "../../ext/interrupt.pxi"

cdef enum:
    pdUNIVAR
    pdVAR
    pdCOEFF
    pdSQR
    pdADD
    pdMUL
    pdPOW
    pdABC


cdef class generic_pd

cdef class CompiledPolynomialFunction:
    cdef generic_pd _dag
    cdef object _coeffs

    cdef object _parse_structure(CompiledPolynomialFunction)
    cdef generic_pd _get_gap(CompiledPolynomialFunction, BinaryTree, int)
    cdef void _fill_gaps_quick(CompiledPolynomialFunction, BinaryTree)
    cdef object eval(CompiledPolynomialFunction, object)

cdef class generic_pd:
    cdef object value
    cdef int refs, hits
    cdef int label
    cdef void eval(generic_pd, object, object)
    cdef generic_pd nodummies(generic_pd)

cdef class dummy_pd(generic_pd):
    cdef generic_pd link
    cdef void set_mul(dummy_pd, generic_pd, generic_pd)
    cdef void set_sqr(dummy_pd, generic_pd)
    cdef generic_pd nodummies(dummy_pd)

cdef class var_pd(generic_pd):
    cdef int index
    cdef void eval(var_pd, object, object)

cdef class univar_pd(generic_pd):
    cdef void eval(univar_pd, object, object)

cdef class coeff_pd(generic_pd):
    cdef int index
    cdef void eval(con_pd, object, object)


cdef class unary_pd(generic_pd):
    cdef generic_pd operand
    cdef generic_pd nodummies(unary_pd)

cdef class sqr_pd(unary_pd):
    cdef void eval(sqr_pd, object, object)

cdef class pow_pd(unary_pd):
    cdef object exponent
    cdef void eval(pow_pd, object, object)


cdef class binary_pd(generic_pd):
    cdef generic_pd left, right
    cdef generic_pd nodummies(binary_pd)

cdef class add_pd(binary_pd):
    cdef void eval(add_pd, object, object)

cdef class mul_pd(binary_pd):
    cdef void eval(mul_pd, object, object)

cdef class abc_pd(binary_pd):
    cdef int index
    cdef void eval(abc_pd, object, object)






